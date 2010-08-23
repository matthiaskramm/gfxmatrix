/* txtsynth.h

   Texture synthesis.

   Copyright (c) 2007,2008,2009,2010 Matthias Kramm <kramm@quiss.org>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */

#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <math.h>
#include <assert.h>
#include "matrix.h"
#include "image.h"
#include "stat.h"
#include "corr.h"
#include "txtsynth.h"

/*
   Features in S&P, but not in this model so far:

   * Skew and Kurtosis across scales
   * E(Re(x*x*x)) for subbands between octaves (=correlation between subbands)
   * E(Re(x*x*x)) for subbands between directions (=correlation between scales)

*/

static inline double clamp01(double lr)
{
    if(lr>1) lr = 1;
    if(lr<0) lr = 0;
    return lr;
}

static inline double sqr(double lr)
{
    return lr*lr;
}

static int config_usepowerspectrum=1;
static int config_useautocorrelation=0;
static int config_spyr=0;

bytearrayset_t* make_masks2(int scales, int width, int height)
{
    char filename[80];
    bytearrayset_t*masks = bytearrayset_new_alloc((scales-1)*3+1,width,height);
    int t;
    int pos = 0;
    for(t=0;t<scales;t++) {
        int dir=0;
        int s=0;
        if(t==scales-1)
            s = 2; //highband
        for(dir=s;dir<=2;dir++) {
            int x,y;
            int mx=width/2,my=height/2;
            
            //double r = (1<<(t+1+highband-(t==scales-1 && highband)))/(double)(2<<scales);
            double r = (t+1-(t==scales-1))/(double)(scales-1);
            double pr = t/(double)(scales-1);

            memset(masks->m[pos]->data, 0, width*height);
            if(t==scales-1 && dir==2) {
                memset(masks->m[pos]->data, 255, width*height);
            }
            for(y=0;y<my;y++)
            for(x=0;x<mx;x++) {
                double xx = x/(double)mx;
                double yy = y/(double)my;
                double rr = sqrt(xx*xx+yy*yy);
                double lr = log(rr) / log(2) + (scales-t);
                double m = cos(clamp01(lr) * M_PI/2);
                unsigned char dot = clamp01(m)*255;
                dot = (rr*1.1)<r?255:0;
                if(dir<2) {
                    if(rr*1.1<pr)
                        dot = 0;
                }
                if(t==scales-1 && dir==2) {
                    dot = (rr*1.1)>=r?255:0;
                }
                if(dir == 0 || dir==2) 
                    masks->m[pos]->data[(height-y)%height*masks->m[t]->width+(width-x)%width] = dot;
                if(dir == 0 || dir==2) 
                    masks->m[pos]->data[(y)*masks->m[t]->width+(x)] = dot;
                if(dir == 1 || dir==2) 
                    masks->m[pos]->data[(y)*masks->m[t]->width+(width-x)%width] = dot;
                if(dir == 1 || dir==2) 
                    masks->m[pos]->data[(height-y)%height*masks->m[t]->width+(x)] = dot;
            }
            pos++;
        }
    }
    return masks;
}

static int mask_sizediv(int t, int scales)
{
    int size[] = {16,8,4,2,1,1};
    int sz = (t==scales-1)?1:(1<<(scales-2-t));

    //printf("scales=%d %d %d\n", scales, sz, size[t]);
    //assert(scales=6);
    //assert(sz == size[t]); // for testing, remove for scales != 6

    return sz;
}

bytearrayset_t* make_masks(int scales, int width, int height)
{
    char highband=1;
    char filename[80];
    bytearrayset_t*masks = bytearrayset_new_alloc(scales,width,height);
    int t;
    highband = !!highband;

    for(t=0;t<scales;t++) {
        int x,y;
        int mx=width/2,my=height/2;
        
        double r = (1<<(t+1+highband-(t==scales-1 && highband)))/(double)(1<<scales);
        //double r = (t+1-(t==scales-1 && highband))/(double)(scales-1);

        memset(masks->m[t]->data, 0, width*height);
        if(t==scales-1 && highband)
            memset(masks->m[t]->data, 255, width*height);
        int highx = 0;
        for(y=0;y<my;y++)
        for(x=0;x<mx;x++) {
            double xx = x/(double)mx;
            double yy = y/(double)my;
            double rr = sqrt(xx*xx+yy*yy);
            unsigned char dot = 0;
            if(config_spyr) { 
                double lr = log(rr) / log(2) + (scales-t-highband+(t==scales-1 && highband));
                double m = cos(clamp01(lr) * M_PI/2);
                dot = clamp01(m)*255;
                if(t==scales-1 && highband) {
                    dot = 255-dot;
                }
            } else {
                dot = (rr*1.1)<r?255:0;
                if(t==scales-1 && highband) {
                    dot = (rr*1.1)>=r?255:0;
                }
            }
            if(dot && x>highx)
                highx = x;
            masks->m[t]->data[(height-y)%height*masks->m[t]->width+(width-x)%width] = dot;
            masks->m[t]->data[(y)*masks->m[t]->width+(x)] = dot;
            masks->m[t]->data[(y)*masks->m[t]->width+(width-x)%width] = dot;
            masks->m[t]->data[(height-y)%height*masks->m[t]->width+(x)] = dot;
        }
        //printf("%d %d/%d %d\n", t, highx, mx);
    }
    return masks;
}

#define CACHE_WIDTH 513
#define CACHE_SCALES 8
static bytearrayset_t**mymasks = 0;

static bytearrayset_t* get_masks(int scales, int width, int height)
{
    if(!mymasks) {
        mymasks = malloc(sizeof(bytearrayset_t*)*CACHE_WIDTH*CACHE_SCALES);
        memset(mymasks, 0, sizeof(bytearrayset_t*)*CACHE_WIDTH*CACHE_SCALES);
    }
    if(width<CACHE_WIDTH && scales<CACHE_SCALES && mymasks[scales*CACHE_WIDTH+width])
        return mymasks[scales*CACHE_WIDTH+width];
    bytearrayset_t*masks = make_masks(scales, width, height);
    if(width<CACHE_WIDTH && scales<CACHE_SCALES) 
        mymasks[scales*CACHE_WIDTH+width] = masks;
    return masks;
}

static void dealloc_all()
{
    if(mymasks) {
        int t;
        for(t=0;t<CACHE_WIDTH*CACHE_SCALES;t++) {
            if(mymasks[t])
                bytearrayset_delete(mymasks[t]);
            mymasks=0;
        }
        free(mymasks);
    }
}

int getmean(image_t*img,matrix_t*m)
{
    // statistics_t*st_u = statistics_new_frommatrix(dst_u);clamp00ff(st_u->mean);
    // statistics_t*st_v = statistics_new_frommatrix(dst_v);clamp00ff(st_v->mean);

    int size = img->width*img->height;
    assert(img->width == m->width);
    assert(img->height == m->height);
   
    int t;
    int num = 0;
    double mean = 0;
    for(t=0;t<size;t++) {
        int lum = img->data[t].r+img->data[t].g+img->data[t].b;
        if(lum > 32*3 && lum < 224*3) {
            // only use color samples which are neither to black or too white
            num++;
            mean += m->data[t];
        }
    }
    if(num)
        mean /= num;
    if(mean<0) mean = 0;
    if(mean>255) mean = 255;
    return mean;
}

image_t*synth_simple(image_t*alpha, image_t*img)
{
    matrix_t*dst = image_extractchannel(img, IMAGE_YUV_Y);
    statistics_t*st = statistics_new_frommatrix(dst);
    complex_matrix_t*c = matrix_fft(dst);
    matrix_t*m = matrix_new_gaussrandom(dst->width, dst->height, 0, 1);
    int it;
    for(it=0;it<10;it++) {
        complex_matrix_t*s = matrix_fft(m);
        int t;
        for(t=0;t<s->width*s->height;t++) {
            double r = sqrt(c->data[t].real*c->data[t].real + c->data[t].imag*c->data[t].imag);
            double r2 = sqrt(s->data[t].real*s->data[t].real + s->data[t].imag*s->data[t].imag);
            if(r2) {
                s->data[t].real *= r / r2 / (s->width*s->height);
                s->data[t].imag *= r / r2 / (s->width*s->height);
            }
            //s->data[t].imag = 0;
        }
        matrix_delete(m);
        m = complex_matrix_ifft_real(s);
        matrix_adjust_statistics(m, st);
        matrix_adjust_minmax(m, st);
        complex_matrix_delete(s);
    }
    
    //image_t*result = image_from_matrix(m,IMAGE_CLAMP);

    matrix_t*dst_u = image_extractchannel(img, IMAGE_YUV_U);
    matrix_t*dst_v = image_extractchannel(img, IMAGE_YUV_V);
    int u = getmean(img, dst_u);
    int v = getmean(img, dst_v);

    image_t*result = image_new(img->width, img->height);
    int t;
    for(t=0;t<img->width*img->height;t++) {
	int yy = m->data[t];
        //result->data[t].r = clamp00ff(yy);
        //result->data[t].g = clamp00ff(yy);
        //result->data[t].b = clamp00ff(yy);
        //u = 128;
        //v = 128;
        result->data[t].r = clamp00ff(yy + ((360*(v-128))>>8));
        result->data[t].g = clamp00ff(yy - ((88*(u-128)+183*(v-128))>>8));
        result->data[t].b = clamp00ff(yy + ((455 * (u-128))>>8));
        result->data[t].a = 255;
    }
    matrix_delete(m);
    matrix_delete(dst);
    matrix_delete(dst_u);
    matrix_delete(dst_v);
    complex_matrix_delete(c);
    statistics_delete(st);
    return result;
}

typedef struct _texturedata_internal1 {
    int u;
    int v;
    statistics_t*st;
    matrix_t*powerspectrum;
} texturedata_internal1_t;

typedef struct _texturedata_internal2 {
    int u;
    int v;
    statistics_t*st;
    int scales;
    int num_masks;
    
    statistics_t**scale_st;
    matrix_t**scale_autocorr;

    matrix_t*powerspectrum;
} texturedata_internal2_t;

static texturedata_internal2_t* analyze_texture12(image_t*image)
{
    texturedata_internal2_t*data = malloc(sizeof(texturedata_internal2_t));
    memset(data, 0, sizeof(texturedata_internal2_t));
    
    matrix_t*texture = image_extractchannel(image, IMAGE_YUV_Y);

    /* get mean u/v values */
    int r=0,g=0,b=0;
    int size = image->width*image->height;
    int t;
    for(t=0;t<size;t++) {
        r += image->data[t].r;
        g += image->data[t].g;
        b += image->data[t].b;
    }
    r/=size;
    g/=size;
    b/=size;
    data->u = (r*((int)(-0.169*256)) + g*((int)(-0.332*256)) + b*((int)( 0.500 *256))+ 128*256)>>8;
    data->v = (r*((int)( 0.500*256)) + g*((int)(-0.419*256)) + b*((int)(-0.0813*256))+ 128*256)>>8;

    /* calculate number of scales needed */
    int wx=image->width,wy=image->height;
    int scales = 1;
    while(1<<(scales+1) < wx && 1<<(scales+1) < wy)
        scales++;
    data->scales = scales;

    /* get first-order statistic of full texture */
    data->st = statistics_new_frommatrix(texture);
  
    /* get power spectrum */
    complex_matrix_t*c = matrix_fft(texture);
    data->powerspectrum = complex_matrix_abs(c);
   
    /* make subband statistics */
    bytearrayset_t*txtmasks = get_masks(scales,texture->width,texture->height);
    
    data->scale_st = malloc(sizeof(statistics_t*)*txtmasks->num);
    memset(data->scale_st, 0, sizeof(statistics_t*)*txtmasks->num);
    data->scale_autocorr = malloc(sizeof(matrix_t*)*txtmasks->num);
    memset(data->scale_autocorr, 0, sizeof(matrix_t*)*txtmasks->num);
    char filename[128];
    for(t=0;t<txtmasks->num;t++) {
        complex_matrix_t*part = complex_matrix_clone(c);
        bytearray_t*mask = txtmasks->m[t];
        assert(mask->width == part->width);
        assert(mask->height == part->height);
        int x,y;
        double d = 1.0 / (255*part->width*part->height);
        for(x=0;x<part->width;x++)
        for(y=0;y<part->height;y++) {
            part->data[y*part->width+x].real *= mask->data[y*mask->width+x]*d;
            part->data[y*part->width+x].imag *= mask->data[y*mask->width+x]*d;
        }
        //sprintf(filename, "txtmask%02d.png", t);
        //bytearray_save(txtmasks->m[t], filename);
        //sprintf(filename, "mask%02d.png", t);
        //bytearray_save(masks->m[t], filename);

        matrix_t*r = complex_matrix_ifft_real(part);
        //sprintf(filename, "txtscale%d.png", t);
        //matrix_save(r, filename);
        //sprintf(filename, "part%d.png", t);
        //complex_matrix_save(part, filename);
        data->scale_st[t] = statistics_new_frommatrix(r);
        data->scale_st[t]->mean = 0;

        matrix_adjust_mean(r, 0);

        int wx = r->width/mask_sizediv(t, scales);
        int wy = r->height/mask_sizediv(t, scales);
        int zx = wx/2>9?9:wx/2;
        int zy = wy/2>9?9:wy/2;
        matrix_t*r2 = matrix_scaledown(r, r->width / mask_sizediv(t, data->scales), 
                                          r->height/ mask_sizediv(t, data->scales));
        data->scale_autocorr[t] = matrix_getautocorrelation(r2, zx, zy);
        matrix_delete(r2);
        //printf("%d %dx%d, %dx%d autocorrelation\n", t, wx, wy, zx, zy);

        //printf("scale %d:", t);
        //statistics_print(scale_st[t]);
        complex_matrix_delete(part);part=0;
        matrix_delete(r);r=0;
    }
    matrix_delete(texture);texture=0;
    return data;
}

static void texturedata_free2(texturedata_t*_data)
{
    assert(_data->type == 2);
    texturedata_internal2_t*data = (texturedata_internal2_t*)_data->internal;
    int t;
    for(t=0;t<data->num_masks;t++) {
        if(data->scale_st[t]) {
            free(data->scale_st[t]);data->scale_st[t]=0;
        }
        if(data->scale_autocorr[t]) {
            free(data->scale_autocorr[t]);data->scale_autocorr[t]=0;
        }
    }
    free(data->scale_st);data->scale_st=0;
}

void texturedata_set_parameter(char*key, char*value)
{
    if(!strcmp(key, "powerspectrum")) {
        config_usepowerspectrum = atoi(value);
    } else if(!strcmp(key, "autocorrelation")) {
        config_useautocorrelation = atoi(value);
    }
}


image_t*synthesize1(image_t*alpha, texturedata_t*_data)
{
    texturedata_internal2_t*data = (texturedata_internal2_t*)_data->internal;

    matrix_t*m = matrix_new_gaussrandom(data->powerspectrum->width, data->powerspectrum->height, 0, 1);
    int it;
    for(it=0;it<10;it++) {
        complex_matrix_t*s = matrix_fft(m);
        int t;
        assert(s->width == data->powerspectrum->width);
        assert(s->height == data->powerspectrum->height);
        for(t=0;t<s->width*s->height;t++) {
            double r = data->powerspectrum->data[t];
            double r2 = sqrt(s->data[t].real*s->data[t].real + s->data[t].imag*s->data[t].imag);
            if(r2) {
                double f = r / r2 / (s->width*s->height);
                s->data[t].real *= f;
                s->data[t].imag *= f;
            }
            //s->data[t].imag = 0;
        }
        matrix_delete(m);
        m = complex_matrix_ifft_real(s);
        matrix_adjust_statistics(m, data->st);
        matrix_adjust_minmax(m, data->st);
        complex_matrix_delete(s);
    }
    
    image_t*result = image_new(m->width, m->height);
    int t;
    for(t=0;t<result->width*result->height;t++) {
	int yy = m->data[t];
        result->data[t].r = clamp00ff(yy + ((360*(data->v-128))>>8));
        result->data[t].g = clamp00ff(yy - ((88*(data->u-128)+183*(data->v-128))>>8));
        result->data[t].b = clamp00ff(yy + ((455 * (data->u-128))>>8));
        result->data[t].a = 255;
    }
    matrix_delete(m);
    return result;
}

/* fill holes in image with elements from texture */
static image_t* synthesize2(image_t*image, texturedata_t*_data)
{
    assert(_data->type == 2);
    texturedata_internal2_t*data = (texturedata_internal2_t*)_data->internal;
/*    image_t*srcimg = image_load("texture_grass.png");
    matrix_t*m = image_extractchannel(srcimg, IMAGE_GRAY);*/
    image_t*image_tofree=0;
    if(!image) {
        image = image_tofree = image_new(256,256);
    }

    matrix_t*orig = image_extractchannel(image, IMAGE_YUV_Y);
    char filename[80];

    //image_save(image_from_matrix(m,IMAGE_CLAMP), "src.png");
    //image_save(image_from_matrix(texture,IMAGE_CLAMP), "texture.png");
    //image_save_and_free(image_from_complex_matrix(c,IMAGE_CLAMP), "dest_fft.png");
    
    matrix_t*m = matrix_new_gaussrandom(image->width, image->height, 0, 255);
   
    /* adapt radius information to our output size */
    matrix_t*destradius = matrix_new(m->width, m->height);
    int x,y;
    
    assert(destradius->width >= data->powerspectrum->width);
    assert(destradius->height >= data->powerspectrum->height);

    for(y=0;y<destradius->height;y++)
    for(x=0;x<destradius->width;x++) {
        int x2 = data->powerspectrum->width*x / destradius->width;
        int y2 = data->powerspectrum->height*y / destradius->height;
        destradius->data[y*destradius->width+x] = data->powerspectrum->data[y2*data->powerspectrum->width+x2];
    }
    
    //image_save(image_from_matrix(destradius,IMAGE_CLAMP), "radius.png");
    
    bytearrayset_t*masks = get_masks(data->scales,m->width,m->height);
    bytearray_t*highmask = masks->m[masks->num-1];
    
    int it;
    for(it=0;it<10;it++) {
        //printf("iteration %d\n", it);
        //sprintf(filename, "image%03d.png", it);
        //matrix_save(m, filename);
        complex_matrix_t*src_fft = matrix_fft(m);
        /** adjust skewness and kurtosis across scales **/
        int t;
        //for(t=0;t<txtmasks->num;t++) {
        int j;
        int len = src_fft->width*src_fft->height;
        for(t=0;t<masks->num;t++) {
            complex_matrix_t*part = complex_matrix_clone(src_fft);
            bytearray_t*mask = masks->m[t];
            int x,y,s;
            assert(mask->width == part->width);
            assert(mask->height == part->height);
            double d = 1.0 / (255*part->width*part->height);
            //sprintf(filename, "scale%d_fft_beforemask.png", t);
            //complex_matrix_save(part, filename);
            for(s=0;s<len;s++) {
                part->data[s].real *= mask->data[s]*d;
                part->data[s].imag *= mask->data[s]*d;
            }

            matrix_t*r = complex_matrix_ifft_real(part);
        
            //statistics_print(statistics_new_frommatrix(r));

            //printf("before%d:",t);statistics_print(statistics_new_frommatrix(r));
            //printf("goal%d :",t);statistics_print(scale_st[t]);
            
            //sprintf(filename, "before_acorr_adjust_%d.png", t);
            //matrix_save(r, filename);

            if(t==masks->num-1) {
                // highband
                matrix_adjust_variance(r, data->scale_st[t]->var);
                //matrix_adjust_variance(r, 0);
                matrix_adjust_mean(r, 0);
                memset(r->data, 0, len*sizeof(double));
            } else {
                if(config_useautocorrelation) {
                    matrix_adjust_mean(r, 0);
                    matrix_t*r2 = matrix_scaledown(r, m->width / mask_sizediv(t, data->scales), 
                                                      m->height/ mask_sizediv(t, data->scales));
                    matrix_adjustautocorrelation3(r2, data->scale_autocorr[t]);
                    matrix_delete(r);
                    r = matrix_scaleup(r2, m->width, m->height);
                }
                matrix_adjust_statistics(r, data->scale_st[t]);
            }
            //printf("after%d :",t);statistics_print(statistics_new_frommatrix(r));

            complex_matrix_t*adjusted = matrix_fft(r);
            adjusted->data[0].real = 0; //set mean to 0
            adjusted->data[0].imag = 0; 

            for(s=0;s<len;s++) {
                if(config_spyr) {
                    src_fft->data[s].real *= 1-mask->data[s]/(255.0);
                    src_fft->data[s].imag *= 1-mask->data[s]/(255.0);
                    src_fft->data[s].real += adjusted->data[s].real;
                    src_fft->data[s].imag += adjusted->data[s].imag;
                } else {
                    if(mask->data[s]) {
                        src_fft->data[s].real = adjusted->data[s].real;
                        src_fft->data[s].imag = adjusted->data[s].imag;
                    }
                }
            }

            matrix_delete(r);r=0;
            complex_matrix_delete(part);part=0;
            src_fft->data[0].real = 0;
            src_fft->data[0].imag = 0;
            //sprintf(filename, "image%03d_fft.png",t);
            //complex_matrix_save(src_fft, filename);
            //sprintf(filename, "image%03d.png", t);
            //matrix_save(complex_matrix_ifft_real(src_fft), filename);
        }
        int s;
        for(s=0;s<len;s++) {
            if(highmask->data[s]) {
                src_fft->data[s].real = 0;
                src_fft->data[s].imag = 0;
            }
            src_fft->data[s].real /= len;
            src_fft->data[s].imag /= len;
        }

        matrix_t*mnew = complex_matrix_ifft_real(src_fft);
        matrix_adjust_statistics(mnew, data->st);
        double fade = 1.0;
        if(config_usepowerspectrum)
            fade = 0.3;
        for(t=0;t<m->width*m->height;t++) {
            m->data[t] = m->data[t]*(1-fade) + mnew->data[t]*fade;// + (lrand48()/2147483648.0);
        }
        matrix_delete(mnew);mnew=0;
        complex_matrix_delete(src_fft);src_fft=0;

        if(config_usepowerspectrum) {
            /** adjust power spectrum **/
            src_fft = matrix_fft(m);

            assert(highmask->width == destradius->width);
            assert(highmask->height == destradius->height);
            assert(src_fft->width == destradius->width);
            assert(src_fft->height == destradius->height);
            for(y=0;y<destradius->height;y++)
            for(x=0;x<destradius->width;x++) {
                double r = destradius->data[y*destradius->width+x];
                complex_t*ss = &src_fft->data[y*src_fft->width+x];
                double r2 = sqrt(ss->real*ss->real + ss->imag*ss->imag);
                if(!r) {
                    r = 0;//r2/32;
                }
                //r = (r+r2)/2;

                if(highmask->data[y*highmask->width+x]) {
                    r = 0; //remove high mask
                }

                if(r2) {
                    ss->real *= r / r2 / (src_fft->width*src_fft->height);
                    ss->imag *= r / r2 / (src_fft->width*src_fft->height);
                }
                //src_fft->data[t].imag = 0;
            }
            mnew = complex_matrix_ifft_real(src_fft);
            matrix_adjust_variance(mnew, data->st->var);
            matrix_adjust_mean(mnew, data->st->mean);
            complex_matrix_delete(src_fft);src_fft=0;
            fade = 0.3;
            for(t=0;t<m->width*m->height;t++) {
                m->data[t] = m->data[t]*(1-fade) + mnew->data[t]*fade;// + (lrand48()/2147483648.0);
            }
            matrix_delete(mnew);mnew=0;
        }

        matrix_adjust_statistics(m, data->st);
        matrix_adjust_minmax(m, data->st);
#if 1
        /* copy image into the center of the synthesized texture */
        int x1=(m->width-orig->width)/2;
        int y1=(m->height-orig->height)/2;
        for(y=0;y<orig->height;y++)
        for(x=0;x<orig->width;x++) {
            double*from = &orig->data[y*orig->width+x];
            double*to = &m->data[(y+y1)*m->width+(x+x1)];
            if(image->data[y*image->width+x].a) 
                *to = *from;//*to*0.25 + *from*0.75;
        }
#endif
    }

    //complex_matrix_t*fft = matrix_fft(m);
    //fft->data[0].real = fft->data[0].imag = 0;
    //complex_matrix_save(fft, "result_fft.png");
    //printf("goal  :");statistics_print(statistics_new_frommatrix(texture));
    //printf("result:");statistics_print(statistics_new_frommatrix(m));

    image_t*result = image_new(m->width, m->height);
    int t;
    for(t=0;t<m->width*m->height;t++) {
	int yy = m->data[t];
        result->data[t].r = clamp00ff(yy + ((360*(data->v-128))>>8));
        result->data[t].g = clamp00ff(yy - ((88*(data->u-128)+183*(data->v-128))>>8));
        result->data[t].b = clamp00ff(yy + ((455 * (data->u-128))>>8));
        result->data[t].a = 255;
    }
    matrix_delete(destradius);destradius=0;
    matrix_delete(m);m=0;
    matrix_delete(orig);orig=0;
    if(image_tofree) {
        image_delete(image_tofree);
    }
    return result;
}

texturedata_t* texturedata_fromimage(image_t*image, int method)
{
    texturedata_t*data = malloc(sizeof(texturedata_t));
    memset(data, 0, sizeof(texturedata_t));
    data->type = method;

    if(method <= 1) {
        data->internal = analyze_texture12(image);
    } else if(method == 2) {
        data->internal = analyze_texture12(image);
    } else {
        fprintf(stderr, "Internal error: invalid texture analyze method\n");
    }
    return data;
}

image_t* texturedata_synthesize_image(texturedata_t*data, image_t*alpha)
{
    if(data->type <= 1) {
        return synthesize1(alpha, data);
    } else if(data->type == 2) {
        return synthesize2(alpha, data);
    } else {
        fprintf(stderr, "Internal error: invalid texture synth method\n");
    }
    return 0;
}

void texturedata_delete(texturedata_t*data) 
{
    if(data->type == 2) {
        texturedata_free2(data);
    }
    memset(data, 0, sizeof(data));
    free(data);
}

double texturedata_dist(texturedata_t*_txt1, texturedata_t*_txt2)
{
    texturedata_internal2_t*txt1 = _txt1->internal;
    texturedata_internal2_t*txt2 = _txt2->internal;
    double dist = 0;
    dist += fabs(txt1->u - txt2->u)*32;
    dist += fabs(txt1->v - txt2->v)*32;
    dist += fabs(txt1->st->mean - txt2->st->mean)*32;
    dist += fabs(txt1->st->dev - txt2->st->dev);
    dist += fabs(txt1->st->skew - txt2->st->skew)*256;
    dist += fabs(txt1->st->kurt - txt2->st->kurt)*256;
    dist += fabs(txt1->st->min - txt2->st->min);
    dist += fabs(txt1->st->max - txt2->st->max);
    if(txt1->scales!=txt2->scales) {
        dist += 5000; // two textures with *very* different resolutions, add a static total
    } else {
        int ms = txt1->scales;
        int t;
        for(t=0;t<ms;t++) {
            dist += fabs(txt1->scale_st[t]->var - txt2->scale_st[t]->var) / ms;
            dist += fabs(txt1->scale_st[t]->skew - txt2->scale_st[t]->skew)*256 / ms;
            dist += fabs(txt1->scale_st[t]->kurt - txt2->scale_st[t]->kurt)*256 / ms;
        }
    }
    return dist;
}
