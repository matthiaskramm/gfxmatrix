/* image.c 

   Image loading + saving + conversion.

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

#include <math.h>
#include <assert.h>
#include "../png.c"
#include "image.h"
#include "stat.h"

inline static void*malloc16(int len)
{
    void*ptr=0;
    posix_memalign(&ptr, 16, len);
    return ptr;
}

static double lab_f(double t)
{
    return t;
    if(t > 0.008856) {
        return pow(t, 1/3.0);
    } else {
        return 7.787*t + 16.0/116.0;
    }
}

static void rgb2xyz(RGBA* rgba, double*dest)
{
    dest[0] = rgba->r*0.4124 + rgba->g*0.3576 + rgba->b*0.1805;
    dest[1] = rgba->r*0.2126 + rgba->g*0.7152 + rgba->b*0.0722;
    dest[2] = rgba->r*0.0193 + rgba->g*0.1191 + rgba->b*0.9505;
}

static RGBA white = {255,255,255,255};

matrix_t* image_extractchannel(image_t*img, int mode)
{
    int t;
    int l = img->width*img->height;
    matrix_t* m = matrix_new(img->width, img->height);
    int channel = 0;

    if(mode&IMAGE_GRAY) {
        for(t=0;t<l;t++) {
            m->data[t] = ((double)img->data[t].r + (double)img->data[t].g + (double)img->data[t].b)/3.0;
        }
    } else {
        if(mode&(IMAGE_LAB_L|IMAGE_LAB_A|IMAGE_LAB_B)) {
            double xyzn[3]; 
            rgb2xyz(&white, xyzn);
            for(t=0;t<l;t++) {
                double xyz[3];
                rgb2xyz(&img->data[t], xyz);
                if(mode&IMAGE_LAB_L)  {
                    m->data[t] = 116 * lab_f(xyz[1] / xyzn[1]) - 16;
                }
                if(mode&IMAGE_LAB_A) {
                    m->data[t] = 500 * lab_f(xyz[0] / xyzn[0]) - lab_f(xyz[1] / xyzn[1]);
                }
                if(mode&IMAGE_LAB_B) {
                    m->data[t] = 200 * lab_f(xyz[1] / xyzn[1]) - lab_f(xyz[2] / xyzn[2]);
                }
            }
        } else if(mode&(IMAGE_YUV_Y|IMAGE_YUV_U|IMAGE_YUV_V)) {
	    if(mode&IMAGE_YUV_Y) {
                for(t=0;t<l;t++) 
                    m->data[t] = (img->data[t].r*((int)( 0.299*256)) + img->data[t].g*((int)( 0.587*256)) + img->data[t].b*((int)( 0.114 *256)))>>8;
            } else if(mode&IMAGE_YUV_U) {
                for(t=0;t<l;t++)
	            m->data[t] = (img->data[t].r*((int)(-0.169*256)) + img->data[t].g*((int)(-0.332*256)) + img->data[t].b*((int)( 0.500 *256))+ 128*256)>>8;
            } else if(mode&IMAGE_YUV_V) {
                for(t=0;t<l;t++)
                    m->data[t] = (img->data[t].r*((int)( 0.500*256)) + img->data[t].g*((int)(-0.419*256)) + img->data[t].b*((int)(-0.0813*256))+ 128*256)>>8;
            }
        } else {
            //if(mode&IMAGE_ALPHA)
            //    channel = 0;
            if(mode&IMAGE_RED)
                channel = 1;
            if(mode&IMAGE_GREEN)
                channel = 2;
            if(mode&IMAGE_BLUE)
                channel = 3;
            for(t=0;t<l;t++) {
                m->data[t] = ((unsigned char*)&img->data[t])[channel];
            }
        }
    }
   
    if(mode&IMAGE_NORMALIZE) {
        double sum = 0;
        for(t=0;t<l;t++) {
            sum += m->data[t];
        }
        sum /= l;
        for(t=0;t<l;t++) {
            m->data[t] -= sum;
        }
    }

    return m;
}

bytearray_t* image_getchannel(image_t*img, int mode)
{
    int t;
    int l = img->width*img->height;
    bytearray_t* m = bytearray_new(img->width, img->height);

    if(mode&IMAGE_GRAY) {
        for(t=0;t<l;t++) {
            m->data[t] = (img->data[t].r + img->data[t].g + img->data[t].b) / 3;
        }
    } else {
        if(mode&(IMAGE_LAB_L|IMAGE_LAB_A|IMAGE_LAB_B)) {
            fprintf(stderr, "Error: LAB not supported for bytearrays\n");
        } else if(mode&(IMAGE_YUV_Y|IMAGE_YUV_U|IMAGE_YUV_V)) {
	    if(mode&IMAGE_YUV_Y) {
                for(t=0;t<l;t++) 
                    m->data[t] = (img->data[t].r*((int)( 0.299*256)) + img->data[t].g*((int)( 0.587*256)) + img->data[t].b*((int)( 0.114 *256)))>>8;
            } else if(mode&IMAGE_YUV_U) {
                for(t=0;t<l;t++)
	            m->data[t] = (img->data[t].r*((int)(-0.169*256)) + img->data[t].g*((int)(-0.332*256)) + img->data[t].b*((int)( 0.500 *256))+ 128*256)>>8;
            } else if(mode&IMAGE_YUV_V) {
                for(t=0;t<l;t++)
                    m->data[t] = (img->data[t].r*((int)( 0.500*256)) + img->data[t].g*((int)(-0.419*256)) + img->data[t].b*((int)(-0.0813*256))+ 128*256)>>8;
            }
        } else {
            int channel = 0;
            if(mode&IMAGE_RED) channel = 1;
            else if(mode&IMAGE_GREEN) channel = 2;
            else if(mode&IMAGE_BLUE) channel = 3;
            for(t=0;t<l;t++) {
                m->data[t] = ((unsigned char*)&img->data[t])[channel];
            }
        }
    }
    if(mode&IMAGE_DIV2) {
        for(t=0;t<l;t++) 
            m->data[t] = m->data[t]>>1;
    }
    return m;
}

image_t*image_load(char*filename)
{
    image_t*img = malloc(sizeof(image_t));
    getPNG(filename, &img->width, &img->height, (unsigned char**)&img->data);
    return img;
}

int isfinite(double x)
{
    __finite(x);
}

void image_update_from_matrix2(image_t*img, matrix_t*m, int flags, int xpos, int ypos, double min, double max)
{
    double*pic = m->data;
    assert(m->width+xpos <= img->width && m->height+ypos <= img->height);

    int t;
    if(min == max) {
        max = -HUGE_VAL;
        min = HUGE_VAL;

        for(t=0;t<m->width*m->height;t++) {
            if(__finite(pic[t])) {
                if(pic[t]>max)
                    max = pic[t];
                if(pic[t]<min)
                    min = pic[t];
            }
        }
        if(max==min) {
            min = 0;
            max = 1.0;
        }
    }
    if(flags&IMAGE_ABSMINMAX) {
        if(max < 0) ;
        else if(-min > max) max=-min;
        else min = -max;
    }

    int x,y;
    for(y=0;y<m->height;y++) {
        double*line = &m->data[m->width*y];
        RGBA*out = &img->data[img->width*(y+ypos)+xpos];
        for(x=0;x<m->width;x++) {
            if(!__finite(line[x])) {
                out[x].a = 255;
                out[x].r = 255;
                out[x].g = 0;
                out[x].b = 0;
            } else if(flags&(IMAGE_MINMAX|IMAGE_ABSMINMAX)) {
                int val = (line[x]-min)*255/(max-min);
                int rr=255,gg=255,bb=255;
                if(val>255) 
                    val=255;
                if(val<0) 
                    val=0;
                if(flags&IMAGE_MARKNEGATIVE) {
                    if(line[x]<0) {
                        gg = 0;
                        bb = 0;
                    } else {
                        rr = 0;
                    }
                }
                out[x].a = 255;
                out[x].r = ((unsigned char)val)&rr;
                out[x].g = ((unsigned char)val)&gg;
                out[x].b = ((unsigned char)val)&bb;
            } else if(flags&IMAGE_CLAMP) {
                int val = line[x];
                char overflow = 0;
                if(val>255) {
                    overflow=1;
                    val=255;
                }
                if(val<0) {
                    val=0;overflow=1;
                }
                out[x].a = 255;
                out[x].r = (unsigned char)val;
                out[x].g = (unsigned char)val;
                out[x].b = (unsigned char)val;
                if(overflow && (flags&IMAGE_MARKOVERFLOW)) {
                    out[x].b = 0;
                }
            } else {
                int val = (line[x]+1.0)*255/2;
                if(val>255) 
                    val=255;
                if(val<0) 
                    val=0;
                out[x].a = 255;
                out[x].r = (unsigned char)val;
                out[x].g = (unsigned char)val;
                out[x].b = (unsigned char)val;
            }
            if(flags&IMAGE_OPENGL) {
                unsigned char a = out[x].a;
                unsigned char r = out[x].r;
                unsigned char g = out[x].g;
                unsigned char b = out[x].b;
                out[x].a = r;
                out[x].r = g;
                out[x].g = b;
                out[x].b = a;
            }
        }
    }
}
void image_update_from_matrix(image_t*img, matrix_t*m, int flags, int xpos, int ypos)
{
    image_update_from_matrix2(img, m, flags, xpos, ypos, 0,0);
}


image_t*image_from_complex_matrix(complex_matrix_t*m, int flags)
{
    image_t*img = image_new(m->width*2 + 1, m->height);
    memset(img->data, 255, img->width*img->height*sizeof(RGBA));
    matrix_t*m1,*m2;
    m1 = complex_matrix_realpart(m);
    m2 = complex_matrix_imagpart(m);
    statistics_t*stat1 = statistics_new_frommatrix(m1);
    statistics_t*stat2 = statistics_new_frommatrix(m2);
    
    if(stat2->min < stat1->min)
        stat1->min = stat2->min;
    if(stat2->max > stat1->max)
        stat1->max = stat2->max;

    image_update_from_matrix2(img, m1, flags, 0, 0, stat1->min, stat1->max);
    image_update_from_matrix2(img, m2, flags, m->width + 1, 0, stat1->min, stat1->max);

    statistics_delete(stat1);
    statistics_delete(stat2);

    matrix_delete(m1);
    matrix_delete(m2);
    return img;
}


image_t*image_from_matrix(matrix_t*m, int flags)
{
    image_t*img = malloc(sizeof(image_t));
    img->width = m->width;
    img->height = m->height;
    img->data = (RGBA*)malloc16(m->width*m->height*4);
    image_update_from_matrix(img, m, flags, 0, 0);
    return img;
}

void image_save(image_t*img, char*filename) 
{
    writePNG(filename, (unsigned char*)img->data, img->width, img->height);
}
void image_save_and_free(image_t*img, char*filename)
{
    image_save(img, filename);
    image_delete(img);
}
image_t*image_new(int width, int height)
{
    image_t*img = malloc(sizeof(image_t));
    img->width = width;
    img->height = height;
    img->data = malloc16(4*width*height);
    memset(img->data, 0, 4*width*height);
    return img;
}
image_t*image_clone(image_t*img)
{
    image_t*n= image_new(img->width, img->height);
    memcpy(n->data, img->data, n->width*n->height*4);
    return n;
}
void image_crop(image_t*img, int width, int height)
{
    assert(width>0);
    assert(height>0);
    assert(width<=img->width);
    assert(height<=img->height);
    RGBA*data = malloc16(width*height*sizeof(RGBA));
    int y;
    for(y=0;y<height;y++)
        memcpy(&data[y*width], &img->data[y*img->width], width*sizeof(RGBA));
    free(img->data);
    img->data = data;
    img->width = width;
    img->height = height;
}

void image_delete(image_t*img)
{
    free(img->data);
    memset(img, 0, sizeof(image_t));
    free(img);
}
void rgb2rgba(unsigned char*src, unsigned char*dest, int len)
{
    int t;
    for(t=0;t<len;t++) {
        dest[t*4+0] = src[t*3+0];
        dest[t*4+1] = src[t*3+1];
        dest[t*4+2] = src[t*3+2];
        dest[t*4+3] = 255;
    }
}
void image_paste(image_t*dest, int x, int y, image_t*src)
{
    int width=src->width;
    int height=src->height;
    if(x+width > dest->width)
        width = dest->width - x;
    if(y+height > dest->height)
        height = dest->height - y;
    if(x<0 || y<0) {
        fprintf(stderr, "negative coordinates not allowed in paste");
        return;
    }
    if(width<=0 || height<=0) {
        fprintf(stderr, "Area (%d,%d,%d,%d) is outside of (0,0,%d,%d)\n", 
                x, y, x+src->width, y+src->height, dest->width, dest->height);
        return;
    }
    int startx = x;
    int starty = y;
    for(y=0;y<height;y++) {
        RGBA*dline = &dest->data[(y+starty)*dest->width+startx];
        RGBA*sline = &src->data[y*src->width];
        for(x=0;x<width;x++) {
            dline[x] = sline[x];
        }
    }
}
image_t* image_cut(image_t*src, int x, int y, int width, int height)
{
    if(x+width > src->width)
        width = src->width - x;
    if(y+height > src->height)
        height = src->height - y;
    if(x<0 || y<0) {
        fprintf(stderr, "negative coordinates not allowed in cut");
        return;
    }
    if(width<=0 || height<=0) {
        fprintf(stderr, "Area (%d,%d,%d,%d) is outside of (0,0,%d,%d)\n", 
                x, y, x+width, y+height, width, height);
        return;
    }
    image_t*dest = image_new(width, height);
    int startx = x;
    int starty = y;
    for(y=0;y<height;y++) {
        RGBA*sline = &src->data[(y+starty)*src->width+startx];
        RGBA*dline = &dest->data[y*dest->width];
        memcpy(dline, sline, sizeof(dline[0])*width);
    }
    return dest;
}

void bytearray_save_map(bytearray_t*b, char*filename)
{
    image_t*img = image_new(b->width, b->height);

    int t;
    for(t=0;t<b->width*b->height;t++) {
        img->data[t].r = (b->data[t]<<0)&0xe0;
        img->data[t].g = (b->data[t]<<3)&0xe0;
        img->data[t].b = (b->data[t]<<6)&0xc0;
        img->data[t].a = 255;
    }
    image_save(img, filename);
    image_delete(img);
}

static int*find_n_maximas(int*data, int len, int num, int cutoff)
{
    int* max = malloc(sizeof(int)*num);
    int* maxp = malloc(sizeof(int)*num);
    memset(max, 0, sizeof(int)*num);
    memset(maxp, -1, sizeof(int)*num);
    int t;
    int lowestmax = 0;
    for(t=0;t<len;t++) {
        int s;
        if(data[t]>lowestmax) {
            for(s=0;s<num;s++) {
                if(data[t] > max[s]) {
                    int r;
                    for(r=num-1;r>s;r--) {
                        max[r] = max[r-1];
                        maxp[r] = maxp[r-1];
                    }
                    max[s] = data[t];
                    maxp[s] = t;
                    lowestmax = max[num-1];
                    break;
                }
            }
        }
    }
    for(t=0;t<num;t++) {
        if(max[t]<cutoff)
            maxp[t] = -1;
    }
    free(max);
    return maxp;
}

static inline char colcmp(RGBA c1, RGBA c2)
{
    return c1.r == c2.r && c1.g == c2.g && c1.b == c2.b;
}


bytearray_t* bytearray_load_map(char*filename)
{
    image_t*mask = image_load(filename);
    
    int*count = malloc(256*256*256*sizeof(int));
    memset(count, 0, 256*256*256*sizeof(int));
   
    int x,y;
    int t;

    char ismap=1;

    for(y=1;y<mask->height-1;y++) {
        RGBA* m0= &mask->data[(y-1)*mask->width];
        RGBA* m1= &mask->data[y*mask->width];
        RGBA* m2= &mask->data[(y+1)*mask->width];
        for(x=1;x<mask->width-1;x++) {
            RGBA c = m1[x];
            if(colcmp(m0[x-1],c) && colcmp(m0[x],c) && colcmp(m0[x+1],c) &&
               colcmp(m1[x-1],c) && colcmp(m1[x],c) && colcmp(m1[x+1],c) &&
               colcmp(m2[x-1],c) && colcmp(m2[x],c) && colcmp(m2[x+1],c)) {
                count[c.r*65536 + c.g*256 +c.b]++;
            }
        }
    }
    int* max = find_n_maximas(count, 256*256*256, 64, 0);
    free(count);
    
    RGBA cols[64];
    int num=0;
    for(t=0;t<64;t++) {
        if(max[t]>=0) {
            cols[num].r = max[t]>>16;
            cols[num].g = max[t]>>8;
            cols[num].b = max[t];
            cols[num].a = 255;
            num++;
        }
    }
    bytearray_t*b = bytearray_new(mask->width, mask->height);
    RGBA white = {255,255,255,255};
    
    int s;
    for(s=0;s<mask->width*mask->height;s++) {
        if((mask->data[s].r&~0xe0) ||
           (mask->data[s].g&~0xe0) ||
           (mask->data[s].b&~0xc0))
            ismap = 0;
    }
    if(ismap) {
        printf("image %s is in map format\n", filename);
    }

    for(s=0;s<mask->width*mask->height;s++) {
        int t;
        if(colcmp(white, mask->data[s])) {
            b->data[s] = 255;
            continue;
        }
        RGBA col = mask->data[s];
        for(t=0;t<num;t++) {
            if(colcmp(cols[t], col))
                break;
        }
        if(t<num) {
            if(ismap) {
                t = col.r|col.g>>3|col.b>>6;
            }
            b->data[s] = t;
        } else {
            b->data[s] = 255;
        }
    }
    free(max);
    image_delete(mask);
    return b;
}

image_t* image_extract(image_t*img, window_t*w)
{
    image_t*i = image_new(w->x2 - w->x1, w->y2 - w->y1);
    int x,y;
    for(y=0;y<i->height;y++) {
        int x;
        for(x=0;x<i->width;x++) {
            i->data[y*i->width+x] = img->data[(y+w->y1)*img->width+w->x1+x];
        }
        //memcpy(&i->data[y*i->width], &img->data[(y+w->y1)*img->width+w->x1], i->width*sizeof(img->data[0]));
    }
    return i;
}

void image_insert_tiling(image_t*img, bytearray_t*map, image_t*replacement, int t)
{
    int x,y;
    for(y=0;y<img->height;y++) {
        RGBA*ysrc = &replacement->data[(y%replacement->height)*replacement->width];
        RGBA*ydest = &img->data[y*img->width];
        unsigned char*ymap = 0;
        if(map) {
            ymap = &map->data[y*img->width];
        }
        for(x=0;x<img->width;x++) {
            if(!ymap || ymap[x]==t) {
                ydest[x] = ysrc[x%replacement->width];
            }
        }
    }
}

static inline double sqr(double x)
{
    return x*x;
}

double image_compare(image_t*src, image_t*dest)
{
    int width = src->width<dest->width? src->width : dest->width;
    int height = src->height<dest->height? src->height : dest->height;
    
    double corrdiff = 0;
    double corrabs = 0;
    int x,y;
    for(y=0;y<height;y++) {
        RGBA*dline = &dest->data[y*dest->width];
        RGBA*sline = &src->data[y*src->width];
        for(x=0;x<width;x++) {
            corrdiff += sqr(dline[x].r - sline[x].r);
            corrdiff += sqr(dline[x].g - sline[x].g);
            corrdiff += sqr(dline[x].b - sline[x].b);
            corrabs += sqr(dline[x].r);
            corrabs += sqr(dline[x].g);
            corrabs += sqr(dline[x].b);
        }
    }

    return 10 * log(corrabs / corrdiff) / log(10);
}

static double inline Cut_4th(double to_cut)
{
    return ((int)(0.5 + to_cut * 1000))/1000.0;
}

static double inline max(double d1, double d2) {
    if(d1>d2) return d1;
    else return d2;
}

static double inline min(double d1, double d2) {
    if(d1<d2) return d1;
    else return d2;
}
HLS rgb2hls(RGBA c)
{
    HLS color_hls;

    double Wrgb[] = {0.299, 0.587, 0.144};
    double sum = Wrgb[0] + Wrgb[1] + Wrgb[2];
    Wrgb[0] /= sum; Wrgb[1] /= sum; Wrgb[2] /= sum;

    double R = c.r / 255.0;
    double G = c.g / 255.0;
    double B = c.b / 255.0;

    double L = R*Wrgb[0]+G*Wrgb[1]+B*Wrgb[2];
    double H,S;

    if(L==1 || L==0) {
        S = 0;
        H = 0;
    } else {
        double RGBMin = R<B?(R<G?R:G):(B<G?B:G);
        double RGBMax = R>B?(R>G?R:G):(B>G?B:G);

        assert(RGBMax >= R && RGBMax >= G && RGBMax >= B);
        assert(RGBMin <= R && RGBMin <= G && RGBMin <= B);

        S = max(1 - RGBMin/L, (RGBMax-L)/(1-L));

        if (R >= G && G > B) { // r max, b min
            H = ( 1 - ( R - G ) / ( R - B ) ) / 6.0;
        } else if (G >= R && R > B) { // g max, b min
            H = ( 1 + ( G - R ) / ( G - B ) ) / 6.0;
        } else if (G >= B && B > R) { // g max, r min
            H = ( 3 - ( G - B ) / ( G - R ) ) / 6.0;
        } else if (B >= G && G > R) { // b max, r min
            H = ( 3 + ( B - G ) / ( B - R ) ) / 6.0;
        } else if (B >= R && R > G) { // b max, g min
            H = ( 5 - ( B - R ) / ( B - G ) ) / 6.0;
        } else if (R >= B && B > G) { // r max, g min
            H = ( 5 + ( R - B ) / ( R - G ) ) / 6.0;    
        } else {
            H = 0;
            S = 0;
        }
    }
    
    double Hnew = 0;
    if (0<=H && H<1/6.0)
        Hnew = H * 2.0;
    else if (1/6.0<=H && H<1/3.0)
        Hnew = H+1/6.0;
    else if (1/3.0<=H && H<2/3.0)
        Hnew = H/2.0 + 1/3.0;
    else if (2/3.0<=H && H<1.0)
        Hnew = H;

    H = Hnew;

    color_hls.h = max(min(Cut_4th(H), 1.0), 0.0);
    color_hls.l = max(min(Cut_4th(L), 1.0), 0.0);
    color_hls.s = max(min(Cut_4th(S), 1.0), 0.0);

    return color_hls;
}
