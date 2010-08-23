/* matrix.c 

   Matrix + linear algebra.

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
#include <math.h>
#include <memory.h>
#include <assert.h>
#include "matrix.h"
#include "image.h"
#include "stat.h"
#include "gradopt.h"
#include <fftw3.h>

inline static void*malloc16(int len)
{
    void*ptr=0;
    posix_memalign(&ptr, 16, len);
    return ptr;
}

matrix_t*matrix_clone(matrix_t*m)
{
    matrix_t*n = (matrix_t*)malloc(sizeof(matrix_t));
    n->width = m->width;
    n->height = m->height;
    n->data = (double*)malloc16(m->width*m->height*sizeof(m->data[0]));
    memcpy(n->data, m->data, m->width*m->height*sizeof(m->data[0]));
    return n;
}
complex_matrix_t*complex_matrix_clone(complex_matrix_t*m)
{
    complex_matrix_t*n = (complex_matrix_t*)malloc(sizeof(complex_matrix_t));
    n->width = m->width;
    n->height = m->height;
    n->data = (complex_t*)malloc16(m->width*m->height*sizeof(m->data[0]));
    memcpy(n->data, m->data, m->width*m->height*sizeof(m->data[0]));
    return n;
}
void matrix_clear(matrix_t*m)
{
    free(m->data);m->data = 0;
    m->width = 0;
    m->height = 0;
}
void matrix_print_stats(matrix_t*m)
{
    statistics_t*stat = statistics_new_frommatrix(m);
    statistics_print(stat);
    statistics_delete(stat);
}
void complex_matrix_clear(complex_matrix_t*m)
{
    free(m->data);m->data = 0;
    m->width = 0;
    m->height = 0;
}
void matrix_delete(matrix_t*m)
{
    matrix_clear(m);free(m);
}
void complex_matrix_delete(complex_matrix_t*m)
{
    complex_matrix_clear(m);free(m);
}
matrix_t* matrix_new(int width, int height)
{
    matrix_t*m = (matrix_t*)malloc(sizeof(matrix_t));
    m->data = (double*)malloc16(sizeof(double)*width*height);
    memset(m->data, 0, sizeof(double)*width*height);
    m->width = width;
    m->height = height;
    return m;
}
complex_matrix_t* complex_matrix_new(int width, int height)
{
    complex_matrix_t*m = (complex_matrix_t*)malloc(sizeof(complex_matrix_t));
    m->data = (complex_t*)malloc16(width*height*sizeof(complex_t));
    memset(m->data, 0, sizeof(complex_t)*width*height);
    m->width = width;
    m->height = height;
    return m;
}
matrix_t*matrix_new_random(int width, int height, double min, double max)
{
    matrix_t*m = matrix_new(width, height);
    int x,y;
    for(y=0;y<height;y++)
    for(x=0;x<width;x++) {
        m->data[y*width+x] = drand48() * (max-min) + min;
    }
    return m;
}
matrix_t*matrix_new_gaussrandom(int width, int height, double mean, double var)
{
    matrix_t*m = matrix_new(width, height);

    int num_iterations = 50;

    var = 2 * sqrt(var*3.0 / (double)num_iterations);

    int x,y;
    for(y=0;y<height;y++)
    for(x=0;x<width;x++) {
	int i;
	double r = mean;
	for(i=0;i<num_iterations;i++) {
	    r += (drand48()-0.5)*var;
	}
        m->data[y*width+x] = r;
    }
    return m;
}
matrix_t*matrix_new_unit(int width, int height)
{
    matrix_t*m = matrix_new(width, height);
    int x,y;
    for(y=0;y<height;y++)
    for(x=0;x<width;x++) {
        m->data[y*width+x] = x==y;
    }
    return m;
}
matrix_t* matrix_new_fromdata(double*data, int width, int height)
{
    matrix_t*m = (matrix_t*)malloc(sizeof(matrix_t));
    m->data = data;
    m->width = width;
    m->height = height;
    return m;
}
void matrix_add_inplace(matrix_t*m,matrix_t*add)
{
    int t,l=m->width*m->height;
    if(m->width != add->width || m->height != add->height) {
        fprintf(stderr, "Matrix sizes mismatch: %d,%d <-> %d,%d\n", m->width, m->height, add->width, add->height);
        return;
    }
    for(t=0;t<l;t++) {
        m->data[t] += add->data[t];
    }
}
void matrix_update(matrix_t*dest, matrix_t*src)
{
    assert(src->width == dest->width && src->height == dest->height);
    memcpy(dest->data, src->data, dest->width*dest->height*sizeof(dest->data[0]));
}
void matrix_store(matrix_t*m, char*filename)
{
    FILE*fi = fopen(filename, "wb");
    int x,y;
    for(y=0;y<m->height;y++) {
        for(x=0;x<m->width;x++) {
            fprintf(fi, "%f\t", m->data[y*m->width+x]);
        }
        fprintf(fi, "\n");
    }
    fclose(fi);
}
void bytearray_store(bytearray_t*m, char*filename)
{
    FILE*fi = fopen(filename, "wb");
    int x,y;
    for(y=0;y<m->height;y++) {
        for(x=0;x<m->width;x++) {
            fprintf(fi, "%d\t", m->data[y*m->width+x]);
        }
        fprintf(fi, "\n");
    }
    fclose(fi);
}
void matrix_save(matrix_t*b, char*filename)
{
    image_save_and_free(image_from_matrix(b, IMAGE_MINMAX), filename);
}
void complex_matrix_save(complex_matrix_t*b, char*filename)
{
    image_save_and_free(image_from_complex_matrix(b, IMAGE_MINMAX), filename);
}


void matrix_print(matrix_t*m)
{
    if(!m)
        return;
    FILE*fi = stdout;
    int x,y;
    for(y=0;y<m->height;y++) {
        if(y==0) fprintf(fi, "/ ");
        else if(y<m->height-1) fprintf(fi, "| ");
        else fprintf(fi, "\\ ");
        for(x=0;x<m->width;x++) {
            fprintf(fi, "%10.4f ", m->data[y*m->width+x]);
        }
        if(y==0) fprintf(fi, " \\\n");
        else if(y<m->height-1) fprintf(fi, " |\n");
        else fprintf(fi, " /\n");
    }
    fprintf(fi, "\n");
}

void complex_matrix_print(complex_matrix_t*m)
{
    if(!m)
        return;
    FILE*fi = stdout;
    int x,y;
    for(y=0;y<m->height;y++) {
        if(y==0) fprintf(fi, "/ ");
        else if(y<m->height-1) fprintf(fi, "| ");
        else fprintf(fi, "\\ ");
        for(x=0;x<m->width;x++) {
            fprintf(fi, "%8.3f+%8.3fj ", m->data[y*m->width+x].real, m->data[y*m->width+x].imag);
        }
        if(y==0) fprintf(fi, " \\\n");
        else if(y<m->height-1) fprintf(fi, " |\n");
        else fprintf(fi, " /\n");
    }
    fprintf(fi, "\n");
}

double*matrix_multiply_with_vector(matrix_t*m, double*vector)
{
    double*result = (double*)malloc16(m->height*sizeof(double));
    int x,y;
    for(y=0;y<m->height;y++) {
        double sum = 0;
        double*row = &m->data[y*m->width];
        for(x=0;x<m->width;x++) {
            sum += vector[x] * row[x];
        }
        result[y] = sum;
    }
    return result;
}

double*matrix_multiply_with_vector_t(matrix_t*m, double*vector)
{
    double*result = (double*)malloc16(m->width*sizeof(double));
    int x,y;
    for(x=0;x<m->width;x++) {
        double sum = 0;
        double*column = &m->data[x];
        for(y=0;y<m->height;y++) {
            sum += vector[y] * column[y*m->width];
        }
        result[x] = sum;
    }
    return result;
}

matrix_t*matrix_transpose(matrix_t*m)
{
    matrix_t*n = matrix_new(m->height, m->width);
    int x,y;
    for(y=0;y<m->height;y++) {
        for(x=0;x<m->width;x++) {
            n->data[x*n->width+y] = m->data[y*m->width+x];
        }
    }
    return n;
}

matrix_t*matrix_multiply(matrix_t*m1, matrix_t*m2)
{
    matrix_t*r = matrix_new(m2->width, m1->height);
    int v,x,y;
    for(v=0;v<m2->width;v++) {
        for(y=0;y<m1->height;y++) {
            double sum = 0;
            double*row = &m1->data[y*m1->width];
            double*column = &m2->data[v];
            for(x=0;x<m1->width;x++) {
                sum += *column * *row;
                column += m2->width;row++;
            }
            r->data[y*r->width+v] = sum;
        }
    }
    return r;
}

matrix_t*matrix_multiply_t(matrix_t*m1, char t1, matrix_t*m2, char t2)
{
    if((t1?m1->height:m1->width) != (t2?m2->width:m2->height)) {
        fprintf(stderr, "Matrix sizes mismatch: %dx%d%s * %dx%d%s\n",
                m1->width, m1->height, t1?" (transposed)":"",
                m2->width, m2->height, t2?" (transposed)":"");
        exit(1);
    }
    matrix_t*r = matrix_new(t2?m2->height:m2->width, t1?m1->width:m1->height);
    int v,x,y;
    int add1x = t1?1:m1->width;
    int add1y = t1?m1->width:1;
    int add2x = t2?m2->width:1;
    int add2y = t2?1:m2->width;
    int len1 = t1?m1->height:m1->width;

    double*column = m2->data;
    double*dest = r->data;
    for(v=0;v<r->width;v++) {
        double*destcolumn = dest++;
        double*row = m1->data;
        for(y=0;y<r->height;y++) {
            double sum = 0;
            //double*row = &m1->data[y*add1x];
            //double*column = &m2->data[v*add2x];
            double*cc = column;
            double*rr = row;
            int x=len1;
            do {
                sum += *cc * *rr;
                rr += add1y;
                cc += add2y;
            } while(--x);
            
            *destcolumn = sum;

            destcolumn+=r->width;
            row += add1x;
        }
        column += add2x;
    }
    return r;
}

matrix_t*matrix_add(matrix_t*m1, double v1, matrix_t*m2, double v2)
{
    if(m1->width != m2->width || m1->height != m2->height) {
        fprintf(stderr, "Matrix sizes mismatch: %dx%d + %dx%d\n",
                m1->width, m1->height, m2->width, m2->height);
        exit(1);
    }

    int t, l = m1->width*m1->height;
    matrix_t*r = matrix_new(m1->width, m1->height);
    for(t=0;t<l;t++) {
        r->data[t] = m1->data[t]*v1 + m2->data[t]*v2;
    }
    return r;
}

matrix_t*matrix_scaledown(matrix_t*m, int width, int height)
{
    assert(m->width >= width);
    assert(m->height >= height);
    
    //assert(!(m->width%width));
    //assert(!(m->height%height));

    int x,y;
    int fx = m->width/width;
    int fy = m->height/height;
    matrix_t*dm = matrix_new(width,height);
    for(y=0;y<height;y++) {
        double*from = &m->data[y*fy*m->width];
        double*to = &dm->data[y*dm->width];
        for(x=0;x<width;x++) {
            to[x] = from[x*fx];
        }
    }
    return dm;
}

matrix_t*matrix_scaleup(matrix_t*m, int width, int height)
{
    assert(!(width%m->width));
    assert(!(height%m->height));
    assert(width >= m->width);
    assert(height >= m->height);

    assert(!(m->width&1|m->height&1));
   
    int x,y;
        
    complex_matrix_t*fft = matrix_fft(m);
    complex_matrix_t*fft2 = complex_matrix_new(width, height);

    for(y=0;y<=m->height/2;y++) {
        for(x=0;x<=m->width/2;x++) {
            int x1 = x;
            int y1 = y;
            int x2 = fft2->width-x;
            int y2 = fft2->height-y;
            int ax1 = x;
            int ay1 = y;
            int ax2 = fft->width-x;
            int ay2 = fft->height-y;
            fft2->data[y1*fft2->width+x1] = fft->data[ay1*fft->width+ax1];
            if(y) {
                fft2->data[y2*fft2->width+x1] = fft->data[ay2*fft->width+ax1];
            }
            if(x) {
                fft2->data[y1*fft2->width+x2] = fft->data[ay1*fft->width+ax2];
                if(y) {
                    fft2->data[y2*fft2->width+x2] = fft->data[ay2*fft->width+ax2];
                }
            }
        }
    }
    
    double f = 1.0/(m->width*m->height);

    int size = fft2->width*fft2->height;
    int t;
    for(t=0;t<size;t++) {
        fft2->data[t].real *= f;
        fft2->data[t].imag *= f;
    }
   
    complex_matrix_delete(fft);
    matrix_t*inv = complex_matrix_ifft_real(fft2);
    complex_matrix_delete(fft2);
    
    return inv;
}

matrix_t* matrix_getautocorrelation_old(matrix_t*img, int csizex, int csizey)
{
    matrix_t*m = matrix_new(csizex,csizey);
    int midx = (1 + (csizex / img->width))*img->width - csizex/2;
    int midy = (1 + (csizey / img->height))*img->height - csizey/2;
    int x,y,xx,yy;
    double div = img->width*img->height;
    for(yy=0;yy<csizey;yy++) {
        double*out = &m->data[yy*csizex];
        for(xx=0;xx<csizex;xx++) {
            double c = 0;
            for(y=0;y<img->height;y++) {
                int yy1 = y;
                int yy2 = (y + midy + yy) % img->height;
                double*line1 = &img->data[yy1*img->width];
                double*line2 = &img->data[yy2*img->width];
                for(x=0;x<img->width;x++) {
                    int xx1 = x;
                    int xx2 = (x + midx + xx) % img->width;
                    c += line1[xx1]*line2[xx2];
                }
            }
            out[xx] = c / div;
        }
    }
    return m;
}

double matrix_getcrosscorrelation(matrix_t*m1, matrix_t*m2)
{
    if(m2->width < m1->width) {
        matrix_t*i = m1;m1=m2;m2=i;
    }
    if(m1->width*m2->height != m2->width*m1->height) {
        fprintf(stderr, "Bad ratio: %dx%d <-> %dx%d\n", 
                m1->width, m1->height, m2->width, m2->height);
        return 0.0;
    }
    int stepx = m2->width / m1->width;
    int stepy = m2->height / m1->height;
    int y1,y2;
    double sum = 0;
    for(y1=0,y2=0;y1<m1->height;y1++,y2+=stepy) {
        double*line1 = &m1->data[y1*m1->width];
        double*line2 = &m2->data[y2*m2->width];
        if(stepx==1) {
            int x,l = m1->width;
            for(x=0;x<l;x++) {
                sum += line1[x]*line2[x];
            }
        } else {
            int x1,x2;
            for(x1=0,x2=0;x1<m1->width;x1++,x2+=stepx) {
                sum += line1[x1]*line2[x2];
            }
        }
    }
    return sum / (m1->width * m1->height);
}

matrix_t* matrixset_getmean(matrixset_t*set)
{
    if(!set->num)
        return 0;
    int width = set->m[0]->width;
    int height = set->m[0]->height;
    int l = width * height;
    matrix_t* mean = matrix_new(width,height);
    int s,t;
    for(t=0;t<l;t++) {
        double sum = 0.0;
        for(s=0;s<set->num;s++) {
             sum += set->m[s]->data[t];
        }
        mean->data[t] = sum / s;
    }
    return mean;
}

matrix_t* matrix_getcrosscorrelations(matrixset_t*set1, matrixset_t*set2, char subtractmean)
{
    matrix_t*result = matrix_new(set2->num, set1->num);
    matrix_t*mean1 = matrixset_getmean(set1);
    matrix_t*mean2 = matrixset_getmean(set2);
    int n1,n2;
    for(n1=0;n1<set1->num;n1++)
    for(n2=0;n2<set2->num;n2++) {
        if(subtractmean) {
            matrix_t*m1 = matrix_add(set1->m[n1], 1.0, mean1, -1.0);
            matrix_t*m2 = matrix_add(set2->m[n2], 1.0, mean2, -1.0);
            result->data[n1*result->width+n2] = matrix_getcrosscorrelation(m1, m2);
            matrix_delete(m1);
            matrix_delete(m2);
        } else {
            result->data[n1*result->width+n2] = matrix_getcrosscorrelation(set1->m[n1], set2->m[n2]);
        }
    }
    matrix_delete(mean1);
    matrix_delete(mean2);
    return result;
}

matrix_t*matrix_convolve(matrix_t*src_image, matrix_t*m, int xres, int yres, char edgetype)
{
    matrix_t *dst_image = matrix_new(src_image->width/xres,src_image->height/yres);
    int x_src,y_src,x_dst,y_dst;
    int x_mid = m->width / 2;
    int y_mid = m->height / 2;
    for(y_src=0,y_dst=0; y_dst<dst_image->height; y_src+=yres,y_dst++) {
        double*dst = &dst_image->data[y_dst*dst_image->width];
        for(x_src=0,x_dst=0; x_dst<dst_image->width; x_src+=xres,x_dst++) {
            int x,y;
            double val = 0;
            if(y_src-y_mid<0 || y_src-y_mid+m->height>=src_image->height ||
               x_src-x_mid<0 || x_src-x_mid+m->width>=src_image->width) {
                // special mode, for reflections
                for(y=0;y<m->height;y++) {
                    int yy = y_src + y - y_mid;
                    if(edgetype==EDGE_REFLECT) {
                        if(yy<0) {yy=(-yy-1)%src_image->height;}
                        if(yy>=src_image->height) {yy=src_image->height-1-(yy%src_image->height);}
                    } else if(edgetype==EDGE_WRAP) {
                        while(yy<0) yy+=src_image->height;
                        yy%=src_image->height;
                    } else if(edgetype==EDGE_ZERO) {
                        continue;
                    }
                    double*f = &m->data[y*m->height];
                    for(x=0;x<m->width;x++) {
                        char xclip = 0;
                        int xx = x_src + x - x_mid;
                        if(edgetype==EDGE_REFLECT) {
                            if(xx<0) {xx=(-xx-1)%src_image->width;}
                            if(xx>=src_image->width) {xx=src_image->width-1-(xx%src_image->width);}
                        } else if(edgetype==EDGE_WRAP) {
                            while(xx<0) xx+=src_image->width;
                            xx%=src_image->width;
                        } else if(edgetype==EDGE_ZERO) {
                            continue;
                        }
                        val += src_image->data[yy*src_image->width + xx]*f[x];
                    }
                }
            } else {
                for(y=0;y<m->height;y++) {
                    double*src2 = &src_image->data[x_src-x_mid+(y_src+y-y_mid)*src_image->width];
                    double*f = &m->data[y*m->height];
                    for(x=0;x<m->width;x++) {
                        val += src2[x]*f[x];
                    }
                }
            }
            dst[x_dst] = val / xres; // THIS USED TO BE just "val"
        }
    }
    return dst_image;
}

double matrix_diff(matrix_t*m1, matrix_t*m2)
{
    if(m1->width != m2->width) return m1->width*m1->height;
    if(m1->height != m2->height) return m1->height*m1->height;
    int t;
    double d = 0;
    for(t=0;t<m1->width*m1->height;t++) {
        double r = m1->data[t]-m2->data[t];
        d += r*r;
    }
    return d;
}

complex_matrix_t*complex_matrix_expand(complex_matrix_t*m, int width, int height, int flags)
{
    complex_matrix_t*newmatrix = complex_matrix_new(width, height);
    assert(width >= m->width && height >= m->height);
    if(!(flags&EXPAND_RIGHTDOWN)) {
        int x,y;
        for(y=0;y<m->height;y++) {
            for(x=0;x<m->width;x++) {
                int xx = (x + newmatrix->width - m->width/2) % newmatrix->width;
                int yy = (y + newmatrix->height - m->height/2) % newmatrix->height;
                newmatrix->data[yy*newmatrix->width+xx] = m->data[y*m->width+x];
            }
        }
    } else {
        int x,y;
        for(y=0;y<m->height;y++) {
            memcpy(&newmatrix->data[y*newmatrix->width], &m->data[y*m->width], m->width*sizeof(m->data[0]));
        }
    }
    return newmatrix;
}

matrix_t*matrix_convolve_fft(matrix_t*img, matrix_t*m)
{
    matrix_t*filter2 = matrix_new(img->width, img->height);
    assert(img->width >= m->width*2);
    int x,y;
    for(y=0;y<m->height;y++) {
        for(x=0;x<m->width;x++) {
            int xx = (x + filter2->width - m->width/2) % filter2->width;
            int yy = (y + filter2->height - m->height/2) % filter2->height;
            //int xx = (x + filter2->height/2 - m->width/2);
            //int yy = (y + filter2->width/2 - m->height/2);
            //memcpy(&filter2->data[yy*filter2->width+xx], &m->data[y*m->width+x], sizeof(double)*m->width);
            filter2->data[yy*filter2->width+xx] = m->data[y*m->width+x];
        }
    }
    
    //image_save(image_from_matrix(filter2,IMAGE_MINMAX), "filter.png");

    complex_matrix_t*ff_img = matrix_fft(img);
    complex_matrix_t*ff_filter = matrix_fft(filter2);
    matrix_delete(filter2);

    assert(ff_img->width == ff_filter->width);
    assert(ff_img->height == ff_filter->height);
    
//    image_save(image_from_complex_matrix(ff_filter,IMAGE_MINMAX,0), "filter_fft1.png");
//    image_save(image_from_complex_matrix(ff_filter,IMAGE_MINMAX,1), "filter_fft2.png");
//    image_save(image_from_complex_matrix(ff_img,IMAGE_MINMAX,0), "image_fft1.png");
//    image_save(image_from_complex_matrix(ff_img,IMAGE_MINMAX,1), "image_fft2.png");
//
    double size = sqrt(2) / (img->width * img->width);
    for(y=0;y<ff_img->height;y++) 
    for(x=0;x<ff_img->width;x++) {
        complex_t i = ff_img->data[y*ff_img->width+x];
        complex_t f = ff_filter->data[y*ff_img->width+x];
        f.imag = -f.imag; //transpose
        ff_img->data[y*ff_img->width+x].real = (i.real * f.real - i.imag * f.imag) * size;
        ff_img->data[y*ff_img->width+x].imag = (i.real * f.imag + i.imag * f.real) * size;
    }
    matrix_t*result = complex_matrix_ifft_real(ff_img);
    complex_matrix_delete(ff_img);
    complex_matrix_delete(ff_filter);
    return result;
}

static void complex_matrix_multiply_with2(complex_matrix_t*ff_img, complex_matrix_t*ff_filter, double factor, complex_matrix_t*dest)
{
    int x,y;
    for(y=0;y<ff_img->height;y++) {
        complex_t*line1 = &ff_img->data[y*ff_img->width];
        complex_t*line2 = &ff_filter->data[y*ff_img->width];
        complex_t*line3 = &dest->data[y*dest->width];
        for(x=0;x<ff_img->width;x++) {
            complex_t i = line1[x];
            complex_t f = line2[x];
            line3[x].real = (  i.real * f.real + i.imag * f.imag) * factor;
            line3[x].imag = (- i.real * f.imag + i.imag * f.real) * factor;
        }
    }
}

static void complex_matrix_multiply_with(complex_matrix_t*ff_img, complex_matrix_t*ff_filter, complex_matrix_t*dest)
{
    assert(ff_img->width == ff_filter->width && ff_img->width == dest->width);
    assert(ff_img->height == ff_filter->height && ff_img->height == dest->height);
    int t;
    int size = ff_img->width * ff_img->height;
    complex_t*line1 = ff_img->data;
    complex_t*line2 = ff_filter->data;
    complex_t*line3 = dest->data;
    for(t=0;t<size;t++) {
        line3[t].real =   line1[t].real * line2[t].real + line1[t].imag * line2[t].imag;
        line3[t].imag = - line1[t].real * line2[t].imag + line1[t].imag * line2[t].real;
    }
}


typedef struct _convolve_internal {
    complex_matrix_t*ff_filter;
    complex_matrix_t*dest;
    complex_matrix_t*src;
    double size;
    fftw_plan plan;
} convolve_internal_t;

typedef struct _filterbank_internal {
    fftw_plan plan;
    complex_matrix_t*src;
    complex_matrix_t*dest;
    int flags;
} filterbank_internal_t;


convolve_t*convolve_new(matrix_t*filter, int width, int height)
{
    convolve_t*c = malloc(sizeof(convolve_t));
    convolve_internal_t*i = c->internal = malloc(sizeof(convolve_internal_t));
    
    i->dest = complex_matrix_new(width, height);
    i->src = complex_matrix_new(width, height);

    matrix_t*filter2 = matrix_new(width, height);
    assert(width >= filter->width*2);
    assert(height >= filter->height*2);
    int x,y;
    for(y=0;y<filter->height;y++) {
        for(x=0;x<filter->width;x++) {
            int xx = (x + filter2->width - filter->width/2) % filter2->width;
            int yy = (y + filter2->height - filter->height/2) % filter2->height;
            filter2->data[yy*filter2->width+xx] = filter->data[y*filter->width+x];
        }
    }
    i->ff_filter = matrix_fft(filter2);
    i->size = sqrt(2) / (width * width);
    
    matrix_delete(filter2);
    
    i->plan = fftw_plan_dft_2d(i->src->height, i->src->width, (fftw_complex*)i->src->data, 
                               (fftw_complex*)i->dest->data, FFTW_BACKWARD, FFTW_MEASURE);
    return c;
}

convolve_t*convolve_new_complex(complex_matrix_t*filter, int width, int height)
{
    convolve_t*c = malloc(sizeof(convolve_t));
    convolve_internal_t*i = c->internal = malloc(sizeof(convolve_internal_t));

    i->dest = complex_matrix_new(width, height);
    i->src = complex_matrix_new(width, height);

    complex_matrix_t*filter2 = complex_matrix_new(width, height);
    assert(width >= filter->width*2);
    assert(height >= filter->height*2);
    int x,y;
    for(y=0;y<filter->height;y++) {
        for(x=0;x<filter->width;x++) {
            int xx = (x + filter2->width - filter->width/2) % filter2->width;
            int yy = (y + filter2->height - filter->height/2) % filter2->height;
            filter2->data[yy*filter2->width+xx] = filter->data[y*filter->width+x];
        }
    }
    i->ff_filter = complex_matrix_fft(filter2);
    i->size = sqrt(2) / (width * width);

    complex_matrix_delete(filter2);
    
    i->plan = fftw_plan_dft_2d(i->src->height, i->src->width, (fftw_complex*)i->src->data, 
                               (fftw_complex*)i->dest->data, FFTW_BACKWARD, FFTW_MEASURE);
    return c;
}

filterbank_t*filterbank_new(complex_matrixset_t*set, int width, int height, int flags)
{
    filterbank_t*f = malloc(sizeof(filterbank_t));
    filterbank_internal_t*i = malloc(sizeof(filterbank_internal_t));
    f->internal = i;
    f->num = set->num;
    f->width = width;
    f->height = height;
    f->c = malloc(sizeof(complex_matrix_t*)*set->num);

    int size = width*height;
    double factor = sqrt(2) / (width * width);
    int t;
    for(t=0;t<f->num;t++) {
        complex_matrix_t*c = complex_matrix_expand(set->m[t], width, height, flags&EXPAND_RIGHTDOWN);
        f->c[t] = complex_matrix_fft(c);
        int s;
        for(s=0;s<size;s++) {
            f->c[t]->data[s].real *= factor;
            f->c[t]->data[s].imag *= factor;
        }
        complex_matrix_delete(c);
    }
   
    i->flags = flags;
    i->src = complex_matrix_new(width, height);
    i->dest = complex_matrix_new(width, height);
    i->plan = fftw_plan_dft_2d(i->src->height, i->src->width, (fftw_complex*)i->src->data, 
                               (fftw_complex*)i->dest->data, FFTW_BACKWARD, FFTW_MEASURE);
    return f;
}
void filterbank_delete(filterbank_t*f)
{
    filterbank_internal_t*i = (filterbank_internal_t*)f->internal;
    int t;
    for(t=0;t<f->num;t++) {
        complex_matrix_delete(f->c[t]);f->c[t] = 0;
    }
    complex_matrix_delete(i->src);i->src = 0;
    complex_matrix_delete(i->dest);i->dest = 0;
    // todo: free plan
    free(f->c);f->c = 0;
    free(i);f->internal = 0;
}

bytearrayset_t*filterbank_apply_tobytearray(filterbank_t*f, complex_matrix_t*ff_img)
{
    int t;
    filterbank_internal_t*i = (filterbank_internal_t*)f->internal;

    assert(ff_img->width == f->width);
    assert(ff_img->height == f->height);

    bytearrayset_t*set = bytearrayset_new(f->num);
    for(t=0;t<f->num;t++) {
        complex_matrix_multiply_with(ff_img,f->c[t],i->src);
        fftw_execute(i->plan);

        int size = i->dest->width * i->dest->height;
        int s;
        complex_t*data = i->dest->data;
        double min,max;
        min = max = data[0].real*data[0].real + data[0].imag*data[0].imag;

        for(s=1;s<size;s++) {
            double c = data[s].real*data[s].real + data[s].imag*data[s].imag;
            if(i->flags&FILTERBANK_SQRT) 
                c = sqrt(c);
            data[s].real = c;
            if(c<min) min=c;
            if(c>max) max=c;
        }
        if(min==max) max=min+1;
        double f = 127/(max - min);
        unsigned char*dest = malloc16(size);
        for(s=0;s<size;s++) {
            dest[s] = ((int)((data[s].real - min)*f))&127;
        }
        set->m[t] = bytearray_new_fromdata(dest, i->dest->width, i->dest->height);

    }
    return set;
}

matrix_t*convolve_apply(convolve_t*c, complex_matrix_t*ff_img)
{
    convolve_internal_t*i = c->internal;
    assert(ff_img->width == i->ff_filter->width);
    assert(ff_img->height == i->ff_filter->height);
    complex_matrix_multiply_with2(ff_img,i->ff_filter,i->size,i->src);
    fftw_execute(i->plan);
    return complex_matrix_realpart(i->dest);
}
complex_matrix_t*convolve_apply_complex(convolve_t*c, complex_matrix_t*ff_img)
{
    convolve_internal_t*i = c->internal;
    assert(ff_img->width == i->ff_filter->width);
    assert(ff_img->height == i->ff_filter->height);
    complex_matrix_multiply_with2(ff_img,i->ff_filter,i->size,i->src);
    fftw_execute(i->plan);
    return complex_matrix_clone(i->dest);
}

bytearray_t*convolve_apply_complex_tobytearray(convolve_t*c, complex_matrix_t*ff_img)
{
    convolve_internal_t*i = c->internal;
    assert(ff_img->width == i->ff_filter->width);
    assert(ff_img->height == i->ff_filter->height);
    complex_matrix_multiply_with2(ff_img,i->ff_filter,i->size,i->src);
    fftw_execute(i->plan);

    int size = i->dest->width * i->dest->height;
    bytearray_t*b = bytearray_new(i->dest->width, i->dest->height);
    int t;
    complex_t*data = i->dest->data;
    double min,max;
    min=max=data[0].real*data[0].real + data[0].imag*data[0].imag;

    for(t=0;t<size;t++) {
        double c = data[t].real*data[t].real + data[t].imag*data[t].imag;
        data[t].real = c;
        if(c<min) min=c;
        if(c>max) max=c;
    }
    if(min==max) max=min+1;
    double f = 127/(max - min);
    unsigned char*dest = b->data;
    for(t=0;t<size;t++) {
        //dest[t] = (int)(data[t].real - min)*f;
        dest[t] = ((int)((data[t].real - min)*f))&127;
    }
    return b;
}

void convolve_delete(convolve_t*c)
{
    convolve_internal_t*i = c->internal;
    complex_matrix_delete(i->ff_filter);
    complex_matrix_delete(i->dest);
    fftw_destroy_plan(i->plan);
    memset(c->internal, 0, sizeof(convolve_internal_t));
    free(c->internal);c->internal = 0;
    free(c);
}

matrix_t* matrix_inverse_convolve(matrix_t*src_image, matrix_t*m, int xres, int yres, char edgetype)
{
    matrix_t *dst_image = matrix_new(src_image->width*xres,src_image->height*yres);

    int x_src,y_src,x_dst,y_dst;
    int x_mid = m->width / 2;
    int y_mid = m->height / 2;
    for(y_src=0,y_dst=0; y_src<src_image->height; y_src++,y_dst+=yres) {
        for(x_src=0,x_dst=0; x_src<src_image->width; x_src++,x_dst+=xres) {
            int x,y;
            double val = src_image->data[y_src*src_image->width+x_src];
            if(y_dst-y_mid<0 || y_dst-y_mid+m->height>=dst_image->height ||
               x_dst-x_mid<0 || x_dst-x_mid+m->width>=dst_image->width) {
                /* special mode, for reflections */
                for(y=0;y<m->width;y++) {
                    int yy = y_dst + y - y_mid;
                    if(edgetype == EDGE_REFLECT) {
                        if(yy<0) {yy=(-yy-1)%dst_image->height;}
                        if(yy>=dst_image->height) {yy=dst_image->height-1-(yy%dst_image->height);}
                    } else if(edgetype == EDGE_WRAP) {
                        while(yy<0) yy+=dst_image->height;
                        yy%=dst_image->height;
                    } else if(edgetype == EDGE_ZERO) {
                        if(yy<0) continue;
                        if(yy>=dst_image->height) continue;
                    }
                    double*f = &m->data[y*m->height];
                    for(x=0;x<m->width;x++) {
                        int xx = x_dst + x - x_mid;
                        if(edgetype == EDGE_REFLECT) {
                            if(xx<0) {xx=(-xx-1)%dst_image->width;}
                            if(xx>=dst_image->width) {xx=dst_image->width-1-(xx%dst_image->width);}
                        } else if(edgetype == EDGE_WRAP) {
                            while(xx<0) xx+=dst_image->width;
                            xx%=dst_image->width;
                        } else if(edgetype == EDGE_ZERO) {
                            if(xx<0) continue;
                            if(xx>=dst_image->width) continue;
                        }
                        //dst_image->data[yy*dst_image->width + xx] += val*f[m->width-x-1];
                        dst_image->data[yy*dst_image->width + xx] += val*f[x]*xres;
                    }
                }
            } else {
                for(y=0;y<m->width;y++) {
                    double*dst2 = &dst_image->data[x_dst-x_mid+(y_dst+y-y_mid)*dst_image->width];
                    double*f = &m->data[y*m->height];
                    for(x=0;x<m->width;x++) {
                        dst2[x] += val*f[x]*xres;
                    }
                }
            }
        }
    }
    return dst_image;
}

static void matrix_split_xx_init(void*data, double*pos)
{
    matrix_t*_m = (matrix_t*)data;
    int dim = _m->width;
    int t;
    for(t=0;t<dim;t++) {
	pos[t] = 1.0;
    }
}
static void matrix_split_xx_dir(void*data, double*pos, double*dir)
{
    matrix_t*_m = (matrix_t*)data;
    double*m = _m->data;
    int dim = _m->width;
    memset(dir, 0, dim*sizeof(double));
    int s,t;
    /* calculate the gradient direction */
    for(t=0;t<dim;t++) {
        for(s=0;s<dim;s++) {
            dir[t] += 2*pos[s]*(m[s*dim+t] + m[t*dim+s] - 2*pos[t]*pos[s]);
        }
    }
}
static double matrix_split_xx_fitness(void*data, double*v1, double*v2, double f)
{
    matrix_t*_m = (matrix_t*)data;
    double*m = _m->data;
    int dim = _m->width;
    int s,t;
    double diff = 0;
    if(!v2) {
        v2 = v1;f=0;
    }
    for(t=0;t<dim;t++) {
        for(s=0;s<dim;s++) {
            double d = (v1[t]+v2[t]*f)*(v1[s]+v2[s]*f) - m[t*dim+s];
            diff += d*d;
        }
    }
    return diff;
}

double* matrix_split_xx(matrix_t*m)
{
    assert(m->width == m->height);
    return gradient_approximation(m->width, matrix_split_xx_init, matrix_split_xx_dir, matrix_split_xx_fitness, (void*)m, 0);
}

#ifdef HAVE_GSL
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_cblas.h>
#include "gsltools.h"

void matrix_geteigenvectorfactorization(matrix_t*m, matrix_t**dest1, matrix_t**dest2, matrix_t**dest3)
{
    assert(m->width == m->height);
    double*data = (double*)malloc16(sizeof(double)*m->height*m->width);
    memcpy(data, m->data, sizeof(double)*m->height*m->width);

    gsl_matrix_view mm = gsl_matrix_view_array(data, m->height, m->height);
    gsl_vector *s = gsl_vector_alloc(m->width);
    gsl_matrix *v = gsl_matrix_alloc(m->height, m->width);
    gsl_vector *tmp= gsl_vector_alloc(m->width);
    gsl_linalg_SV_decomp(&mm.matrix, v, s, tmp);

    matrix_t*e1 = matrix_new(m->width, m->height);
    matrix_t*e2 = matrix_new(m->width, m->height);
    matrix_t*e3 = matrix_new(m->width, m->height);

    int x,y;
    for(x=0;x<m->width;x++)
    for(y=0;y<m->height;y++) {
        e1->data[y*e1->width+x] = gsl_matrix_get(&mm.matrix, y, x);
        e2->data[y*e2->width+x] = x==y?gsl_vector_get(s, x):0;
        e3->data[y*e3->width+x] = gsl_matrix_get(v, y, x);
    }
    *dest1=e1;
    *dest2=e2;
    *dest3=e3;

    gsl_matrix_free(v);
    gsl_vector_free(s);
    gsl_vector_free(tmp);
}
#define sign(x) ((x)<0?-1:((x)>0?1:0))
void matrix_symm_geteigenvectorfactorization(matrix_t*m, matrix_t**dest1, matrix_t**dest2)
{
    assert(m->width == m->height);
    int n=m->width;
    int x,y;
    for(x=0;x<n;x++)
    for(y=0;y<n;y++) {
        if(fabs(m->data[x*n+y] - m->data[y*n+x]) > 0.001) {
            fprintf(stderr, "Matrix supposed to be symmetric is not symmetric\n");
            exit(1);
        }
    }
   
    matrix_t*dest3;
    matrix_geteigenvectorfactorization(m, dest1, dest2, &dest3);

    /*matrix_print(matrix_multiply_t(matrix_multiply_t(*dest1,1,m,0),0,*dest1,0));
    exit(0);*/
    matrix_t*signmatrix = matrix_multiply_t(*dest1,1,dest3,0);
  
    for(y=0;y<n;y++) {
        for(x=0;x<n;x++) {
            if(fabs(fabs((*dest1)->data[y*n+x]) - fabs(dest3->data[y*n+x])) > 0.001) {
                fprintf(stderr, "Internal error: U!=V^T in symmetric eigenvector factorization (%d,%d: %f<->%f)\n", x, y,
                    (*dest1)->data[x*n+y], dest3->data[x*n+y]);
                exit(1);
            }
            if(x==y) {
                if(fabs(fabs(signmatrix->data[y*n+y])-1) > 0.001) {
                    fprintf(stderr, "Internal error: U^T*V has non -1/+1 diagonal entry (%f)\n", signmatrix[y*n+y]);
                    exit(1);
                }
            } else {
                if(fabs(signmatrix->data[y*n+x]) > 0.001) {
                    fprintf(stderr, "Internal error: U^T*V has off-diagonal entries (%f)\n", signmatrix[y*n+x]);
                    exit(1);
                }
            }
        }
        (*dest2)->data[y*n+y] *= signmatrix->data[y*n+y];
    }
    matrix_delete(signmatrix);
    matrix_delete(dest3);
}
void matrix_symm_geteigenvectorfactorization2(matrix_t*m, matrix_t**dest1, matrix_t**dest2)
{
    assert(m->width == m->height);
    int dim = m->width;
    double*v = malloc16(dim*sizeof(double));
    double*data = (double*)malloc16(sizeof(double)*dim*dim);
    memcpy(data, m->data, sizeof(double)*dim*dim);

    gsl_matrix_view mm = gsl_matrix_view_array(data, dim, dim);
    gsl_vector *vals = gsl_vector_alloc(dim);
    gsl_matrix *vecs = gsl_matrix_alloc(dim, dim);
    
    gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(dim);
    gsl_eigen_symmv(&mm.matrix, vals, vecs, w);
    gsl_eigen_symmv_free(w);
    gsl_eigen_symmv_sort(vals, vecs, GSL_EIGEN_SORT_ABS_DESC);

    matrix_t*diag = matrix_new(dim,dim);
    matrix_t*u = matrix_new(dim,dim);
    int x,y;
    for(y=0;y<dim;y++) {
        diag->data[y*dim+y] = gsl_vector_get(vals, y);
        for(x=0;x<dim;x++) {
            u->data[y*dim+x] = gsl_matrix_get(vecs, y, x);
        }
    }
    *dest1 = u;
    *dest2 = diag;
    free(data);
}

double* matrix_split_xx_2(matrix_t*_m)
{
    assert(_m->width == _m->height);
    int dim = _m->width;
    double*v = malloc16(_m->width*sizeof(double));
    int x;
    for(x=0;x<dim;x++) 
        v[x] = 1.0;
    double maxdiff = matrix_split_xx_fitness((void*)_m, v, 0, 0);
    gsl_matrix_view m = gsl_matrix_view_array(_m->data, dim, dim);
    gsl_vector *vals = gsl_vector_alloc(dim);
    gsl_matrix *vecs = gsl_matrix_alloc(dim, dim);
    
    gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(dim);
    gsl_eigen_symmv(&m.matrix, vals, vecs, w);
    gsl_eigen_symmv_free(w);
    gsl_eigen_symmv_sort(vals, vecs, GSL_EIGEN_SORT_ABS_DESC);

    double val0 = gsl_vector_get(vals, 0);
    double r = sign(val0)*sqrt(fabs(val0));
    for(x=0;x<dim;x++) {
        v[x] = gsl_matrix_get(vecs, x, 0);// * r;
    }
    double currentdiff = matrix_split_xx_fitness((void*)_m, v, 0, 0);
    printf("%4.2f%% precision\n", 100 - (currentdiff * 100 / maxdiff));
    return v;
}

double*matrix_solve(matrix_t*_A, double*_b)
{
    int width = _A->width;
    int height= _A->height;
    double*data = (double*)malloc16(sizeof(double)*height*width);
    memcpy(data, _A->data, sizeof(double)*height*width);
    int t,x,y;
    for(y=0;y<height;y++)
    for(x=0;x<width;x++) {
        if(!__finite(data[y*width+x])) {
            data[y*width+x] = 0;
            fprintf(stderr, "Error: NaN values in solve() input data\n");
        }
/*        _A->data[y*width+x] = 0;
        if(x==y)
            _A->data[y*width+x] = 1;
        data[y*width+x] = _A->data[y*width+x];*/
        //data[y*width+x] += drand48()/100000;
    }
    for(t=0;t<height;t++) {
        if(!__finite(_b[t])) {
            _b[t] = 0;
            fprintf(stderr, "Error: NaN values in solve() input data\n");
        }
    }
   
    gsl_matrix_view A = gsl_matrix_view_array(data, height, width);
    gsl_vector_view b = gsl_vector_view_array(_b, height);
    gsl_vector *xx = gsl_vector_alloc(width);

#if 0
    int s=0,ret;
    gsl_permutation * p = gsl_permutation_alloc(width);
    ret = my_gsl_linalg_LU_decomp (&A.matrix, p, &s);
    if(ret!=0) 
        fprintf(stderr, "Error: gsl_linalg_LU_decomp returned non-zero value\n");
    gsl_linalg_LU_solve (&A.matrix, p, &b.vector, xx);
    if(ret!=0) 
        fprintf(stderr, "Error: gsl_linalg_LU_solve returned non-zero value\n");
    gsl_permutation_free (p);
#else
    gsl_matrix*V = gsl_matrix_alloc(height, width);
    gsl_vector*S = gsl_vector_alloc(width);
    gsl_vector*work = gsl_vector_alloc(width);
    gsl_linalg_SV_decomp(&A.matrix, V, S, work);
    gsl_linalg_SV_solve(&A.matrix, V, S, &b.vector, xx);
    gsl_vector_free(work);
    gsl_vector_free(S);
#endif

    free(data);

    double*result = (double*)malloc16(sizeof(double)*width);
    int error = 0;
    for(t=0;t<width;t++) {
        double r = gsl_vector_get(xx,t);
        if(!__finite(r)) {
            result[t] = 0;
            error++;
        } else {
            result[t] = r;
        }
    }
    if(error) {
        printf("Error: Encountered NaN (%d values, %2d%%) during matrix inversion\n", error, (error*100)/width);
    }
    
    /* checks out!
    double*tmp = (double*)malloc(sizeof(double)*height);
    int xx,yy;
    for(yy=0;yy<height;yy++) {
        double sum = 0;
        for(xx=0;xx<width;xx++) {
            sum += gsl_vector_get(x, xx) * _A->data[yy*width+xx];
        }
        tmp[yy] = sum;
        printf("%f %f\n", _b[yy], tmp[yy]);
    }
    //matrix_t*check = matrix_new_fromdata(tmp, 9, 9);
    //image_save(image_from_matrix(check,IMAGE_MINMAX), "check2.png");
    */
    
    gsl_vector_free (xx);

    return result;

}
double*matrix_solve_approx(matrix_t*_A, double*_b)
{
    if(_A->width > _A->height) {
        return matrix_solve_underdetermined(_A, _b);
    }

    double*data = (double*)malloc16(sizeof(double)*_A->width*_A->height);
    memcpy(data, _A->data, sizeof(double)*_A->width*_A->height);
    gsl_matrix_view Av = gsl_matrix_view_array(data, _A->height, _A->width);
    gsl_vector_view bv = gsl_vector_view_array(_b, _A->height);

    gsl_matrix*A = &Av.matrix;
    gsl_vector*b = &bv.vector;

    gsl_vector *x = gsl_vector_alloc(A->size2);

    gsl_matrix* cov = gsl_matrix_alloc(A->size2, A->size2);
    double chi;
    gsl_multifit_linear_workspace*work = gsl_multifit_linear_alloc(A->size1,A->size2);
    gsl_multifit_linear(A, b, x, cov, &chi, work);
    gsl_multifit_linear_free(work);
    free(data);
    gsl_matrix_free(cov);
    
    double*result = (double*)malloc16(sizeof(double)*A->size1);
    int t;
    for(t=0;t<_A->width;t++) {
        result[t] = gsl_vector_get(x,t);
    }
    gsl_vector_free(x);
    return result;
}

matrix_t* matrix_invert_gsl(matrix_t*A)
{
    int n = A->width;
    assert(A->width == A->height);

    double*data = (double*)malloc16(sizeof(double)*A->width*A->height);
    memcpy(data, A->data, sizeof(double)*A->width*A->height);
    
    int s = 0;
    gsl_matrix_view m = gsl_matrix_view_array(data, A->height, A->width);
    gsl_matrix*inverse = gsl_matrix_alloc(n, n);
    gsl_permutation*perm = gsl_permutation_alloc(n);
    gsl_permutation_init (perm);
    //gsl_linalg_LU_decomp(&m.matrix, perm, &s);
    //gsl_linalg_LU_invert(&m.matrix, perm, inverse);
    gsl_eigen_invert_jacobi(&m.matrix, inverse, 10000);
   
    matrix_t*result = matrix_new(A->width, A->height);
    int x,y;
    for(x=0;x<A->width;x++)
    for(y=0;y<A->height;y++) {
        result->data[y*A->width+x] = gsl_matrix_get(inverse, y, x);
    }
    gsl_matrix_free(inverse);
    free(data);
    return result;
}

matrix_t* matrix_invert(matrix_t*A)
{
    assert(A->width == A->height);
    int n = A->width;
    double*r = malloc16(sizeof(double)*n*n);
    memset(r, 0, sizeof(double)*n*n);
    double*m = malloc16(sizeof(double)*n*n);
    memcpy(m, A->data, sizeof(double)*n*n);
    int* row = malloc16(sizeof(int)*n);
    int x,y,s;
    for(s=0;s<n;s++) {
        row[s] = s*n;
        r[s*n+s] = 1.0;
    }
    int rr;
    for(rr=0;rr<n;rr++) {
        int pivot_row=-1;
        double max = 0.0;
        /* find maximum row entry */
        for(y=rr;y<n;y++) {
            double c = m[row[y]+rr],cabs=fabs(c);
            if(cabs > max) {
                max = cabs;
                pivot_row = y;
            }
        }
        if(pivot_row < 0) {
            free(r);
            free(m);
            free(row);
            return 0;
        }
        if(pivot_row!=rr) {
            int p = row[rr];row[rr] = row[pivot_row];row[pivot_row]=p;
        }
        double pivot = m[row[rr]+rr];

        /*int xx,yy;
        for(yy=0;yy<n;yy++) {
            for(xx=0;xx<n;xx++) {
                if(xx==yy && yy==rr)
                    printf("[%6.3f] ", m[row[yy]+xx]);
                else
                    printf(" %6.3f  ", m[row[yy]+xx]);
            }
            for(xx=0;xx<n;xx++) {
                printf("%6.3f ", r[row[yy]+xx]);
            }
            printf("\n");
        }
        printf("\n");*/
        for(x=rr;x<n;x++) {
            m[row[rr]+x] /= pivot;
        }
        for(x=0;x<n;x++) {
            r[row[rr]+x] /= pivot;
        }
        for(y=rr+1;y<n;y++) {
            double d = m[row[y]+rr];
            m[row[y]+rr] = 0;
            for(x=rr+1;x<n;x++) {
                m[row[y]+x] -= d*m[row[rr]+x];
            }
            for(x=0;x<n;x++) {
                r[row[y]+x] -= d*r[row[rr]+x];
            }
        }
    }
    for(rr=n-1;rr>=0;rr--) {
        double pivot = m[row[rr]+rr]; // should be 1.0
        /*int xx,yy;
        for(yy=0;yy<n;yy++) {
            for(xx=0;xx<n;xx++) {
                if(xx==yy && yy==rr)
                    printf("[%6.3f] ", m[row[yy]+xx]);
                else
                    printf(" %6.3f  ", m[row[yy]+xx]);
            }
            for(xx=0;xx<n;xx++) {
                printf("%6.3f ", r[row[yy]+xx]);
            }
            printf("\n");
        }
        printf("\n");*/
        for(y=rr-1;y>=0;y--) {
            double d = m[row[y]+rr];
            m[row[y]+rr] = 0;
            for(x=rr+1;x<n;x++) {
                m[row[y]+x] -= d*m[row[rr]+x];
            }
            for(x=0;x<n;x++) {
                r[row[y]+x] -= d*r[row[rr]+x];
            }
        }
    }
    for(y=0;y<n;y++) {
        for(x=0;x<n;x++) {
            m[y*n+x] = r[row[y]+x];
        }
    }
    free(r);
    return matrix_new_fromdata(m, n, n);
}

gsl_matrix*matrix_to_gsl(matrix_t*m)
{
    gsl_matrix*g = malloc(sizeof(gsl_matrix));
    g->size1 = m->height;
    g->size2 = m->width;
    g->tda = g->size2;
    g->owner = 0;
    g->block = 0;
    return g;
}

#endif

bytearray_t*bytearray_clone(bytearray_t*m)
{
    bytearray_t*n = (bytearray_t*)malloc(sizeof(bytearray_t));
    n->width = m->width;
    n->height = m->height;
    n->data = (unsigned char*)malloc16(m->width*m->height*sizeof(m->data[0]));
    memcpy(n->data, m->data, m->width*m->height*sizeof(m->data[0]));
    return n;
}
bytearrayset_t* bytearrayset_new(int num)
{
    bytearrayset_t*s = malloc(sizeof(bytearrayset_t));
    s->m = malloc(sizeof(bytearray_t*)*num);
    s->num = num;
    return s;
}
void bytearray_delete(bytearray_t*m)
{
    free(m->data);m->data=0;m->width = 0;m->height = 0;
    free(m);
}
bytearray_t* bytearray_new(int width, int height)
{
    bytearray_t*m = (bytearray_t*)malloc(sizeof(bytearray_t));
    m->data = (unsigned char*)malloc16(width*height);
    memset(m->data, 0, width*height);
    m->width = width; m->height = height;
    return m;
}
bytearray_t* bytearray_new_fromdata(unsigned char*data, int width, int height)
{
    bytearray_t*m = (bytearray_t*)malloc(sizeof(bytearray_t));
    m->data = data;
    m->width = width; m->height = height;
    return m;
}

bytearrayset_t* bytearrayset_new_alloc(int num, int width, int height)
{
    bytearrayset_t*s = malloc(sizeof(bytearrayset_t));
    s->m = malloc(sizeof(bytearray_t*)*num);
    s->num = num;
    int t;
    for(t=0;t<num;t++)
        s->m[t] = bytearray_new(width, height);
    return s;
}
void bytearrayset_extend(bytearrayset_t*s, int num)
{
    s->num+=num;
    s->m = realloc(s->m, sizeof(bytearray_t*)*(s->num));
}
void bytearrayset_append(bytearrayset_t*s, bytearrayset_t*a)
{
    s->num+=a->num;
    s->m = realloc(s->m, sizeof(bytearray_t*)*(s->num));
    int t;
    for(t=0;t<a->num;t++) {
        s->m[s->num-a->num+t] = a->m[t];
    }
}
void bytearrayset_extend_alloc(bytearrayset_t*s, int num, int width, int height)
{
    int oldnum = s->num;
    bytearrayset_extend(s, num);
    int t;
    for(t=oldnum;t<s->num;t++)
        s->m[t] = bytearray_new(width, height);
}

void bytearrayset_delete(bytearrayset_t*b)
{
    int t;
    for(t=0;t<b->num;t++) {
        if(b->m[t]) {
            bytearray_delete(b->m[t]);
            b->m[t] = 0;
        }
    }
    free(b->m);b->m = 0;
}

void bytearray_save(bytearray_t*b, char*filename)
{
    image_t*img = image_new(b->width, b->height);
    int t;
    for(t=0;t<b->width*b->height;t++) {
        img->data[t].r = b->data[t];
        img->data[t].g = b->data[t];
        img->data[t].b = b->data[t];
        img->data[t].a = 255;
    }
    image_save(img, filename);
    image_delete(img);
}


bytearray_t*bytearray_maxfilter(bytearray_t*b)
{
    bytearray_t*tmp = bytearray_clone(b);
    
    int x,y;
    for(x=0;x<b->width;x++) {
        unsigned char*l0 = &tmp->data[x-tmp->width*2];
        unsigned char*l1 = &tmp->data[x-tmp->width];
        unsigned char*l2 = &tmp->data[x];
        unsigned char*l3 = &tmp->data[x+tmp->width];
        unsigned char*l4 = &tmp->data[x+tmp->width*2];
        int end = tmp->width*(tmp->height-2);
        for(y=tmp->width*2;y<end;y+=tmp->width) {
            if(l0[y] >= l2[y] || l1[y] >= l2[y] || l3[y] >= l2[y] || l4[y] >= l2[y])
                l2[y] = 0;
        }
    }

    for(y=0;y<b->height;y++) {
        unsigned char*l1 = &b->data[y*b->width];
        for(x=2;x<b->width-2;x++) {
            if(l1[x-2]>=l1[x] || l1[x-1]>=l1[x] || l1[x+1]>=l1[x] || l1[x+2]>=l1[x])
                l1[x] = 0;
        }
        unsigned char*l2 = &tmp->data[y*b->width];
        for(x=0;x<b->width;x++) {
            if(l2[x]>l1[x])
                l1[x] = l2[x];
        }
    }
    bytearray_delete(tmp);
}

void bytearray_average_rect(bytearray_t*in, int width, int height)
{
    bytearray_t*b = bytearray_clone(in);
    int x,y,t;
    for(y=0;y<b->height;y++) {
        unsigned char*sumdata = &b->data[y*b->width];
        unsigned char*data = &in->data[y*in->width];
        int x=b->width-1;
        int t;
        int sum = 0;
        for(t=1;t<=width;t++) {
            sum += data[x];
            sumdata[x] = sum/t;
            x--;
        }
        do {
            sum += data[x];
            sum -= data[x+width];
            sumdata[x] = sum/width;
        } while(--x>=0);
    }
    memcpy(in->data, b->data, b->width*b->height);
    for(x=0;x<b->width;x++) {
        unsigned char*sumdata = &b->data[x];
        unsigned char*data = &in->data[x];
        int y=b->width*(b->height-1);
        int t;
        int sum = 0;
        int step = height*b->width;
        for(t=1;t<=height;t++) {
            sum += data[y];
            sumdata[y] = sum/t;
            y-=b->width;
        }
        do {
            sum += data[y];
            sum -= data[y+step];
            sumdata[y] = sum/height;
            y-=b->width;
        } while(y>=0);
    }
    memcpy(in->data, b->data, b->width*b->height);
    bytearray_delete(b);
}

void bytearray_normalize(bytearray_t*b)
{
    int size = b->width*b->height;
    unsigned char min = 255;
    unsigned char max = 0;
    int t;
    for(t=0;t<size;t++) {
        if(b->data[t]<min)
            min = b->data[t];
        if(b->data[t]>max)
            max = b->data[t];
    }
    if(min==max) {
        min = min&0xfe;
        max = min|0x01;
    }

    int mul = 0x7f00 / (max-min);
    for(t=0;t<size;t++) {
        b->data[t] = ((b->data[t]-min)*mul)>>8;
    }
}

bytearrayset_t*bytearrayset_select_n(bytearrayset_t*b, int newwidth, int newheight)
{
    int width = b->m[0]->width;
    int height = b->m[0]->height;
    int size = width*height;
    int newsize = newwidth*newheight;
    unsigned char*occupied = malloc16(size);
    int*pos = malloc16(sizeof(int)*newsize);
    memset(occupied, 0, size);
    int t;
    for(t=0;t<newsize;t++) {
        int x,y,p;
        do {
            x = lrand48()%width;
            y = lrand48()%height;
            p = y*width+x;
        } while(occupied[p]);
        pos[t] = p;
        occupied[p] = 1;
    }
    bytearrayset_t*subset = bytearrayset_new_alloc(b->num, newwidth, newheight);
    for(t=0;t<b->num;t++) {
        unsigned char*src=b->m[t]->data;
        unsigned char*dest=subset->m[t]->data;
        int s;
        for(s=0;s<newsize;s++) {
            dest[s] = src[pos[s]];
        }
    }
    free(occupied);
    free(pos);
    return subset;
}

unsigned long* bytearray_integral(bytearray_t*b)
{
    unsigned long *sum = malloc16(sizeof(unsigned long)*b->width*b->height);
    unsigned long *lastline = sum;
    int width = b->width;
    memset(lastline, 0, sizeof(lastline[0])*width);

    int x,y;
    for(y=0;y<b->height;y++)  {
        double v=0;
        unsigned long *line = &sum[y*width];
        unsigned char *bline = &b->data[y*width];

        for(x=0;x<width;x++) {
            v += bline[x];
            line[x] = lastline[x] + v;
        }
        lastline = line;
    }
    return sum;
}

bytearray_t* bytearray_cut(bytearray_t*src, int x, int y, int width, int height)
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
    bytearray_t*dest = bytearray_new(width, height);
    int startx = x;
    int starty = y;
    for(y=0;y<height;y++) {
        unsigned char*sline = &src->data[(y+starty)*src->width+startx];
        unsigned char*dline = &dest->data[y*dest->width];
        memcpy(dline, sline, sizeof(dline[0])*width);
    }
    return dest;
}

bytearray_t* bytearray_combine_maps(bytearray_t*b1, bytearray_t*b2)
{
    int width = b1->width<b2->width?b1->width:b2->width;
    int height = b1->height<b2->height?b1->height:b2->height;

    unsigned char*used = malloc(256 * 256);
    memset(used, 0, 256 * 256);

    int x,y;
    for(y=0;y<height;y++) {
        unsigned char*line1 = &b1->data[y*b1->width];
        unsigned char*line2 = &b2->data[y*b2->width];
        for(x=0;x<width;x++) {
            int p= (line1[x] << 8 | line2[x]);
            used[p] = 1;
        }
    }

    int s,t;

    int pos = 0;

    for(s=0;s<256;s++)
    for(t=0;t<256;t++) {
        int p = s << 8 | t;
        if(used[p]) {
            if(pos>=255) {
                free(used);
                fprintf(stderr, "Can't combine maps %dx%d and %dx%d (more than 255 intersecting values)\n", b1->width, b1->height, b2->width, b2->height);
                return 0;
            }
            used[p] = pos;
            pos++;
        }
    }


    bytearray_t*b = bytearray_new(width, height);

    for(y=0;y<height;y++) {
        unsigned char*line1 = &b1->data[y*b1->width];
        unsigned char*line2 = &b2->data[y*b2->width];
        unsigned char*line = &b->data[y*b->width];
        for(x=0;x<width;x++) {
            if(line1[x]==255 || line2[x]==255)
                line[x] = 255;
            else {
                int p= (line1[x] << 8 | line2[x]);
                line[x] = used[p];
            }
        }
    }
    free(used);
    return b;
}

complex_matrixset_t*complex_matrixset_new(int num)
{
    complex_matrixset_t* set = malloc(sizeof(complex_matrixset_t));
    set->num = num;
    set->m = malloc(sizeof(complex_matrix_t)*num);
    memset(set->m, 0, sizeof(complex_matrix_t)*num);
    return set;
}

complex_matrixset_t*complex_matrixset_clone(complex_matrixset_t*src)
{
    complex_matrixset_t* dest = malloc(sizeof(complex_matrixset_t));
    dest->num = src->num;
    dest->m = malloc(sizeof(complex_matrix_t)*src->num);
    memset(dest->m, 0, sizeof(complex_matrix_t)*src->num);
    int t;
    for(t=0;t<dest->num;t++) {
        dest->m[t] = complex_matrix_clone(src->m[t]);
    }
    return dest;
}

void complex_matrixset_delete(complex_matrixset_t*src)
{
    int t;
    for(t=0;t<src->num;t++) {
        if(src->m[t]) {
            complex_matrix_delete(src->m[t]);
            src->m[t] = 0;
        }
    }
    free(src->m);src->m = 0;
}
    
matrixset_t*matrixset_new(int num)
{
    matrixset_t* set = malloc(sizeof(matrixset_t));
    set->num = num;
    set->m = malloc(sizeof(matrix_t)*num);
    memset(set->m, 0, sizeof(matrix_t)*num);
    return set;
}
matrixset_t* matrixset_new_alloc(int num, int width, int height)
{
    matrixset_t*s = malloc(sizeof(matrixset_t));
    s->m = malloc(sizeof(matrix_t*)*num);
    s->num = num;
    int t;
    for(t=0;t<num;t++)
        s->m[t] = matrix_new(width, height);
    return s;
}


matrixset_t*matrixset_clone(matrixset_t*src)
{
    matrixset_t* dest = malloc(sizeof(matrixset_t));
    dest->num = src->num;
    dest->m = malloc(sizeof(matrix_t)*src->num);
    memset(dest->m, 0, sizeof(matrix_t)*src->num);
    int t;
    for(t=0;t<dest->num;t++) {
        dest->m[t] = matrix_clone(src->m[t]);
    }
    return dest;
}

void matrixset_delete(matrixset_t*src)
{
    int t;
    for(t=0;t<src->num;t++) {
        if(src->m[t]) {
            matrix_delete(src->m[t]);
            src->m[t] = 0;
        }
    }
    free(src->m);src->m = 0;
}

matrix_t* complex_matrix_realpart(complex_matrix_t*src)
{
    matrix_t* dest = matrix_new(src->width, src->height); 
    int t;
    int size = src->width*src->height;
    for(t=0;t<size;t++) {
        dest->data[t] = src->data[t].real;
    }
    return dest;
}

matrix_t* complex_matrix_imagpart(complex_matrix_t*src)
{
    matrix_t* dest = matrix_new(src->width, src->height); 
    int t;
    int size = src->width*src->height;
    for(t=0;t<size;t++) {
        dest->data[t] = src->data[t].imag;
    }
    return dest;
}

complex_matrix_t* matrix_to_complex(matrix_t*src)
{
    complex_matrix_t* dest = complex_matrix_new(src->width, src->height); 
    int t;
    int size = src->width*src->height;
    for(t=0;t<size;t++) {
        dest->data[t].real = src->data[t];
    }
    return dest;
}

matrix_t* matrix_getautocorrelation(matrix_t*src, int csizex, int csizey)
{
    int datalen = (src->width/2+1)*src->height;
    complex_t* data = malloc16(sizeof(complex_t)*datalen);

    fftw_plan plan = fftw_plan_dft_r2c_2d(src->height, src->width, src->data, (fftw_complex*)data, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    
    int t;
    int size = src->width*src->height;
    for(t=0;t<datalen;t++) {
        double r = (data[t].imag * data[t].imag + data[t].real * data[t].real) / size;
        data[t].real = r;
        data[t].imag = 0;
    }

    matrix_t* cov = matrix_new(src->width, src->height);
    plan = fftw_plan_dft_c2r_2d(cov->height, cov->width, (fftw_complex*)data, cov->data, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    free(data);

    matrix_t*cov2 = matrix_new(csizex, csizey);
    int x,y;
    for(y=0;y<csizex;y++) {
        int yy = (y+cov->height-csizey/2)%cov->height;
        for(x=0;x<csizey;x++) {
            int xx = (x+cov->width-csizex/2)%cov->width;
            cov2->data[y*cov2->width+x] = cov->data[yy*cov->width+xx] / size;
        }
    }
    matrix_delete(cov);
    return cov2;
}

complex_matrix_t* complex_matrix_fft(complex_matrix_t*src)
{
    assert(sizeof(fftw_complex) == sizeof(complex_t));
    complex_matrix_t* dest = complex_matrix_new(src->width, src->height);

    fftw_plan plan = fftw_plan_dft_2d(src->height, src->width, (fftw_complex*)src->data, 
                                                               (fftw_complex*)dest->data, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    return dest;
}

complex_matrix_t* complex_matrix_ifft(complex_matrix_t*src)
{
    assert(sizeof(fftw_complex) == sizeof(complex_t));
    complex_matrix_t* dest = complex_matrix_new(src->width, src->height);

    fftw_plan plan = fftw_plan_dft_2d(src->height, src->width, (fftw_complex*)src->data, 
                                                               (fftw_complex*)dest->data, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    return dest;
}

complex_matrix_t* matrix_fft_old(matrix_t*src)
{
    complex_matrix_t* src2 = matrix_to_complex(src);
    complex_matrix_t* result = complex_matrix_fft(src2);
    complex_matrix_delete(src2);
    return result;
}

void fix_fft_complex_matrix(complex_matrix_t*dest)
{
    int y;
    int xs = dest->width/2+1;
    int rxs = dest->width-xs-1;
    int xs4 = xs*4;
    int bx = 2*xs-1-((dest->width&1)^1);
    /* create lower half */
    for(y=dest->height-1;y>=0;y--) {
        complex_t*ind = &dest->data[y*xs];
        complex_t*outd = &dest->data[y*dest->width];
        int x;
        for(x=xs-1;x>=0;x--) {
            outd[x] = ind[x];
        }
    }
    int yr = 0;
    for(y=0;y<dest->height;y++) {
        complex_t*ind = &dest->data[yr*dest->width];
        complex_t*outd = &dest->data[y*dest->width];
        int x;
        for(x=xs;x<dest->width;x++) {
            outd[x].real = ind[bx-x].real;
            outd[x].imag = -ind[bx-x].imag;
            //unsigned long*in = (unsigned long*)&ind[bx-x];
            //unsigned long*out = (unsigned long*)&outd[x];
            //out[0] = in[0];
            //out[1] = in[1];
            //out[2] = in[2];
            //out[3] = in[3]^0x80000000;
        }
        yr=dest->height-y-1;
    }
}

complex_matrix_t* matrix_fft(matrix_t*src)
{
    assert(sizeof(fftw_complex) == sizeof(complex_t));
   // complex_matrix_t* dest = complex_matrix_new((src->width+1)/2, (src->height+1)/2);
    complex_matrix_t* dest = complex_matrix_new(src->width, src->height);

    fftw_plan plan = fftw_plan_dft_r2c_2d(src->height, src->width, src->data, 
                                                               (fftw_complex*)dest->data, 
                                                               FFTW_ESTIMATE);
    fftw_execute(plan);

    fix_fft_complex_matrix(dest);
    fftw_destroy_plan(plan);
    return dest;
}
matrix_t* complex_matrix_ifft_real(complex_matrix_t*src)
{
    complex_matrix_t* c = complex_matrix_ifft(src);
    matrix_t* m = complex_matrix_realpart(c);
    complex_matrix_delete(c);
    return m;
}
/*complex_matrix_t* matrix_ifft(complex_matrix_t*src)
{
    complex_matrix_t* dest = complex_matrix_new(src->width, src->height);
    fftw_plan plan = fftw_plan_dft_2d(src->width, src->height, (fftw_complex*)src->data, (fftw_complex*)dest->data, FFTW_BACKWARD, 0);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    return dest;
}*/

complex_matrix_t* matrix_fft_real(matrix_t*src)
{
    assert(sizeof(fftw_complex) == sizeof(complex_t));
    complex_matrix_t* dest = complex_matrix_new(src->width / 2 + 1, src->height);
    fftw_plan plan = fftw_plan_dft_r2c_2d(src->height, src->width, src->data, (fftw_complex*)dest->data, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    return dest;
}
matrix_t* matrix_ifft_real2(complex_matrix_t*src)
{
    assert(sizeof(fftw_complex) == sizeof(complex_t));
    matrix_t* dest = matrix_new((src->width - 1) *2, src->height);
    fftw_plan plan = fftw_plan_dft_c2r_2d(dest->height, dest->width, (fftw_complex*)src->data, dest->data, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    return dest;
}
matrix_t* complex_matrix_abs(complex_matrix_t*src)
{
    int t;
    int size = src->width*src->height;
    matrix_t*m = matrix_new(src->width, src->height);
    for(t=0;t<size;t++) {
        m->data[t] = sqrt(src->data[t].real*src->data[t].real + src->data[t].imag*src->data[t].imag);
    }
    return m;
}

double* matrix_solve_rectangular(matrix_t*m, double*b)
{
    int tries;
    matrix_t*initial = matrix_new(m->height, m->height);
    double*solution = malloc16(m->width*sizeof(double));
    memset(solution, 0, m->width*sizeof(double));
    char*occupied = malloc16(m->width);
    int*column = malloc(m->height*sizeof(int));
    int*goodcolumns = malloc(sizeof(int)*m->width);
    memset(goodcolumns, 0, sizeof(int)*m->width);
    int num_goodcolumns = 0;

    int x;
    for(x=0;x<m->width;x++) {
        int y;
        double sum = 0;
        for(y=0;y<m->height;y++) {
            double r = m->data[m->width*y+x];
            sum += r*r;
        }
        if(sum>0.0001) 
            goodcolumns[num_goodcolumns++] = x;
    }
    if(num_goodcolumns < m->height) {
        fprintf(stderr, "matrix_solve_rectangular: Couldn't find any base for near-zero matrix\n");
        return solution;
    }

    assert(m->width >= m->height);
    while(1) {
        tries++;
        if(tries==10) {
            fprintf(stderr, "matrix_solve_rectangular: Couldn't find a valid solution for %dx%d matrix\n", m->width, m->height);
            matrix_delete(initial);
            free(occupied);
            free(column);
            return solution;
        }
        int x;
        memset(occupied, 0, m->width);
        for(x=0;x<m->height;x++) {
            int c;
            do {
                c = goodcolumns[lrand48()%num_goodcolumns];
            } while(occupied[c]);
            occupied[c] = 1;
            column[x] = c;
            int y;
            for(y=0;y<m->height;y++) {
                initial->data[y*initial->width+x] = m->data[y*m->width+c];
            }
        }
        matrix_t*inv = matrix_invert(initial);
        if(inv) {
            double*r = matrix_multiply_with_vector(inv, b);
            matrix_delete(inv);
            for(x=0;x<m->height;x++) {
                solution[column[x]] = r[x];
            }
            free(r);
            matrix_delete(initial);
            free(goodcolumns);
            free(occupied);
            free(column);
            return solution;
        }
    }
    matrix_delete(initial);
    free(goodcolumns);
    free(occupied);
    free(column);
    return solution;
}

void vector_print(double*v, int len)
{
    int t;
    for(t=0;t<len;t++) {
        printf("%f ", v[t]);
    }
    printf("\n");
}


matrix_t* matrix_getzerospace(matrix_t*A)
{
    int n = A->width;
    int h = A->height;
    assert(h<=n);
    double*m = malloc(sizeof(double)*n*h);
    memcpy(m, A->data, sizeof(double)*n*h);
    int* row = malloc(sizeof(int)*n);
    int x,y,s;
    int*is_nonzero_column=malloc(sizeof(int)*n);
    memset(is_nonzero_column, 0, sizeof(int)*n);
    int*zero_columns=malloc(sizeof(int)*n);
    int rank = 0;
    int freedom = 0;

    for(s=0;s<h;s++) {
        row[s] = s*n;
    }

    int rr=0;
    int column;
    for(column=0;column<n&&rr<h;column++) {
        int pivot_row=-1;
        double max = 0.0;

        /* find maximum row entry */
        for(y=rr;y<h;y++) {
            double c = m[row[y]+column],cabs=fabs(c);
            if(cabs > max && cabs>0.00001) {
                max = cabs;
                pivot_row = y;
            }
        }
        //printf("pivoting using row %d\n", pivot_row);
        if(pivot_row < 0) {
            zero_columns[freedom++] = column;
            continue;
        } else {
            is_nonzero_column[column] = 1+rank;
            rank++;
        }
        if(pivot_row!=rr) {
            int p = row[rr];row[rr] = row[pivot_row];row[pivot_row]=p;
        }
        double pivot = m[row[rr]+column];

        for(x=column;x<n;x++) {
            m[row[rr]+x] /= pivot;
        }

        for(y=0;y<h;y++) {
            if(row[y]==row[rr])
                continue;
            double d = m[row[y]+column];
            m[row[y]+column] = 0;
            for(x=rr+1;x<n;x++) {
                m[row[y]+x] -= d*m[row[rr]+x];
            }
        }
        
        rr++;
    }
    for(;column<n;column++) {
        zero_columns[freedom++] = column;
    }

    /*{
        printf("rank: %d, freedom: %d\n", rank, freedom);
        int x,y;
        for(x=0;x<n;x++) {
            printf("[%2d ] ", is_nonzero_column[x]);
        }
        printf("\n");
        for(y=0;y<h;y++) {
            for(x=0;x<n;x++) {
                printf("%5.3f ", m[row[y]+x]);
            }
            printf("\n");
        }
        exit(0);
    }*/

    matrix_t*result = matrix_new(freedom, n);
    
    /* now take a look at all zero columns */
    for(x=0;x<freedom;x++) {
        double sum=0;
        for(y=0;y<n;y++) {
            double r = 0;
            if(is_nonzero_column[y]) {
                int pos = is_nonzero_column[y]-1;
                r = -m[row[pos]+zero_columns[x]];
            } else if(y==zero_columns[x]) {
                r = 1;
            } else {
                r = 0;
            }
            sum += r*r;
            result->data[y*result->width+x] = r;
        }
        /* normalize */
        sum = sqrt(sum);
        if(sum > 0.00001) {
            for(y=0;y<n;y++) {
                result->data[y*result->width+x] /= sum;
            }
        }
    }
    free(m);
    free(zero_columns);
    free(is_nonzero_column);
    free(row);
    return result;
}

double* matrix_solve2(matrix_t*m, double*v)
{
    matrix_t*i = matrix_invert(m);
    if(!i)
        return 0;
    double*r = matrix_multiply_with_vector(i, v);
    matrix_delete(i);
    return r;
}

double*matrix_solve_underdetermined(matrix_t*m, double*b)
{
    assert(m->height < m->width);

    double*x0 = matrix_solve_rectangular(m, b);

    matrix_t*z = matrix_getzerospace(m);
    //double*zx0 = matrix_multiply_with_vector(z,x0);
    matrix_t*ztz = matrix_multiply_t(z,1,z,0);
    assert(m->width == z->height);
    assert(z->width == ztz->width);
    //double*ztx0 = matrix_multiply_with_vector_t(z, zx0);
    double*ztx0 = matrix_multiply_with_vector_t(z, x0);
    int t;
    for(t=0;t<z->width;t++) ztx0[t] = -ztx0[t];
    double* w = matrix_solve(ztz, ztx0);
    double* result = matrix_multiply_with_vector(z, w);

    for(t=0;t<z->height;t++) {
        result[t] += x0[t];
    }
    //free(zx0);
    matrix_delete(ztz);
    matrix_delete(z);

    free(w);
    free(ztx0);
    free(x0);
    return result;
}



#ifdef MAIN

int test_eigenvectors()
{
    matrix_t*e1,*e2,*e3;
    int t;
    for(t=2;t<50;t++) {
        matrix_t*m = matrix_new_random(t,t,-1,1);
        int x,y;
        for(y=0;y<t;y++)
        for(x=0;x<y;x++) {
            m->data[y*m->width+x] = m->data[x*m->width+y];
        }
        matrix_symm_geteigenvectorfactorization2(m, &e1, &e2);
        //matrix_symm_geteigenvectorfactorization(m, &e1,&e2);
        matrix_t*tmp0,*tmp1,*tmp2;
        tmp0 = matrix_multiply(e1,e2);
        tmp1 = matrix_multiply_t(tmp0, 0, e1,1);
        tmp2 = matrix_add(tmp1, 1, m, -1);
        for(y=0;y<tmp2->height;y++)
        for(x=0;x<tmp2->width;x++) {
            if(fabs(tmp2->data[y*tmp2->width+x])>0.001) {
                matrix_print(tmp2);
                return 0;
            }
        }
        matrix_delete(tmp0);
        matrix_delete(tmp1);
        matrix_delete(tmp2);
    }
    return 1;
}

int test_inverse()
{
    int t;
    for(t=2;t<50;t++) {
        matrix_t*m = matrix_new_random(t,t,-1,1);
        //matrix_t*i = matrix_invert(m);
        matrix_t*i = matrix_invert(m);
        matrix_t*r = matrix_multiply(m,i);
        int x,y;

        if(r) // result is nonsingular
        for(y=0;y<r->height;y++)
        for(x=0;x<r->width;x++) {
            if(x!=y && fabs(r->data[y*r->width+x])>0.01) {
                printf("Bad Inversion result\n");
                matrix_print(m);
                matrix_print(i);
                matrix_print(r);
                return 0;
            }
        }
        matrix_delete(m);
        matrix_delete(i);
        matrix_delete(r);
    }
    return 1;
}
void test_autocorr()
{
    matrix_t*m = matrix_new_gaussrandom(200,200,0,1);
    double size = (m->width*m->height);

    matrix_t* m2 = matrix_getautocorrelation_old(m, 9, 9);
    matrix_print(m2);

    complex_matrix_t*c = matrix_fft(m);
    int t;
    for(t=0;t<c->width*c->height;t++) {
        double r = (c->data[t].imag * c->data[t].imag + c->data[t].real * c->data[t].real) / size;
        c->data[t].real = r;
        c->data[t].imag = 0;
    }
    matrix_t*cov = complex_matrix_ifft_real(c);

    int x,y;
    for(y=0;y<9;y++) {
        int yy = (y+cov->height-4)%cov->height;
        for(x=0;x<9;x++) {
            int xx = (x+cov->width-4)%cov->width;
            printf("%10.4f ", cov->data[yy*cov->width+xx] / size);
        }
        printf("\n");
    }
    printf("\n");
}


#include <unistd.h>
#include <time.h>
#include "image.h"
#include "stat.h"
#include <sys/times.h>

void test_convolve()
{
    matrix_t*m = image_extractchannel(image_load("../images/flower.png"), IMAGE_GRAY);
    double blabla_data[] = {
    -8.1125000725E-4, 3.9103501476E-3, 1.3462699717E-3, 7.4700999539E-4, 0.0000000000, -7.4700999539E-4, -1.3462699717E-3, -3.9103501476E-3, 8.1125000725E-4, 
    4.4451598078E-3, 4.4565401040E-3, -3.7740699481E-3, -3.6522001028E-4, 0.0000000000, 3.6522001028E-4, 3.7740699481E-3, -4.4565401040E-3, -4.4451598078E-3, 
    1.2316980399E-2, -5.8724298142E-3, 8.2581602037E-3, -2.2522680461E-2, 0.0000000000, 2.2522680461E-2, -8.2581602037E-3, 5.8724298142E-3, -1.2316980399E-2, 
    1.3955879956E-2, -2.8760801069E-3, 3.9442278445E-2, -0.1105690673, 0.0000000000, 0.1105690673, -3.9442278445E-2, 2.8760801069E-3, -1.3955879956E-2, 
    1.4179450460E-2, 8.5267601535E-3, 5.3605638444E-2, -0.1768419296, 0.0000000000, 0.1768419296, -5.3605638444E-2, -8.5267601535E-3, -1.4179450460E-2, 
    1.3955879956E-2, -2.8760801069E-3, 3.9442278445E-2, -0.1105690673, 0.0000000000, 0.1105690673, -3.9442278445E-2, 2.8760801069E-3, -1.3955879956E-2, 
    1.2316980399E-2, -5.8724298142E-3, 8.2581602037E-3, -2.2522680461E-2, 0.0000000000, 2.2522680461E-2, -8.2581602037E-3, 5.8724298142E-3, -1.2316980399E-2, 
    4.4451598078E-3, 4.4565401040E-3, -3.7740699481E-3, -3.6522001028E-4, 0.0000000000, 3.6522001028E-4, 3.7740699481E-3, -4.4565401040E-3, -4.4451598078E-3, 
    -8.1125000725E-4, 3.9103501476E-3, 1.3462699717E-3, 7.4700999539E-4, 0.0000000000, -7.4700999539E-4, -1.3462699717E-3, -3.9103501476E-3, 8.1125000725E-4, 
    };
    matrix_t filter = {blabla_data, 9,9};

    matrix_t*m1 = matrix_convolve(m, &filter, 1, 1, EDGE_ZERO);
    matrix_t*m2 = matrix_convolve_fft(m, &filter);

    struct tms a1,b1,a2,b2;
    times(&a1);
    int t;
    for(t=0;t<20;t++) {
        matrix_delete(matrix_convolve(m, &filter, 1, 1, EDGE_ZERO));
    }
    times(&b1);
    times(&a2);
    for(t=0;t<20;t++) {
        matrix_delete(matrix_convolve_fft(m, &filter));
    }
    times(&b2);
    printf("%d %d\n", b1.tms_utime - a1.tms_utime, b2.tms_utime - a2.tms_utime);

    if(m1) {
        image_save(image_from_matrix(m1,IMAGE_MINMAX), "test1.png");
        statistics_t*stat = statistics_new_frommatrix(m1);
        statistics_print(stat);
    }
    if(m2) {
        image_save(image_from_matrix(m2,IMAGE_MINMAX), "test2.png");
        statistics_t*stat = statistics_new_frommatrix(m2);
        statistics_print(stat);
    }
}

void test_fft()
{
    matrix_t*m = image_extractchannel(image_load("../images/flower.png"), IMAGE_GRAY);
    complex_matrix_t*cm = matrix_to_complex(m);
    complex_matrix_t*fft1 = complex_matrix_fft(cm);
    complex_matrix_t*fft2 = matrix_fft(m);
        
    image_save(image_from_complex_matrix(fft1,IMAGE_MINMAX), "test1.png");
    image_save(image_from_complex_matrix(fft2,IMAGE_MINMAX), "test2.png");
}

int test_zeromatrix()
{
    int t;
    for(t=0;t<1000;t++) {
        matrix_t*m = matrix_new_random(8,4,-10,10);
        if(t>500) {
            int s;
            for(s=0;s<m->width*m->height;s++) {
                if(m->data[s] < 0) m->data[s] = 0;
                else               m->data[s] = 1;
            }
        }

/*        double v[15*15*3] = {
         0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000001, 0.000028, 0.000076, 0.000028, 0.000001, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000001, 0.000208, 0.004169, 0.011332, 0.004169, 0.000208, 0.000001, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000028, 0.004169, 0.083731, 0.227605, 0.083731, 0.004169, 0.000028, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000076, 0.011332, 0.227605, 0.618693, 0.227605, 0.011332, 0.000076, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000028, 0.004169, 0.083731, 0.227605, 0.083731, 0.004169, 0.000028, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000001, 0.000208, 0.004169, 0.011332, 0.004169, 0.000208, 0.000001, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000001, 0.000028, 0.000076, 0.000028, 0.000001, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
         0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000001, 0.000001, 0.000001, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000002, 0.000007, 0.000015, 0.000020, 0.000015, 0.000007, 0.000002, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000001, 0.000006, 0.000032, 0.000113, 0.000239, 0.000307, 0.000239, 0.000113, 0.000032, 0.000006, 0.000001, 0.000000, 0.000000, 0.000000, 0.000000, 0.000006, 0.000053, 0.000307, 0.001072, 0.002270, 0.002915, 0.002270, 0.001072, 0.000307, 0.000053, 0.000006, 0.000000, 0.000000, 0.000000, 0.000002, 0.000032, 0.000307, 0.001768, 0.006171, 0.013064, 0.016775, 0.013064, 0.006171, 0.001768, 0.000307, 0.000032, 0.000002, 0.000000, 0.000000, 0.000007, 0.000113, 0.001072, 0.006171, 0.021539, 0.045599, 0.058550, 0.045599, 0.021539, 0.006171, 0.001072, 0.000113, 0.000007, 0.000000, 0.000001, 0.000015, 0.000239, 0.002270, 0.013064, 0.045599, 0.096532, 0.123950, 0.096532, 0.045599, 0.013064, 0.002270, 0.000239, 0.000015, 0.000001, 0.000001, 0.000020, 0.000307, 0.002915, 0.016775, 0.058550, 0.123950, 0.159155, 0.123950, 0.058550, 0.016775, 0.002915, 0.000307, 0.000020, 0.000001, 0.000001, 0.000015, 0.000239, 0.002270, 0.013064, 0.045599, 0.096532, 0.123950, 0.096532, 0.045599, 0.013064, 0.002270, 0.000239, 0.000015, 0.000001, 0.000000, 0.000007, 0.000113, 0.001072, 0.006171, 0.021539, 0.045599, 0.058550, 0.045599, 0.021539, 0.006171, 0.001072, 0.000113, 0.000007, 0.000000, 0.000000, 0.000002, 0.000032, 0.000307, 0.001768, 0.006171, 0.013064, 0.016775, 0.013064, 0.006171, 0.001768, 0.000307, 0.000032, 0.000002, 0.000000, 0.000000, 0.000000, 0.000006, 0.000053, 0.000307, 0.001072, 0.002270, 0.002915, 0.002270, 0.001072, 0.000307, 0.000053, 0.000006, 0.000000, 0.000000, 0.000000, 0.000000, 0.000001, 0.000006, 0.000032, 0.000113, 0.000239, 0.000307, 0.000239, 0.000113, 0.000032, 0.000006, 0.000001, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000002, 0.000007, 0.000015, 0.000020, 0.000015, 0.000007, 0.000002, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000001, 0.000001, 0.000001, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
         0.000008, 0.000027, 0.000073, 0.000165, 0.000310, 0.000486, 0.000637, 0.000696, 0.000637, 0.000486, 0.000310, 0.000165, 0.000073, 0.000027, 0.000008, 0.000027, 0.000088, 0.000237, 0.000532, 0.000998, 0.001566, 0.002051, 0.002244, 0.002051, 0.001566, 0.000998, 0.000532, 0.000237, 0.000088, 0.000027, 0.000073, 0.000237, 0.000637, 0.001431, 0.002686, 0.004213, 0.005519, 0.006039, 0.005519, 0.004213, 0.002686, 0.001431, 0.000637, 0.000237, 0.000073, 0.000165, 0.000532, 0.001431, 0.003216, 0.006039, 0.009471, 0.012407, 0.013575, 0.012407, 0.009471, 0.006039, 0.003216, 0.001431, 0.000532, 0.000165, 0.000310, 0.000998, 0.002686, 0.006039, 0.011339, 0.017783, 0.023295, 0.025489, 0.023295, 0.017783, 0.011339, 0.006039, 0.002686, 0.000998, 0.000310, 0.000486, 0.001566, 0.004213, 0.009471, 0.017783, 0.027889, 0.036534, 0.039974, 0.036534, 0.027889, 0.017783, 0.009471, 0.004213, 0.001566, 0.000486, 0.000637, 0.002051, 0.005519, 0.012407, 0.023295, 0.036534, 0.047858, 0.052365, 0.047858, 0.036534, 0.023295, 0.012407, 0.005519, 0.002051, 0.000637, 0.000696, 0.002244, 0.006039, 0.013575, 0.025489, 0.039974, 0.052365, 0.057296, 0.052365, 0.039974, 0.025489, 0.013575, 0.006039, 0.002244, 0.000696, 0.000637, 0.002051, 0.005519, 0.012407, 0.023295, 0.036534, 0.047858, 0.052365, 0.047858, 0.036534, 0.023295, 0.012407, 0.005519, 0.002051, 0.000637, 0.000486, 0.001566, 0.004213, 0.009471, 0.017783, 0.027889, 0.036534, 0.039974, 0.036534, 0.027889, 0.017783, 0.009471, 0.004213, 0.001566, 0.000486, 0.000310, 0.000998, 0.002686, 0.006039, 0.011339, 0.017783, 0.023295, 0.025489, 0.023295, 0.017783, 0.011339, 0.006039, 0.002686, 0.000998, 0.000310, 0.000165, 0.000532, 0.001431, 0.003216, 0.006039, 0.009471, 0.012407, 0.013575, 0.012407, 0.009471, 0.006039, 0.003216, 0.001431, 0.000532, 0.000165, 0.000073, 0.000237, 0.000637, 0.001431, 0.002686, 0.004213, 0.005519, 0.006039, 0.005519, 0.004213, 0.002686, 0.001431, 0.000637, 0.000237, 0.000073, 0.000027, 0.000088, 0.000237, 0.000532, 0.000998, 0.001566, 0.002051, 0.002244, 0.002051, 0.001566, 0.000998, 0.000532, 0.000237, 0.000088, 0.000027, 0.000008, 0.000027, 0.000073, 0.000165, 0.000310, 0.000486, 0.000637, 0.000696, 0.000637, 0.000486, 0.000310, 0.000165, 0.000073, 0.000027, 0.000008, 
        };
        m = matrix_new_fromdata(v, 15*15, 3);*/

        matrix_t*z = matrix_getzerospace(m);
        matrix_t* tst = matrix_multiply(m, z);
        int x,y;
        for(y=0;y<tst->height;y++)
        for(x=0;x<tst->width;x++) {
            if(fabs(tst->data[y*tst->width+x])>0.001) {
                matrix_print(tst);
                return 0;
            }
        }
        matrix_delete(m);
        matrix_delete(z);
        matrix_delete(tst);
    }

    return 1;
}

int test_solve2()
{
    int step;
    for(step=0;step<1000;step++) {
        matrix_t*m = matrix_new_random(50,10,-10,10);
        matrix_t*bb = matrix_new_random(1,10,-10,10);
        double*b = bb->data;
        double*initial = matrix_solve_rectangular(m, b);
        double*rr = matrix_multiply_with_vector(m, initial);
        int t;
        for(t=0;t<m->height;t++) {
            if(fabs(rr[t]-b[t]) > 0.001) {
                matrix_print(m);
                return 0;
            }
        }
        free(b);
        free(initial);
        free(rr);
        matrix_delete(m);
        matrix_delete(bb);
    }
    return 1;
}

int test_mult()
{
    int step;
    for(step=0;step<1000;step++) {
        matrix_t*m1 = matrix_new_random(50,50,-10,10);
        matrix_t*m2= matrix_new_random(50,50,-10,10);

        matrix_t*m1t = matrix_transpose(m1);
        matrix_t*m2t = matrix_transpose(m2);

        int t1,t2;
        for(t1=0;t1<2;t1++)
        for(t2=0;t2<2;t2++) {
            matrix_t*c1 = matrix_multiply(t1?m1t:m1, t2?m2t:m2);
            matrix_t*c2 = matrix_multiply_t(m1, t1, m2, t2);
            double x = matrix_diff(c1,c2);
            if(x > 0.001) {
                printf("test_mult failed %f\n", x);
                return 0;
            }
            matrix_delete(c1);
            matrix_delete(c2);
        }
        matrix_delete(m1t);
        matrix_delete(m2t);
        matrix_delete(m1);
        matrix_delete(m2);
    }
    return 1;
}

int test_quadopt()
{
    int step;
    for(step=0;step<10;step++) {
        matrix_t*m = matrix_new_random(50,10,-10,10);

        /*double v[15*15*3] = {
         0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000001, 0.000028, 0.000076, 0.000028, 0.000001, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000001, 0.000208, 0.004169, 0.011332, 0.004169, 0.000208, 0.000001, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000028, 0.004169, 0.083731, 0.227605, 0.083731, 0.004169, 0.000028, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000076, 0.011332, 0.227605, 0.618693, 0.227605, 0.011332, 0.000076, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000028, 0.004169, 0.083731, 0.227605, 0.083731, 0.004169, 0.000028, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000001, 0.000208, 0.004169, 0.011332, 0.004169, 0.000208, 0.000001, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000001, 0.000028, 0.000076, 0.000028, 0.000001, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
         0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000001, 0.000001, 0.000001, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000002, 0.000007, 0.000015, 0.000020, 0.000015, 0.000007, 0.000002, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000001, 0.000006, 0.000032, 0.000113, 0.000239, 0.000307, 0.000239, 0.000113, 0.000032, 0.000006, 0.000001, 0.000000, 0.000000, 0.000000, 0.000000, 0.000006, 0.000053, 0.000307, 0.001072, 0.002270, 0.002915, 0.002270, 0.001072, 0.000307, 0.000053, 0.000006, 0.000000, 0.000000, 0.000000, 0.000002, 0.000032, 0.000307, 0.001768, 0.006171, 0.013064, 0.016775, 0.013064, 0.006171, 0.001768, 0.000307, 0.000032, 0.000002, 0.000000, 0.000000, 0.000007, 0.000113, 0.001072, 0.006171, 0.021539, 0.045599, 0.058550, 0.045599, 0.021539, 0.006171, 0.001072, 0.000113, 0.000007, 0.000000, 0.000001, 0.000015, 0.000239, 0.002270, 0.013064, 0.045599, 0.096532, 0.123950, 0.096532, 0.045599, 0.013064, 0.002270, 0.000239, 0.000015, 0.000001, 0.000001, 0.000020, 0.000307, 0.002915, 0.016775, 0.058550, 0.123950, 0.159155, 0.123950, 0.058550, 0.016775, 0.002915, 0.000307, 0.000020, 0.000001, 0.000001, 0.000015, 0.000239, 0.002270, 0.013064, 0.045599, 0.096532, 0.123950, 0.096532, 0.045599, 0.013064, 0.002270, 0.000239, 0.000015, 0.000001, 0.000000, 0.000007, 0.000113, 0.001072, 0.006171, 0.021539, 0.045599, 0.058550, 0.045599, 0.021539, 0.006171, 0.001072, 0.000113, 0.000007, 0.000000, 0.000000, 0.000002, 0.000032, 0.000307, 0.001768, 0.006171, 0.013064, 0.016775, 0.013064, 0.006171, 0.001768, 0.000307, 0.000032, 0.000002, 0.000000, 0.000000, 0.000000, 0.000006, 0.000053, 0.000307, 0.001072, 0.002270, 0.002915, 0.002270, 0.001072, 0.000307, 0.000053, 0.000006, 0.000000, 0.000000, 0.000000, 0.000000, 0.000001, 0.000006, 0.000032, 0.000113, 0.000239, 0.000307, 0.000239, 0.000113, 0.000032, 0.000006, 0.000001, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000002, 0.000007, 0.000015, 0.000020, 0.000015, 0.000007, 0.000002, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000001, 0.000001, 0.000001, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
         0.000008, 0.000027, 0.000073, 0.000165, 0.000310, 0.000486, 0.000637, 0.000696, 0.000637, 0.000486, 0.000310, 0.000165, 0.000073, 0.000027, 0.000008, 0.000027, 0.000088, 0.000237, 0.000532, 0.000998, 0.001566, 0.002051, 0.002244, 0.002051, 0.001566, 0.000998, 0.000532, 0.000237, 0.000088, 0.000027, 0.000073, 0.000237, 0.000637, 0.001431, 0.002686, 0.004213, 0.005519, 0.006039, 0.005519, 0.004213, 0.002686, 0.001431, 0.000637, 0.000237, 0.000073, 0.000165, 0.000532, 0.001431, 0.003216, 0.006039, 0.009471, 0.012407, 0.013575, 0.012407, 0.009471, 0.006039, 0.003216, 0.001431, 0.000532, 0.000165, 0.000310, 0.000998, 0.002686, 0.006039, 0.011339, 0.017783, 0.023295, 0.025489, 0.023295, 0.017783, 0.011339, 0.006039, 0.002686, 0.000998, 0.000310, 0.000486, 0.001566, 0.004213, 0.009471, 0.017783, 0.027889, 0.036534, 0.039974, 0.036534, 0.027889, 0.017783, 0.009471, 0.004213, 0.001566, 0.000486, 0.000637, 0.002051, 0.005519, 0.012407, 0.023295, 0.036534, 0.047858, 0.052365, 0.047858, 0.036534, 0.023295, 0.012407, 0.005519, 0.002051, 0.000637, 0.000696, 0.002244, 0.006039, 0.013575, 0.025489, 0.039974, 0.052365, 0.057296, 0.052365, 0.039974, 0.025489, 0.013575, 0.006039, 0.002244, 0.000696, 0.000637, 0.002051, 0.005519, 0.012407, 0.023295, 0.036534, 0.047858, 0.052365, 0.047858, 0.036534, 0.023295, 0.012407, 0.005519, 0.002051, 0.000637, 0.000486, 0.001566, 0.004213, 0.009471, 0.017783, 0.027889, 0.036534, 0.039974, 0.036534, 0.027889, 0.017783, 0.009471, 0.004213, 0.001566, 0.000486, 0.000310, 0.000998, 0.002686, 0.006039, 0.011339, 0.017783, 0.023295, 0.025489, 0.023295, 0.017783, 0.011339, 0.006039, 0.002686, 0.000998, 0.000310, 0.000165, 0.000532, 0.001431, 0.003216, 0.006039, 0.009471, 0.012407, 0.013575, 0.012407, 0.009471, 0.006039, 0.003216, 0.001431, 0.000532, 0.000165, 0.000073, 0.000237, 0.000637, 0.001431, 0.002686, 0.004213, 0.005519, 0.006039, 0.005519, 0.004213, 0.002686, 0.001431, 0.000637, 0.000237, 0.000073, 0.000027, 0.000088, 0.000237, 0.000532, 0.000998, 0.001566, 0.002051, 0.002244, 0.002051, 0.001566, 0.000998, 0.000532, 0.000237, 0.000088, 0.000027, 0.000008, 0.000027, 0.000073, 0.000165, 0.000310, 0.000486, 0.000637, 0.000696, 0.000637, 0.000486, 0.000310, 0.000165, 0.000073, 0.000027, 0.000008, 
        };
        m = matrix_new_fromdata(v, 15*15, 3);*/

        matrix_t*bb = matrix_new_random(1,10,-10,10);
        double*b = bb->data;

        double*x = matrix_solve_underdetermined(m, b);

        double* test = matrix_multiply_with_vector(m, x);
        int t;
        for(t=0;t<m->height;t++) {
            if(fabs(test[t]-b[t]) > 0.001) {
                printf("quadopt test failed, matrix: \n");
                matrix_print(m);
                return 0;
            }
        }
        matrix_delete(m);
        matrix_delete(bb);
        free(b);
    }
    return 1;
}
int test_newfft()
{
    if(0) {
        matrix_t* m = matrix_new_gaussrandom(4,3,0,10);
        complex_matrix_print(matrix_fft(m));
        printf("---------------\n");
        complex_matrix_print(matrix_fft_old(m));
        exit(0);
    }
       
    int t;
    matrix_t* m = matrix_new_gaussrandom(512,512,0,1);
    struct tms a1,b1,a2,b2;
    times(&a1);
    for(t=0;t<100;t++) {
        complex_matrix_delete(matrix_fft(m));
    }
    times(&b1);
    times(&a2);
    for(t=0;t<100;t++) {
        complex_matrix_delete(matrix_fft_old(m));
    }
    times(&b2);
    printf("%d %d\n", b1.tms_utime - a1.tms_utime, b2.tms_utime - a2.tms_utime);

    for(t=0;t<40;t++) {
        int width=(lrand48()%64)+1;
        int height=(lrand48()%64)+1;
        matrix_t* m = matrix_new_gaussrandom(width,height,0,1);
        complex_matrix_t*m1 = matrix_fft(m);
        complex_matrix_t*m2 = matrix_fft_old(m);
        int t;
        for(t=0;t<m1->width*m1->height;t++) {
            if(fabs(m1->data[t].real-m2->data[t].real) > 0.0001)
                return 0;
            if(fabs(m1->data[t].real-m2->data[t].real) > 0.0001)
                return 0;
        }
        complex_matrix_delete(m1);
        complex_matrix_delete(m2);
        matrix_delete(m);
    }
    return 1;
}

int test_average()
{
    bytearray_t*b = bytearray_new(320,199);

    int t;
    for(t=0;t<100;t++) {
        int distx=(lrand48()%10)+2,disty=(lrand48()%10)+2;
        bytearray_t*b2 = bytearray_clone(b);
        bytearray_average_rect(b2, distx,disty);
        int xx,yy;
        for(yy=0;yy<b->height;yy+=disty)
        for(xx=0;xx<b->width;xx+=distx) {
            int x,y;
            int dx = b->width-xx;
            int dy = b->height-yy;
            if(dx>distx)
                dx = distx;
            if(dy>disty)
                dy = disty;
            int sum = 0;
            int count = 0;
            for(y=0;y<disty;y++)
            for(x=0;x<distx;x++) {
                sum += b->data[(yy+y)*b2->width+xx+x];
                count++;
            }
            int diff = abs(b2->data[yy*b2->width+xx] - sum/count);
            if(diff>16) {
                printf("test_average: discrepancy between averaged values: %d!=%d, during %dx%d average, at position (%d,%d) of (%d,%d)\n", 
                        b2->data[yy*b2->width+xx], sum/count, distx, disty, xx,yy, b->width, b->height);
                return 0;
            }
        }
        bytearray_delete(b2);
    }
    return 1;
}
    
int run_tests()
{
    srand48(time(0));
    printf("test_average\n");
    if(!test_average()) return 0;
    printf("test_newfft\n");
    if(!test_newfft()) return 0;
    printf("test_mult\n");
    if(!test_mult()) return 0;
    printf("test_quadopt\n");
    if(!test_quadopt()) return 0;
    printf("test_zeromatrix\n");
    if(!test_zeromatrix()) return 0;
    printf("test_solve2\n");
    if(!test_solve2()) return 0;
    printf("test_eigenvectors\n");
    if(!test_eigenvectors()) return 0;
    printf("test_inverse\n");
    if(!test_inverse()) return 0;
    return 1;
}

int main(int argn, char*argv[])
{
    matrix_t*m = image_extractchannel(image_load("../images/gazelle.png"), IMAGE_YUV_Y);
    matrix_t*m1 = matrix_scaleup(m, m->width*2, m->height*2);
    return;

    if(!run_tests()) printf("one or more tests failed\n");
    else printf("all tests ok\n");

    //test_autocorr();
    //test_convolve();
    //test_fft();

}

#endif
