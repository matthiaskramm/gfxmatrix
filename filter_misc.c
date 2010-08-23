/* filter_misc.c

   Various filter kernels

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
#include <assert.h>
#include <math.h>
#include "matrix.h"
#include "filter.h"
#include "image.h"

static void normalize(matrix_t*m)
{
    double sum = 0;
    int t;
    int l = m->width*m->height;
    for(t=0;t<l;t++) {
        sum += (m->data[t]*m->data[t]);
    }
    if(fabs(sum<0.0001))
        return;
    //sum = sqrt(sum);

    for(t=0;t<l;t++) {
        m->data[t] /= sum;
    }
}

matrix_t* gabor_filter(int width, int height, float ax, float ay, float phase_shift, float ifreq, float angle, char complex_counterpart)
{
    matrix_t*m = matrix_new(width, height);
    int x,y;
    float ccos = cos(angle*M_PI/180);
    float csin = sin(angle*M_PI/180);
    int x0 = width/2;
    int y0 = height/2;
    for(y=0;y<height;y++)
    for(x=0;x<width;x++) {
        float xx = (x-x0)*ccos + (y-y0)*csin;
        float yy = (x-x0)*csin - (y-y0)*ccos;
        float r = ax*ax*xx*xx + ay*ay*yy*yy;
        if(complex_counterpart) 
            m->data[y*width+x] = exp(-r)*cos(2*M_PI*ifreq*xx+phase_shift);
        else
            m->data[y*width+x] = exp(-r)*sin(2*M_PI*ifreq*xx+phase_shift);
    }
    normalize(m);
    return m;
}

static void normalize_complex(complex_matrix_t*m)
{
    double sum = 0;
    int t;
    int l = m->width*m->height;
    for(t=0;t<l;t++) {
        sum += (m->data[t].real*m->data[t].real + m->data[t].imag*m->data[t].imag);
    }
    if(fabs(sum<0.0001))
        return;
    //sum = sqrt(sum);

    for(t=0;t<l;t++) {
        m->data[t].real /= sum;
        m->data[t].imag /= sum;
    }
}

complex_matrix_t* gabor2_filter(int width, int height, float a, float f1, float f2)
{
    complex_matrix_t*m = complex_matrix_new(width, height);
    int x,y;
    int x0 = width/2;
    int y0 = height/2;
    for(y=0;y<height;y++)
    for(x=0;x<width;x++) {
        float xx = (x-x0);
        float yy = (y-y0);
        float r = a*a*xx*xx + a*a*yy*yy;
        m->data[y*width+x].real = exp(-r)*cos(2*M_PI*(f1*xx+f2*yy) / width);
        m->data[y*width+x].imag = exp(-r)*sin(2*M_PI*(f1*xx+f2*yy) / width);
    }
    normalize_complex(m);
    return m;
}


matrix_t* gauss_diff(int width, int height, float a, float angle, float len)
{
    matrix_t*m = matrix_new(width, height);
    int x,y;
    float ccos = cos(angle*M_PI/180);
    float csin = sin(angle*M_PI/180);
    float x0 = width/2 - ccos*len;
    float y0 = height/2 - csin*len;
    float x1 = width/2 + ccos*len;
    float y1 = height/2 + csin*len;
    for(y=0;y<height;y++)
    for(x=0;x<width;x++) {
        float r0 = a*a*(x-x0)*(x-x0) + a*a*(y-y0)*(y-y0);
        float r1 = a*a*(x-x1)*(x-x1) + a*a*(y-y1)*(y-y1);
        float e0 = exp(-r0);
        float e1 = exp(-r1);
        m->data[y*width+x] = e0-e1;
    }
    normalize(m);
    return m;
}

/* derivative of a gaussian */
matrix_t* gauss_diff2(int width, int height, float d, float angle)
{
    matrix_t*m = matrix_new(width, height);
    int x,y;
    float x0 = width/2;
    float y0 = height/2;
    float ccos = cos(angle*M_PI/180);
    float csin = sin(angle*M_PI/180);
    for(y=0;y<height;y++)
    for(x=0;x<width;x++) {
        float xx = x-x0;
        float yy = y-y0;
        float r = (xx*xx+yy*yy)/(2*d*d);
        float l = ccos*xx + csin*yy;
        float e = l*exp(-r);
        m->data[y*width+x] = e;
    }
    normalize(m);
    return m;
}

matrix_t* gauss_highpass(int width, int height, float r0, float a)
{
    matrix_t*m = matrix_new(width, height);
    int x,y;
    float x0 = width/2;
    float y0 = height/2;
    for(y=0;y<height;y++)
    for(x=0;x<width;x++) {
        float xx = x-x0;
        float yy = y-y0;
        float r = (sqrt(xx*xx+yy*yy) - r0)*a;
        float e = exp(-r*r);
        m->data[y*width+x] = e;
    }
    normalize(m);
    return m;
}

matrix_t* gauss_filter(int width, int height, float d)
{
    matrix_t*m = matrix_new(width, height);
    int x,y;
    float x0 = width/2;
    float y0 = height/2;
    for(y=0;y<height;y++)
    for(x=0;x<width;x++) {
        float xx = x-x0;
        float yy = y-y0;
        float r = (xx*xx+yy*yy) / (2*d*d);
        float e = exp(-r);
        m->data[y*width+x] = e;
    }
    normalize(m);
    return m;
}

matrix_t* lgauss_filter(int width, int height, float d)
{
    matrix_t*m = matrix_new(width, height);
    int x,y;
    float x0 = width/2;
    float y0 = height/2;
    for(y=0;y<height;y++)
    for(x=0;x<width;x++) {
        float xx = x-x0;
        float yy = y-y0;
        float r = (xx*xx+yy*yy) / (2*d*d);
        float e = (1.0/(M_PI*d*d*d*d))*(1-r)*exp(-r);
        m->data[y*width+x] = e;
    }
    normalize(m);
    return m;
}

#define ADD_MATRIX(m) {if(pass==1) {filter_t*f = filter_new(filtertype_convolve, FILTER_DIRECTIONAL, m, 1, 1);base->children[num] = filtertree_new(0, f);} num++;}
#define ADD_COMPLEX_MATRIX(m) {if(pass==1) {filter_t*f = filter_new(filtertype_convolve, FILTER_DIRECTIONAL|FILTER_COMPLEX, (matrix_t*)m, 1, 1);base->children[num] = filtertree_new(0, f);} num++;}

static void check_filtertree_uniqueness(filtertree_t*base)
{
    int s,t;
    for(t=0;t<base->num_children;t++) {
        for(s=0;s<base->num_children;s++) {
            if(s!=t && matrix_diff(base->children[t]->filter->matrix, base->children[s]->filter->matrix) < 0.001)
                printf("matrix %d and %d are equal\n", s, t);
        }
    }
}

filtertree_t* makeGaborTree(int width, int height, int xx, int yy)
{
    int pass;
    int num;
    filtertree_t*base = 0;

    for(pass=0;pass<2;pass++) {
        if(pass==1)
            base = filtertree_new(num, 0);
        num = 0;

        int x,y;
        for(x=0;x<xx;x++) {
            for(y=0;y<yy;y++) {
                if(x||y) {
                    ADD_COMPLEX_MATRIX(gabor2_filter(width,height, 80.0 / (width*height), x, y));
                }
            }
        }
    }
    check_filtertree_uniqueness(base);
    return base;
}

complex_matrix_t*gauss_diff_complex(int width, int height, double a)
{
    matrix_t*g1 = gauss_diff2(width, height, a, 0);
    matrix_t*g2 = gauss_diff2(width, height, a, 90);
    complex_matrix_t*g = complex_matrix_new(g1->width, g1->height);
    int t;
    for(t=0;t<g->width*g->height;t++) {
        g->data[t].real = g1->data[t];
        g->data[t].imag = g2->data[t];
    }
    matrix_delete(g1);
    matrix_delete(g2);
    return g;
}

complex_matrixset_t* filterset_new_dct(int width, int height, int stepx, int stepy, char complexversion)
{
    complex_matrixset_t* f = complex_matrixset_new((width/stepx)*(height/stepy));

    int x1,y1,x2,y2;

    int size = width*height;

    int pos = 0;
    for(y1=0;y1<height;y1+=stepy)
    for(x1=0;x1<width;x1+=stepx) {
        f->m[pos] = complex_matrix_new(width,height);
        complex_t*data = f->m[pos]->data;
        double ax = 3.14159265358979*x1/(double)(2*width);
        double ay = 3.14159265358979*y1/(double)(2*height);
        double sum = 0;

        double c = 0.25;
        if(!x1) c *= 0.70710678118654746;
        if(!y1) c *= 0.70710678118654746;

        for(x2=0;x2<width;x2++)
        for(y2=0;y2<height;y2++) {
            double xx2 = (x2*2+1);
            double yy2 = (y2*2+1);
            if(complexversion) {
                data[y2*width+x2].real = c*(cos(ax*xx2)*cos(ay*yy2) - sin(ax*xx2)*sin(ay*yy2));
                data[y2*width+x2].imag = c*(cos(ax*xx2)*sin(ay*yy2) + sin(ax*xx2)*cos(ay*yy2));
            } else {
                data[y2*width+x2].real = c*(cos(ax*xx2)*cos(ay*yy2));
                data[y2*width+x2].imag = 0;
            }
        }
        pos++;
    }
    return f;
}
