/* gradopt.h 

   Autocorrelation and crosscorrelation adjustment.

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
#include "corr.h"
#include "gradopt.h"
#include "stat.h"

#define sign(x) ((x)<0?-1:((x)>0?1:0))

static int debug = 0;

static inline double sqr(double x)
{
    return x*x;
}

void matrix_adjustautocorrelation(matrix_t*img, matrix_t*corr)
{
    /* size of the autocorrelation area to adjust */
    int csizex = corr->width;
    int csizey = corr->height;
    int csize = csizex*csizey;
    int midx = csizex/2;
    int midy = csizey/2;

    /* size of our filter */
    int fsizex = csizex;
    int fsizey = csizey;
    int fsize = fsizex*fsizey;
    int fmidx = fsizex/2;
    int fmidy = fsizey/2;

    /* size of the correlation area we need from the image meant to be adjusted */
    int rx = csizex*3;
    int ry = csizey*3;
    int rmidx = rx/2;
    int rmidy = ry/2;

    if(debug)
        printf("getting %dx%d autocorrelation from source image (%dx%d)\n", rx, ry, img->width, img->height);
    matrix_t*r = matrix_getautocorrelation(img, rx, ry);
    if(debug) {
        image_save(image_from_matrix(r,IMAGE_MINMAX), "autocorrelation_source.png");
    image_save(image_from_matrix(corr,IMAGE_MINMAX), "autocorrelation_dest.png");
    }

    //matrix_t*A = matrix_new(fsize*fsize, csize);
    matrix_t*A = matrix_new(fsize*fsize, csize*csize);

    double*b = malloc(A->height*sizeof(double));
    memset(b, 0, A->height*sizeof(double));
    memcpy(b, corr->data, csize*sizeof(double));

    int x1,y1,x2,y2,xx,yy;
    int row = 0;
    for(yy=0;yy<csizey;yy++)
    for(xx=0;xx<csizex;xx++) {
        if(row<A->height) {
            for(y2=0;y2<fsizey;y2++)
            for(x2=0;x2<fsizex;x2++) {
                for(y1=0;y1<fsizey;y1++)
                for(x1=0;x1<fsizex;x1++) {
                    int column = (y2*fsizex+x2)*fsize + (y1*fsizex+x1); // variable nr
                    int u1 = (x1-fmidx)+(xx-midx)-(x2-fmidx)+rmidx;
                    int v1 = (y1-fmidy)+(yy-midy)-(y2-fmidy)+rmidy;
                    assert(u1>=0 && u1<r->width);
                    assert(v1>=0 && v1<r->height);
                    A->data[row*A->width + column] = r->data[v1*r->width+u1];
                }
            }
        }
        row++;
    }
    if(debug) {
        image_save(image_from_matrix(A,IMAGE_MINMAX), "circular_correlation.png");
    }

#if 0
    printf("checking matrix for duplicate rows\n");
    matrix_t*sim = matrix_new(A->height, A->height);
    for(x1=0;x1<A->height;x1++)
    for(x2=0;x2<A->height;x2++) {
        int r=0;
        double*line1 = &A->data[x1*A->width];
        double*line2 = &A->data[x2*A->width];
        int count = 0;
        for(r=0;r<A->width;r++) {
            if(fabs(line1[r]-line2[r])<0.1)
                count++;
        }
        if(count == A->width)
            sim->data[x1*sim->width+x2] = 1.0;
        else
            sim->data[x1*sim->width+x2] = 0.0;
    }
    image_save(image_from_matrix(sim,IMAGE_MINMAX), "sim.png");
#endif

    matrix_delete(r);
    matrix_t*l = (matrix_t*)malloc(sizeof(matrix_t));
    if(debug) {
        printf("solving %dx%d matrix for correlation filter matrix of size %d\n", A->width, A->height, A->width);
    }
    //l->data = matrix_solve_approx(A, b);
    l->data = matrix_solve(A, b);
    l->width = fsize;
    l->height = fsize;

    if(debug) {
        image_save(image_from_matrix(l,IMAGE_MINMAX), "basefilter.png");
        matrix_t*check = matrix_new_fromdata(matrix_multiply_with_vector(A, l->data), csizex, csizey);
        image_save(image_from_matrix(check,IMAGE_MINMAX), "autocorrelation_dest_check.png");
    }

    matrix_delete(A);

    if(debug) {
        printf("resolving %dx%d correlation filter matrix into filter product of size %dx%d\n", l->width, l->height, fsizex, fsizey);
    }
    matrix_t*filter  = (matrix_t*)malloc(sizeof(matrix_t));
    filter->data = matrix_split_xx(l);
    //filter->data = matrix_split_xx_2(l);
    filter->width = fsizex;
    filter->height = fsizey;
    if(debug) {
        image_save(image_from_matrix(filter,IMAGE_MINMAX), "filter.png");
        matrix_t* tst = matrix_new(fsize,fsize);
        for(y1=0;y1<fsize;y1++) 
        for(x1=0;x1<fsize;x1++) 
        tst->data[y1*fsize+x1]=filter->data[x1]*filter->data[y1];
        image_save(image_from_matrix(tst,IMAGE_MINMAX), "basefilter_approx.png");
    }

    matrix_delete(l);
    matrix_t*newimage = matrix_convolve(img, filter, 1, 1, EDGE_WRAP);
    matrix_delete(filter);
    if(debug) {
        matrix_t*newcorr = matrix_getautocorrelation(newimage, csizex,csizey);
        image_save(image_from_matrix(newcorr,IMAGE_MINMAX), "autocorrelation_dest_approx.png");
    }

    assert(newimage->width == img->width && newimage->height == img->height);
    memcpy(img->data, newimage->data, img->width*img->height*sizeof(img->data[0]));
}

typedef struct _mopt
{
    matrix_t*m;
    double*b;
    int len;
    double centerforce;
} mopt_t;

//#define FIX0

static void matrix_find_xx_dir(void*data, double*pos, double*dir)
{
    mopt_t*d = (mopt_t*)data;
    matrix_t*m = d->m;
    int len = d->len;
    memset(dir, 0, d->len*sizeof(double));
    int i,n,j,k;
    int start = 0;
    
    /* calculate the gradient direction */
    for(i=0;i<d->len;i++) {
        for(n=0;n<m->height;n++) {
            double x1 = 0;
            double*row = &m->data[n*m->width];
            double subderivate = 0;
            double val = 0;
            for(j=0;j<len;j++)
            for(k=0;k<len;k++) {
                val += pos[j]*pos[k]*row[len*k+j];
            }
            for(j=0;j<len;j++) {
                subderivate += (row[len*i+j]+row[len*j+i])*pos[j];
            }
            dir[i] += 2*(val - d->b[n])*subderivate;
        }
        dir[i] += d->centerforce*2*pos[i];
        dir[i] = -dir[i];
    }
#ifdef FIX0
    dir[0] = 0;
#endif
}
static void matrix_find_xx_init(void*data, double*pos)
{
    mopt_t*d = (mopt_t*)data;
    int t;
    for(t=0;t<d->len;t++) {
	pos[t] = 1-(lrand48()&2);
    }
#ifdef FIX0
    pos[0] = 1;
#endif
}
static double matrix_find_xx_fitness(void*data, double*v1, double*v2, double f)
{
    mopt_t*d = (mopt_t*)data;
    matrix_t*m = d->m;
    int len = d->len;
    int dim = m->width;
    double diff = 0;
    int n,i,j;
    double*v3;
    if(!v2) {
        v3 = v1;
    } else {
        v3 = malloc(sizeof(double)*d->len);
        int t;
        for(t=0;t<d->len;t++)
            v3[t] = v1[t]+f*v2[t];
    }
    for(n=0;n<m->height;n++) {
        double*row = &m->data[n*m->width];
        double val = 0;
        for(i=0;i<len;i++)
        for(j=0;j<len;j++) {
            val += v3[i]*v3[j]*row[len*j+i];
        }
        diff += sqr(val - d->b[n]);
    }
    for(i=0;i<len;i++) {
        diff += d->centerforce*sqr(v3[i]);
    }
    if(v3!=v1)
        free(v3);
    return diff;
}
double* matrix_find_xx_b_opt(matrix_t*m, double*b, int len)
{
    mopt_t data;
    data.m = m;
    data.b = b;
    data.len = len;
    data.centerforce = 0;
    double*best = 0;
    double bestdiff = 0;
    int t;
    for(t=0;t<3;t++) {
        double*v = gradient_approximation(len, matrix_find_xx_init, matrix_find_xx_dir, matrix_find_xx_fitness, (void*)&data, 400);
        double d = matrix_find_xx_fitness(&data, v, 0, 0); 
        if(!best || d<bestdiff) {
            bestdiff = d;
            if(best)
                free(best);
            best = v;
        } else {
            free(v);
        }
    }
    /*for(t=0;t<len;t++) {
        printf("%7.3f ", best[t]);
    }
    printf("Difference: %f\n", bestdiff);*/
    return best;
}

void showautocorrdiff(matrix_t*r, matrix_t*corr)
{
    double corrdiff = 0;
    double corrabs = 0;
    int xx,yy;
    for(yy=0;yy<corr->height;yy++) {
        int x1 = (r->width-corr->width)/2;
        int y1 = yy+(r->height-corr->height)/2;
        double*src = &r->data[y1*r->width + x1];
        double*dest = &corr->data[yy*corr->width];
        for(xx=0;xx<corr->width;xx++) {
            corrdiff += sqr(dest[xx] - src[xx]);
            corrabs += sqr(dest[xx]);
        }
    }
    printf("%2.2f%% (%2.2f SNR) of goal correlation (diff: %f/%f)\n", 
            100.0 - 100*corrdiff / corrabs,
            10 * (log(corrabs / corrdiff) / log(10)),
            corrdiff, corrabs);
}

void matrix_adjustautocorrelation2(matrix_t*img, matrix_t*corr)
{
    /* size of the autocorrelation area to adjust */
    int csizex = corr->width;
    int csizey = corr->height;
    int csize = csizex*csizey;
    int midx = csizex/2;
    int midy = csizey/2;

    /* size of our filter */
    int fsizex = csizex;
    int fsizey = csizey;
    int fsize = fsizex*fsizey;
    int fmidx = fsizex/2;
    int fmidy = fsizey/2;

    int f2sizex = csizex/2+1;
    int f2sizey = csizey/2+1;
    int f2size = f2sizex*f2sizey;

    /* size of the correlation area we need from the image meant to be adjusted */
    int rx = csizex*3;
    int ry = csizey*3;
    int rmidx = rx/2;
    int rmidy = ry/2;

    if(debug)
        printf("getting %dx%d autocorrelation from source image (%dx%d)\n", rx, ry, img->width, img->height);

    matrix_t*r = matrix_getautocorrelation(img, rx, ry);

    int x1,y1,x2,y2,xx,yy;

    if(debug) {
        image_save(image_from_matrix(r,IMAGE_MINMAX), "autocorrelation_source.png");
        image_save(image_from_matrix(corr,IMAGE_MINMAX), "autocorrelation_dest.png");
    }

    if(debug) {
        showautocorrdiff(r, corr);
    }

    matrix_t*A = matrix_new(f2size*f2size, csize);

    double*b = malloc(A->height*sizeof(double));
    memset(b, 0, A->height*sizeof(double));
    memcpy(b, corr->data, csize*sizeof(double));

    int row = 0;
    for(yy=0;yy<csizey;yy++)
    for(xx=0;xx<csizex;xx++) {
        if(row<A->height) {
            for(y2=0;y2<f2sizey;y2++)
            for(x2=0;x2<f2sizex;x2++) {
                for(y1=0;y1<f2sizey;y1++)
                for(x1=0;x1<f2sizex;x1++) {
                    int column = (y2*f2sizex+x2)*f2size + (y1*f2sizex+x1); // variable nr
                    int u1 = (xx-midx)+(x1)-(x2)+rmidx;
                    int v1 = (yy-midy)+(y1)-(y2)+rmidy;
                    int u2 = (xx-midx)-(x1)-(x2)+rmidx;
                    int v2 = (yy-midy)-(y1)-(y2)+rmidy;
                    int u3 = (xx-midx)+(x1)+(x2)+rmidx;
                    int v3 = (yy-midy)+(y1)+(y2)+rmidy;
                    int u4 = (xx-midx)-(x1)+(x2)+rmidx;
                    int v4 = (yy-midy)-(y1)+(y2)+rmidy;
                    assert(u1>=0 && u1<r->width);
                    assert(v1>=0 && v1<r->height);
                    assert(u2>=0 && u2<r->width);
                    assert(v2>=0 && v2<r->height);
                    assert(u3>=0 && u3<r->width);
                    assert(v3>=0 && v3<r->height);
                    assert(u4>=0 && u4<r->width);
                    assert(v4>=0 && v4<r->height);
                    A->data[row*A->width + column] = 
                        r->data[v1*r->width+u1]+
                        r->data[v2*r->width+u2]+
                        r->data[v3*r->width+u3]+
                        r->data[v4*r->width+u4];
                }
            }
        }
        row++;
    }
    if(debug) {
        image_save(image_from_matrix(A,IMAGE_MINMAX), "circular_correlation.png");
    }

    /*matrix_delete(r);
    matrix_t*l = (matrix_t*)malloc(sizeof(matrix_t));
    if(debug) {
        printf("solving %dx%d matrix for correlation filter matrix of size %d\n", A->width, A->height, A->width);
    }
    l->data = matrix_solve_approx(A, b);
    //l->data = matrix_solve(A, b);
    l->width = f2size;
    l->height = f2size;

    if(debug) {
        image_save(image_from_matrix(l,IMAGE_MINMAX), "basefilter.png");
        matrix_t*check = matrix_new_fromdata(matrix_multiply_with_vector(A, l->data), csizex, csizey);
        image_save(image_from_matrix(check,IMAGE_MINMAX), "autocorrelation_dest_check.png");
    }

    matrix_delete(A);

    if(debug) {
        printf("resolving %dx%d correlation filter matrix into filter product of size %dx%d\n", l->width, l->height, f2sizex, f2sizey);
    }*/
    
    matrix_t*halffilter  = (matrix_t*)malloc(sizeof(matrix_t));
    //halffilter->data = matrix_split_xx(l);
    //filter->data = matrix_split_xx_2(l);
    halffilter->data = matrix_find_xx_b_opt(A, b, f2size);
    halffilter->width = f2sizex;
    halffilter->height = f2sizey;
    if(debug) {
        image_save(image_from_matrix(halffilter,IMAGE_MINMAX), "halffilter.png");
        matrix_t* tst = matrix_new(f2size,f2size);
        for(y1=0;y1<f2size;y1++) 
        for(x1=0;x1<f2size;x1++) 
        tst->data[y1*f2size+x1]=halffilter->data[x1]*halffilter->data[y1];
        image_save(image_from_matrix(tst,IMAGE_MINMAX), "basefilter_approx.png");
    }

    matrix_t*filter = matrix_new(fsizex,fsizey);
    for(y1=0;y1<halffilter->height;y1++) 
    for(x1=0;x1<halffilter->width;x1++) {
        filter->data[(y1+fmidy)*filter->width+(x1+fmidx)] = halffilter->data[y1*halffilter->width+x1];
        filter->data[(fmidy-y1)*filter->width+(fmidx-x1)] = halffilter->data[y1*halffilter->width+x1];
    }
    image_save(image_from_matrix(filter,IMAGE_MINMAX), "filter.png");

    matrix_t*newimage = matrix_convolve(img, filter, 1, 1, EDGE_WRAP);
    matrix_delete(filter);
        
    matrix_t*newcorr = matrix_getautocorrelation(newimage, csizex,csizey);
    if(debug) {
        image_save(image_from_matrix(newcorr,IMAGE_MINMAX), "autocorrelation_dest_approx.png");
    }
    //showautocorrdiff(newcorr, corr);
    matrix_delete(newcorr);

    assert(newimage->width == img->width && newimage->height == img->height);
    memcpy(img->data, newimage->data, img->width*img->height*sizeof(img->data[0]));
}

void matrix_adjustcrosscorrelation(matrixset_t*varset, matrix_t*varcorr, 
                                   matrixset_t*staticset, matrix_t*staticcorr)
{
    matrix_t*xx = matrix_getcrosscorrelations(varset, varset, 0);
    matrix_t*xy=0,*yy=0,*yyi=0;

    if(staticset && staticcorr) {
        xy = matrix_getcrosscorrelations(varset, staticset, 0);
        yy = matrix_getcrosscorrelations(staticset, staticset, 0);
        yyi = matrix_invert(yy); // symmetric
        if(!yyi) {
            /* matrix is singular- we could now reduce the staticset step by step
               to determine a unique basis.
               For now, just leave this alone */
            return;
        }
    }

    if(debug) {
        printf("xx / current crosscorrelation:\n");
        matrix_print(xx);
        printf("xy / current crossband crosscorrelation:\n");
        matrix_print(xy);
        printf("yy / current crosscorrelation of other band:\n");
        matrix_print(yy);
        printf("yyi / current crosscorrelation of other band (inverted):\n");
        matrix_print(yyi);
    }

    /* vv = xx - xy yyi xy^T  (inter+intra case)
       vv = xx                (intra case)
     */
    matrix_t*vv = 0;
    if(staticset && staticcorr) {
        matrix_t*t1,*t2;
        t1 = matrix_multiply_t(xy, 0, yyi, 0);
        t2 = matrix_multiply_t(t1, 0, xy, 1);
        vv = matrix_add(xx, 1, t2, -1);
        matrix_delete(t1);
        matrix_delete(t2);
    } else {
        vv = xx;
    }
   
    /* rh = varcorr - staticcorr * yyi * staticcorr^T  (inter+intra case)
       rh = varcorr                                    (intra case)
     */
    matrix_t*rh = 0;
    if(staticset && staticcorr) {
        matrix_t*t1 = matrix_multiply_t(staticcorr, 0, yyi, 0);
        matrix_t*t2 = matrix_multiply_t(t1, 0, staticcorr, 1);
        rh = matrix_add(varcorr, 1, t2, -1);
        matrix_delete(t1);
        matrix_delete(t2);
    } else {
        rh = varcorr;
    }
    
    /* now solve
           m*vv*m^T = rh
       by the formula
           m = e1 * d1^1/2 * d2^(-1/2) * e2^T
           with 
            rh = e1*d1*e1^T
            vv = e2*d2*e2^T
    */
    if(debug) {
        printf("rh matrix:\n");
        matrix_print(rh);
        printf("vv matrix:\n");
        matrix_print(vv);
    }

    matrix_t*e1,*d1,*e2,*d2;
    matrix_symm_geteigenvectorfactorization2(rh, &e1,&d1);
    matrix_symm_geteigenvectorfactorization2(vv, &e2,&d2);
    if(debug) {
        printf("vv eigenvectors:\n");
        matrix_print(e1);
        matrix_print(d1);
        printf("rh eigenvectors:\n");
        matrix_print(e2);
        matrix_print(d2);
    }

    char not_positive_semidefinite = 0;
    int t;
    for(t=0;t<d1->width;t++) {
        double val = d1->data[t*d1->width+t];
        if(val < 0)
            not_positive_semidefinite |= 1;
        d1->data[t*d1->width+t] = sqrt(fabs(val));
    }
    for(t=0;t<d2->width;t++) {
        double val = d2->data[t*d2->width+t];
        if(val < 0) {
            not_positive_semidefinite |= 2;
            val = 0;
        }
        if(val < 0.00001)
            d2->data[t*d2->width+t] = 0;
        else
            d2->data[t*d2->width+t] = 1.0 / sqrt(fabs(val));
    }
    if(debug) {
        printf("diagonals:\n");
        matrix_print(d1);
        matrix_print(d2);
    }
   
    matrix_t*m = 0;
    if(not_positive_semidefinite) {
        m = matrix_new(varcorr->width, varcorr->height);
        int x,y;
        fprintf(stderr, "matrix %d not positive-semidefinite- can't adjust cross-correlation\n", not_positive_semidefinite);
        if(staticcorr && staticset) {
            fprintf(stderr, "retrying without cross-band adjustment\n");
            matrix_adjustcrosscorrelation(varset, varcorr, 0, 0);
        }
        return;
        // store unit matrix
        for(y=0;y<m->height;y++) for(x=0;x<m->width;x++) m->data[y*m->width+x] = x==y;
    } else {
        matrix_t*tmp0 = matrix_multiply(e1, d1);
#if 0
        matrix_t*tmp01 = matrix_multiply_t(tmp0, 0, e1, 1);
        tmp0 = matrix_multiply_t(tmp01, 0, e2, 0);
#endif
        matrix_t*tmp1 = matrix_multiply(tmp0, d2);
        m = matrix_multiply_t(tmp1, 0, e2, 1);
        matrix_delete(tmp0);
        matrix_delete(tmp1);
    }

#if 0
    // check that m*vv*m^T == rh
    printf("-----\n");
    matrix_print(matrix_multiply(m, matrix_multiply_t(vv, 0, m, 1)));
    printf("-----\n");
    matrix_print(rh);
    printf("-----\n");
#endif

    /* k = (staticcorr - m*xy)*yyi */
    matrix_t*k = 0;
    if(staticset && staticcorr) {
        matrix_t*t1 = matrix_multiply_t(m, 0, xy, 0);
        matrix_t*t2 = matrix_add(staticcorr, 1, t1, -1);
        k = matrix_multiply_t(t2, 0, yyi, 0);
        matrix_delete(t1);
        matrix_delete(t2);
    }

    if(debug) {
        printf("Adjustment matrix for interband crosscorrelation:\n");
        matrix_print(m);
        printf("Adjustment matrix for crossband crosscorrelation:\n");
        matrix_print(k);
    }

#define MATRIX_VALUE(m,x,y) ((m)->data[(((y)*(m)->height)/height)*(m)->width+((x)*(m)->width)/width])

    int x,y;
    int width = varset->m[0]->width;
    int height = varset->m[0]->height;
    double*v = (double*)malloc(sizeof(double)*varset->num);
    for(y=0;y<height;y++)
    for(x=0;x<width;x++) {
        int i;
        memset(v, 0, sizeof(double)*varset->num);
        for(i=0;i<varset->num;i++) {
            int j;
            for(j=0;j<varset->num;j++) {
                v[i] += MATRIX_VALUE(varset->m[j],x,y)*m->data[i*m->width+j];
            }
            if(k) {
                for(j=0;j<staticset->num;j++) {
                    assert(i<m->height && j<m->width);
                    v[i] += MATRIX_VALUE(staticset->m[j],x,y)*k->data[i*k->width+j];
                }
            }
        }
        for(i=0;i<varset->num;i++) {
            MATRIX_VALUE(varset->m[i],x,y) = v[i];
        }
    }

    if(debug) {
        printf("Dest correlation:\n");
        matrix_print(varcorr);
        if(staticset) {
            matrix_print(staticcorr);
        }
    }
    if(debug) {
        printf("New source correlation:\n");
        matrix_t*crosscorr1 = matrix_getcrosscorrelations(varset,varset, 0);
        matrix_print(crosscorr1);
        if(staticset) {
           matrix_t*crosscorr2 = matrix_getcrosscorrelations(varset,staticset, 0);
           matrix_print(crosscorr2);
        }
    }
}

void matrix_adjustautocorrelation3(matrix_t*img, matrix_t*corr)
{
    int size = img->width*img->height;
    int imidx = img->width/2;
    int imidy = img->height/2;

    /* size of the autocorrelation area to adjust */
    int csizex = corr->width;
    int csizey = corr->height;
    int csize = csizex*csizey;
    int cmidx = csizex/2;
    int cmidy = csizey/2;

    /* size of our filter */
    int fsizex = csizex;
    int fsizey = csizey;
    int fsize = fsizex*fsizey;
    int fmidx = fsizex/2;
    int fmidy = fsizey/2;

    /* size of the correlation area we need from the image meant to be adjusted */
    int rx = csizex+fsizex;
    int ry = csizey+fsizey;
    int rmidx = rx/2;
    int rmidy = ry/2;
    if(debug)
        printf("getting %dx%d autocorrelation from source image (%dx%d)\n", rx, ry, img->width, img->height);
    matrix_t*r = matrix_getautocorrelation(img, rx, ry);

    int x1,y1,x2,y2,xx,yy;

    if(debug) {
        image_save(image_from_matrix(r,IMAGE_MINMAX), "autocorrelation_source.png");
        image_save(image_from_matrix(corr,IMAGE_MINMAX), "autocorrelation_dest.png");
    }

    if(debug) {
        showautocorrdiff(r, corr);
    }

    matrix_t*A = matrix_new(fsize, csize);
    double*b = malloc(A->height*sizeof(double));
    memset(b, 0, A->height*sizeof(double));
    memcpy(b, corr->data, csize*sizeof(double));

    int row = 0;
    int x,y;
    for(yy=0;yy<csizey;yy++)
    for(xx=0;xx<csizex;xx++) {
        double*line = &A->data[(yy*csizex+xx)*A->width];
        for(y=0;y<fsizey;y++)
        for(x=0;x<fsizex;x++) {
            assert(y-fmidy+yy-cmidy+rmidy >= 0 && y-fmidy+yy-cmidy+rmidx<r->height);
            assert(x-fmidx+xx-cmidx+rmidx >= 0 && y-fmidx+xx-cmidx+rmidx<r->width);
            line[y*fsizex+x] = r->data[(y-fmidy+yy-cmidy+rmidy)*r->width + (x-fmidx+xx-cmidx+rmidx)];
        }
        row++;
    }
    matrix_delete(r);

    if(debug) {
        image_save(image_from_matrix(A,IMAGE_MINMAX), "circular_correlation.png");
    }
    
    matrix_t*l = (matrix_t*)malloc(sizeof(matrix_t));
    if(debug) {
        printf("solving %dx%d matrix for correlation filter matrix of size %d\n", A->width, A->height, A->width);
    }
    l->data = matrix_solve(A, b);
    l->width = fsizex;
    l->height = fsizey;
    
    if(debug) {
        image_save(image_from_matrix(l,IMAGE_MINMAX), "basefilter.png");
        matrix_t*check = matrix_new_fromdata(matrix_multiply_with_vector(A, l->data), csizex, csizey);
        image_save(image_from_matrix(check,IMAGE_MINMAX), "autocorrelation_dest_check.png");
    }
    matrix_delete(A);

    matrix_t*filter2 = matrix_new(img->width, img->height);
    for(y=0;y<fsizey;y++) 
    for(x=0;x<fsizex;x++) {
        int xx = (x+filter2->height-fmidx)%filter2->height;
        int yy = (y+filter2->width-fmidy)%filter2->width;
        assert(xx >= 0 && xx < filter2->width);
        assert(yy >= 0 && yy < filter2->height);
        filter2->data[yy*filter2->width+xx] = l->data[y*l->width+x];
    }
    if(debug) {
        image_save(image_from_matrix(filter2,IMAGE_MINMAX), "basefilter_extended.png");
    }

    complex_matrix_t*ff_img = matrix_fft(img);
    complex_matrix_t*ff_filter = matrix_fft(filter2);
    
    matrix_delete(filter2);

    assert(ff_img->width == ff_filter->width);
    assert(ff_img->height  == ff_filter->height);
    
    if(debug) {
        image_save(image_from_complex_matrix(ff_filter,IMAGE_MINMAX), "basefilter_fft.png");
    }
    
    for(y=0;y<ff_img->height;y++) 
    for(x=0;x<ff_img->width;x++) {
        complex_t*f = &ff_filter->data[y*ff_img->width+x];
        f->real = sqrt(fabs(f->real)) / size;
        f->imag = 0;
    }

//#if USE_CONVOLVE
    matrix_t*filter = complex_matrix_ifft_real(ff_filter);
    matrix_t*smallfilter = matrix_new(21,21);
    for(y=0;y<smallfilter->height;y++) 
    for(x=0;x<smallfilter->width;x++) {
        smallfilter->data[y*smallfilter->width+x] =  
            filter->data[((y+filter->height-smallfilter->height/2)%filter->height)*filter->width+((x+filter->width-smallfilter->width/2)%filter->width)];
    }
    static int counter = 1;
    char buf[256];
    //sprintf(buf, "filter%d.png", counter++);
    //image_save(image_from_matrix(smallfilter,IMAGE_MINMAX), buf);

//#endif

    if(debug) {
        image_save(image_from_complex_matrix(ff_filter,IMAGE_MINMAX), "filter_fft.png");
#if USE_CONVOLVE
        image_save(image_from_matrix(filter,IMAGE_MINMAX), "filter.png");
        image_save(image_from_matrix(smallfilter,IMAGE_MINMAX), "filter_extract.png");
#endif
    }
    
    if(debug) {
        image_save(image_from_complex_matrix(ff_img,IMAGE_MINMAX), "image_fft.png");
    }

    for(y=0;y<ff_img->height;y++) 
    for(x=0;x<ff_img->width;x++) {
        complex_t i = ff_img->data[y*ff_img->width+x];
        complex_t f = ff_filter->data[y*ff_img->width+x];
        f.imag = -f.imag; //transpose
        ff_img->data[y*ff_img->width+x].real = (i.real * f.real - i.imag * f.imag);
        ff_img->data[y*ff_img->width+x].imag = (i.real * f.imag + i.imag * f.real);
    }

    if(debug) {
        image_save(image_from_complex_matrix(ff_img,IMAGE_MINMAX), "image_filtered_fft.png");
    }


    matrix_t*r1=0,*r2 = 0;
#if USE_CONVOLVE
    matrix_t* newimg = matrix_convolve(img, smallfilter, 1, 1, EDGE_WRAP);
    statistics_print(statistics_new_frommatrix(newimg));
    r1 = matrix_getautocorrelation(newimg, corr->width, corr->height);
    showautocorrdiff(r1, corr);
    matrix_delete(smallfilter);
#endif
    
    matrix_t* newimg = complex_matrix_ifft_real(ff_img);

    assert(newimg->width == img->width);
    assert(newimg->height  == img->height);
    memcpy(img->data, newimg->data, img->width*img->width*sizeof(img->data[0]));
    
    if(debug) {
        image_save(image_from_matrix(r1,IMAGE_MINMAX), "autocorrelation_dest_approx1.png");
        image_save(image_from_matrix(newimg,IMAGE_MINMAX), "newimg1.png");
#if USE_CONVOLVE
        image_save(image_from_matrix(r2,IMAGE_MINMAX), "autocorrelation_dest_approx2.png");
        image_save(image_from_matrix(newimg2,IMAGE_MINMAX), "newimg2.png");
#endif
    }

    //r2 = matrix_getautocorrelation(img, corr->width, corr->height);
    //showautocorrdiff(r2, corr);
    //matrix_delete(r2);

#if USE_CONVOLVE
    matrix_delete(filter);
#endif
    matrix_delete(newimg);
}

