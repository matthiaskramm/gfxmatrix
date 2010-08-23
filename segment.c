/* segment.c 

   Texton based segmenting.

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
#include <assert.h>
#include <memory.h>
#include "image.h"
#include "filter.h"
#include "filter_misc.h"
#include "segment.h"
#include "stat.h"
#include "kmeans.h"

#define ADD_MATRIX(m) {if(pass==1) {filter_t*f = filter_new(filtertype_convolve, FILTER_DIRECTIONAL, m, 1, 1);assert(num<base->num_children);base->children[num] = filtertree_new(0, f);} num++;}
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

filtertree_t* makeFilterTree(int flags)
{
    int pass;
    int num=0;
    filtertree_t*base = 0;

    for(pass=0;pass<2;pass++) {
        if(pass==1)
            base = filtertree_new(num, 0);
        num = 0;

        if(flags&FILTER_GAUSS) {
            ADD_MATRIX(gauss_highpass(15, 15, 0.0, 1.0));
            ADD_MATRIX(gauss_highpass(15, 15, 0.0, 0.5));
            ADD_MATRIX(gauss_highpass(15, 15, 0.0, 0.3));
        }

        if(flags&FILTER_GAUSS_HIGHPASS) {
            ADD_MATRIX(gauss_highpass(15, 15, 5.0, 0.5));
            //ADD_MATRIX(gauss_highpass(15, 15, 5.0, 1.0));
            //ADD_MATRIX(gauss_highpass(15, 15, 5.0, 2.0));
            //ADD_MATRIX(gauss_highpass(15, 15, 3.0, 0.5));
            //ADD_MATRIX(gauss_highpass(15, 15, 3.0, 1.0));
            //ADD_MATRIX(gauss_highpass(15, 15, 3.0, 2.0));
        }

        if(flags&FILTER_GAUSS_DIFF) {
            int t;
            for(t=0;t<4;t++) {
                //ADD_MATRIX(gauss_diff(15, 15, 1/4.0, t*45, 2.0));
                ADD_MATRIX(gauss_diff(15, 15, 1/4.0, t*45, 1.0));
                //ADD_MATRIX(gauss_diff(15, 15, 1/2.0, t*45, 2.0));
                ADD_MATRIX(gauss_diff(15, 15, 1/2.0, t*45, 1.0));
            }
        }

        if(flags&FILTER_GABOR) {
            int s;
            ADD_MATRIX(gabor_filter(15,15, 0.25,0.25, 0, 1/8.0, 0, 0));
            ADD_MATRIX(gabor_filter(15,15, 0.25,0.25, 0, 1/8.0, 90, 0));

            ADD_MATRIX(gabor_filter(15,15, 0.25,0.25, 0, 1/8.0, 0, 1));
            ADD_MATRIX(gabor_filter(15,15, 0.25,0.25, 0, 1/8.0, 45, 1));
            ADD_MATRIX(gabor_filter(15,15, 0.25,0.25, 0, 1/8.0, 90, 1));
            ADD_MATRIX(gabor_filter(15,15, 0.25,0.25, 0, 1/8.0, 135, 1));
                    
            ADD_MATRIX(gabor_filter(15,15, 0.25,0.25, 0, 1/16.0, 0, 0));
            ADD_MATRIX(gabor_filter(15,15, 0.25,0.25, 0, 1/16.0, 45, 0));
            ADD_MATRIX(gabor_filter(15,15, 0.25,0.25, 0, 1/16.0, 90, 0));
            ADD_MATRIX(gabor_filter(15,15, 0.25,0.25, 0, 1/16.0, 135, 0));

            ADD_MATRIX(gabor_filter(15,15, 0.25,0.75, 0, 0.0, 45, 1));
            ADD_MATRIX(gabor_filter(15,15, 0.25,0.75, 0, 0.0, -45, 1));
            //ADD_MATRIX(gabor_filter(15,15, 0.25,0.75, 0, 0.0, 45, 0));
            //ADD_MATRIX(gabor_filter(15,15, 0.25,0.75, 0, 1, 90, 0));
            //ADD_MATRIX(gabor_filter(15,15, 0.25,0.75, 0, 1, 135, 0));

            ADD_MATRIX(gabor_filter(15,15, 0.25,0.75, 0, 1/16.0, 0, 0));
            ADD_MATRIX(gabor_filter(15,15, 0.25,0.75, 0, 1/16.0, 45, 0));
            ADD_MATRIX(gabor_filter(15,15, 0.25,0.75, 0, 1/16.0, 90, 0));
            ADD_MATRIX(gabor_filter(15,15, 0.25,0.75, 0, 1/16.0, 135, 0));
            
            ADD_MATRIX(gabor_filter(15,15, 0.25,0.25, 0, 1/4.0, 0, 1));
            ADD_MATRIX(gabor_filter(15,15, 0.25,0.25, 0, 1/4.0, 45, 1));
            ADD_MATRIX(gabor_filter(15,15, 0.25,0.25, 0, 1/4.0, 90, 1));
            ADD_MATRIX(gabor_filter(15,15, 0.25,0.25, 0, 1/4.0, 135, 1));
            
            for(s=1;s<4;s+=2) {
                int t;
                float scale = 1/(1.0+s);

                for(t=0;t<4;t++) {
                    //ADD_MATRIX(gabor_filter(15,15, scale,scale*3, 0, 0, t*45, 1));
                    
                    //ADD_MATRIX(gabor_filter(15,15, scale,scale, 0, 1/8.0, t*45, 1));
                    
                    //ADD_MATRIX(gabor_filter(15,15, scale,scale, 0, 1/4.0, t*45, 0));
                    //ADD_MATRIX(gabor_filter(15,15, scale,scale, 0, 1/4.0, t*45, 1));
                    //ADD_MATRIX(gabor_filter(15,15, scale,scale*3, 0, 1/8.0, t*45, 0));
                    //ADD_MATRIX(gabor_filter(15,15, scale,scale*3, 0, 1/8.0, t*45, 1));
                    //ADD_MATRIX(gabor_filter(15,15, scale,scale*3, 0, 1/4.0, t*45, 0));
                    //ADD_MATRIX(gabor_filter(15,15, scale,scale*3, 0, 1/4.0, t*45, 1));

                    //ADD_MATRIX(gabor_filter(15,15, scale,scale, 0, 1/2.0, t*45, 1));
                    //ADD_MATRIX(gabor_filter(15,15, scale,scale*3, 0, 1/2.0, t*45, 1));
                }
            }
        }
    }
    check_filtertree_uniqueness(base);
    return base;
}

filtertree_t* makeEdgeFilter()
{
    int pass;
    int num;
    filtertree_t*base = 0;

    for(pass=0;pass<2;pass++) {
        if(pass==1)
            base = filtertree_new(num, 0);
        num = 0;
        
        ADD_MATRIX(gauss_diff2(15, 15, 1.0, 0));
        ADD_MATRIX(gauss_diff2(15, 15, 1.0, 90));

        //ADD_MATRIX(gabor_filter(15,15, 0.25,0.25, 0, 1/8.0, 0, 0));
        //ADD_MATRIX(gabor_filter(15,15, 0.25,0.25, 0, 1/8.0, 90, 0));
    }
    check_filtertree_uniqueness(base);
    return base;
}

filtertree_t* makeSimpleFilterTree(char gauss, char lgauss, char gauss_diff)
{
    int pass;
    int num;
    filtertree_t*base = 0;

    for(pass=0;pass<2;pass++) {
        if(pass==1)
            base = filtertree_new(num, 0);
        num = 0;
        if(gauss) {
            ADD_MATRIX(gauss_filter(15, 15, 1.0));
            ADD_MATRIX(gauss_filter(23, 23, 2.0));
            ADD_MATRIX(gauss_filter(31, 31, 4.0));
        }
       
        if(lgauss) {
            ADD_MATRIX(lgauss_filter(15, 15, 1.0));
            ADD_MATRIX(lgauss_filter(23, 23, 2.0));
            ADD_MATRIX(lgauss_filter(31, 31, 4.0));
            ADD_MATRIX(lgauss_filter(39, 39, 8.0));
        }
        
        if(gauss_diff) {
            ADD_MATRIX(gauss_diff2(15, 15, 2.0, 0));
            ADD_MATRIX(gauss_diff2(15, 15, 2.0, 90));
            ADD_MATRIX(gauss_diff2(31, 31, 4.0, 0));
            ADD_MATRIX(gauss_diff2(31, 31, 4.0, 90));
        }
    }
    check_filtertree_uniqueness(base);
    return base;
}

bytearray_t* kmeans(matrixset_t*set, int num_centers, matrix_t**centersptr, int flags)
{
    int size = set->m[0]->width*set->m[0]->height;
    int t;
    int depth = set->num;
    if(num_centers>255) {
        fprintf(stderr, "kmeans(): Can only cluster up to 255 centers\n");
        return 0;
    }
    matrix_t*cov=0,*icov=0;
    matrixset_t*iset = 0;

    cov = matrix_getcrosscorrelations(set,set,0);
    if(flags&CLUSTER_MAHALANOBIS) {
        matrix_t*mean = matrixset_getmean(set);
        int r;
        for(r=0;r<depth*depth;r++) {
            if(fabs(cov->data[r])<0.0001)
                cov->data[r] = 0;
        }
        icov = matrix_invert(cov);
        if(!icov) {
            fprintf(stderr, "Bad covariance matrix (singular)- check your filterbank\n");
            exit(1);
        }
        iset = matrixset_clone(set);

        for(t=0;t<size;t++) {
            int i,j;
            for(i=0;i<depth;i++) {
                double sum = 0;
                double*row = &icov->data[i*icov->width];
                for(j=0;j<depth;j++) {
                    sum += row[j]*set->m[j]->data[t];
                }
                iset->m[i]->data[t] = sum;
            }
        }
    }
    unsigned char* map = (unsigned char*)malloc(size);
    memset(map, 255, size);

    matrix_t*centers = matrix_new(num_centers, depth); // each column a center
    int minchanges = size;
    int minchanges_count = 0;
    int iteration = 0;

#define RANDOM_RANGE
#if defined(RANDOM_NODES)
    /* randomly assign all centers to different points */
    for(t=0;t<num_centers;t++) {
        int pos = lrand48()%size;
        int r;
        for(r=0;r<depth;r++) {
            centers->data[r*centers->width+t] = /*set->m[r]->data[pos] +*/ (lrand48()%4096)/4096.0;
        }
    }
#elif defined(RANDOM_RANGE)
    int r;
    for(r=0;r<depth;r++) {
        double min = set->m[r]->data[0];
        double max = set->m[r]->data[0];
        int t;
        for(t=1;t<size;t++) {
            if(set->m[r]->data[t] > max) max = set->m[r]->data[t];
            if(set->m[r]->data[t] < min) min = set->m[r]->data[t];
        }
        for(t=0;t<num_centers;t++) {
            centers->data[r*centers->width+t] = ((lrand48()&0x7fffffff)*(max-min)/((double)0x7fffffff))+min;
        }
    }
#endif

    while(1) {
        matrix_t*icenters = 0;
        if(icov) 
            icenters = matrix_multiply(icov, centers);
    
        int t;
        int changes = 0;
        for(t=0;t<size;t++) {
            int s;
            int best;
            double bestv;
            for(s=0;s<num_centers;s++) {
                int r;
                float dist = 0.0;
                if(icenters && iset) {
                    for(r=0;r<depth;r++) {
                        float c1 = centers->data[r*centers->width+s] - set->m[r]->data[t];
                        float c2 = icenters->data[r*icenters->width+s] - iset->m[r]->data[t];
                        dist += c1*c2;
                    }
                } else {
                    for(r=0;r<depth;r++) {
                        float c1 = centers->data[r*centers->width+s] - set->m[r]->data[t];
                        dist += c1*c1 / cov->data[r*cov->width+r];
                    }
                }
                //if(dist < bestv || (dist==bestv && map[t]==s)) 
                if(!s || dist < bestv) {
                    best = s;
                    bestv = dist;
                }
            }
            if(map[t]!=best) {
                changes++;
                map[t] = best;
            }
        }
        memset(centers->data, 0, sizeof(centers->data[0])*centers->width*centers->height);
        for(t=0;t<num_centers;t++) {
            int count = 0;
            int s;
            for(s=0;s<size;s++) {
                if(map[s] == t) {
                    int r;
                    for(r=0;r<depth;r++) {
                        centers->data[centers->width*r+t] += set->m[r]->data[t];
                    }
                    count ++;
                }
            }
            if(count) {
                int r;
                //printf("cluster %d: %d data points\n", t, count);
                for(r=0;r<depth;r++) {
                    centers->data[centers->width*r+t] /= count;
                    //printf("%f\n", centers->data[centers->width*r+t]);
                }
            }
            //printf("%d ", count);
        }
        //printf(" (%d changes, min=%d, count=%d)\n", changes, minchanges, minchanges_count);
        if(icenters)
            matrix_delete(icenters);

        if(changes < minchanges) {
            minchanges = changes;
            minchanges_count = 0;
        } else {
            minchanges_count++;
            if(minchanges_count >= num_centers*2)
                break;
        }
        if(!changes)
            break;
        iteration++;
    }

    bytearray_t*b = malloc(sizeof(bytearray_t));
    b->width = set->m[0]->width;
    b->height = set->m[0]->height;
    b->data = map;

    if(centersptr)
        *centersptr = centers;
    else
        matrix_delete(centers);

    if(cov)
        matrix_delete(cov);
    if(icov)
        matrix_delete(icov);
    
    //for(t=0;t<iset->num;t++) set->m[t] = iset->m[t];
    if(iset)
        matrixset_delete(iset);

    return b;
}

matrix_t* mkTexton(filtertree_t*tree, double*values)
{
    int t;
    int width=0,height=0;
    for(t=0;t<tree->num_children;t++) {
        matrix_t*m = tree->children[t]->filter->matrix;
        if(m->width>width) width=m->width;
        if(m->height>height) height=m->height;
    }
    /*
    matrix_t*texton = matrix_new(width,height);
    for(t=0;t<tree->num_children;t++) {
        matrix_t*m = tree->children[t]->filter->matrix;
        if(tree->children[t]->filter->flags&FILTER_COMPLEX) 
            m = complex_matrix_realpart((complex_matrix_t*)m);

        int midx = (width - m->width)/2;
        int midy = (height - m->height)/2;
        int x,y;
        for(y=0;y<m->height;y++)
        for(x=0;x<m->width;x++) {
            texton->data[(y+midy)*texton->width+(x+midx)] += m->data[y*m->width + x]*values[t];
        }
    }
    return texton;*/

    //width = 11;
    //height = 11;

    matrix_t*m2 = matrix_new(width*height, tree->num_children);
    for(t=0;t<tree->num_children;t++) {
        matrix_t*m = tree->children[t]->filter->matrix;
        if(tree->children[t]->filter->flags&FILTER_COMPLEX) 
            m = complex_matrix_realpart((complex_matrix_t*)m);
        int midx = (width - m->width)/2;
        int midy = (height - m->height)/2;
        int x,y;
        for(y=0;y<m->height;y++)
        for(x=0;x<m->width;x++) {
            if(x+midx >=0 && x+midx < width &&
               y+midy >=0 && y+midy < height) {
                m2->data[t*m2->width + (y+midy)*width+(x+midx)] = m->data[y*m->width + x];
            }
        }
    }
    //image_save_and_free(image_from_matrix(m2, IMAGE_MINMAX), "equation.png");
    //matrix_save(m2, "test.matrix");

    /*double*data = malloc(sizeof(double)*width*height);
    double*values2 = malloc(sizeof(double)*tree->num_children);
    for(t=0;t<1;t++) {
        int s;
        for(s=0;s<tree->num_children;s++) {
            values2[s] = values[s];// + ((lrand48()%255)-128) / 64.0;
        }
        double*data2 = matrix_solve_underdetermined(m2, values2);
        for(s=0;s<width*height;s++) {
            data[s] += data2[s];
        }
        free(data2);
    }*/
    double*data = matrix_solve_underdetermined(m2, values);

    matrix_t*m = malloc(sizeof(matrix_t));
    m->data = data;
    m->width = width;
    m->height = height;

    /*double*r = matrix_multiply_with_vector(m2, data);
    for(t=0;t<width*height;t++) {
        printf("%f ", data[t]);
    }
    printf("\n");
    for(t=0;t<tree->num_children;t++) {
        printf("%f %f %f\n", r[t], values[t], r[t] - values[t]);
    }*/

    matrix_delete(m2);
    return m;
}

int update_kmeans_map(bytearray_t*map, bytearrayset_t*set, bytearray_t*centers, int maxdist)
{
    int width = set->m[0]->width;
    int height = set->m[0]->height;
    int size = width*height;
    int depth = set->num;
    int num_centers = centers->height;
    int t;
    int changes = 0;
    unsigned char vector[depth];

    for(t=0;t<size;t++) {
        int s;
        int best = 0;
        unsigned long bestv = 256*256*depth;
        int r;
        for(r=0;r<depth;r++) {
            vector[r] = set->m[r]->data[t];
        }
        unsigned char*c = centers->data;
        for(s=0;s<num_centers;s++) {
            int r;
            unsigned long dist = 0;
            for(r=0;r<depth;r++) {
                int c1 = c[r] - vector[r];
#define NORM_MAX
#if defined(NORM_MAX)
                if(abs(c1) > dist) 
                    dist=abs(c1);
#elif defined(NORM_SQR)
                dist+=c1*c1;
#else //defined(NORM_ABS)
                dist += abs(c1);
#endif
            }
            c+=centers->width;
            if(dist < bestv) {
                best = s;
                bestv = dist;
            }
        }
        if(map->data[t]!=best && bestv<=maxdist) {
            map->data[t] = best;
            changes++;
        }
    }
    return changes;
}

static inline void*malloc16(int len)
{
    void*ptr=0;
    posix_memalign(&ptr, 16, len);
    return ptr;
}

int update_kmeans_map_sse(bytearray_t*map, bytearrayset_t*set, bytearray_t*centers, int dist)
{
    int width = set->m[0]->width;
    int height = set->m[0]->height;
    int size = width*height;
    int depth = set->num;
    int num_centers = centers->height;
    if((int)map->data&15 || (int)set->m[0]->data&15 || (int)centers->data&15) {
        fprintf(stderr, "Warning (kmeans): bad data alignment. Switching back to (slow) c-function\n");
        return update_kmeans_map(map,set,centers,0xffffffff);
    }
    if((width*height) & 15) {
        fprintf(stderr, "Warning (kmeans): bad data size. Switching back to (slow) c-function\n");
        return update_kmeans_map(map,set,centers,0xffffffff);
    }
    unsigned char**clusters = malloc(sizeof(unsigned char*)*num_centers);
    int t;
    unsigned char* val = centers->data;
    for(t=0;t<num_centers;t++) {
        clusters[t] = malloc16(16*depth);
        int s;
        int r;
        unsigned char*b = clusters[t];
        for(r=0;r<depth;r++) {
            for(s=0;s<16;s++) {
                *b = *val;
                b++;
            }
            val++;
        }
    }
    unsigned char**vectors = malloc(sizeof(unsigned char*)*depth);
    for(t=0;t<depth;t++) {
        vectors[t] = set->m[t]->data;
    }
    if(dist>127) dist = 127;
    if(dist<0) dist = 0;
    
    int changes = distance8(vectors, clusters, depth, num_centers, map->data, size/16, dist);

    for(t=0;t<num_centers;t++) {
        free(clusters[t]);clusters[t] = 0;
    }
    free(vectors);
    free(clusters);
    return changes;
}

void calculate_kmeans_centers_2(bytearray_t*map, bytearrayset_t*set, bytearray_t*centers)
{
    int width = set->m[0]->width;
    int height = set->m[0]->height;
    int size = width*height;
    int depth = set->num;
    int num_centers = centers->height;
    int r,t;
    
    int sizes[num_centers];
    int count[num_centers];

    for(r=0;r<depth;r++) {
        memset(sizes, 0, sizeof(int)*num_centers);
        memset(count, 0, sizeof(int)*num_centers);
        unsigned char*data = set->m[r]->data;
        int s;
        for(s=0;s<size;s++) {
            int p = map->data[s];
            if(p==255)
                continue;
            sizes[p] += data[s];
            count[p] ++;
        }
        //printf("count: ");
        for(t=0;t<num_centers;t++) {
            int v = sizes[t];
            //printf("%d/%d ", sizes[t], count[t]);
            if(count[t]) {
                v /= count[t];
            } else {
                //v = lrand48()&255;
                v = lrand48()&127; //?
            }
            centers->data[centers->width*t+r] = v;
        }
        //printf("\n");
    }
}

bytearray_t* kmeans_bytearray(bytearrayset_t*origset, int num_centers, bytearray_t**centersptr, int flags, int dist)
{
    int owidth = origset->m[0]->width;
    int oheight = origset->m[0]->height;
    int osize = owidth*oheight;
    int t;
    int depth = origset->num;

    int width = (int)(sqrt(owidth)*2);
    int height = (int)(sqrt(oheight)*2);
    width=width&~15;
    height=height&~15;
    bytearrayset_t*set = bytearrayset_select_n(origset, width, height);
    int size = width*height;
    char sse = 1;

    /*int s;
    int zeroes = 0;
    for(s=0;s<osize;s++) {
        int r;
        for(r=0;r<depth;r++) {
            if(origset->m[r]->data[s]!=0)
                break;
        }
        if(r==depth)
            zeroes++;
    }
    printf("%d\n", zeroes);*/

    if(num_centers>255) {
        fprintf(stderr, "kmeans(): Can only cluster up to 255 centers\n");
        return 0;
    }
    bytearray_t*map = bytearray_new(width, height);

    int x,y;
    unsigned char highbit=0;
    for(t=0;t<set->num;t++) {
        int s;
        for(s=0;s<size;s++) {
            highbit |= set->m[t]->data[s];
        }
    }
    if(highbit&0x80) {
        fprintf(stderr, "********************************************************************************\n");
        fprintf(stderr, "Warning: high bit (0x80) set for input data, reverting to non-sse implementation\n");
        fprintf(stderr, "********************************************************************************\n");
        sse = 0;
    }

    bytearray_t*centers = bytearray_new(depth, num_centers); // each row (!) a center
    int minchanges = width;
    int minchanges_count = 0;
    int iteration = 0;
    unsigned char*vector = malloc(depth);
    unsigned long*wvector = malloc(depth*sizeof(unsigned long));

    int r;
    for(r=0;r<size;r++) {
        map->data[r] = r<num_centers?r:255;//lrand48()%num_centers;
    }

    while(1) {
        char allfull=1;
        if(flags&VERBOSE) {
            int t;
            printf("[");
            for(t=0;t<=num_centers;t++) {
                int count = 0;
                int s;
                if(t==num_centers) {
                    t = 255;
                    printf("|");
                }
                for(s=0;s<size;s++) {
                    if(map->data[s] == t) {
                        count++;
                    }
                }
                if(!count)
                    allfull=0;
                printf("%d ", count);
            }
            printf("]\n");
        }
        
        calculate_kmeans_centers_2(map, set, centers);

        int changes = 0;
        if(!sse) {
            changes = update_kmeans_map(map, set, centers,dist);
        } else {
            changes = update_kmeans_map_sse(map, set, centers, dist);
        }

        if(flags&VERBOSE) {
            printf(" (changes: %d)\n", changes);
        }
       
        if(changes < minchanges) {
            minchanges = changes;
            minchanges_count = 0;
        } else {
            minchanges_count++;
            if(minchanges_count >= num_centers*2)
                break;
        }
        if(!changes && allfull)
            break;
        iteration++;
    }
    calculate_kmeans_centers_2(map, set, centers);
       
    if(centersptr)
        *centersptr = centers;
    else
        bytearray_delete(centers);
     
    bytearray_delete(map);
    
    map = bytearray_new(owidth, oheight);
    memset(map->data, 255, owidth*oheight);
    if(!sse) {
        update_kmeans_map(map, origset, centers, dist);
    } else {
        update_kmeans_map_sse(map, origset, centers, dist);
    }

    if(flags&VERBOSE) {
        printf("Final run: \n");
        printf("[");
        for(t=0;t<=num_centers;t++) {
            int count = 0;
            int s;
            if(t==num_centers) {
                t = 255;
                printf("|");
            }
            for(s=0;s<map->width*map->height;s++) {
                if(map->data[s] == t) {
                    count++;
                }
            }
            printf("%d ", count);
        }
        printf("]\n");
    }


    if(vector)
        free(vector);

    return map;
}

int kmeans_add_rects_fast(bytearray_t*map, bytearrayset_t*set, bytearray_t*centers, int blockx, int blocky, int dist)
{
    memset(map->data, 255, map->width*map->height);
    update_kmeans_map_sse(map, set, centers, dist);
    bytearray_save_map(map, "map_before.png");
    int missing = 0;
    unsigned char*line = malloc(map->width);
    unsigned int*distx = malloc(map->width*sizeof(int));
    unsigned int*disty = malloc(map->width*sizeof(int));
    memcpy(line, map->data, map->width);
    memset(distx, 0, map->width*sizeof(int));
    memset(disty, 0, map->width*sizeof(int));
    
    int x,y;
    for(y=0;y<map->height;y++) {
        unsigned char last = 255;
        int lastdistx = 0;
        int lastdisty = 0;
        unsigned char*data = &map->data[y*map->width];
        for(x=0;x<map->width;x++) {
            if(data[x]!=255) {
                last = data[x];
                lastdistx = 0;
                lastdisty = 0;
            }
            if(distx[x]+disty[x]+1 < lastdistx+lastdisty && disty[x]<blocky) {
                last = line[x];
                lastdistx = distx[x];
                lastdisty = disty[x]+1;
            }
            line[x] = last;
            distx[x] = lastdistx;
            disty[x] = lastdisty;
            lastdistx++;
        }
        for(x=0;x<map->width;x++) {
            if(data[x]==255) {
                if(distx[x]<blockx && disty[x]<blocky) {
                    data[x] = line[x];
                } else {
                    missing++;
                }
            }
        }
    }
    free(line);
    free(distx);
    free(disty);
    return missing;
}

int kmeans_add_rects(bytearray_t*map, bytearrayset_t*set, bytearray_t*centers, int blockx, int blocky, int dist)
{
    memset(map->data, 255, map->width*map->height);
    update_kmeans_map_sse(map, set, centers, dist);
    bytearray_save_map(map, "map_before.png");

    bytearray_t*map2 = bytearray_clone(map);
    int size = map->width*map->height;
    
    int x,y,t;
    for(y=0;y<map->height;y++) {
        unsigned char*data = &map->data[y*map->width];
        unsigned char*data2 = &map2->data[y*map2->width];
        int lefty = map->height - y;
        if(lefty > blocky) {
            lefty = blocky;
        }
        for(x=0;x<map->width;x++) {
            if(data[x]!=255) {
                unsigned char b = data[x];
                int leftx = map->width - x;
                if(leftx > blockx) {
                    leftx = blockx;
                }
                int yy;
                unsigned char*d = &data2[x];
                for(yy=0;yy<lefty;yy++) {
                    memset(d, b, leftx);
                    d+=map->width;
                }
            }
        }
    }
    memcpy(map->data, map2->data, map->width*map->height);
    int missing = 0;
    for(t=0;t<size;t++) {
        missing += (map->data[t]==255);
    }
    bytearray_delete(map2);

    return missing;
}

int kmeans_add_rects_ho(bytearray_t*map, bytearrayset_t*set, bytearray_t*centers, int blockx, int blocky, int dist)
{
    memset(map->data, 255, map->width*map->height);
    update_kmeans_map_sse(map, set, centers, dist);

    bytearray_t*map2 = bytearray_clone(map);
    int size = map->width*map->height;
    
    int x,y,t;
    for(y=0;y<map->height;y++) {
        unsigned char*data = &map->data[y*map->width];
        unsigned char*data2 = &map2->data[y*map2->width];
        int lefty = y;
        if(lefty > blocky) {
            lefty = blocky;
        }
        for(x=0;x<map->width;x++) {
            int leftx = x;
            if(leftx > blockx) {
                leftx = blockx;
            }
            unsigned char*d = &data[x-(leftx-1)-(lefty-1)*map->width];
            int best = 0;
            int bestt = 255;
            unsigned char count[256];
            memset(count, 0, sizeof(count));
            int yy;
            for(yy=0;yy<lefty;yy++) {
                int xx;
                for(xx=0;xx<leftx;xx++) {
                    count[d[xx]]++;
                    if(count[d[xx]]>best && d[xx]!=255) {
                        bestt = d[xx];
                        best = count[bestt];
                    }
                }
                d+=map->width;
            }
            data2[x] = bestt;
        }
    }

    memcpy(map->data, map2->data, map->width*map->height);
    int missing = 0;
    for(t=0;t<size;t++) {
        missing += (map->data[t]==255);
    }
    bytearray_delete(map2);

    return missing;
}

void image_mark(image_t*img, bytearray_t*map)
{
    int t;
    int size = img->width*img->height;
    for(t=img->width+1;t<size;t++) {
        if(map->data[t] != map->data[t-1] || 
           map->data[t] != map->data[t-img->width]) {
            img->data[t].r = 255;
            img->data[t].g = 255;
            img->data[t].b = 0;
            img->data[t].a = 255;
        }
    }
}

window_t window_new(int x1,int y1,int x2,int y2)
{
    window_t c;
    c.x1 = x1;
    c.y1 = y1;
    c.x2 = x2;
    c.y2 = y2;
    c.width = x2-x1;
    c.height = y2-y1;
    return c;
}

window_t find_max_window(bytearray_t*map, int marker)
{
    int*count = malloc(sizeof(int)*map->width*map->height);
    memset(count, 0, sizeof(int)*map->width*map->height);
    
    int x,y;
    for(x=0;x<map->width;x++) {
        int length = 0;
        int yy = map->width*(map->height-1);
        for(y=0;y<map->height;y++) {
            if(map->data[yy+x] == marker)
                length++;
            else
                length=0;
            count[yy+x] = length;
            yy -= map->width;
        }
    }

    window_t bestw = window_new(0,0,0,0);
    int bestsize = 0;

    for(y=0;y<map->height;y++) {
        int*c = &count[y*map->width];
        for(x=0;x<map->width;x++) {
            int t;
            int left = map->width - x;
            if(left<bestsize)
                break;
            int min = c[x];
            for(t=0;t<left;t++) {
                if(c[x+t]<min) 
                    min = c[x+t];
                if(min<=t)
                    break;
            }
            if(t>bestsize) {
                bestw = window_new(x,y,x+t,y+t);
                bestsize = t;
            }
        }
    }
    free(count);
    //if(bestsize)
    //    printf("Found window of size %dx%d for marker %d\n", bestsize, bestsize, marker);

    for(y=bestw.y1;y<bestw.y2;y++) 
    for(x=bestw.x1;x<bestw.x2;x++) {
        if(map->data[y*map->width+x] != marker) {
            /*printf("%d,%d\n", x-bestw.x1,y-bestw.y1);
            for(y=bestw.y1;y<bestw.y2;y++) {
                for(x=bestw.x1;x<bestw.x2;x++) {
                    if(map->data[y*map->width+x]==marker)
                        printf(" %02x/%02dl ", map->data[y*map->width+x], count[y*map->width+x]);
                    else
                        printf("[%02x/%02d] ", map->data[y*map->width+x], count[y*map->width+x]);
                }
                printf("\n");
            }*/
            fprintf(stderr, "Internal error in find_max_window()\n");
            exit(0);
        }
    }
    return bestw;
}

typedef struct _filterbanklist
{
    filterbank_t*filterbank;
    struct _filterbanklist*next;
} filterbanklist_t;

static filterbanklist_t* filterbanklist = 0;

bytearray_t* segment(image_t*img, int num_centers, int bx, int by, int dist, bytearray_t**centersdest, int flags)
{
    matrix_t*m = image_extractchannel(img, IMAGE_GRAY);
    bytearray_t*u = image_getchannel(img, IMAGE_YUV_U);
    bytearray_t*v = image_getchannel(img, IMAGE_YUV_V);
    complex_matrix_t*img_dft = matrix_fft(m);

    filterbanklist_t*list = filterbanklist;
    filterbank_t*filterbank = 0;
    while(list) {
        if(list->filterbank->width == img->width && list->filterbank->height == img->height) {
            filterbank = list->filterbank;
        }
        list = list->next;
    }
    if(!filterbank)  {
        printf("Creating %dx%d filterbank...\n", img->width, img->height);
        complex_matrixset_t* dct = filterset_new_dct(bx,by,2,2,1);
        filterbank = filterbank_new(dct, img->width, img->height, EXPAND_RIGHTDOWN);
        filterbanklist_t*nl = malloc(sizeof(filterbanklist_t));
        memset(nl, 0, sizeof(filterbanklist_t));
        nl->filterbank = filterbank;
        nl->next = filterbanklist;
        filterbanklist = nl;
    }
    
    bytearrayset_t*set = filterbank_apply_tobytearray(filterbank, img_dft);

    if(flags&SEGMENT_SPATIAL) {
        bytearrayset_extend_alloc(set, 2, img->width, img->height);
        int x,y;
        for(y=0;y<img->height;y++)
        for(x=0;x<img->width;x++) {
            set->m[set->num-2]->data[y*img->width+x] = (x^((x>>7&1)*127))&127;
            set->m[set->num-1]->data[y*img->width+x] = (y^((y>>7&1)*127))&127;
        }
    }
    if(flags&SEGMENT_UV) {
        bytearrayset_extend(set, 2);
        bytearray_average_rect(u,bx,by);
        bytearray_normalize(u);
        set->m[set->num-2] = u;
        bytearray_average_rect(v,bx,by);
        bytearray_normalize(v);
        set->m[set->num-1] = v;
    }
    if(flags&SEGMENT_MARGINAL) {
        bytearrayset_append(set, matrix_blockstats(m, bx, by));
    }

    bytearray_t*centers = 0;
    bytearray_t*map = kmeans_bytearray(set, num_centers, &centers, /*VERBOSE*/0, dist);
    //bytearray_save_map(map, "map_before.png");
    
    kmeans_add_rects_ho(map, set, centers, bx, by, dist);
    
    bytearrayset_delete(set);

    if(centersdest) {
        *centersdest = centers;
    } else {
        bytearray_delete(centers);
    }
    return map;
}

