/* filter.c 

   Filter pyramid functions.

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

#include <memory.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "filter.h"
#include "image.h"
#include "matrix.h"

static matrix_t*filter_apply_subsample(filter_t*filter, matrix_t*src_image)
{
    matrix_t *dst_image = matrix_new(src_image->width/filter->xres,src_image->height/filter->yres);
    int x_src,y_src,x_dst,y_dst;
    double divide = 1 / (filter->xres*filter->yres);
    for(y_src=0,y_dst=0; y_src<src_image->height; y_src+=filter->yres,y_dst++) {
        double*src = &src_image->data[y_src*src_image->width];
        double*dst = &dst_image->data[y_dst*dst_image->width];
        for(x_src=0,x_dst=0; x_src<src_image->width; x_src+=filter->xres,x_dst++) {
            int x,y;
            double val = 0;
            for(y=0;y<filter->yres;y++) {
                double*src2 = &src[x_src+y*src_image->height];
                for(x=0;x<filter->xres;x++) {
                    val += src2[x];
                }
            }
            dst[x_src] = val * divide;
        }
    }
    return dst_image;
}



static matrix_t*filter_apply_convolution(filter_t*filter, matrix_t*src_image)
{
    if(!(filter->flags&FILTER_COMPLEX)) {
        matrix_t*r = matrix_convolve(src_image, filter->matrix, filter->xres, filter->yres, EDGE_REFLECT);
        return r;
    } else {
        matrix_t*real = complex_matrix_realpart((complex_matrix_t*)filter->matrix);
        matrix_t*imag = complex_matrix_imagpart((complex_matrix_t*)filter->matrix);
        matrix_t*r = matrix_convolve(src_image, real, filter->xres, filter->yres, EDGE_REFLECT);
        matrix_t*i = matrix_convolve(src_image, imag, filter->xres, filter->yres, EDGE_REFLECT);
        int t;
        int size = r->width*r->height;
        matrix_t*result = matrix_new(r->width, r->height);
        for(t=0;t<size;t++) {
            //result->data[t] = sqrt(r->data[t]*r->data[t] + i->data[t]*i->data[t]);
            result->data[t] = (r->data[t]*r->data[t] + i->data[t]*i->data[t]);
        }
        matrix_delete(real);
        matrix_delete(imag);
        matrix_delete(r);
        matrix_delete(i);
        return result;
        //return result;
    }
    //return matrix_convolve(src_image, filter->matrix, filter->xres, filter->yres, EDGE_WRAP);
}

static matrix_t*filter_reverse_subsample(filter_t*filter, matrix_t*src_image)
{
    matrix_t *dst_image = matrix_new(src_image->width*filter->xres,src_image->height*filter->yres);
    int x_src,y_src,x_dst,y_dst;
    for(y_src=0,y_dst=0; y_src<src_image->height; y_src++,y_dst+=filter->yres) {
        double*dst = &dst_image->data[y_dst*dst_image->width];
        for(x_src=0,x_dst=0; x_src<src_image->width; x_src++,x_dst+=filter->xres) {
            int x,y;
            double val = src_image->data[y_src*src_image->width+x_src];
            for(y=0;y<filter->yres;y++) {
                double*dst2 = &dst[y*dst_image->width+x_dst];
                for(x=0;x<filter->xres;x++) {
                    dst2[x] = val;
                }
            }
        }
    }
    return dst_image;
}

static matrix_t*filter_reverse_convolution(filter_t*filter, matrix_t*src_image)
{
    return matrix_inverse_convolve(src_image, filter->matrix, filter->xres, filter->yres, EDGE_ZERO);
    //return matrix_inverse_convolve(src_image, filter->matrix, filter->xres, filter->yres, EDGE_WRAP);
}

matrix_t*filter_apply(filter_t*filter, matrix_t*image)
{
    if(filter->type == filtertype_nop) {
        return matrix_clone(image);
    } else if(filter->type == filtertype_subsample) {
        return filter_apply_subsample(filter, image);
    } else if(filter->type == filtertype_convolve) {
        return filter_apply_convolution(filter, image);
    } else {
        fprintf(stderr, "unknown filter type %02x\n", filter->type);
        *(int*)0 = 0xdead;
    }
    return 0;
}
matrix_t*filter_reverse(filter_t*filter, matrix_t*image)
{
    if(filter->type == filtertype_nop) {
        return matrix_clone(image);
    } else if(filter->type == filtertype_subsample) {
        return filter_reverse_subsample(filter, image);
    } else if(filter->type == filtertype_convolve) {
        return filter_reverse_convolution(filter, image);
    } else {
        fprintf(stderr, "unknown filter type %02x\n", filter->type);
        *(int*)0 = 0xdead;
    }
    return 0;
}

static int count_nodes(filtertree_t*tree)
{
    if(tree->num_children) {
        int count = 0;
        /*if(tree->filter) {
            count++;
        }*/
        int t;
        for(t=0;t<tree->num_children;t++) {
            count += count_nodes(tree->children[t]);
        }
        return count;
    } else {
        return 1;
    }
}

static void apply_nodes(filtertree_t*tree, matrix_t*image, matrix_t***dest)
{
    if(!tree->num_children) {
        if(tree->filter) {
            **dest = filter_apply(tree->filter, image);
        } else {
            **dest = matrix_clone(image);
        }
        (*dest)++;
    } else {
        if(tree->filter) {
            matrix_t*new_image = filter_apply(tree->filter, image);
            //**dest = new_image;
            //(*dest)++;
            int t;
            for(t=0;t<tree->num_children;t++) {
                apply_nodes(tree->children[t], new_image, dest);
            }
        } else {
            int t;
            for(t=0;t<tree->num_children;t++) {
                apply_nodes(tree->children[t], image, dest);
            }
        }
    }
}

static void fill_nodes(filtertree_t*tree, matrix_t*image)
{
    if(!tree->num_children) {
        if(tree->filter) {
            tree->response = filter_apply(tree->filter, image);
        } else {
            tree->response = matrix_clone(image);
        }
    } else {
        if(tree->filter) {
            matrix_t*new_image = filter_apply(tree->filter, image);
            tree->response = new_image;
            int t;
            for(t=0;t<tree->num_children;t++) {
                fill_nodes(tree->children[t], new_image);
            }
        } else {
            /* dummy node, which just passes the image on to children 
               without modifying it */
            int t;
            for(t=0;t<tree->num_children;t++) {
                fill_nodes(tree->children[t], image);
            }
            /* duplicate matrix- will also have been stored in the parent */
            tree->response = matrix_clone(image);
        }
    }
}

matrix_t*reverse_nodes(filtertree_t*tree, matrix_t***matrixlist)
{
    int t;
    matrix_t*m = 0;
    if(tree->num_children) {
        if(tree->filter)
            (*matrixlist)++;
        for(t=0;t<tree->num_children;t++) {
            matrix_t*g = reverse_nodes(tree->children[t], matrixlist);
            if(!m) m = g;
            else matrix_add_inplace(m,g);
        }
    } else {
        m = **matrixlist;
        (*matrixlist)++;
    }
    if(tree->filter) {
        m = filter_reverse(tree->filter, m);
    }
    /*
    if(tree->filter && tree->num_children) {
        static int counter=1;
        char filename[80];
        //sprintf(filename, "rev%d.dat", counter);
        //matrix_save(m, filename);
        sprintf(filename, "rev%02d.png", counter);
        image_save(image_from_matrix(m,IMAGE_MINMAX), filename);
        if(tree->filter && tree->filter->matrix) {
            sprintf(filename, "filter%d.png", counter);
            image_save(image_from_matrix(tree->filter->matrix,IMAGE_MINMAX), filename);
        }
        counter++;
    }
    */
    return m;
}

matrixset_t* filtertree_apply(filtertree_t*tree, matrix_t*image)
{
    int num = count_nodes(tree);
    matrix_t**matrices = (matrix_t**)malloc(sizeof(matrix_t*)*num);
    matrixset_t*set = malloc(sizeof(matrixset_t));
    set->m = matrices;
    set->num = num;
    matrix_t**mpos = matrices;
    apply_nodes(tree, image, &mpos);
    return set;
}

void filtertree_fill(filtertree_t*tree, matrix_t*image)
{
    fill_nodes(tree, image);
}

matrix_t*filtertree_reverse(filtertree_t*tree, matrixset_t*set)
{
    int num = count_nodes(tree);
    matrix_t**mpos = set->m;
    return reverse_nodes(tree, &mpos);
}

filtertree_t*filtertree_new(int num_children, filter_t*filter)
{
    filtertree_t*tree = malloc(sizeof(filtertree_t));
    memset(tree, 0, sizeof(filtertree_t));
    tree->num_children = num_children;
    tree->filter = filter;
    if(tree->num_children) {
        tree->children = malloc(sizeof(filtertree_t)*num_children);
        memset(tree->children,0,sizeof(filtertree_t)*num_children);
    } else {
        tree->children = 0;
    }
    return tree;
}
matrixset_t*filtertree_getfilters(filtertree_t*tree)
{
    matrixset_t*filters = matrixset_new(tree->num_children);
    int t;
    for(t=0;t<tree->num_children;t++) {
        if(!tree->children[t]->filter || tree->children[t]->num_children) {
            fprintf(stderr, "Warning: recursive filtertree collapsing not implemented yet\n");
            return 0;
        }
        filters->m[t] = matrix_clone(tree->children[t]->filter->matrix);
    }
    return filters;
}
filter_t*filter_new(char type, char flags, matrix_t*matrix, int xres, int yres)
{
    filter_t*f = malloc(sizeof(filter_t));
    memset(f, 0, sizeof(filter_t));
    f->type = type;
    f->flags = flags;
    f->matrix = matrix;
    f->xres = xres;
    f->yres = yres;
    return f;
}
filter_t*filter_new_convolution(matrix_t*matrix, int xres, int yres)
{
    filter_t*f = malloc(sizeof(filter_t));
    memset(f, 0, sizeof(filter_t));
    f->type = filtertype_convolve;
    f->flags = 0;
    f->matrix = matrix;
    f->xres = xres;
    f->yres = yres;
    return f;
}
