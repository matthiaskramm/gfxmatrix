/* filter.h 

   Filters and pyramids of filters.

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

#ifndef __filter_h__
#define __filter_h__

#include "common.h"

enum {filtertype_nop, filtertype_subsample, filtertype_convolve} filtertype;

#define FILTER_HIGHBAND 1
#define FILTER_LOWBAND 2
#define FILTER_DIRECTIONAL 4
#define FILTER_INTERMEDIATE 8
#define FILTER_COMPLEX 64

typedef struct _filter
{
    char type;
    char flags;
    matrix_t* matrix;
    int xres,yres;
    void*internal;
} filter_t;

typedef struct _filterset
{
    filter_t**f;
    int num;
} filterset_t;

typedef struct _filtertree
{
    filter_t* filter;
    int num_children;
    struct _filtertree ** children;
    char*description;

    matrix_t*response;
} filtertree_t;

matrix_t*filter_apply(filter_t*filter, matrix_t*image);
matrix_t*filter_reverse(filter_t*filter, matrix_t*image);
matrixset_t*filtertree_apply(filtertree_t*tree, matrix_t*image);
void filtertree_fill(filtertree_t*tree, matrix_t*image);
matrix_t*filtertree_reverse(filtertree_t*tree, matrixset_t*set);

filter_t*filter_new(char type, char flags, matrix_t*matrix, int xres, int yres);
filter_t*filter_new_convolution(matrix_t*matrix, int xres, int yres);
filtertree_t*filtertree_new(int num_children, filter_t*filter);
matrixset_t*filtertree_getfilters(filtertree_t*tree);

#endif
