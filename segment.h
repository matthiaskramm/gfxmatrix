/* segment.h 

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

#ifndef __textonutils_c__
#define __textonutils_c__

#include "filter.h"

#define FILTER_GABOR 1
#define FILTER_GAUSS_DIFF 2 
#define FILTER_GAUSS_HIGHPASS 4
#define FILTER_GAUSS 8
#define FILTER_ALL 255
#define CLUSTER_MAHALANOBIS 1
#define VERBOSE 2

filtertree_t* makeFilterTree(int flags);
filtertree_t* makeEdgeFilter();
filtertree_t* makeSimpleFilterTree(char gauss, char lgauss, char gauss_diff);
bytearray_t* kmeans(matrixset_t*set, int num_centers, matrix_t**centersptr, int flags);
matrix_t* mkTexton(filtertree_t*tree, double*values);
bytearray_t* kmeans_bytearray(bytearrayset_t*set, int num_centers, bytearray_t**centersptr, int flags, int dist);
int update_kmeans_map(bytearray_t*map, bytearrayset_t*set, bytearray_t*centers, int maxdist);
int kmeans_add_rects(bytearray_t*map, bytearrayset_t*set, bytearray_t*centers, int blockx, int blocky, int dist);
int kmeans_add_rects_fast(bytearray_t*map, bytearrayset_t*set, bytearray_t*centers, int blockx, int blocky, int dist);
int kmeans_add_rects_ho(bytearray_t*map, bytearrayset_t*set, bytearray_t*centers, int blockx, int blocky, int dist);
void image_mark(image_t*img, bytearray_t*map);

window_t window_new(int x1,int y1,int x2,int y2);
window_t find_max_window(bytearray_t*map, int marker);

#define SEGMENT_SPATIAL 1
#define SEGMENT_UV 2
#define SEGMENT_MARGINAL 4
#define SEGMENT_DEFAULTS (SEGMENT_UV|SEGMENT_MARGINAL)
bytearray_t* segment(image_t*img, int num_centers, int bx, int by, int dist, bytearray_t**centersdest, int flags);
#endif
