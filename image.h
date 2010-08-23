/* image.h 

   Image loading, saving and to/from matrix conversion.

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

#ifndef __image_h__
#define __image_h__

#include "matrix.h"

#define IMAGE_GRAY 1
#define IMAGE_RED 2
#define IMAGE_GREEN 4
#define IMAGE_BLUE 8
#define IMAGE_NORMALIZE 16
#define IMAGE_OPENGL 32
#define IMAGE_MARKNEGATIVE 64
#define IMAGE_LAB_L 1024
#define IMAGE_LAB_A 2048
#define IMAGE_LAB_B 4096
#define IMAGE_YUV_U 8192
#define IMAGE_YUV_V 16384
#define IMAGE_YUV_Y 32768
#define IMAGE_DIV2  65536

matrix_t* image_extractchannel(image_t*img, int mode);
bytearray_t* image_getchannel(image_t*img, int mode);
image_t*image_load(char*filename);
image_t*image_new(int width, int height);
image_t*image_clone(image_t*img);
void image_crop(image_t*img, int width, int height);
void image_delete(image_t*img);

#define IMAGE_CLAMP 1
#define IMAGE_MINMAX 2
#define IMAGE_MARKOVERFLOW 4
#define IMAGE_ABSMINMAX 8
image_t*image_from_matrix(matrix_t*m, int flags);
image_t*image_from_complex_matrix(complex_matrix_t*m, int flags);
void image_update_from_matrix(image_t*img, matrix_t*m, int flags, int xpos, int ypos);
image_t* image_cut(image_t*src, int x, int y, int width, int height);
void image_paste(image_t*dest, int x, int y, image_t*src);
void image_save(image_t*img, char*filename);
void image_save_and_free(image_t*img, char*filename);

void rgb2rgba(unsigned char*src, unsigned char*dest, int len);

image_t* image_extract(image_t*img, window_t*w);

void bytearray_save_map(bytearray_t*b, char*filename);
bytearray_t* bytearray_load_map(char*filename);

void image_insert_tiling(image_t*img, bytearray_t*map, image_t*replacement, int t);

double image_compare(image_t*src, image_t*dest);

HLS rgb2hls(RGBA c);

#define clamp00ff(x) ((x)>=0?((x)<256?(x):255):0)

#endif
