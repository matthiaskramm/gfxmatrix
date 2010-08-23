/* graphcut.h 

   2D graph cut library.

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

#ifndef __graphcut_h__
#define __graphcut_h__

#include "../common.h"
#include "matrix.h"
#include "image.h"

typedef signed int weight_t;
#define MAX_WEIGHT 0x07ffffff

#define RIGHT 0
#define UP 1
#define LEFT 2
#define DOWN 3

typedef struct _map {
    weight_t*(weight[4]);
    int width;
    int height;
} map_t;

map_t* map_new(int width, int height);
map_t* map_from_matrix(matrix_t*m);
map_t* map_from_image(image_t*img);
void inverse_map(map_t*map);
void map_delete(map_t*);

unsigned char* find_cut(map_t*map, int x1, int y1, int x2, int y2);

#endif
