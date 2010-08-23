/* filter_misc.h 

   Filter kernels.
   
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

#ifndef __filter_misc_h__
#define __filter_misc_h__

#include "filter.h"

matrix_t* gabor_filter(int width, int height, float ax, float ay, float phase_shift, float ifreq, float angle, char complex_counterpart);
complex_matrix_t* gabor2_filter(int width, int height, float a, float f1, float f2);
matrix_t* gauss_diff(int width, int height, float a, float angle, float len);
matrix_t* gauss_diff2(int width, int height, float d, float angle);
matrix_t* gauss_highpass(int width, int height, float r0, float a);
matrix_t* gauss_filter(int width, int height, float d);
matrix_t* lgauss_filter(int width, int height, float d);
complex_matrix_t*gauss_diff_complex(int width, int height, double a);
filtertree_t* makeGaborTree(int width, int height, int xx, int yy);
complex_matrixset_t* filterset_new_dct(int width, int height, int stepx, int stepy, char complexversion);
#endif
