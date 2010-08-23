/* corr.h 

   Matrix autocorrelation and crosscorrelation adjustment.

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

#ifndef __corr_h__
#define __corr_h__
void matrix_adjustautocorrelation(matrix_t*img, matrix_t*corr);
void matrix_adjustautocorrelation2(matrix_t*img, matrix_t*corr);
void matrix_adjustcrosscorrelation(matrixset_t*varset, matrix_t*varcorr, 
                                   matrixset_t*staticset, matrix_t*staticcorr);
void matrix_adjustautocorrelation3(matrix_t*img, matrix_t*corr);

#endif
