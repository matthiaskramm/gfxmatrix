/* texture.h 

   Texture synthesis.

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

#ifndef __texture_h__
#define __texture_h__
#include "stat.h"
#include "filter.h"

typedef struct _textureparameters
{
    int bands;
    filtertree_t*filter;
    matrixset_t*interband_crosscorr;
    matrixset_t*intraband_crosscorr;
    int autox,autoy;
    matrixset_t*image_autocorr;

    statistics_t**stat;
    statistics_t*imgstat;

    matrixset_t*original;
} textureparameters_t;

textureparameters_t*textureparameters_new_fromimage(image_t*img);
void textureparameters_adjust_crosscorrelation(textureparameters_t*p, matrixset_t*s, int band);
void textureparameters_adjust_statistics(textureparameters_t*p, matrixset_t*s, int image);
void textureparameters_adjust_autocorrelation(textureparameters_t*p, matrixset_t*s, int image);
void textureparameters_adjust_finalimage(textureparameters_t*p, matrix_t*image);
void textureparameters_copy_from_original(textureparameters_t*p, matrixset_t*s, int image);
void textureparameters_iterate(textureparameters_t*p, matrix_t*m);
#endif
