/* txtsynth.h

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

#ifndef __txtsynth_h__
#define __txtsynth_h__

#include "image.h"

typedef struct _texturedata {
    int type;
    void*internal;
} texturedata_t;

texturedata_t* texturedata_fromimage(image_t*image, int method);
image_t* texturedata_synthesize_image(texturedata_t*data, image_t*alpha);
void texturedata_delete(texturedata_t*data);

void texturedata_set_parameter(char*key, char*value);

image_t* synth_simple(image_t*alpha, image_t*img);

double texturedata_dist(texturedata_t*txt1, texturedata_t*txt2);

#endif //__txtsynth_h__
