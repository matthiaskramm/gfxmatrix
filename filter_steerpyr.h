/* filter_steerpyr.h 

   Steerable pyramid filters.

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

#ifndef __filter_steerpyr_h__
#define __filter_steerpyr_h__
filtertree_t*filtertree_new_steerpyr(int depth);
filtertree_t*filtertree_new_pyr(int depth);
complex_matrixset_t* filter_steerpyr(matrix_t*img, int bands, int dirs);
#endif
