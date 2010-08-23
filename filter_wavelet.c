/* filter_wavelet.c

   Wavelet filters.

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

#include "common.h"
#include "filter.h"

//--------------- Daubechies 9/7 Wavelet Transform

static double unity_data[1*1] = {1};
static matrix_t m_unity = {unity_data, 1,1};
static filter_t unity = {filtertype_convolve, FILTER_INTERMEDIATE, &m_unity, 1, 1};

double band97_h_data_x[9*9] = {
    0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,
    0.037828, -0.023849, -0.110624, 0.377402, 0.852699, 0.377402, -0.110624, -0.023849, 0.037828,
    0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,
};
double band97_h_data_y[9*9] = {
    0,0,0,0, 0.037828,0,0,0,0,
    0,0,0,0,-0.023849,0,0,0,0,
    0,0,0,0,-0.110624,0,0,0,0,
    0,0,0,0, 0.377402,0,0,0,0,
    0,0,0,0, 0.852699,0,0,0,0,
    0,0,0,0, 0.377402,0,0,0,0,
    0,0,0,0,-0.110624,0,0,0,0,
    0,0,0,0,-0.023849,0,0,0,0,
    0,0,0,0, 0.037828,0,0,0,0,
};
matrix_t m_band97_hx = {band97_h_data_x, 9,9};
filter_t band97_hx = {filtertype_convolve, FILTER_HIGHBAND, &m_band97_hx, 1, 1};
matrix_t m_band97_hy = {band97_h_data_y, 9,9};
filter_t band97_hy = {filtertype_convolve, FILTER_HIGHBAND, &m_band97_hy, 1, 1};

double band97_g_data[7*7] = {
    -0.0645390 * -0.064539 , -0.040689 * -0.064539, 0.418092 * -0.064539, 0.788486 * -0.064539, 0.418092 * -0.064539, -0.040689 * -0.064539, -0.0645390 * -0.064539, 
    -0.0645390 * -0.040689 , -0.040689 * -0.040689, 0.418092 * -0.040689, 0.788486 * -0.040689, 0.418092 * -0.040689, -0.040689 * -0.040689, -0.0645390 * -0.040689, 
    -0.0645390 *  0.418092 , -0.040689 *  0.418092, 0.418092 *  0.418092, 0.788486 *  0.418092, 0.418092 *  0.418092, -0.040689 *  0.418092, -0.0645390 *  0.418092, 
    -0.0645390 *  0.788486 , -0.040689 *  0.788486, 0.418092 *  0.788486, 0.788486 *  0.788486, 0.418092 *  0.788486, -0.040689 *  0.788486, -0.0645390 *  0.788486, 
    -0.0645390 *  0.418092 , -0.040689 *  0.418092, 0.418092 *  0.418092, 0.788486 *  0.418092, 0.418092 *  0.418092, -0.040689 *  0.418092, -0.0645390 *  0.418092, 
    -0.0645390 * -0.040689 , -0.040689 * -0.040689, 0.418092 * -0.040689, 0.788486 * -0.040689, 0.418092 * -0.040689, -0.040689 * -0.040689, -0.0645390 * -0.040689, 
    -0.0645390 * -0.064539 , -0.040689 * -0.064539, 0.418092 * -0.064539, 0.788486 * -0.064539, 0.418092 * -0.064539, -0.040689 * -0.064539, -0.0645390 * -0.064539, 
};
matrix_t m_band97_g = {band97_g_data, 7,7};
filter_t band97_g = {filtertype_convolve, FILTER_LOWBAND, &m_band97_g, 2, 2};

//--------------- LeGall 5/3 Wavelet Transform

double band53_h_data_x[3*3] = {
    0,0,0,
    -0.5,1.0,-0.5,
    0,0,0,
    //-0.176777, 0.353553, 1.060660, 0.353553, -0.176777,
};
double band53_h_data_y[3*3] = {
    0, -0.5, 0,
    0,  1.0, 0,
    0, -0.5, 0,
};
matrix_t m_band53_hx = {band53_h_data_x, 3,3};
filter_t band53_hx = {filtertype_convolve, FILTER_HIGHBAND, &m_band53_hx, 1, 1};
matrix_t m_band53_hy = {band53_h_data_y, 3,3};
filter_t band53_hy = {filtertype_convolve, FILTER_HIGHBAND, &m_band53_hy, 1, 1};

double band53_g_data[5*5] = {
    /*-1/8.0*-1/8.0, -2/8.0*-1/8.0, 6/8.0*-1/8.0, -2/8.0*-1/8.0, -1/8.0*-1/8.0,
    -1/8.0*-2/8.0, -2/8.0*-2/8.0, 6/8.0*-2/8.0, -2/8.0*-2/8.0, -1/8.0*-2/8.0,
    -1/8.0* 6/8.0, -2/8.0* 6/8.0, 6/8.0* 6/8.0, -2/8.0* 6/8.0, -1/8.0* 6/8.0,
    -1/8.0*-2/8.0, -2/8.0*-2/8.0, 6/8.0*-2/8.0, -2/8.0*-2/8.0, -1/8.0*-2/8.0,
    -1/8.0*-1/8.0, -2/8.0*-1/8.0, 6/8.0*-1/8.0, -2/8.0*-1/8.0, -1/8.0*-1/8.0,*/
    -0.176777*-0.176777, 0.353553*-0.176777, 1.060660*-0.176777, 0.353553*-0.176777, -0.176777*-0.176777,
    -0.176777* 0.353553, 0.353553* 0.353553, 1.060660* 0.353553, 0.353553* 0.353553, -0.176777* 0.353553,
    -0.176777* 1.060660, 0.353553* 1.060660, 1.060660* 1.060660, 0.353553* 1.060660, -0.176777* 1.060660,
    -0.176777* 0.353553, 0.353553* 0.353553, 1.060660* 0.353553, 0.353553* 0.353553, -0.176777* 0.353553,
    -0.176777*-0.176777, 0.353553*-0.176777, 1.060660*-0.176777, 0.353553*-0.176777, -0.176777*-0.176777,
};
matrix_t m_band53_g = {band53_g_data, 5,5};
filter_t band53_g = {filtertype_convolve, FILTER_LOWBAND, &m_band53_g, 2, 2};

//---------------

filtertree_t*filtertree_new_wavelet53(int depth)
{
    filtertree_t*r,*n;
    filtertree_t*base = n = filtertree_new(3, 0);
    int t;
    for(t=0;t<depth;t++) {
        n->children[0] =     filtertree_new(0, &band53_hx);r=n;n->description="highbandx";
        n->children[1] =     filtertree_new(0, &band53_hy);r=n;n->description="highbandy";
        n->children[2] = n = filtertree_new(t<depth-1?3:0, &band53_g);r=n;n->description="lowband";
    }
    return base;
}

filtertree_t*filtertree_new_wavelet97(int depth)
{
    filtertree_t*r,*n;
    filtertree_t*base = n = filtertree_new(3, 0);
    int t;
    for(t=0;t<depth;t++) {
        n->children[0] =     filtertree_new(0, &band97_hx);r=n;n->description="highbandx";
        n->children[1] =     filtertree_new(0, &band97_hy);r=n;n->description="highbandy";
        n->children[2] = n = filtertree_new(t<depth-1?3:0, &band97_g);r=n;n->description="lowband";
    }
    return base;
}
