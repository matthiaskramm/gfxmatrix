/* matrix.h 

   Matrix and linalg functions.

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

#ifndef __matrix_h__
#define __matrix_h__

#include "common.h"

void matrix_clear(matrix_t*m);
void matrix_delete(matrix_t*m);
matrix_t*matrix_clone(matrix_t*m);
matrix_t*matrix_new(int width, int height);
matrix_t*matrix_new_random(int width, int height, double min, double max);
matrix_t*matrix_new_gaussrandom(int width, int height, double mean, double var);
matrix_t*matrix_new_unit(int width, int height);
matrix_t*matrix_new_fromdata(double*data, int width, int height);
double*matrix_multiply_with_vector(matrix_t*m, double*vector);
void matrix_add_inplace(matrix_t*m,matrix_t*add);
void matrix_save(matrix_t*m, char*filename);
void matrix_print(matrix_t*m);
void complex_matrix_print(complex_matrix_t*m);

complex_matrix_t*complex_matrix_new(int width, int height);
complex_matrix_t*complex_matrix_clone(complex_matrix_t*m);
void complex_matrix_delete(complex_matrix_t*m);

double matrix_diff(matrix_t*m1, matrix_t*m2);

matrix_t* matrix_invert(matrix_t*m);

matrix_t* matrix_getautocorrelation(matrix_t*img, int csizex, int csizey);
double matrix_getcrosscorrelation(matrix_t*m1, matrix_t*m2);
matrix_t* matrix_getcrosscorrelations(matrixset_t*set1, matrixset_t*set2, char subtract_mean);

void matrix_symm_geteigenvectorfactorization(matrix_t*m, matrix_t**dest1, matrix_t**dest2);
void matrix_symm_geteigenvectorfactorization2(matrix_t*m, matrix_t**dest1, matrix_t**dest2);

/* splits the matrix into a vector product xx^T */
double*matrix_split_xx(matrix_t*_m);

#ifdef HAVE_GSL
/* alternative implementation, based on biggest eigenvalue/eigenvector */
double* matrix_split_xx_2(matrix_t*_m);
//gsl_matrix*matrix_to_gsl(matrix_t*m);
#endif

double*matrix_solve(matrix_t*_A, double*_b);
double*matrix_solve_approx(matrix_t*_A, double*_b);
double*matrix_solve_underdetermined(matrix_t*m, double*b);

void matrix_store(matrix_t*m, char*filename);
void bytearray_store(bytearray_t*m, char*filename);

#define EDGE_ZERO 0
#define EDGE_REFLECT 1
#define EDGE_WRAP 2

matrix_t*matrix_convolve(matrix_t*src_image, matrix_t*m, int xres, int yres, char edgetype);
matrix_t*matrix_inverse_convolve(matrix_t*src_image, matrix_t*m, int xres, int yres, char edgetype);
matrix_t*matrix_convolve_fft(matrix_t*img, matrix_t*m);

typedef struct _convolve {
    void*internal;
} convolve_t;

typedef struct _filterbank {
    complex_matrix_t**c;
    int width;
    int height;
    int num;
    void*internal;
} filterbank_t;

filterbank_t*filterbank_new(complex_matrixset_t*set, int width, int height, int flags);
bytearrayset_t*filterbank_apply_tobytearray(filterbank_t*f, complex_matrix_t*cm);
void filterbank_delete(filterbank_t*fbank);

convolve_t*convolve_new(matrix_t*filter, int width, int height);
convolve_t*convolve_new_complex(complex_matrix_t*filter, int width, int height);
matrix_t*convolve_apply(convolve_t*m, complex_matrix_t*fft_image);
complex_matrix_t*convolve_apply_complex(convolve_t*c, complex_matrix_t*ff_img);
bytearray_t*convolve_apply_complex_tobytearray(convolve_t*c, complex_matrix_t*ff_img);
void convolve_delete(convolve_t*m);

matrixset_t*matrixset_new(int num);
matrixset_t*matrixset_clone(matrixset_t*src);
void matrixset_delete(matrixset_t*src);
matrix_t* matrixset_getmean(matrixset_t*set);
matrixset_t* matrixset_new_alloc(int num, int width, int height);

complex_matrixset_t*complex_matrixset_new(int num);
complex_matrixset_t*complex_matrixset_clone(complex_matrixset_t*src);
void complex_matrixset_delete(complex_matrixset_t*src);

matrix_t*matrix_multiply(matrix_t*m1, matrix_t*m2);
matrix_t*matrix_multiply_t(matrix_t*m1, char t1, matrix_t*m2, char t2);
matrix_t*matrix_add(matrix_t*m1, double v1, matrix_t*m2, double v2);
matrix_t*matrix_transpose(matrix_t*m);
void matrix_update(matrix_t*dest, matrix_t*src);

matrix_t* complex_matrix_abs(complex_matrix_t*src);

complex_matrix_t* matrix_fft(matrix_t*src);
complex_matrix_t* matrix_ifft(complex_matrix_t*src);
matrix_t* matrix_ifft_real(complex_matrix_t*src);
complex_matrix_t* complex_matrix_fft(complex_matrix_t*src);
void complex_matrix_save(complex_matrix_t*b, char*filename);

/* only on matrix of even sizes, only by common divisors */
matrix_t*matrix_scaledown(matrix_t*m, int width, int height);
matrix_t*matrix_scaleup(matrix_t*m, int width, int height);


#define EXPAND_CENTRAL 0
#define EXPAND_RIGHTDOWN 1
#define FILTERBANK_SQRT 4
complex_matrix_t*complex_matrix_expand(complex_matrix_t*m, int width, int height, int flags);

matrix_t* complex_matrix_realpart(complex_matrix_t*src);
matrix_t* complex_matrix_imagpart(complex_matrix_t*src);
matrix_t* complex_matrix_ifft_real(complex_matrix_t*src);
complex_matrix_t* complex_matrix_ifft(complex_matrix_t*src);

void matrix_print_stats(matrix_t*m);

bytearrayset_t* bytearrayset_new(int num);
void bytearrayset_delete(bytearrayset_t*b);
bytearray_t*bytearray_maxfilter(bytearray_t*b);
void bytearray_delete(bytearray_t*m);
bytearray_t* bytearray_cut(bytearray_t*src, int x, int y, int width, int height);
bytearray_t* bytearray_new(int width, int height);
bytearray_t* bytearray_new_fromdata(unsigned char*data, int width, int height);
bytearrayset_t* bytearrayset_new_alloc(int num, int width, int height);
void bytearray_save(bytearray_t*b, char*filename);
unsigned long* bytearray_integral(bytearray_t*b);
bytearray_t*bytearray_clone(bytearray_t*b);
bytearrayset_t*bytearrayset_select_n(bytearrayset_t*b, int newwidth, int newheight);
void bytearrayset_extend(bytearrayset_t*s, int num);
void bytearrayset_extend_alloc(bytearrayset_t*s, int num, int width, int height);
void bytearrayset_append(bytearrayset_t*s, bytearrayset_t*a);
void bytearray_normalize(bytearray_t*b);
bytearrayset_t*bytearrayset_select_n(bytearrayset_t*b, int newwidth, int newheight);
void bytearray_average_rect(bytearray_t*in, int width, int height);
bytearray_t* bytearray_combine_maps(bytearray_t*b1, bytearray_t*b2);


#endif
