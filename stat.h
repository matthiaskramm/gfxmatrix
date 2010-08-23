/* stat.h 

   2nd, 3rd and 4th order statistics.

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

#ifndef __stat_h__
#define __stat_h__

#define STAT_MEAN 1
#define STAT_VARIANCE 2
#define STAT_SKEW 3
#define STAT_KURTOSIS 4

typedef struct _statistics
{
    double v[16];
    double mean,absmean,dev,var,skew,kurt;
    double min,max;
} statistics_t;

void statistics_delete(statistics_t*s);
statistics_t*statistics_merge(statistics_t*s1, statistics_t*s2);
statistics_t*statistics_new_frommatrix(matrix_t*m);
void matrix_normalize(matrix_t*m);
void matrix_adjust_statistics(matrix_t*m, statistics_t*s);
void matrix_adjust_minmax(matrix_t*m, statistics_t*dest);
void statistics_print(statistics_t*stat);
void matrix_adjust_variance(matrix_t*msrc, double var);
void matrix_adjust_mean(matrix_t*msrc, double mean);

bytearrayset_t*matrix_blockstats(matrix_t*m, int w, int h);
#endif

