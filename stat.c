/* stat.c 

   2nd, 3rd and 4th order function adjustment.

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

#include <gsl/gsl_poly.h>
#include <gsl/gsl_errno.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#define __USE_ISOC99
#include <math.h>
#include "matrix.h"
#include "stat.h"

void statistics_delete(statistics_t*s)
{
    free(s);
}
statistics_t*statistics_new()
{
    statistics_t*st = (statistics_t*)malloc(sizeof(statistics_t));
    memset(st, 0, sizeof(statistics_t));
    return st;
}
statistics_t*statistics_new_frommatrix(matrix_t*m)
{
    int s,t;
    int l = m->width*m->height;
    statistics_t*st = (statistics_t*)malloc(sizeof(statistics_t));
    memset(st, 0, sizeof(statistics_t));
    double min=HUGE_VAL,max=-HUGE_VAL;
    double mean = 0;
    double absmean = 0;
    for(t=0;t<l;t++) {
        mean += m->data[t];
        absmean += fabs(m->data[t]);
    }
    mean /= l;
    absmean /= l;
    for(t=0;t<l;t++) {
        //double p = m->data[t]-mean, p2 = 1;
        double p = m->data[t], p2 = 1;
        if(p > max) max = p;
        if(p < min) min = p;
        p-=mean;
        for(s=0;s<16;s++) {
            st->v[s] += p2;
            p2 *= p;
        }
    }
    for(s=0;s<16;s++) {
        st->v[s] /= l;
    }
    st->mean = st->v[1] = mean;
    st->absmean = absmean;
    st->var = st->v[2];
    st->dev = sqrt(st->var);
    if(st->var < 1e-4) {
        st->skew = 0;
        st->kurt = 0;
    } else {
        st->skew = st->v[3] /** (l*(l-1))*/ / (st->dev*st->dev*st->dev);
        st->kurt = st->v[4] /** (l*(l-1))*/ / (st->var*st->var) /* -3 */;
    }

    st->min = min;
    st->max = max;
    return st;
}
void statistics_print(statistics_t*stat)
{
    printf("mean: %f dev/var: %f/%f skew: %f curt: %f min: %f max: %f\n", stat->mean, stat->dev, stat->var, stat->skew, stat->kurt, stat->min, stat->max);
    //printf("v0: %f v1: %f v2: %f v3: %f v4: %f\n", stat->v[0], stat->v[1], stat->v[2], stat->v[3], stat->v[4]);
}
void statistics_print_comparison(statistics_t*a, statistics_t*b)
{
    printf("mean: %f/%f var: %f/%f skew: %f/%f curt: %f/%f min: %f/%f max: %f/%f\n", 
	    a->mean, b->mean,
	    a->var, b->var,
	    a->skew, b->skew,
	    a->kurt, b->kurt, 
	    a->min, b->min, a->max, b->max);
}
void matrix_normalize(matrix_t*m)
{
    statistics_t*st = statistics_new_frommatrix(m);
    int s,t;
    int l = m->width*m->height;
    
    double dev = sqrt(st->v[STAT_VARIANCE]);
    double mean = st->v[STAT_MEAN];
    for(t=0;t<l;t++) {
        m->data[t] = (m->data[t] - mean) / dev;
    }
    statistics_delete(st);
}

static inline double sqr(const double a) {return a*a;}

void matrix_adjust_skewness_old(matrix_t*msrc, double skew)
{
    int l = msrc->width*msrc->height;
    /* adjust skewness */
    statistics_t*src = statistics_new_frommatrix(msrc);
    gsl_poly_complex_workspace* w = gsl_poly_complex_workspace_alloc(4);
    double z[4*2];
    memset(z, 0, sizeof(z));
    double a[4] = {src->v[3]-(skew*src->dev*src->dev*src->dev), 3*src->v[4], 3*src->v[5], src->v[6]};
    gsl_poly_complex_solve(a, 4, w, z);
    gsl_poly_complex_workspace_free(w);
    double shift2 = 0;
    int t;
    for(t=0;t<3;t++) {
	if(fabs(z[t*2+1]) < 0.001) {
	    shift2 = z[t*2];
	}
    }
    for(t=0;t<l;t++) {
	double v = msrc->data[t];
	v = v + shift2*v*v;
	msrc->data[t] = v;
    }
    statistics_delete(src);
}

void matrix_adjust_kurtosis_old(matrix_t*msrc, double k)
{
#if 0

        /* adjust kurtosis */
        src = statistics_new_frommatrix(m);
        gsl_poly_complex_workspace* w5 = gsl_poly_complex_workspace_alloc(5);
        double z5[5*2];
        memset(z5, 0, sizeof(z5));
        //* c[0] + c[1] x + c[2] x^2 + ... + c[len-1] x^(len-1)
        //double a[4] = {-2,-1,2,1}; /* should have -1,1,-2 as zeroes */
        double a5[5] = {src->v[4]-dest->v[4], 4*src->v[6], 6*src->v[8], 4*src->v[10], src->v[12]};
        /*for(t=0;t<5;t++) {
            printf("%f*x^%d +", a5[t], t);
        } printf("\n");*/
        gsl_poly_complex_solve(a5, 5, w5, z5);
        gsl_poly_complex_workspace_free(w5);
        double shift3 = 0;
        for(t=0;t<4;t++) {
            printf("%f %f\n", z5[t*2+0], z5[t*2+1]);
            if(fabs(z5[t*2+1]) < 0.001) {
                shift3 = z5[t*2];
            }
        }
        //printf("using %f as shift3\n", shift3);
        for(t=0;t<l;t++) {
            double v = m->data[t];// -  src->v[STAT_MEAN];
            //v = v + shift*v + shift2*v*v;
            v = v + shift3*v*v*v;
            m->data[t] = v;// + mean_dest;
        }
        statistics_delete(src);

#endif 
}

static int gsl_error_occured = 0;
static char currentaction[4000];
void handleerror(const char * reason, const char * file, int line, int gsl_errno)
{
    printf("FAILED(%d): %s\n", gsl_errno, currentaction);
    printf("%s\n", reason);
    gsl_error_occured = 1;
#ifdef TOUCHY
    exit(1);
#endif
}

void matrix_adjust_skewness(matrix_t*msrc, double sk)
{
    int l = msrc->width*msrc->height;
    /* adjust skewness*/
    statistics_t* src = statistics_new_frommatrix(msrc);

    if(src->var < 1e-4) {
        printf("not adjusting skewness of matrix with variance %f\n", src->var);
        return;
    }
    if(!isfinite(sk)) {
        printf("not adjusting to skewness %f\n", sk);
        return;
    }

    double s = src->skew;
    double*m = src->v;
    double sd = src->dev;
    //printf("sd: %f s: %f\n", sd, s);

    double A = m[6]-3*sd*s*m[5]+3*sd*sd*(s*s-1)*m[4]+sd*sd*sd*sd*sd*sd*(2+3*s*s-s*s*s*s);
    double B = 3*(m[5]-2*sd*s*m[4]+sd*sd*sd*sd*sd*s*s*s);
    double C = 3*(m[4]-sd*sd*sd*sd*(1+s*s));
    double D = s*sd*sd*sd;

    //printf("A: %f B: %f C: %f D: %f\n", A, B, C, D);

    double a7=A*A;
    double a6=2*A*B;
    double a5=B*B+2*A*C;
    double a4=2*(A*D+B*C);
    double a3=C*C+2*B*D;
    double a2=2*C*D;
    double a1=D*D;

    double A2=sd*sd;
    double B2=m[4]-(1+s*s)*sd*sd*sd*sd;

    double b7=B2*B2*B2;
    double b5=3*A2*B2*B2;
    double b3=3*A2*A2*B2;
    double b1=A2*A2*A2;

    double d[8];
    d[7] = B*b7;
    d[6] = 2*C*b7 - A*b5;
    d[5] = 3*D*b7;
    d[4] = C*b5 - 2*A*b3;
    d[3] = 2*D*b5 - B*b3;
    d[2] = -3*A*b1;
    d[1] = D*b3 - 2*B*b1;
    d[0] = -C*b1;

    gsl_poly_complex_workspace* w8 = gsl_poly_complex_workspace_alloc(8);
    double z[8*2];
    memset(z, 0, sizeof(z));
    int t;
    for(t=0;t<8;t++) {
        if(!isfinite(d[t])) {
            printf("Image contains NaN or Inf values, not adjusting skew\n");
#ifdef TOUCHY
            exit(1);
#endif
            return;
        }
    }
    
    gsl_error_handler_t*oldhandler = gsl_set_error_handler(handleerror);
    gsl_error_occured = 0;
    sprintf(currentaction, "skew1 %f %f %f %f %f %f %f %f\n", d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]);

    gsl_poly_complex_solve(d, 8, w8, z);
    gsl_set_error_handler(oldhandler);
    if(gsl_error_occured)
        return;

    double lmi = -1e33;
    double lma = 1e33;

    for(t=0;t<7;t++) {
	//printf("%f+%fi\n", z[t*2], z[t*2+1]);
	if(fabs(z[t*2+1]) > fabs(z[t*2])*1e-6)
	    continue;
	if(z[t*2] < 0 && z[t*2] > lmi)
	    lmi = z[t*2];
	if(z[t*2] >= 0 && z[t*2] < lma)
	    lma = z[t*2];
    }
    if(lmi < -1e32) {
	lmi = -1e16;
    }
    if(lma > 1e32) {
	lma = 1e16;
    }
    //printf("lmi: %f lma: %f\n", lmi, lma);
    double skmin = (D + lmi*(C + lmi*(B + lmi*A))) / sqrt(b7*lmi*lmi*lmi*lmi*lmi*lmi + b5*lmi*lmi*lmi*lmi + b3*lmi*lmi + b1);
    double skmax = (D + lma*(C + lma*(B + lma*A))) / sqrt(b7*lma*lma*lma*lma*lma*lma + b5*lma*lma*lma*lma + b3*lma*lma + b1);
    if(skmin > skmax) {double r = skmin;skmin = skmax;skmax = r;}
    //printf("skmin: %f skmax: %f\n", skmin, skmax);
    

    double lambda = 0;
    if(sk<skmin || sk>skmax) {
	if(sk<skmin)
	    lambda = lmi;
	else
	    lambda = lma;
	//fprintf(stderr, "Warning: Target skew (%f) not in range [%f %f]\n", sk, skmin, skmax);
    } else {
	double c7[7];
	c7[0] = a1-b1*sk*sk;
	c7[1] = a2;
	c7[2] = a3-b3*sk*sk;
	c7[3] = a4;
	c7[4] = a5-b5*sk*sk;
	c7[5] = a6;
	c7[6] = a7-b7*sk*sk;

	gsl_poly_complex_workspace* w7 = gsl_poly_complex_workspace_alloc(7);
        gsl_error_handler_t*oldhandler = gsl_set_error_handler(handleerror);
        gsl_error_occured = 0;
        sprintf(currentaction, "skew2 %f %f %f %f %f %f %f\n", c7[0], c7[1], c7[2], c7[3], c7[4], c7[5], c7[6]);
	gsl_poly_complex_solve(c7, 7, w7, z);
        gsl_set_error_handler(oldhandler);
        if(gsl_error_occured)
            return;
	gsl_poly_complex_workspace_free(w7);
       
	lambda = 1e32;
	char haszero = 0;
	for(t=0;t<6;t++) {
	    if(z[t*2+1] > 0.00001 || (z[t*2]*(sk-s)) < 0)
		continue;
	    double x = z[t*2];
	    double foo = A*x*x*x + B*x*x + C*x + D;
	    if(foo==0.0)
		haszero = 1;
	}
	for(t=0;t<6;t++) {
	    if(z[t*2+1] > 0.00001 || (z[t*2]*(sk-s)) < 0) {
		//printf("(%f+%fi)\n", z[t*2], z[t*2+1]);
		continue;
	    }
	    double x = z[t*2];
	    double foo = A*x*x*x + B*x*x + C*x + D;
	    //printf("%f+%fi, %f\n", z[t*2], z[t*2+1], foo);
	    if(haszero) {
		if(foo!=0.0)
		    continue;
	    } else {
		if(foo*sk < 0)
		    continue;
	    }

	    if(fabs(z[t*2]) < fabs(lambda)) {
		lambda = z[t*2];
	    }
	}
	if(fabs(lambda) > 1e31) {
	    lambda = 0;
	}
    }
    gsl_poly_complex_workspace_free(w8);
    //printf("lambda=%f\n", lambda);
    
    for(t=0;t<l;t++) {
	double v = msrc->data[t] - src->mean;
	v = v + lambda*(v*v - sd*sd - sd*s*v);
	// we don't adjust variance or mean, relying on the later explicit routines
	// to do this for us
	msrc->data[t] = v;
    }

    //printf("%f %f\n", v1, v2);
}


void matrix_adjust_kurtosis(matrix_t*msrc, double k)
{
    int l = msrc->width*msrc->height;
    /* adjust kurtosis */
    statistics_t* src = statistics_new_frommatrix(msrc);

    if(src->var < 1e-4) {
        printf("not adjusting kurtosis of matrix with variance %f\n", src->v);
        return;
    }
    if(!isfinite(k)) {
        printf("not adjusting to kurtosis %f\n", k);
        return;
    }

    double k0 = src->kurt;
    double a = src->v[4] / src->v[2];
    double*m = src->v;

    double A = m[12]-4*a*m[10]-4*m[3]*m[9]+6*a*a*m[8]+12*a*m[3]*m[7]+6*m[3]*m[3]*m[6]-
	       4*a*a*a*m[6]-12*a*a*m[3]*m[5]+a*a*a*a*m[4]-12*a*m[3]*m[3]*m[4]+
	       4*a*a*a*m[3]*m[3]+6*a*a*m[3]*m[3]*m[2]-3*m[3]*m[3]*m[3]*m[3];
    double B = 4*(m[10]-3*a*m[8]-3*m[3]*m[7]+3*a*a*m[6]+6*a*m[3]*m[5]+3*m[3]*m[3]*m[4]-
	       a*a*a*m[4]-3*a*a*m[3]*m[3]-3*m[4]*m[3]*m[3]);
    double C = 6*(m[8]-2*a*m[6]-2*m[3]*m[5]+a*a*m[4]+2*a*m[3]*m[3]+m[3]*m[3]*m[2]);
    double D = 4*(m[6]-a*a*m[2]-m[3]*m[3]);
    double E = m[4];
    double F = D/4;
    double G = m[2];

    double D1 = B*F;
    double D2 = 2*C*F - 4*A*G;
    double D3 = 4*F*D - 3*B*G - D*F;
    double D4 = 4*F*E - 2*C*G;
    double D5 = -D*G;

    //printf("A:%f B:%f C:%f D:%f E:%f\n", A,B,C,D,E);

    double a5[5] = {D5,D4,D3,D2,D1};
    gsl_poly_complex_workspace* w5 = gsl_poly_complex_workspace_alloc(5);
    double z5[5*2];
    memset(z5, 0, sizeof(z5));
    int t;
    for(t=0;t<5;t++) {
        if(!isfinite(a5[t])) {
            printf("Image contains NaN or Inf values, not adjusting kurtosis\n");
#ifdef TOUCHY
            exit(1);
#endif
            return;
        }
    }
    gsl_error_handler_t*oldhandler = gsl_set_error_handler(handleerror);
    gsl_error_occured = 0;
    sprintf(currentaction, "kurtosis1 %f %f %f %f %f\n", a5[0], a5[1], a5[2], a5[3], a5[4]);
    gsl_poly_complex_solve(a5, 5, w5, z5);
    gsl_set_error_handler(oldhandler);
    if(gsl_error_occured)
        return;

    double lmi = -1e33;
    double lma = 1e33;

    for(t=0;t<4;t++) {
	//printf("%f+%fi\n", z5[t*2], z5[t*2+1]);
	if(fabs(z5[t*2+1]) > fabs(z5[t*2])*1e-6)
	    continue;
	if(z5[t*2] < 0 && z5[t*2] > lmi)
	    lmi = z5[t*2];
	if(z5[t*2] >= 0 && z5[t*2] < lma)
	    lma = z5[t*2];
    }
    if(lmi < -1e32) {
	lmi = -1e16;
    }
    if(lma > 1e32) {
	lma = 1e16;
    }
    //printf("lmi: %f lma: %f\n", lmi, lma);
    double kmin = (E + lmi*(D + lmi*(C + lmi*(B + lmi*A)))) / sqr(G+F*lmi*lmi);
    double kmax = (E + lma*(D + lma*(C + lma*(B + lma*A)))) / sqr(G+F*lma*lma);
    if(kmin > kmax) {double r = kmin;kmin = kmax;kmax = r;}
    //printf("kmin: %f kmax: %f\n", kmin, kmax);

    double lambda = 0;
    if(k<kmin || k>kmax) {
	if(k<kmin)
	    lambda = lmi;
	else
	    lambda = lma;

	//fprintf(stderr, "Warning: Target kurtosis (%f) not in range [%f %f]\n", k, kmin, kmax);
    } else {
	a5[0] = E - k*G*G;
	a5[1] = D;
	a5[2] = C - 2*k*F*G;
	a5[3] = B;
	a5[4] = A - k*F*F;
        gsl_error_handler_t*oldhandler = gsl_set_error_handler(handleerror);
        gsl_error_occured = 0;
        sprintf(currentaction, "kurtosis2 %f %f %f %f %f\n", a5[0], a5[1], a5[2], a5[3], a5[4]);
	gsl_poly_complex_solve(a5, 5, w5, z5);
        gsl_set_error_handler(oldhandler);
        if(gsl_error_occured)
            return;
       
	lambda = 1e32;
	for(t=0;t<4;t++) {
	    if(z5[t*2+1] > 0.00001)
		continue;
	    if(fabs(z5[t*2]) < fabs(lambda)) {
		//printf("%f+%fi\n", z5[t*2], z5[t*2+1]);
		lambda = z5[t*2];
	    }
	}
	if(fabs(lambda) > 1e31) {
	    lambda = 0;
	}
    }
    gsl_poly_complex_workspace_free(w5);
    //printf("lambda=%f\n", lambda);
    
    for(t=0;t<l;t++) {
	double v = msrc->data[t] - src->mean;
	v = v + lambda*(v*v*v - a*v - m[3]);
	// we don't adjust variance or mean, relying on the later explicit routines
	// to do this for us
	msrc->data[t] = v;
    }

    //printf("%f %f\n", v1, v2);
}

void matrix_adjust_variance(matrix_t*msrc, double var)
{
    // adjust variance
    int l = msrc->width*msrc->height;

    int t;
    double msrc_mean = 0;
    for(t=0;t<l;t++) {
        msrc_mean += msrc->data[t];
    }
    msrc_mean /= l;
    double msrc_var = 0;
    for(t=0;t<l;t++) {
        double v = msrc->data[t] - msrc_mean;
        msrc_var += v*v;
    }
    msrc_var /= l;

    double mul = sqrt(var / msrc_var);

    if(msrc_mean > 1e-4) {
        for(t=0;t<l;t++) {
            double v = msrc->data[t];
            v = mul*v;
            msrc->data[t] = v;
        }
    } else {
        for(t=0;t<l;t++) {
            double v = msrc->data[t] - msrc_mean;
            v = mul*v;
            msrc->data[t] = v + msrc_mean;
        }
    }
}
        
void matrix_adjust_mean(matrix_t*msrc, double mean)
{
    // adjust mean
    int l = msrc->width*msrc->height;
    int t;
    double msrc_mean = 0;
    for(t=0;t<l;t++) {
        msrc_mean += msrc->data[t];
    }
    msrc_mean /= l;

    double add = mean - msrc_mean;

    for(t=0;t<l;t++) {
	msrc->data[t] += add;
    }
}

void matrix_adjust_statistics(matrix_t*msrc, statistics_t*dest)
{
    int t;
    int l = msrc->width*msrc->height;
    //statistics_print(dest);
    //statistics_print(statistics_new_frommatrix(msrc));
  
    matrix_adjust_skewness(msrc, dest->skew);
    matrix_adjust_kurtosis(msrc, dest->kurt);
    matrix_adjust_variance(msrc, dest->var);
    matrix_adjust_mean(msrc, dest->mean);

    /*printf("after: ");
    statistics_print(statistics_new_frommatrix(msrc));*/
    //statistics_print_comparison(dest, statistics_new_frommatrix(msrc));
}
void matrix_adjust_minmax(matrix_t*m, statistics_t*dest)
{
    int t;
    int l = m->width*m->height;
    for(t=0;t<l;t++) {
        if(m->data[t] < dest->min)
            m->data[t] = dest->min;
        if(m->data[t] > dest->max)
            m->data[t] = dest->max;
    }
}
statistics_t*statistics_merge(statistics_t*s1, statistics_t*s2)
{
    statistics_t*s = statistics_new();
    int t;
    for(t=0;t<16;t++) {
        s->v[t] = (s1->v[t] + s2->v[t])/2;
    }
    s->mean = (s1->mean + s2->mean) / 2;
    s->absmean = (s1->absmean + s2->absmean) / 2;
    s->dev = (s1->dev + s2->dev) / 2;
    s->var = (s1->var + s2->var) / 2;
    s->skew = (s1->skew + s2->skew) / 2;
    s->kurt = (s1->kurt  + s2->kurt) / 2;
    s->min = (s1->min + s2->min)/2;
    s->max = (s1->max + s2->max)/2;
}

bytearrayset_t*matrix_blockstats(matrix_t*m, int w, int h)
{
    matrixset_t*mset = matrixset_new_alloc(6, m->width, m->height);
    int x,y;
    int pos=0;
    int size = m->width*m->height;
    for(y=0;y<m->height;y++) {
        int ly = m->height-y;
        if(ly>h)
            ly=h;
        for(x=0;x<m->width;x++) {
            int xx,yy;
            int lx = m->width-x;
            if(lx>h)
                lx=h;
            double mean=0;
            double*d = &m->data[pos];
            double min=d[0],max=d[0];
            int n = lx*ly;
            /* the following two loops could be optimized by integral images */
            for(yy=0;yy<ly;yy++) {
                for(xx=0;xx<lx;xx++) {
                    mean+=d[xx];
                    if(d[xx]<min) min = d[xx];
                    if(d[xx]>max) max = d[xx];
                }
                d += m->width;
            }
            mean/=n;
            double sum2=0;
            double sum3=0;
            double sum4=0;
            d = &m->data[pos];
            for(yy=0;yy<ly;yy++) {
                for(xx=0;xx<lx;xx++) {
                    double v = d[xx]- mean;
                    sum2+=v*v;
                    sum3+=v*v*v;
                    sum4+=v*v*v*v;
                }
                d+=m->width;
            }
            sum2/=n;
            sum3/=n;
            sum4/=n;
            double dev = sqrt(sum2/n);
            mset->m[0]->data[pos] = mean;
            mset->m[1]->data[pos] = dev;
            mset->m[2]->data[pos] = 0;//sum3/(dev*dev*dev);
            mset->m[3]->data[pos] = sum4/(sum2*sum2);
            mset->m[4]->data[pos] = min;
            mset->m[5]->data[pos] = max;
            pos++;
        }
    }
    bytearrayset_t*set = bytearrayset_new_alloc(mset->num, m->width, m->height);
    int t;
    for(t=0;t<mset->num;t++) {
        double*data = mset->m[t]->data;
        double min = data[0];
        double max = data[0];
        int s;
        for(s=1;s<size;s++) {
            if(data[s]<min) min=data[s];
            if(data[s]>max) max=data[s];
        }
        if(min==max) max=min+1;
        double f = 127/(max - min);
        unsigned char*dest = set->m[t]->data;
        for(s=0;s<size;s++) {
            dest[s] = ((int)((data[s] - min)*f))&127;
        }
        free(mset->m[t]);mset->m[t]=0;
    }
    return set;
}

