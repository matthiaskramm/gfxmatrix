/* gradopt.c 

   Gradient optimization utility.

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

#include <stdlib.h>
#include <stdio.h>
#include <memory.h>

typedef unsigned u32;
typedef unsigned long long int u64;
static int cmpdouble(double*f1, double*f2)
{
    if(*f1<*f2)
        return -1;
    else if(*f1>*f2)
        return 1;
    else 
        return 0;

    if(*(u64*)&f1 < *(u64*)&f2)
        return -1;
    else if(*(u64*)&f1 > *(u64*)&f2)
        return 1;
    else
        return 0;
}

/* 1 / (2 ** (PRECISION-128)) */
#define PRECISION 230
#define sign(x) ((x)<0?-1:((x)>0?1:0))

double* gradient_approximation(int dim, void(*init)(void*, double*), void(*getdir)(void*,double*,double*), double(*fitness)(void*,double*,double*,double), void*data, int maxiterations)
{
    int s,t;
    double*dir = (double*)malloc(dim*sizeof(double));
    double*v = (double*)malloc(dim*sizeof(double));

    int stepsize = 127;
    float step = 1.0;
    double currentdiff;
    int steps = 0;
    double lastdiff = 0;
    char debug = 0;
    int iteration = 0;
  
    init(data, v);
    double maxdiff = currentdiff = fitness(data, v, 0, 0);
    while(1) {
        getdir(data, v, dir);

	/* calculate the ammount we need to move into the gradients direction */
        double d = fitness(data, v, dir, step);
	if(cmpdouble(&d, &currentdiff) <= 0) {
            if(debug) printf("current stepsize (%d/%f) makes things smaller by %f\n", stepsize, step, currentdiff-d);
	    // current stepsize makes things smaller
	    while(stepsize<PRECISION) {
		u32 x = ((stepsize+1) << 23);
		float newstep = *(float*)&x;
                double d2 = fitness(data, v, dir, newstep);
                if(cmpdouble(&d,&d2) <= 0) {
                    /* taking a bigger step may make things smaller, *too*, but not
                       as small as a smaller step will make them */
                    break;
                }
                d = d2;
		if(cmpdouble(&d, &currentdiff) > 0) {
                    if(debug) printf("stepsize of %d/%f would make things bigger by %f\n", stepsize, newstep, d-currentdiff);
                    if(debug) printf("upgraded stepsize to %d\n", stepsize);
		    // we successfully upgraded stepsize to the current possible maximum
		    break;
		}
		stepsize++;step=newstep;
	    }
	    if(stepsize==PRECISION) {
		break;
            }
	} else {
	    // current stepsize makes things bigger
            if(debug) printf("current stepsize (%d) makes things bigger\n", stepsize);
	    while(stepsize>0) {
		stepsize--;
		u32 x = (stepsize << 23);
		step = *(float*)&x;
		double d = fitness(data, v, dir, step);
                //printf("step=%d/%f diffnow=%f diff=%f cmp=%d\n", stepsize, step, d, currentdiff, cmpdouble(&d, &currentdiff));
		if(cmpdouble(&d, &currentdiff) < 0) {
                    if(debug) printf("ok, reducing to stepsize %d makes things smaller again\n", stepsize);
		    // ok, reducing to this stepsize makes things smaller again
		    break;
		}
	    }
            if(!stepsize) {
                if(debug) printf("couldn't find adequate stepsize, even stepping with 0.0 makes things bigger\n");
                int t;
                if(debug) for(t=0;t<dim;t++)
                            printf("%f\n", dir[t]);
                /* The new gradient doesn't make things better, no matter how small
                   the value we multiply with it. 
                   We could now try some bigger values, to protect against a possible integer
                   overflow- or we can just leave it as that, the current point is probably
                   the best maximum we're going to get anyway */
                break;
            }
	}
        /*char ok = 0;*/
	for(t=0;t<dim;t++) {
	    v[t] += dir[t]*step;
            /*if(dir[t]*step > maxstep)
                ok = 1;*/
            //printf("%7.4f  ", v[t]);
	}
	lastdiff = currentdiff;
	currentdiff = fitness(data, v, 0, 0);
        //printf("diff:%f\n", currentdiff);

        if(maxiterations>0 && iteration++>=maxiterations) {
            //printf("breaking gradient approach after %d iterations\n", iteration);
            break;
        }
	if(cmpdouble(&lastdiff, &currentdiff) <= 0) {
	    break;
	}
	steps++;
    }
   
    //printf("%4.2f%% precision (fitness=%f)\n", 100 - (currentdiff * 100 / maxdiff), currentdiff);


    free(dir);
    return v;
}

