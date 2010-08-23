/* color.c 

   Color utilities.

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
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <assert.h>
#include <dirent.h>
#include "color.h"

typedef RGBA gfxcolor_t;

static void gfximage_transform(image_t*img, cxform_t*cxform)
{
    int t;
    int size = img->width*img->height;

    int rr,rg,rb,ra, tr;
    int gr,gg,gb,ga, tg;
    int br,bg,bb,ba, tb;
    int ar,ag,ab,aa, ta;
    rr = (int)(cxform->rr*256);gr = (int)(cxform->gr*256);
    rg = (int)(cxform->rg*256);gg = (int)(cxform->gg*256);
    rb = (int)(cxform->rb*256);gb = (int)(cxform->gb*256);
    ra = (int)(cxform->ra*256);ga = (int)(cxform->ga*256);
    br = (int)(cxform->br*256);ar = (int)(cxform->ar*256);tr = (int)(cxform->tr*256);
    bg = (int)(cxform->bg*256);ag = (int)(cxform->ag*256);tg = (int)(cxform->tg*256);
    bb = (int)(cxform->bb*256);ab = (int)(cxform->ab*256);tb = (int)(cxform->tb*256);
    ba = (int)(cxform->ba*256);aa = (int)(cxform->aa*256);ta = (int)(cxform->ta*256);

    for(t=0;t<size;t++) {
        RGBA*pixel = &img->data[t];
        U8 r = (pixel->r * rr + pixel->g * rg + pixel->b * rb + pixel->a * ra + tr) / 256;
        U8 g = (pixel->r * gr + pixel->g * gg + pixel->b * gb + pixel->a * ga + tg) / 256;
        U8 b = (pixel->r * br + pixel->g * bg + pixel->b * bb + pixel->a * ba + tb) / 256;
        U8 a = (pixel->r * ar + pixel->g * ag + pixel->b * ab + pixel->a * aa + ta) / 256;
        pixel->r = r;
        pixel->g = g;
        pixel->b = b;
        pixel->a = a;
    }
}

double clamp256(double x)
{
    if(x<0) {printf("value %f clamped to 0\n",x);x=0;}
    if(x>=256) {printf("value %f clamped to 255\n",x);x=255;}
    return x;
}

void gfxcolor_transform(gfxcolor_t*col, cxform_t*cx)
{
    U8 r = (U8)clamp256(col->r * cx->rr + col->g * cx->rg + col->b * cx->rb + col->a * cx->ra + cx->tr);
    U8 g = (U8)clamp256(col->r * cx->gr + col->g * cx->gg + col->b * cx->gb + col->a * cx->ga + cx->tg);
    U8 b = (U8)clamp256(col->r * cx->br + col->g * cx->bg + col->b * cx->bb + col->a * cx->ba + cx->tb);
    U8 a = (U8)clamp256(col->r * cx->ar + col->g * cx->ag + col->b * cx->ab + col->a * cx->aa + cx->ta);
    col->r = r;
    col->g = g;
    col->b = b;
    col->a = a;
}

static void gfxcxform_show(cxform_t*cx)
{
    printf("%5.2f %5.2f %5.2f %5.2f  %5.2f\n", cx->rr, cx->rg, cx->rb, cx->ra, cx->tr);
    printf("%5.2f %5.2f %5.2f %5.2f  %5.2f\n", cx->gr, cx->gg, cx->gb, cx->ga, cx->tg);
    printf("%5.2f %5.2f %5.2f %5.2f  %5.2f\n", cx->br, cx->bg, cx->bb, cx->ba, cx->tb);
    printf("%5.2f %5.2f %5.2f %5.2f  %5.2f\n", cx->ar, cx->ag, cx->ab, cx->aa, cx->ta);
}

cxform_t* seperateColors(image_t*img)
{
    int t;
    int size = img->width*img->height;
    double corr[4][4];
    double e[4][4];
    double i[4][4];
    memset(&corr, 0, sizeof(corr));
    cxform_t*cxform = (cxform_t*)malloc(sizeof(cxform_t));
    cxform_t inverse;

    for(t=0;t<size;t++) {
        U8*pixel = (U8*)&img->data[t];
        int u,v;
        for(u=0;u<4;u++)
        for(v=0;v<4;v++) {
            corr[u][v] += pixel[u]*pixel[v];
        }
    }
    
    gsl_matrix *v= gsl_matrix_alloc(4, 4);
    gsl_vector *s = gsl_vector_alloc(4);
    gsl_vector *tmp= gsl_vector_alloc(4);
    gsl_matrix_view m = gsl_matrix_view_array((double*)&corr, 4, 4);
    gsl_linalg_SV_decomp(&m.matrix, v, s, tmp);
    int x,y,z;
    double min=500000,max=-500000;
    double rad[4];
    for(y=0;y<4;y++) {
        double val = gsl_vector_get(s, y);
        double min=0,max=0;
        for(x=0;x<4;x++) {
            i[x][y] = gsl_matrix_get(v, x, y);
            e[y][x] = gsl_matrix_get(&m.matrix, x, y);
            if(i[x][y]<0) {
                min += i[x][y];
            } else {
                max += i[x][y];
            }
        }
        rad[y] = ((-min)>max?-min:max) * 2.01;
    }
    /*double tst[4][4];
    for(y=0;y<4;y++) {
        for(x=0;x<4;x++) {
            double sum = 0;
            for(z=0;z<4;z++) {
                sum += e[x][z] * i[z][y] * 256;
            }
            printf("%8.4f  ", sum);
        }
        printf("\n");
    }*/

    gsl_matrix_free(v);
    gsl_vector_free(s);
    gsl_vector_free(tmp);

    cxform_t tr;
    tr.rr = i[0][0]/rad[0]; tr.rg = i[1][0]/rad[0]; tr.rb = i[2][0]/rad[0]; tr.ra = i[3][0]/rad[0]; tr.tr = 128;
    tr.gr = i[0][1]/rad[1]; tr.gg = i[1][1]/rad[1]; tr.gb = i[2][1]/rad[1]; tr.ga = i[3][1]/rad[1]; tr.tg = 128;
    tr.br = i[0][2]/rad[2]; tr.bg = i[1][2]/rad[2]; tr.bb = i[2][2]/rad[2]; tr.ba = i[3][2]/rad[2]; tr.tb = 128;
    tr.ar = i[0][3]/rad[3]; tr.ag = i[1][3]/rad[3]; tr.ab = i[2][3]/rad[3]; tr.aa = i[3][3]/rad[3]; tr.ta = 128;

    cxform->rr = e[0][0]*rad[0]; cxform->rg = e[1][0]*rad[1]; cxform->rb = e[2][0]*rad[2]; cxform->ra = e[3][0]*rad[3];
    cxform->gr = e[0][1]*rad[0]; cxform->gg = e[1][1]*rad[1]; cxform->gb = e[2][1]*rad[2]; cxform->ga = e[3][1]*rad[3];
    cxform->br = e[0][2]*rad[0]; cxform->bg = e[1][2]*rad[1]; cxform->bb = e[2][2]*rad[2]; cxform->ba = e[3][2]*rad[3];
    cxform->ar = e[0][3]*rad[0]; cxform->ag = e[1][3]*rad[1]; cxform->ab = e[2][3]*rad[2]; cxform->aa = e[3][3]*rad[3];
    cxform->tr = - (cxform->rr * tr.tr + cxform->rg * tr.tg + cxform->rb * tr.tb + cxform->ra * tr.ta);
    cxform->tg = - (cxform->gr * tr.tr + cxform->gg * tr.tg + cxform->gb * tr.tb + cxform->ga * tr.ta);
    cxform->tb = - (cxform->br * tr.tr + cxform->bg * tr.tg + cxform->bb * tr.tb + cxform->ba * tr.ta);
    cxform->ta = - (cxform->ar * tr.tr + cxform->ag * tr.tg + cxform->ab * tr.tb + cxform->aa * tr.ta);
    
    gfximage_transform(img, &tr);
    //gfximage_transform(img, cxform);

    /*gfxcolor_t red;
    //red.r = 255; red.g = red.b = red.a = 0;
    red.r = lrand48();
    red.g = lrand48();
    red.b = lrand48();
    red.a = lrand48();
    printf("%02x%02x%02x%02x\n", red.r,red.g,red.b,red.a);
    gfxcolor_transform(&red, &tr);
    printf("%02x%02x%02x%02x\n", red.r,red.g,red.b,red.a);
    gfxcolor_transform(&red, cxform);
    printf("%02x%02x%02x%02x\n", red.r,red.g,red.b,red.a);

    RGBA cmin = {255,255,255,255};
    RGBA cmax = {0,0,0,0};
    for(t=0;t<size;t++) {
        if(img->data[t].r < cmin.r) cmin.r = img->data[t].r;
        if(img->data[t].g < cmin.g) cmin.g = img->data[t].g;
        if(img->data[t].b < cmin.b) cmin.b = img->data[t].b;
        if(img->data[t].a < cmin.a) cmin.a = img->data[t].a;
        if(img->data[t].r > cmax.r) cmax.r = img->data[t].r;
        if(img->data[t].g > cmax.g) cmax.g = img->data[t].g;
        if(img->data[t].b > cmax.b) cmax.b = img->data[t].b;
        if(img->data[t].a > cmax.a) cmax.a = img->data[t].a;
    }*/

    return cxform;
}
