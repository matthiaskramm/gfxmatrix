/* texture.c 

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

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <memory.h>
#include "image.h"
#include "filter.h"
#include "filter_steerpyr.h"
#include "corr.h"
#include "stat.h"
#include "texture.h"

typedef struct _filteriterator
{
    int band2start;
    int band2end;
    int band1start;
    int band1end;
    matrixset_t*s;
} filteriterator_t;

filteriterator_t*filteriterator_new(matrixset_t*s)
{
    filteriterator_t*i = malloc(sizeof(filteriterator_t));
    i->band2start=-1;
    i->band2end=s->num;
    i->band1start=s->num-1;
    i->band1end=s->num;
    i->s = s;
}
#define BANDSTART 1
matrixset_t*filteriterator_getband(filteriterator_t*i)
{
    if(i->band1start<BANDSTART)
        return 0;
    while(i->band1start>BANDSTART && i->s->m[i->band1start]->width == i->s->m[i->band1start-1]->width)
        i->band1start--;
    while(i->band2start>BANDSTART && i->s->m[i->band2start]->width == i->s->m[i->band2start-1]->width)
        i->band2start--;

    /*int t;
    for(t=0;t<i->s->num;t++) {
        if(i->band1start == t)
            printf("[");
        if(i->band2start == t)
            printf("{");
        printf("%dx%d ", i->s->m[t]->width, i->s->m[t]->height);
        if(i->band1end == t+1)
            printf("]");
        if(i->band2end == t+1)
            printf("}");
    }
    printf("\n");*/
    matrixset_t*set = malloc(sizeof(matrixset_t));
    set->num = i->band1end-i->band1start;
    set->m = &i->s->m[i->band1start];

    i->band1end=i->band1start;
    i->band2end=i->band2start;
    i->band1start--;
    i->band2start--;
    if(i->band2start<BANDSTART) 
        i->band2start = i->s->num-1;
    return set;
}

textureparameters_t*textureparameters_new_fromimage(image_t*img)
{
    textureparameters_t*p = malloc(sizeof(textureparameters_t));
    memset(p, 0, sizeof(textureparameters_t));

    matrix_t*imgm = image_extractchannel(img, IMAGE_GRAY);
    p->filter = filtertree_new_steerpyr(3);
    matrixset_t*s = p->original = filtertree_apply(p->filter, imgm);
    int j;

    p->imgstat = statistics_new_frommatrix(imgm);
   
    /* get statistics */
    p->stat = malloc(sizeof(statistics_t*)*s->num);
    //printf("initial statistics:\n");
    for(j=0;j<s->num;j++) {
        p->stat[j] = statistics_new_frommatrix(s->m[j]);
        //statistics_print(p->stat[j]);
    }
    
    /* get autocorrelations */
    p->image_autocorr = matrixset_new(s->num);
    p->autox=9;
    p->autoy=9;
    for(j=0;j<s->num;j++) {
        p->image_autocorr->m[j] = matrix_getautocorrelation(s->m[j], p->autox, p->autoy);
    }

    /* get interband and intraband correlations */
    filteriterator_t*i = filteriterator_new(s);
    matrixset_t*innerset=0,*outerset=0;
    p->bands = 0;
    while((innerset = filteriterator_getband(i))) {
        p->bands++;
        outerset = innerset;
    }
    p->intraband_crosscorr = matrixset_new(p->bands);
    p->interband_crosscorr = matrixset_new(p->bands-1);
    i = filteriterator_new(s);
    innerset=0;outerset=0;
    j = 0;
    while((innerset = filteriterator_getband(i))) {
        if(outerset) {
            p->interband_crosscorr->m[j] = matrix_getcrosscorrelations(innerset, outerset,0);
        } else {
            p->interband_crosscorr->m[j] = 0;
        }
        p->intraband_crosscorr->m[j] = matrix_getcrosscorrelations(innerset, innerset,0);
        j++;
        outerset = innerset;
    }


    return p;
}

void textureparameters_adjust_crosscorrelation(textureparameters_t*p, matrixset_t*s, int image)
{
    filteriterator_t*i = filteriterator_new(s);
    matrixset_t*innerset=0,*outerset=0;
    int j = 0;
    while((innerset = filteriterator_getband(i))) {
        if(image==-1 || (image >= (innerset->m - s->m) && image < (innerset->m - s->m + innerset->num))) {
            if(outerset && p->interband_crosscorr->m[j]) {
                matrix_adjustcrosscorrelation(innerset, p->intraband_crosscorr->m[j], outerset, p->interband_crosscorr->m[j]);
            } else {
                matrix_adjustcrosscorrelation(innerset, p->intraband_crosscorr->m[j], 0, 0);
            }
        }
        j++;outerset = innerset;
    }
}
void textureparameters_adjust_statistics(textureparameters_t*p, matrixset_t*s, int image)
{
    if(image==-1) {
        int t;
        for(t=0;t<s->num;t++) {
            matrix_adjust_statistics(s->m[t],p->stat[t]);
        }
    } else {
        matrix_adjust_statistics(s->m[image],p->stat[image]);
    }
}
void textureparameters_adjust_autocorrelation(textureparameters_t*p, matrixset_t*s, int image)
{
    matrix_adjustautocorrelation2(s->m[image], p->image_autocorr->m[image]);
}
void textureparameters_adjust_finalimage(textureparameters_t*p, matrix_t*image)
{
    matrix_adjust_statistics(image, p->imgstat);
    matrix_adjust_minmax(image, p->imgstat);
}
void textureparameters_copy_from_original(textureparameters_t*p, matrixset_t*s, int image)
{
    assert(s->m[image]->width == p->original->m[image]->width);
    assert(s->m[image]->height == p->original->m[image]->height);

    memcpy(s->m[image]->data, p->original->m[image]->data, s->m[image]->width*s->m[image]->height*sizeof(s->m[image]->data[0]));
}

void dump_matrix_to_file(matrix_t*m)
{
    static int counter=1;
    char filename[80];
    sprintf(filename, "adjusted%02d.png", counter);
    image_save(image_from_matrix(m,IMAGE_MINMAX), filename);
    counter++;
}

matrix_t*txt_reverse_nodes(textureparameters_t*p, filtertree_t*tree, matrixset_t*s, int*listpos)
{
    int t;
    matrix_t*m = 0;
    int autocorr = -1;
    if(tree->filter) {
        autocorr = *listpos;
    }
    if(tree->num_children) {
        if(tree->filter) {
            (*listpos)++;
        }
        for(t=0;t<tree->num_children;t++) {
            int l = *listpos;
            matrix_t*g = txt_reverse_nodes(p, tree->children[t], s, listpos);
            if(!m) m = g;
            else matrix_add_inplace(m,g);
        }
    } else {
        m = s->m[(*listpos)];
        (*listpos)++;
    }
    if(tree->filter) {
        m = filter_reverse(tree->filter, m);
        if(tree->filter->flags&FILTER_HIGHBAND) {
            memset(m->data, 0, sizeof(m->data[0])*m->width*m->height);
        }
    }
    if(autocorr>=0) {
        /* adjust autocorr */
        if(tree->filter && !(tree->filter->flags&FILTER_HIGHBAND)) {
            printf("Adjusting %dx%d subband autocorr (%s) (source:%dx%d, pos %d)\n", m->width, m->height, tree->description,
                    p->original->m[autocorr]->width,
                    p->original->m[autocorr]->height, autocorr);
            matrix_adjustautocorrelation3(m, p->image_autocorr->m[autocorr]);
        }
        matrix_adjust_statistics(m, p->stat[autocorr]);
    }
    return m;
}

void textureparameters_iterate(textureparameters_t*p, matrix_t*m)
{
    matrixset_t*s = filtertree_apply(p->filter, m);
    filtertree_t*tree = p->filter;
   
    int listpos = 0;
    matrix_t*new = txt_reverse_nodes(p, tree, s, &listpos);

    matrix_adjust_statistics(new, p->imgstat);
    dump_matrix_to_file(m);
    matrix_adjust_minmax(new, p->imgstat);
    dump_matrix_to_file(m);
    matrix_update(m, new);
    matrix_delete(new);
        
}

