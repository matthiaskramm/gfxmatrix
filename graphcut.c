/* graphcut.c 

   2D Graphcut utility.

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
#include <math.h>
#include <memory.h>
#include "graphcut.h"

//#define DEBUG

//#define CHECKS
//#define PREPROCESS

#ifdef DEBUG
#define DBG
#else
#define DBG if(0)
#endif

#define ACTIVE 0x10
#define IN_TREE 0x20

#define RIGHT_FLAG 1
#define UP_FLAG 2
#define LEFT_FLAG 4
#define DOWN_FLAG 8

#define TWOTREES

/* [REVERSE-DIR(2)] [IN_TREE] [ACTIVE] [UP][DOWN][LEFT][RIGHT] */

static int inv[4] = {2, 3, 0, 1};
static int bit[4] = {1, 2, 4, 8};
static int xadd[4] = {1,0,-1,0};
static int yadd[4] = {0,-1,0,1};
static char*dirname[] = {"right","up","left","down","stay"};

typedef struct _posqueue_entry {
    int pos;
    struct _posqueue_entry*next;
} posqueue_entry_t;

typedef struct _posqueue {
    posqueue_entry_t*list;
} posqueue_t;

typedef struct _graphcut_workspace {
    unsigned char*flags1;
    unsigned char*flags2;
    map_t*map;
    int width;
    int height;
    int size;
    int add[4];
    int pos1;
    int pos2;
    posqueue_t*queue1;
    posqueue_t*queue2;
    posqueue_t*tmpqueue;
} graphcut_workspace_t;

static posqueue_t*posqueue_new() 
{
    posqueue_t*m = (posqueue_t*)malloc(sizeof(posqueue_t));
    memset(m, 0, sizeof(posqueue_t));
    return m;
}
static void posqueue_delete(posqueue_t*q)
{
    free(q);
}
static inline void posqueue_addpos(posqueue_t*queue, int pos)
{
    posqueue_entry_t*old = queue->list;
    queue->list = malloc(sizeof(posqueue_entry_t));
    queue->list->pos = pos;
    queue->list->next = old;
}
static inline int posqueue_extract(posqueue_t*queue)
{
    posqueue_entry_t*item = queue->list;
    int pos;
    if(!item)
	return -1;
    pos = item->pos;
    queue->list = queue->list->next;
    free(item);
    return pos;
}
static inline int posqueue_notempty(posqueue_t*queue)
{
    return (int)queue->list;
}
static void posqueue_print(posqueue_t*queue)
{
    posqueue_entry_t*e = queue->list;
    while(e) {
	printf("%d ", e->pos);
	e = e->next;
    }
    printf("\n");
}
static void posqueue_purge(posqueue_t*queue)
{
    posqueue_entry_t*e = queue->list;
    while(e) {
	posqueue_entry_t*next = e->next;
	e->next = 0;free(e);
	e = next;
    }
    queue->list = 0;
}

map_t* map_new(int width, int height)
{
    map_t*map = malloc(sizeof(map_t));
    map->width = width;
    map->height = height;
    map->weight[0] = malloc(sizeof(map->weight[0][0])*width*height);
    map->weight[1] = malloc(sizeof(map->weight[1][0])*width*height);
    map->weight[2] = malloc(sizeof(map->weight[2][0])*width*height);
    map->weight[3] = malloc(sizeof(map->weight[3][0])*width*height);
    memset(map->weight[0], 0, sizeof(map->weight[0][0])*width*height);
    memset(map->weight[1], 0, sizeof(map->weight[1][0])*width*height);
    memset(map->weight[2], 0, sizeof(map->weight[2][0])*width*height);
    memset(map->weight[3], 0, sizeof(map->weight[3][0])*width*height);
    return map;
}

void map_delete(map_t*map)
{
    free(map->weight[0]);map->weight[0] = 0;
    free(map->weight[1]);map->weight[1] = 0;
    free(map->weight[2]);map->weight[2] = 0;
    free(map->weight[3]);map->weight[3] = 0;
    free(map);
}

map_t* map_from_matrix(matrix_t*m)
{
    map_t*map = map_new(m->width, m->height);
    int x,y;
    for(y=0;y<m->height;y++) {
        double*prev = y?&m->data[(y-1)*m->width]:&m->data[y*m->width];
        double*line = &m->data[y*m->width];
        double*next= (y<m->height-1)?&m->data[(y+1)*m->width]:&m->data[y*m->width];
        int prevx = 0;
        int nextx = 1;
        int pos = y*map->width;
        for(x=0;x<m->width;x++) {
            map->weight[RIGHT][pos] = (weight_t)fabs(line[x] - line[nextx]); //right
            map->weight[UP   ][pos] = (weight_t)fabs(prev[x] - line[x]); //up
            map->weight[LEFT ][pos] = (weight_t)fabs(line[x] - line[prevx]); //left
            map->weight[DOWN ][pos] = (weight_t)fabs(next[x] - line[x]); //down
            pos++;
            prevx = x;
            if(x<m->width-2) {
                nextx++;
            }
        }
    }
    inverse_map(map);

    return map;
}

int coldiff(RGBA a, RGBA b)
{
    int dr = a.r - b.r;
    int dg = a.g - b.g;
    int db = a.b - b.b;
    return (int)(1+sqrt(dr*dr + dg*dg + db*db));
}

map_t* map_from_image(image_t*img)
{
    map_t*map = map_new(img->width, img->height);
    int x,y;
    for(y=0;y<img->height;y++) {
        RGBA*prev = y?&img->data[(y-1)*img->width]:&img->data[y*img->width];
        RGBA*line = &img->data[y*img->width];
        RGBA*next= (y<img->height-1)?&img->data[(y+1)*img->width]:&img->data[y*img->width];
        int prevx = 0;
        int nextx = 1;
        int pos = y*map->width;
        for(x=0;x<img->width;x++) {
            map->weight[RIGHT][pos] = (weight_t)coldiff(line[x],line[nextx]); //right
            map->weight[UP   ][pos] = (weight_t)coldiff(prev[x],line[x]); //up
            map->weight[LEFT ][pos] = (weight_t)coldiff(line[x],line[prevx]); //left
            map->weight[DOWN ][pos] = (weight_t)coldiff(next[x],line[x]); //down
            pos++;
            prevx = x;
            if(x<img->width-2) {
                nextx++;
            }
        }
    }
    inverse_map(map);

    return map;
}

static void make_map_borders(map_t*map)
{
    int x,y;
    for(x=0;x<map->width;x++) {
	map->weight[DOWN][(map->height-1)*map->width+x] = 0;
	map->weight[UP][0*map->width+x] = 0;
    }
    for(y=0;y<map->width;y++) {
	map->weight[RIGHT][map->width-1+y*map->width] = 0;
	map->weight[LEFT][y*map->width] = 0;
    }
}
void inverse_map(map_t*map)
{
    int t,size=map->width*map->height;
    int d=0;
    weight_t max = 0;
    weight_t min = MAX_WEIGHT;
    for(d=0;d<4;d++) {
        for(t=0;t<size;t++) {
            if(map->weight[d][t] < min)
                min = map->weight[d][t];
            if(map->weight[d][t] > max)
                max = map->weight[d][t];
        }
    }
    for(d=0;d<4;d++) {
        for(t=0;t<size;t++) {
            U32 v = map->weight[d][t];
            //v = max-(v-min);
            //v = 15 * ((v-min))/(max-min) + 1;
            v = ((max-min) - (v-min));// + 1;
            v = (U32)((double)v*(double)v*(double)v / 280);
            map->weight[d][t] = v;
        }
    }
    make_map_borders(map);
}

void map_save(map_t*map, char*filename)
{
    FILE*fi = fopen(filename, "wb");

    fprintf(fi, "width: %d\n", map->width);
    fprintf(fi, "height: %d\n", map->width);
    fprintf(fi, "[left]\n");
    int x,y;
    for(y=0;y<map->height;y++) {
        for(x=0;x<map->width;x++) {
            fprintf(fi, "%d ", map->weight[LEFT][y*map->width+x]);
        }
        fprintf(fi, "\n");
    }
    fprintf(fi, "[up]\n");
    for(y=0;y<map->height;y++) {
        for(x=0;x<map->width;x++) {
            fprintf(fi, "%d ", map->weight[UP][y*map->width+x]);
        }
        fprintf(fi, "\n");
    }
    fprintf(fi, "[right]\n");
    for(y=0;y<map->height;y++) {
        for(x=0;x<map->width;x++) {
            fprintf(fi, "%d ", map->weight[RIGHT][y*map->width+x]);
        }
        fprintf(fi, "\n");
    }
    fprintf(fi, "[down]\n");
    for(y=0;y<map->height;y++) {
        for(x=0;x<map->width;x++) {
            fprintf(fi, "%d ", map->weight[DOWN][y*map->width+x]);
        }
        fprintf(fi, "\n");
    }

    fclose(fi);
}

static graphcut_workspace_t*graphcut_workspace_new(map_t*map, int width, int height, int pos1, int pos2)
{
    graphcut_workspace_t*workspace = malloc(sizeof(graphcut_workspace_t));
    workspace->width = width;
    workspace->height = height;
    workspace->flags1 = malloc(width*height);
    memset(workspace->flags1, 0, width*height);
    workspace->flags2 = malloc(width*height);
    memset(workspace->flags2, 0, width*height);
    workspace->size = width*height;
    workspace->add[0] = 1;
    workspace->add[1] = -width;
    workspace->add[2] = -1;
    workspace->add[3] = width;
    workspace->pos1 = pos1;
    workspace->pos2 = pos2;
    workspace->map = map;
    workspace->queue1 = posqueue_new();
    workspace->queue2 = posqueue_new();
    workspace->tmpqueue = posqueue_new();
    return workspace;
}
static void graphcut_workspace_delete(graphcut_workspace_t*w) 
{
    posqueue_delete(w->queue1);w->queue1=0;
    posqueue_delete(w->queue2);w->queue2=0;
    if(w->flags1) free(w->flags1);w->flags1=0;
    if(w->flags2) free(w->flags2);w->flags2=0;
    free(w);
}


typedef struct _path {
    int*pos;
    unsigned char*dir;
    unsigned char*firsthalf;
    int length;
} path_t;

static path_t*path_new(int len)
{
    path_t*p = malloc(sizeof(path_t));
    p->pos = malloc(sizeof(int)*len);
    p->dir = malloc(sizeof(unsigned char)*len);
    p->firsthalf = malloc(sizeof(unsigned char)*len);
    p->length = len;
    return p;
}
static void path_delete(path_t*path)
{
    free(path->pos);path->pos = 0;
    free(path->dir);path->dir = 0;
    free(path->firsthalf);path->firsthalf = 0;
    free(path);
}

static path_t*extract_path(graphcut_workspace_t*work, unsigned char*mytree, unsigned char*othertree, int pos, int newpos, int dir)
{
    int t;
    int p = pos;
    int len1 = 0;
    /* walk up tree1 */
    while(p != work->pos1) {
	p += work->add[mytree[p]>>6];
	len1++;
    }
    p = newpos;
    int len2 = 0;
    /* walk up tree2 */
    while(p != work->pos2) {
	p += work->add[othertree[p]>>6];
	len2++;
    }
    path_t*path = path_new(len1+len2+2);

    t = len1;
    path->pos[t] = p = pos;
    path->dir[t] = dir;
    path->firsthalf[t] = 1;
    while(p != work->pos1) {
	assert(mytree[p]&IN_TREE);
	char dir = mytree[p]>>6;
	p += work->add[dir];
	t--;
	path->pos[t] = p;
	path->dir[t] = inv[dir];
	path->firsthalf[t] = 1;
    }
    assert(!t);

    t = len1+1;

    p = newpos;
    while(p != work->pos2) {
	assert(othertree[p]&IN_TREE);
	path->pos[t] = p;
	path->dir[t] = othertree[p]>>6;
	path->firsthalf[t] = 0;

	p += work->add[othertree[p]>>6];
	t++;
    }

    /* terminator */
    path->pos[t] = p;
    path->dir[t] = 4; // last node
    path->firsthalf[t] = 0;

    assert(t == len1+len2+1);
    return path;
}

static void map_print(graphcut_workspace_t*w)
{
    int y;
    map_t*map = w->map;
    for(y=0;y<w->height;y++) {
	int c=0;
	for(c=0;c<4;c++) {
	    int x;
	    for(x=0;x<w->width;x++) {
		int p = y*w->width+x;
		unsigned char f = w->flags1[p]|w->flags2[p];
		char cc1='-', cc2='.', cc3='.', cc4='.',cc5='.', cc6='-';
		char mark=' ';
		if(f) {
		    if(f&LEFT_FLAG) cc1='<';
		    if(f&RIGHT_FLAG) cc6='>';
		    if(f&UP_FLAG) cc2='^';
		    if(f&DOWN_FLAG) cc5='v';
		    char*cc = " atA";
		    char*dd = "ruld";
		    cc3=cc[((f>>4)&3)];
		    if(w->flags2[p]) {
			mark = '*';
		    }
		    if(w->flags1[p] && w->flags2[p]) {
			cc3 = '?';
			mark = '?';
		    }
		    cc4=dd[f>>6];
		}

		if(c==0) printf("  %c %2d%c  ", mark, map->weight[UP][p], cc2);
		if(c==1) printf(" %2d%c%c%c%c%-2d", map->weight[LEFT][p], cc1, cc3, cc4, cc6, map->weight[RIGHT][p]);
		if(c==2) printf("    %c%-2d  ", cc5, map->weight[DOWN][p]);
		if(c==3) printf("         ");
	    }
	    printf(" |\n");
	}
    }
}

static void path_print(path_t*path, int width)
{
    int t;
    for(t=0;t<path->length;t++) {
	int p = path->pos[t];
	printf("%d/%d(%d)", p%width, p/width, path->firsthalf[t]);
	if(t<path->length-1)
	    printf(" -%s-> ", dirname[path->dir[t]]);
    }
    printf("\n");
}


static void workspace_print(graphcut_workspace_t*w)
{
    printf("queue1: ");posqueue_print(w->queue1);
    printf("queue2: ");posqueue_print(w->queue2);
    map_print(w);
}

static void myassert(graphcut_workspace_t*w, char assertion, const char*file, int line, const char*func)
{
    if(!assertion) {
	printf("Assertion %s:%d (%s) failed:\n", file, line, func);
	workspace_print(w);
	exit(0);
    }
}

#define ASSERT(w,c) {myassert(w,c,__FILE__,__LINE__,__func__);}

static path_t* expand_pos(graphcut_workspace_t*w, posqueue_t*queue, int pos, char reverse, unsigned char*mytree, unsigned char*othertree)
{
    map_t*map = w->map;
    int dir;
    if(!(mytree[pos]&IN_TREE) || !(mytree[pos]&ACTIVE)) {
	/* this node got deleted or marked inactive in the meantime. ignore it */
	return 0;
    }

    for(dir=0;dir<4;dir++) {
	int newpos = pos + w->add[dir];
	if(newpos<0 || newpos>=w->size || mytree[newpos])
	    continue;
	int weight=0;
	if(!reverse) {
	    weight = map->weight[dir][pos];
	} else {
	    weight = map->weight[inv[dir]][newpos];
	}
	if(weight) {
	    if(othertree[newpos]) {
		
                posqueue_addpos(queue, pos); mytree[pos] |= ACTIVE; // re-add, this vertice might have other connections

		path_t*path;
		if(reverse) {
		    path = extract_path(w, othertree, mytree, newpos, pos, inv[dir]);
		} else {
		    path = extract_path(w, mytree, othertree, pos, newpos, dir);
		}
		return path;
	    } else {
		DBG printf("advance into direction %d to new pos %d\n", dir, newpos);
		mytree[pos] |= bit[dir];
		mytree[newpos] = (inv[dir]<<6)|ACTIVE|IN_TREE;
		posqueue_addpos(queue, newpos);

#if 0
		posqueue_addpos(queue, pos); mytree[pos] |= ACTIVE; // re-add, this vertice might have other connections
		return 0;
#endif
	    }
	}
    }
    /* if we couldn't expand this node anymore, it's now an inactive node */
    mytree[pos] &= ~ACTIVE;
    return 0;
}

static weight_t decrease_weights(map_t*map, path_t*path)
{
    int t;
    int min = 0;

#ifdef CHECKS
    for(t=0;t<path->length-1;t++) {
        assert(path->pos[t]<map->width*map->height);
        int x1 = path->pos[t] % map->width;
        int y1 = path->pos[t] / map->width;
        int x2 = path->pos[t+1] % map->width;
        int y2 = path->pos[t+1] / map->width;
        assert(path->dir[t]<4);
        assert(xadd[path->dir[t]] == x2-x1);
        assert(yadd[path->dir[t]] == y2-y1);
    }
    assert(path->dir[path->length-1] == 4);
#endif

    for(t=0;t<path->length-1;t++) {
	int w = map->weight[path->dir[t]][path->pos[t]];
	if(t == 0 || w < min)
	    min = w;
    }
    if(min<=0) 
        return;

    for(t=0;t<path->length-1;t++) {
	map->weight[path->dir[t]][path->pos[t]] -= min;
	map->weight[inv[path->dir[t]]][path->pos[t+1]] += min;
    }
    return min;
}

static void bool_op(graphcut_workspace_t*w, unsigned char*flags, int pos, unsigned char and, unsigned char or)
{
    posqueue_t*q = w->tmpqueue;
    posqueue_purge(q);
    posqueue_addpos(q, pos);

    while(posqueue_notempty(q)) {
	int p = posqueue_extract(q);
	flags[p] = (flags[p]&and)|or;
	int dir;
	for(dir=0;dir<4;dir++) {
	    if(flags[p] & bit[dir]) {
		posqueue_addpos(q, p+w->add[dir]);
	    }
	}
    }
}

static int reconnect(graphcut_workspace_t*w, unsigned char*flags, int pos, char reverse)
{
    map_t*map = w->map;
    int d;
    for(d=0;d<4;d++) {
	int add = w->add[d];
        int newpos = pos+add;
        if(newpos<0 || newpos>=w->size)
            continue;
	int weight;
	if(!reverse) {
	    weight = map->weight[inv[d]][newpos];
	} else {
	    weight = map->weight[d][pos];
	}
	if(weight && (flags[pos+add]&IN_TREE)) {
	    DBG printf("successfully reconnected node %d/%d to %d/%d\n", pos%w->width, pos/w->width, (pos+add)%w->width, (pos+add)/w->width);
	    flags[pos] = (flags[pos]&0x3f) | d<<6; // reverse
	    flags[pos + add] |= bit[inv[d]]; // forward
	    return 1;
	}
    }
    return 0;
}

//#define SIMPLE

static void destroy_subtree(graphcut_workspace_t*w, unsigned char*flags, int pos, posqueue_t*posqueue, int startpos)
{
    DBG printf("destroying subtree starting with %d/%d\n", pos%w->width, pos/w->width);
   
#ifdef SIMPLE
    posqueue_purge(q);
    posqueue_addpos(q, startpos);
    memset(flags, 0, w->width*w->height);
    return;
#else 
    posqueue_t*q = w->tmpqueue;
    posqueue_purge(q);
    posqueue_addpos(q, pos);

    while(posqueue_notempty(q)) {
	int p = posqueue_extract(q);
	int dir;
	for(dir=0;dir<4;dir++) {
	    int newpos = p+w->add[dir];
	    if(newpos<0 || newpos>w->size)
		continue;
	    if(flags[p] & bit[dir]) {
		posqueue_addpos(q, newpos);
	    } else if((flags[newpos]&(ACTIVE|IN_TREE)) == IN_TREE) {
		// TODO: we should check the weight of that edge right now, right here.
		// if it's zero, we don't need to do anything.
		posqueue_addpos(posqueue, newpos);
		flags[newpos]|=ACTIVE;
	    }
	}
	DBG printf("removed pos %d/%d\n", p%w->width, p/w->width);
	flags[p] = 0;
    }
#endif
}

static void combust_tree(graphcut_workspace_t*w, posqueue_t*q1, posqueue_t*q2, path_t*path)
{
    map_t*map = w->map;
    int t;
    for(t=0;t<path->length-1 && path->firsthalf[t+1];t++) {
	int pos = path->pos[t];
	int dir = path->dir[t];
	int newpos = pos+w->add[dir];
	if(!map->weight[dir][pos]) {
	    /* disconnect node */
	    DBG printf("remove link %d/%d->%d/%d from tree 1\n", pos%w->width,pos/w->width, newpos%w->width,newpos/w->width);

	    w->flags1[pos] &= ~(bit[dir]);
	    w->flags1[newpos] &= 0x0f|ACTIVE;
	    bool_op(w, w->flags1, newpos, ~IN_TREE, 0);

	    /* try to reconnect the path to some other tree part */
	    if(reconnect(w, w->flags1, newpos, 0)) {
		bool_op(w, w->flags1, newpos, ~0, IN_TREE);
	    } else {
		destroy_subtree(w, w->flags1, newpos, q1, w->pos1);
		break;
	    }
	}
    }

    for(t=path->length-1;t>0 && !path->firsthalf[t-1];t--) {
	int pos = path->pos[t];
	int newpos = path->pos[t-1];
	int dir = inv[path->dir[t-1]];
	int newpos2 = pos+w->add[dir];
	assert(newpos == newpos2);
	if(!map->weight[inv[dir]][newpos]) {
	    /* disconnect node */
	    DBG printf("remove link %d/%d->%d/%d from tree 2\n", pos%w->width,pos/w->width, newpos%w->width,newpos/w->width);
	    w->flags2[pos] &= ~(bit[dir]);
	    w->flags2[newpos] &= 0x0f|ACTIVE;
	    bool_op(w, w->flags2, newpos, ~IN_TREE, 0);

	    /* try to reconnect the path to some other tree part */
	    if(reconnect(w, w->flags2, newpos, 1)) {
		bool_op(w, w->flags2, newpos, ~0, IN_TREE);
	    } else {
		destroy_subtree(w, w->flags2, newpos, q2, w->pos2);
		break;
	    }
	}
    }
}

static void check_map(map_t*map)
{
    int x,y,d;
    int val = 0;
    for(x=0;x<map->width;x++) {
	if(map->weight[DOWN][(map->height-1)*map->width+x])
            val |= DOWN_FLAG;
	if(map->weight[UP][0*map->width+x])
            val |= UP_FLAG;
    }
    for(y=0;y<map->height;y++) {
	if(map->weight[RIGHT][map->width-1+y*map->width])
            val |= RIGHT_FLAG;
	if(map->weight[LEFT][y*map->width])
            val |= LEFT_FLAG;
    }
    if(val) {
	fprintf(stderr, "Warning: The map hasn't proper edges (positive weights extend beyond the bitmap, %02x)\n", val);
	fprintf(stderr, "I'm fixing this now, but you should be more careful.\n");
    }
    val = 0;
    for(x=0;x<map->width;x++) 
    for(y=0;y<map->height;y++)
    for(d=0;d<4;d++) {
        if(map->weight[d]<=0) 
            val = 1;
    }
    if(val) {
	fprintf(stderr, "Warning: The map has non-positive values\n");
    }
    make_map_borders(map);
}

static int check3path(graphcut_workspace_t*w, int x1, int y1, int x2, int y2, int x3, int y3, int xfirst1, int xfirst2)
{
    map_t*map = w->map;
    int l = abs(x1-x2) + abs(y1-y2) + abs(x2-x3) + abs(y2-y3) + 1;
    path_t* path = path_new(l);

    int x = x1;
    int y = y1;
    int pos = y1*map->width+x1;
    int i = 0;
    int dir;

    dir = xfirst1?(x1<x2?RIGHT:LEFT):(y1<y2?DOWN:UP);
    if(xfirst1 && (x!=x2) || !xfirst1 && (y!=y2)) {
        if(map->weight[dir][pos] <= 0) {
            path_delete(path);
            return 0;
        }
        assert(i<path->length);
        path->pos[i] = pos;
        path->dir[i] = dir;
        i++;
        x += xadd[dir];
        y += yadd[dir];
        pos += w->add[dir];
    }

    dir = (!xfirst1)?(x1<x2?RIGHT:LEFT):(y1<y2?DOWN:UP);
    while(x!=x2 || y!=y2) {
        if(map->weight[dir][pos] <= 0) {
            path_delete(path);
            return 0;
        }
        assert(i<path->length);
        path->pos[i] = pos;
        path->dir[i] = dir;
        i++;
        x += xadd[dir];
        y += yadd[dir];
        pos += w->add[dir];
    }
    
    dir = xfirst2?(x2<x3?RIGHT:LEFT):(y2<y3?DOWN:UP);
    while(xfirst2 && (x!=x3) || !xfirst1 && (y!=y3)) {
        if(map->weight[dir][pos] <= 0) {
            path_delete(path);
            return 0;
        }
        assert(i<path->length);
        path->pos[i] = pos;
        path->dir[i] = dir;
        i++;
        x += xadd[dir];
        y += yadd[dir];
        pos += w->add[dir];
    }

    dir = (!xfirst2)?(x2<x3?RIGHT:LEFT):(y2<y3?DOWN:UP);
    while(x!=x3 || y!=y3) {
        if(map->weight[dir][pos] <= 0) {
            path_delete(path);
            return 0;
        }
        assert(i<path->length);
        path->pos[i] = pos;
        path->dir[i] = dir;
        i++;
        x += xadd[dir];
        y += yadd[dir];
        pos += w->add[dir];
    }
    path->pos[i] = pos;
    path->dir[i] = 4;
    i++;
    assert(i == path->length);
    return decrease_weights(map, path);
}

static int preprocess(graphcut_workspace_t*w, int x1, int y1, int x2, int y2)
{
    int xmid,ymid;
    int sum = 0;
    for(ymid=0;ymid<w->map->height;ymid++)
    for(xmid=0;xmid<w->map->width;xmid++) {
        sum += check3path(w, x1, y1, xmid, ymid, x2, y2, 0, 0);
        sum += check3path(w, x1, y1, xmid, ymid, x2, y2, 0, 1);
        sum += check3path(w, x1, y1, xmid, ymid, x2, y2, 1, 0);
        sum += check3path(w, x1, y1, xmid, ymid, x2, y2, 1, 1);
    }
    return sum;
}

unsigned char* find_cut(map_t*map, int x1, int y1, int x2, int y2)
{
    int pos1 = y1*map->width+x1;
    int pos2 = y2*map->width+x2;

    int max_flow = 0;
    assert(pos1 != pos2);
    graphcut_workspace_t* w = graphcut_workspace_new(map, map->width, map->height, pos1, pos2);
   
    int preprocess_flow = 0;
#ifdef PREPROCESS 
    preprocess_flow = preprocess(w, x1, y1, x2, y2);
    max_flow += preprocess_flow;
#endif

    posqueue_addpos(w->queue1, pos1); w->flags1[pos1] |= ACTIVE|IN_TREE; 
    posqueue_addpos(w->queue2, pos2); w->flags2[pos2] |= ACTIVE|IN_TREE; 
  
    check_map(map);

    while(1) {
	path_t*path;
	while(1) {
	    char done1=0,done2=0;
	    int p1 = posqueue_extract(w->queue1);
	    DBG printf("extend 1 from %d/%d\n", p1%w->width, p1/w->width);
	    if(p1<0) {
                int t;
                /* convert flags1 to map */
                unsigned char*result = w->flags1;
                for(t=0;t<w->size;t++) {
                    if(w->flags1[t]&IN_TREE) result[t] = 1;
                    else result[t] = 0;
                }
                w->flags1 = 0;
		graphcut_workspace_delete(w);
                printf("max flow: %d (%d (%.2f%%) from preprocessing)\n", max_flow, preprocess_flow, preprocess_flow*100.0 / max_flow);
		return result;
	    }
	    path = expand_pos(w, w->queue1, p1, 0, w->flags1, w->flags2);
	    DBG workspace_print(w);
	    if(path)
		break;
	   
#ifdef TWOTREES
	    int p2 = posqueue_extract(w->queue2);
	    DBG printf("extend 2 from %d/%d\n", p2%w->width, p2/w->width);
	    if(p2<0) {
                int t;
                /* convert flags2 to map */
                unsigned char*result = w->flags2;
                for(t=0;t<w->size;t++) {
                    if(w->flags2[t]&IN_TREE) result[t] = 0;
                    else result[t] = 1;
                }
                w->flags2 = 0;
		graphcut_workspace_delete(w);
                printf("max flow: %d (%d (%.2f%%) from preprocessing)\n", max_flow, preprocess_flow, preprocess_flow*100.0 / max_flow);
		return result;
	    }
	    path = expand_pos(w, w->queue2, p2, 1, w->flags2, w->flags1);
	    DBG workspace_print(w);
	    if(path)
		break;
#endif

	}
	DBG printf("found connection between tree1 and tree2\n");
	DBG path_print(path, w->width);

	DBG printf("decreasing weights\n");
	max_flow += decrease_weights(map, path);

	DBG workspace_print(w);

	DBG printf("destroying trees\n");
	combust_tree(w, w->queue1, w->queue2, path);
	DBG workspace_print(w);

	path_delete(path);
    }
}

#ifdef MAIN
#define DIM 512

static int test_small()
{
    weight_t weight[4][16] = 
    {{3,1,2,0, //right
      2,4,1,0,
      3,3,3,0,
      4,1,4,0},
     {0,0,0,0, //up
      4,2,3,2,
      1,3,7,5,
      1,2,5,1},
     {0,3,3,3, //left
      0,3,3,3,
      0,3,1,3,
      0,1,1,1},
     {1,5,1,1, //down
      1,1,1,1,
      1,3,1,3,
      0,0,0,0},
    };

    map_t map;
    map.width = 4;
    map.height = 4;
    map.weight[RIGHT] = weight[0];
    map.weight[UP] = weight[1];
    map.weight[LEFT] = weight[2];
    map.weight[DOWN] = weight[3];

    map_save(&map, "smallmap.tst");

    DBG printf("RIGHT: %02x %02x\n", RIGHT<<6, bit[RIGHT]);
    DBG printf("UP:    %02x %02x\n", UP<<6, bit[UP]);
    DBG printf("LEFT:  %02x %02x\n", LEFT<<6, bit[LEFT]);
    DBG printf("DOWN:  %02x %02x\n", DOWN<<6, bit[DOWN]);
    DBG printf("\n");

    find_cut(&map, 0, 0, 3, 3);
}

void test_big()
{
    weight_t *(weight[4]);
    weight[0] = malloc(sizeof(int)*DIM*DIM);
    weight[1] = malloc(sizeof(int)*DIM*DIM);
    weight[2] = malloc(sizeof(int)*DIM*DIM);
    weight[3] = malloc(sizeof(int)*DIM*DIM);

    int t;
    for(t=0;t<DIM*DIM;t++) {
	weight[0][t] = lrand48() % 256;
	weight[1][t] = lrand48() % 256;
	weight[2][t] = lrand48() % 256;
	weight[3][t] = lrand48() % 256;
    }
    for(t=0;t<DIM;t++) {
	weight[DOWN][(DIM-1)*DIM+t] = 0;
	weight[UP][t] = 0;
	weight[RIGHT][DIM-1+t*DIM] = 0;
	weight[LEFT][t*DIM] = 0;
    }
    
    map_t map;
    map.width = DIM;
    map.height = DIM;
    map.weight[0] = weight[0];
    map.weight[1] = weight[1];
    map.weight[2] = weight[2];
    map.weight[3] = weight[3];

    find_cut(&map, 0,0, DIM-1, DIM-1);
    find_cut(&map, (DIM*1/3),(DIM*1/3), (DIM*2/3),(DIM*2/3)); // at 1/3,1/3 to 2/3,2/3
}

#include "image.h"

void test_matrix()
{
    matrix_t* m = matrix_new(5,4);
    double data[5*4] = { 1.0, 2.0, 3.0, 4.0, 5.0,
                         3.0, 4.0, 2.0, 1.0, 6.0,
                         4.0, 3.0, 1.0, 4.0, 3.0,
                         4.0, 2.0, 5.0, 3.0, 2.0};
    memcpy(m->data, data, sizeof(data));
    
    map_t*map = map_from_matrix(m);
    unsigned char*cut = find_cut(map, 1,1, 3,2);
}

void test_image_extract()
{
    image_t*img1 = image_load("../images/flower_tiny.png");
    image_t*img2 = image_load("../images/flower_tiny_marked2.png");
    
    map_t*map = map_from_image(img1);
   
    int x1=-1,y1=-1;
    int x,y;
    for(y=1;y<map->height-1;y++) {
        for(x=1;x<map->width-1;x++) {
            if(img2->data[y*map->width+x].r == 0
            && img2->data[y*map->width+x].g == 0
            && img2->data[y*map->width+x].b == 255) { // mark color
                map->weight[DOWN][y*map->width+x] = MAX_WEIGHT;
                map->weight[UP][y*map->width+x] = MAX_WEIGHT;
                map->weight[LEFT][y*map->width+x] = MAX_WEIGHT;
                map->weight[RIGHT][y*map->width+x] = MAX_WEIGHT;
                x1 = x;
                y1 = y;
            }
        }
    }
    for(y=0;y<img1->height;y++) {
        for(x=0;x<img1->width;x++) {
            if(img2->data[y*img1->width+x].r == 0
            && img2->data[y*img1->width+x].g == 0
            && img2->data[y*img1->width+x].b == 255) { // mark color
            } else {
                img2->data[y*img1->width+x].r = img1->data[y*img1->width+x].r;
                img2->data[y*img1->width+x].g = img1->data[y*img1->width+x].g;
                img2->data[y*img1->width+x].b = img1->data[y*img1->width+x].b;
            }
            img2->data[y*img1->width+x].a = 255;
        }
    }
    image_save(img2, "mask.png");

    if(x1 < 0 || y1<0) {
        printf("no mark colors found\n");
        return;
    }
    int x2 = 0;
    int y2 = 0;

    /* draw a frame "around" the map */
    for(x=0;x<map->width-1;x++) {
        map->weight[RIGHT][0*map->width+x] = MAX_WEIGHT;
        map->weight[LEFT][0*map->width+x+1] = MAX_WEIGHT;
        map->weight[RIGHT][(map->height-1)*map->width+x] = MAX_WEIGHT;
        map->weight[LEFT][(map->height-1)*map->width+x+1] = MAX_WEIGHT;
    }
    for(y=0;y<map->height-1;y++) {
        map->weight[DOWN][y*map->width+0] = MAX_WEIGHT;
        map->weight[UP][(y+1)*map->width+0] = MAX_WEIGHT;
        map->weight[DOWN][y*map->width+map->width-1] = MAX_WEIGHT;
        map->weight[UP][(y+1)*map->width+map->width-1] = MAX_WEIGHT;
    }
    
    unsigned char*cut = find_cut(map, x1,y1,x2,y2);
    
    int t;
    for(t=0;t<img1->width*img1->height;t++) {
        if(cut[t]) {
            img1->data[t].r = 128+img1->data[t].r/2;
            img1->data[t].g = 128+img1->data[t].g/2;
            img1->data[t].b = 128+img1->data[t].b/2;
        } else {
            img1->data[t].r = img1->data[t].r/2;
            img1->data[t].g = img1->data[t].g/2;
            img1->data[t].b = img1->data[t].b/2;
        }
        img1->data[t].a = 255;
    }
    image_save(img1, "test.png");

    image_delete(img1);
    image_delete(img2);
}

void test_image()
{
    //image_t*img = image_load("../textures/segtest2.png");
    //image_t*img = image_load("../dominator.png");
    image_t*img = image_load("../textures/peppers128.png");

    matrix_t* m = image_extractchannel(img, IMAGE_GRAY);
    
    int x1 = 0;
    int x2 = m->width-1;
    int y1 = (m->height*1/2);
    int y2 = (m->height*1/2);
    map_t*map = map_from_matrix(m);

    int x,y;
    for(y=1;y<m->height-1;y++) {
        map->weight[DOWN][y*map->width+0] = MAX_WEIGHT;
        map->weight[UP  ][y*map->width+0] = MAX_WEIGHT;
        map->weight[DOWN][y*map->width+map->width-1] = MAX_WEIGHT;
        map->weight[UP  ][y*map->width+map->width-1] = MAX_WEIGHT;
    }

    unsigned char*cut1 = find_cut(map, x1, y1, x2, y2);
   
    map_delete(map);

    map = map_from_matrix(m);
    for(y=0;y<map->width;y++) {
        for(x=0;x<map->width;x++) {
            if(x && cut1[y*map->width+x] != cut1[y*map->width+x-1])
                map->weight[LEFT][y*map->width+x] = 0;
            if(x<map->width-1 && cut1[y*map->width+x] != cut1[y*map->width+x+1])
                map->weight[RIGHT][y*map->width+x] = 0;
            if(y && cut1[y*map->width+x] != cut1[(y-1)*map->width+x])
                map->weight[UP][y*map->width+x] = 0;
            if(y<map->height-1 && cut1[y*map->width+x] != cut1[(y+1)*map->width+x])
                map->weight[DOWN][y*map->width+x] = 0;
        }
    }
    for(x=1;x<m->width-1;x++) {
        map->weight[LEFT][x] = MAX_WEIGHT;
        map->weight[RIGHT][x] = MAX_WEIGHT;
        map->weight[LEFT][(map->height-1)*map->width+x] = MAX_WEIGHT;
        map->weight[RIGHT][(map->height-1)*map->width+x] = MAX_WEIGHT;
    }
    
    unsigned char*cut2 = find_cut(map, map->width/2, 0, map->width/2, map->height-1);
    
    int t;
    int oneside=0,otherside=0;
    for(t=0;t<img->width*img->height;t++) {
        if(!cut1[t]) {
            img->data[t].r = 0;
        } 
        if(!cut2[t]) {
            img->data[t].b = 0;
        }
        img->data[t].a = 255;
    }
    image_save(img, "test.png");

    image_delete(img);
    
    
}

void test_twoimages()
{
    //image_t*img1 = image_load("../textures/DARK75.png");
    //image_t*img2 = image_load("../textures/DARK186.png");
    //image_t*img1 = image_load("../textures/DARK75.png");
    //image_t*img2 = image_load("../textures/DARK75.png");
    image_t*img1 = image_load("../textures/DARK159.png");
    image_t*img2 = image_load("../textures/DARK159.png");
    
    matrix_t* m1 = image_extractchannel(img1, IMAGE_GRAY);
    matrix_t* m2 = image_extractchannel(img2, IMAGE_GRAY);

    int xpos = 23;

    assert(img1->height == img2->height);
    assert(img1->width <= img2->width + xpos);

    int xsize = img2->width + xpos;
    int ysize = m1->height;

    map_t*m = map_new(xsize, m1->height);
    int x,y;
    for(y=0;y<ysize;y++) {
        int yup = y>0? y-1: 0;
        int ydown = y<ysize-1? y+1: 0;
        for(x=0;x<xsize;x++) {
            int xleft = x>0? x-1: 0;
            int xright = x<xsize-1? x+1: 0;
            int p = y*m->width+x;
            if(x<xpos) {
                m->weight[LEFT][p] = MAX_WEIGHT;
                m->weight[RIGHT][p] = MAX_WEIGHT;
                m->weight[UP][p] = MAX_WEIGHT;
                m->weight[DOWN][p] = MAX_WEIGHT;
            } else if(x>=img1->width) {
                m->weight[LEFT][p] = MAX_WEIGHT;
                m->weight[RIGHT][p] = MAX_WEIGHT;
                m->weight[UP][p] = MAX_WEIGHT;
                m->weight[DOWN][p] = MAX_WEIGHT;
            } else {
                /*double left = fabs(m1->data[y*m1->width+x - 1] - m2->data[y*m2->width-xpos+x]) +  
                              fabs(m2->data[y*m2->width+x - 1] - m1->data[y*m1->width-xpos+x]);
                double right = fabs(m1->data[y*m1->width+x + 1] - m2->data[y*m2->width-xpos+x]) +
                               fabs(m2->data[y*m2->width+x + 1] - m1->data[y*m1->width-xpos+x]);
                double up = fabs(m1->data[yup*m1->width+x] - m2->data[y*m2->width-xpos+x]) +
                            fabs(m2->data[yup*m2->width+x] - m1->data[y*m1->width-xpos+x]);
                double down = fabs(m1->data[ydown*m1->width+x] - m2->data[y*m2->width-xpos+x]) +
                              fabs(m2->data[ydown*m2->width+x] - m2->data[y*m1->width-xpos+x]);*/
                
                double right = fabs(m1->data[y*m1->width + x] - m2->data[y*m2->width-xpos + xright]);
                double left = fabs(m1->data[y*m1->width + x] - m2->data[y*m2->width-xpos + xleft]);
                double up = fabs(m1->data[y*m1->width+x] - m2->data[yup*m2->width-xpos+x]);
                double down = fabs(m1->data[y*m1->width+x] - m2->data[ydown*m2->width-xpos+x]);

                //m->weight[LEFT][p] = (int)fabs(left) + 1;
                //m->weight[RIGHT][p] = (int)fabs(right) + 1;
                //m->weight[UP][p] = (int)fabs(up) + 1;
                //m->weight[DOWN][p] = (int)fabs(down) + 1;
                m->weight[LEFT][p] = (int)(left*left) + 1;
                m->weight[RIGHT][p] = (int)(right*right) + 1;
                m->weight[UP][p] = (int)(up*up) + 1;
                m->weight[DOWN][p] = (int)(down*down) + 1;
            }
            m->weight[UP][x] = 0;
            m->weight[DOWN][m->width*(m->height-1) + x] = 0;
        }
        m->weight[LEFT][y*m->width] = 0;
        m->weight[RIGHT][y*m->width+m->width-1] = 0;
    }
    int x1 = 0;
    int y1 = (m->height*1/2);
    int x2 = m->width-1;
    int y2 = (m->height*1/2);
    unsigned char*cut = find_cut(m, x1, y1, x2, y2);
    
    int t;
    int oneside=0,otherside=0;
    image_t*img = image_new(xsize,ysize);
    RGBA red = {255,255,0,0};
    for(y=0;y<ysize;y++) {
        for(x=0;x<xsize;x++) {
            if(cut[y*m->width+x] && x<img1->width) {
                img->data[y*img->width+x] = img1->data[y*img1->width+x];
                unsigned char r = img->data[y*img->width+x].r;
                unsigned char b = img->data[y*img->width+x].b;
                //img->data[y*img->width+x].r = b;
                //img->data[y*img->width+x].b = r;
                oneside++;
            } else if(!cut[y*m->width+x] && x>=xpos) {
                img->data[y*img->width+x] = img2->data[y*img1->width+x-xpos];
                otherside++;
            } else {
                img->data[y*img->width+x] = red;
            }
        }
    }
    printf("%d pixels in area1, %d pixels in area2\n", oneside, otherside);
    image_save(img, "test.png");

    image_delete(img);
}


int main()
{
    test_small();
    //test_big();
    //test_image();
    //test_image_extract();
    ///test_matrix();
    //test_twoimages();
}
#endif
