/* gfxwindow.c 

   Simple GUI abstraction.

   Part of the swftools package.

   Copyright (c) 2005 Matthias Kramm <kramm@quiss.org> 

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

#ifdef WIN32
#include "gfxwindow_win32.c"
#else
#include "gfxwindow_unix.c"
#endif

void gfxwindow_paste_image(gfxwindow_t*dest, image_t*src, int x,int y)
{
    int width=src->width;
    int height=src->height;
    if(x+width > dest->width)
        width = dest->width - x;
    if(y+height > dest->height)
        height = dest->height - y;
    if(x<0 || y<0) {
        fprintf(stderr, "negative coordinates not allowed in paste");
        return;
    }
    if(width<=0 || height<=0) {
        fprintf(stderr, "Area (%d,%d,%d,%d) is outside of (0,0,%d,%d)\n", 
                x, y, x+src->width, y+src->height, dest->width, dest->height);
        return;
    }
    int startx = x;
    int starty = y;
    for(y=0;y<height;y++) {
        unsigned char*dline = &dest->currentscr[(y+starty)*dest->width*4+startx*4];
        RGBA*sline = &src->data[y*src->width];
        for(x=0;x<width;x++) {
            dline[x*4+0] = sline[x].r;
            dline[x*4+1] = sline[x].g;
            dline[x*4+2] = sline[x].b;
        }
    }
}

void gfxwindow_paste_bytearray(gfxwindow_t*dest, bytearray_t*src, int x,int y)
{
    int width=src->width;
    int height=src->height;
    if(x+width > dest->width)
        width = dest->width - x;
    if(y+height > dest->height)
        height = dest->height - y;
    if(x<0 || y<0) {
        fprintf(stderr, "negative coordinates not allowed in paste");
        return;
    }
    if(width<=0 || height<=0) {
        fprintf(stderr, "Area (%d,%d,%d,%d) is outside of (0,0,%d,%d)\n", 
                x, y, x+src->width, y+src->height, dest->width, dest->height);
        return;
    }
    int startx = x;
    int starty = y;
    for(y=0;y<height;y++) {
        unsigned char*dline = &dest->currentscr[(y+starty)*dest->width*4+startx*4];
        unsigned char*sline = &src->data[y*src->width];
        for(x=0;x<width;x++) {
            dline[x*4+0] = (sline[x]<<0)&0xe0;
            dline[x*4+1] = (sline[x]<<3)&0xe0;
            dline[x*4+2] = (sline[x]<<6)&0xc0;
        }
    }
}
