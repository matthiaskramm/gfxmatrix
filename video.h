/* video.h 
   
   video i/o abstraction.

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

#ifndef __video_h__
#define __video_h__

#ifdef __cplusplus
extern "C" {
#endif

#include "image.h"

struct video_device
{
    void*internal;
    int frame;
};

struct video_device* init_video(char*device,int x,int y);
int video_capture(struct video_device * device_wrapper, void*mem);
int video_capture_image(struct video_device * device_wrapper, image_t*img);
void video_close(struct video_device *device);
void video_set_verbose(int level);

#ifdef __cplusplus
}
#endif

#endif
