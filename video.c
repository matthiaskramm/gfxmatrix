/* video.c 
   
   video i/o abstraction- v4l implementation.

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
#include <unistd.h>
#include <memory.h>
#include <assert.h>
#include <fcntl.h>
#include <errno.h>
#include <sys/ioctl.h>
#include <sys/mman.h>
#include <linux/videodev.h>
#include "video.h"
#include "image.h"

#define myioctl(a,b,c) {int i = ioctl(a,b,c);if(i<0) {fprintf(stderr, "ioctrl returned %d\n",i);perror(#b);close(a);exit(1);}}
#define testioctl(a,b,c) {int i = ioctl(a,b,c);if(i<0) {fprintf(stderr, "ioctrl returned %d\n",i);perror(#b);close(a);}}

static int verbose = 1;

struct video_internal
{
    int dev;
    struct video_capability cap;
    struct video_window window;
    struct video_channel chan;
    struct video_tuner tuner;
    struct video_picture picture;
    struct video_buffer buffer;
    struct video_audio audio;
    struct video_mmap map[2];
    struct video_mbuf mbuf;
    struct video_info play;
    struct vbi_format vbi;
    unsigned char * mem; // mmap
    unsigned int picsize;
    unsigned int frame;
};

void info_capabilities(struct video_capability*cap)
{
    printf("Name: %s\n", cap->name);
    printf("Type: 0x%x ", cap->type);
    if(cap->type&VID_TYPE_CAPTURE) printf("(Can capture)");
    if(cap->type&VID_TYPE_TUNER) printf("(Can tune)");
    if(cap->type&VID_TYPE_TELETEXT) printf("(Does teletext)");
    if(cap->type&VID_TYPE_OVERLAY) printf("(Overlay onto frame buffer)");
    if(cap->type&VID_TYPE_CHROMAKEY) printf("(Overlay by chromakey)");
    if(cap->type&VID_TYPE_CLIPPING) printf("(Can clip)");
    if(cap->type&VID_TYPE_FRAMERAM) printf("(Uses the frame buffer memory)");
    if(cap->type&VID_TYPE_SCALES) printf("(Scalable)");
    if(cap->type&VID_TYPE_MONOCHROME) printf("(Monochrome only)");
    if(cap->type&VID_TYPE_SUBCAPTURE) printf("(Can capture subareas of the image)");
    if(cap->type&VID_TYPE_MPEG_DECODER) printf("(Can decode MPEG streams)");
    if(cap->type&VID_TYPE_MPEG_ENCODER) printf("(Can encode MPEG streams)");
    if(cap->type&VID_TYPE_MJPEG_DECODER) printf("(Can decode MJPEG streams)");
    if(cap->type&VID_TYPE_MJPEG_ENCODER) printf("(Can encode MJPEG streams)");
    printf("\n");
    if(cap->channels)
	printf("channels: %d\n", cap->channels);
    if(cap->audios)
	printf("audio devices: %d\n", cap->audios);
    if(cap->minwidth|cap->maxwidth|cap->minheight|cap->maxheight)
	printf("(%d-%d)*(%d-%d)\n", 
		cap->minwidth, cap->maxwidth,
		cap->minheight, cap->maxheight);

}

void info_channel(struct video_channel*chan)
{
    printf("channel: %d (\"%s\") - %s|%s %s|%s norm:%d\n", chan->channel, chan->name,
	    (chan->type&VIDEO_VC_TUNER)?"tuner":"",
	    (chan->type&VIDEO_VC_AUDIO)?"audio":"",
	    (chan->type&VIDEO_TYPE_TV)?"tv":"",
	    (chan->type&VIDEO_TYPE_CAMERA)?"camera":"", chan->norm);
}

void info_tuner(struct video_tuner*tuner)
{
    printf("tuner: %d (\"%s\") %u-%u -",
	    tuner->tuner, tuner->name, tuner->rangelow, tuner->rangehigh);

    if(tuner->flags&VIDEO_TUNER_PAL) printf("PAL-");
    if(tuner->flags&VIDEO_TUNER_NTSC) printf("NTSC-");
    if(tuner->flags&VIDEO_TUNER_SECAM) printf("Secam-");
    if(tuner->flags&VIDEO_TUNER_LOW) printf("Khz-");
    else	    printf("Mhz-");
    if(tuner->flags&VIDEO_TUNER_NORM) printf("Norm-");
    if(tuner->flags&VIDEO_TUNER_STEREO_ON) printf("Stereo-");
    if(tuner->flags&VIDEO_TUNER_RDS_ON) printf("RDS-");
    if(tuner->flags&VIDEO_TUNER_MBS_ON) printf("MBS-");
   
    printf(" mode:[");
    if(tuner->mode&VIDEO_MODE_PAL) printf("PAL");
    if(tuner->mode&VIDEO_MODE_NTSC) printf("NTSC");
    if(tuner->mode&VIDEO_MODE_SECAM) printf("Secam");
    if(tuner->mode&VIDEO_MODE_AUTO) printf("Auto");
    printf("]");

    printf(" signal:%d\n", tuner->signal);
};

void info_picture(struct video_picture*pict)
{
    printf("brightness:%d, hue:%d, colour:%d, contrast:%d, whiteness:%d, depth:%d, palette:%d\n",
	pict->brightness,
	pict->hue,
	pict->colour,
	pict->contrast,
	pict->whiteness,
	pict->depth,
	pict->palette);
}

void info_buffer(struct video_buffer*buffer)
{
    printf("at %08x, %d*%d, %d bpp, %d bytes per line\n",
	    buffer->base,
	    buffer->width,buffer->height,
	    buffer->depth,
	    buffer->bytesperline);
}

void info_window(struct video_window*window)
{
    printf("%d*%d:%d:%d chroma:%d, flags:%x clipcount:%d interlace:%d %s\n",
	    window->width,window->height,window->x,window->y,window->chromakey, window->flags, window->clipcount,
	    window->flags&VIDEO_WINDOW_INTERLACE,
	    (window->flags&VIDEO_WINDOW_CHROMAKEY)?"(overlay by chromakey)":"");
}

void info_audio(struct video_audio*audio)
{
    printf("audio %d (\"%s\") volume:%d bass:%d tremble:%d -", audio->audio,
	    audio->name,
	    audio->volume, audio->bass, audio->treble);
    if(audio->flags&VIDEO_AUDIO_MUTE) printf("MUTE-");
    if(audio->flags&VIDEO_AUDIO_MUTABLE) printf("MUTABLE-");
    if(audio->flags&VIDEO_AUDIO_VOLUME) printf("VOLUME-");
    if(audio->flags&VIDEO_AUDIO_BASS) printf("BASS-");
    if(audio->flags&VIDEO_AUDIO_TREBLE) printf("TREBLE-");
    if(audio->mode & VIDEO_SOUND_MONO) printf("MONO-");
    if(audio->mode & VIDEO_SOUND_STEREO) printf("STEREO-");
    if(audio->mode & VIDEO_SOUND_LANG1) printf("LANG1-");
    if(audio->mode & VIDEO_SOUND_LANG2) ("LANG2-");
    printf(" balance:%d volumestep:%d\n", audio->balance, audio->step);
}

void info_palette(int i)
{
    struct pal_info_t {
	int nr;
	char*name;
    } pal_info[] = {
{VIDEO_PALETTE_GREY,     "GREY"},
{VIDEO_PALETTE_HI240,    "HI240"},
{VIDEO_PALETTE_RGB565,   "RGB565"},
{VIDEO_PALETTE_RGB24,    "RGB24"},
{VIDEO_PALETTE_RGB32,    "RGB32"},
{VIDEO_PALETTE_RGB555,   "RGB555"},
{VIDEO_PALETTE_YUV422,   "YUV422"},
{VIDEO_PALETTE_YUYV,     "YUYV"},
{VIDEO_PALETTE_UYVY,     "UYVY"},
{VIDEO_PALETTE_YUV420,   "YUV420"},
{VIDEO_PALETTE_YUV411,   "YUV411"},
{VIDEO_PALETTE_RAW,      "RAW"},
{VIDEO_PALETTE_YUV422P,  "YUV422P"},
{VIDEO_PALETTE_YUV411P,  "YUV411P"},
{VIDEO_PALETTE_YUV420P,  "YUV420P"},
{VIDEO_PALETTE_YUV410P,  "YUV410P"},
{VIDEO_PALETTE_PLANAR,   "PLANAR"},
{VIDEO_PALETTE_COMPONENT,"COMPONENT"}};
    int t;
    for(t=0;t<sizeof(pal_info)/sizeof(struct pal_info_t);t++)
    {
	if(pal_info[t].nr == i) {
	    printf("%s", pal_info[t].name);
	    return;
	}
    }
    printf("UNKNOWN_PALETTE");
}

void info_mmap(struct video_mmap*map)
{

    printf("frame %d, %d*%d, format:%d (", 
	    map->frame, map->width, map->height, map->format);
    info_palette(map->format);
    printf(")\n");
}

void info_mbuf(struct video_mbuf*mbuf)
{
    int t;
    printf("size: %08x, frames: %d offsets:/", mbuf->size, mbuf->frames);
    for(t=0;t<mbuf->frames;t++)
	printf("%08x/", mbuf->offsets[t]);
    printf("\n");
}
void info_playinfo(struct video_info*play)
{
    printf("framecount:%d, %d*%d, time:%d, pictype:%d, ref:%d, %s\n",
	    play->frame_count, play->h_size, play->v_size, play->smpte_timecode, play->picture_type,
	    play->temporal_reference, play->user_data);

}
void info_vbi(struct vbi_format*vbi)
{
    printf("samplingrate:%d samples_per_line:%d sample_format:%d start:%d:%d count:%d:%d %s %s\n",
	   vbi->sampling_rate, 
	   vbi->samples_per_line, 
	   vbi->sample_format, 
	   vbi->start[0], 
	   vbi->start[1], 
	   vbi->count[0], 
	   vbi->count[1], 
	   vbi->flags&VBI_UNSYNC?"unsync":"",
	   vbi->flags&VBI_INTERLACED?"interlaced":""
	   );
};

    /* myioctl(dev, VIDIOCCAPTURE, &one);*/
    /*myioctl(dev, VIDIOCSWRITEMODE, &two);*/
/*	printf("VBI FMT:\n");
	testioctl(device->vbi_dev, VIDIOCGVBIFMT, &device->vbi);
	info_vbi(&device->vbi);*/

struct video_device* init_video(char*devname,int x,int y)
{
    char* channel = "Composite1";

    const int zero = 0;
    const int one = 1;
    const int two = 2;
    const int three = 3;
    const int four = 4;
    const int five = 5;
    int t;
    struct video_internal* device;
    struct video_device* device_wrapper;

    device_wrapper = (struct video_device*)malloc(sizeof(struct video_device));
    device = (struct video_internal*)malloc(sizeof(struct video_internal));
    device_wrapper->internal = device;

    /* VIDEO */

    printf("Opening %s\n", devname);
    device->dev = open(devname, O_RDWR);

    if(device->dev < 0) {
	char buf[64];
	sprintf(buf, "open %s", devname);
	perror(buf);
	exit(1);
    }

    if(verbose) printf("CAPABILITIES:\n");
    myioctl(device->dev, VIDIOCGCAP, &device->cap);
    if(verbose)
	info_capabilities(&device->cap);

    if(verbose) printf("CHANNELS:\n");
    for(t=0;t<device->cap.channels;t++)
    {
	struct video_channel chan;
	int set = 0;
	chan.channel = t;
	myioctl(device->dev, VIDIOCGCHAN, &chan);

	if(verbose) info_channel(&chan);

	if(t==0)
	    memcpy(&device->chan, &chan, sizeof(chan));

	if(!strcmp(chan.name, channel)) {
	    set = 1;
	    memcpy(&device->chan, &chan, sizeof(chan));
	    if(verbose) {
		printf("^^^^^^^^^^^^^^^^^^^^^^^^ selecting channel %d\n", chan.channel);
	    }
	}

	if(chan.type&VIDEO_VC_TUNER) {
	    unsigned long freq;
	    device->tuner.tuner = 0;
	    myioctl(device->dev, VIDIOCGTUNER, &device->tuner);
	    if(verbose) {
	        info_tuner(&device->tuner);
	    }

	    myioctl(device->dev, VIDIOCGFREQ, &freq);
	    if(verbose) {
		printf("tuner frequency: %d\n", freq);
	    }
	}
    }
    if(verbose) printf("setting channel \"%s\"...\n", device->chan.name);
    myioctl(device->dev, VIDIOCSCHAN, &device->chan);

    if(verbose) printf("device->picture:\n");
    myioctl(device->dev, VIDIOCGPICT, &device->picture);
    if(verbose) info_picture(&device->picture);
    device->picture.brightness = 2;
    device->picture.contrast = 49152;
    myioctl(device->dev, VIDIOCSPICT, &device->picture);

    if(verbose) printf("device->window:\n");
    myioctl(device->dev, VIDIOCGWIN, &device->window);
    if(verbose) info_window(&device->window);
    /* ? framebuffer overlay: 
    myioctl(device->dev, VIDIOCSWIN, &device->window);*/

    if(verbose) printf("device->buffer:\n");
    myioctl(device->dev, VIDIOCGFBUF, &device->buffer);
    if(verbose) info_buffer(&device->buffer);

    if(verbose) printf("device->audio:\n");
    for(t=0;t<device->cap.audios;t++)
    {
	myioctl(device->dev, VIDIOCGAUDIO, &device->audio);
	if(verbose) info_audio(&device->audio);
    }
    
    if(verbose) printf("device->mbuf:\n");
    myioctl(device->dev, VIDIOCGMBUF, &device->mbuf);
    if(verbose) info_mbuf(&device->mbuf);

    device->mem =(unsigned char *)mmap(0,device->mbuf.size,PROT_READ|PROT_WRITE,MAP_SHARED,device->dev,0);
    if((unsigned char*)-1 == device->mem) {
        fprintf(stderr, "Unable to allocate memory!\n");
    }
    if(verbose) printf("Memory: %08x\n", device->mem);

    device->map[0].frame = 0;
    device->map[1].frame = 1;
    device->map[0].width = 
    device->map[1].width = x;
    device->map[0].height = 
    device->map[1].height = y;
    //device->map[0].format = device->map[1].format = VIDEO_PALETTE_GREY;
    device->map[0].format = 
    device->map[1].format = VIDEO_PALETTE_RGB24; //RGB32

    device->picsize = ((3 * x + 3)&~3) * y;

    device_wrapper->frame = 0;
    device->frame = 0;
    return device_wrapper;
}

int video_capture(struct video_device * device_wrapper, void*mem)
{
    struct video_internal*device = (struct video_internal*)(device_wrapper->internal);
    int f,nf;
    if(!device->frame) {
	int n;
	for(n=0;n<device->mbuf.frames;n++) {
	    myioctl(device->dev, VIDIOCMCAPTURE, &device->map[n]);
	}
    }
    f = device->frame % device->mbuf.frames;

again:
    /* returns -1 in case of a sigint interrupting */
    if(ioctl(device->dev, VIDIOCSYNC, &device->map[f].frame)<0) {
	perror("video_capture()");
        if(errno == 4) /* interrupted system call */
            goto again;
	return 0;
    }
    memcpy(mem, device->mem+device->mbuf.offsets[f], device->picsize);
    myioctl(device->dev, VIDIOCMCAPTURE, &device->map[f]);

    device_wrapper->frame ++;
    
    device->frame ++;
    return 1;
}

int video_capture_image(struct video_device * device_wrapper, image_t*img)
{
    struct video_internal*device = (struct video_internal*)(device_wrapper->internal);
    int f,nf;
    if(!device->frame) {
	int n;
	for(n=0;n<device->mbuf.frames;n++) {
	    myioctl(device->dev, VIDIOCMCAPTURE, &device->map[n]);
	}
    }
    f = device->frame % device->mbuf.frames;

again:
    /* returns -1 in case of a sigint interrupting */
    if(ioctl(device->dev, VIDIOCSYNC, &device->map[f].frame)<0) {
	perror("video_capture()");
        if(errno == 4) /* interrupted system call */
            goto again;
	return 0;
    }

    assert(img->width*img->height*3==device->picsize);

    // memcpy(mem, device->mem+device->mbuf.offsets[f], device->picsize);
    unsigned char*src = device->mem+device->mbuf.offsets[f];
    int s=0,t;
    for(t=0;t<device->picsize;t+=3) {
        img->data[s].r = src[t];
        img->data[s].g = src[t+1];
        img->data[s].b = src[t+2];
        img->data[s].a = 0xff;
        s++;
    }

    myioctl(device->dev, VIDIOCMCAPTURE, &device->map[f]);

    device_wrapper->frame ++;
    
    device->frame ++;
    return 1;
}

void video_close(struct video_device *device_wrapper)
{
    struct video_internal*device = (struct video_internal*)(device_wrapper->internal);
   
    munmap(device->mem, device->mbuf.size);
/*    myioctl(device->dev, VIDIOCGPLAYINFO, &device->play);
    printf("PLAYINFO:\n");
    info_playinfo(&device->play);*/

    close(device->dev);
    free(device);
    free(device_wrapper);
}

void video_set_verbose(int level)
{
    verbose = level;
}
