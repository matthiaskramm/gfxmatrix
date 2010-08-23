/* gfximage.h 

   Python image interface.

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

#include <Python.h>
#include "matrix.h"
#include "image.h"
#include "segment.h"
#include "txtsynth.h"

staticforward PyTypeObject ImageClass;
staticforward PyTypeObject BytearrayClass;
staticforward PyTypeObject TexturemodelClass;

typedef struct {
    PyObject_HEAD
    image_t*image;
    PyObject* strrepr;
} ImageObject;

typedef struct {
    PyObject_HEAD
    bytearray_t*bytearray;
    int max;
} BytearrayObject;

typedef struct {
    PyObject_HEAD
    texturedata_t*texturedata;
} TexturemodelObject;

static char* strf(char*format, ...)
{
    char buf[1024];
    int l;
    va_list arglist;
    va_start(arglist, format);
    vsprintf(buf, format, arglist);
    va_end(arglist);
    return strdup(buf);
}
#define PY_ERROR(s,args...) (PyErr_SetString(PyExc_Exception, strf(s, ## args)),NULL)
#define PY_NONE Py_BuildValue("s", 0)

//---------------------------------------------------------------------

staticforward PyObject* py_image_crop(PyObject* _self, PyObject* args, PyObject* kwargs);
staticforward PyObject* py_image_paste(PyObject* _self, PyObject* args, PyObject* kwargs);
staticforward PyObject* py_image_cut(PyObject* _self, PyObject* args, PyObject* kwargs);
staticforward PyObject* py_image_clone(PyObject* _self, PyObject* args, PyObject* kwargs);
staticforward PyObject* py_image_compare(PyObject* _self, PyObject* args, PyObject* kwargs);
staticforward PyObject* py_image_astexture(PyObject* _self, PyObject* args, PyObject* kwargs);
staticforward PyObject* py_image_save(PyObject* _self, PyObject* args, PyObject* kwargs);
staticforward PyObject* py_image_segment(PyObject* _self, PyObject* args, PyObject* kwargs);
staticforward PyObject* py_image_openglbitmap(PyObject* _self, PyObject* args, PyObject* kwargs);
staticforward PyObject* py_image_insert_tiling(PyObject* _self, PyObject* args, PyObject* kwargs);
staticforward PyObject* py_image_linear_transform(PyObject* _self, PyObject* args, PyObject* kwargs);
staticforward PyObject* py_image_resize(PyObject* _self, PyObject* args, PyObject* kwargs);
staticforward PyObject* py_image_psnr(PyObject* _self, PyObject* args, PyObject* kwargs);
staticforward PyObject* py_image_signature(PyObject* _self, PyObject* args, PyObject* kwargs);

static PyMethodDef image_methods[] =
{
    /* Image functions */
    {"crop", (PyCFunction)py_image_crop, METH_KEYWORDS, ""},
    {"paste", (PyCFunction)py_image_paste, METH_KEYWORDS, ""},
    {"cut", (PyCFunction)py_image_cut, METH_VARARGS, ""},
    {"clone", (PyCFunction)py_image_clone, METH_VARARGS, ""},
    {"compare", (PyCFunction)py_image_compare, METH_VARARGS, ""},
    {"astexture", (PyCFunction)py_image_astexture, METH_KEYWORDS, ""},
    {"segment", (PyCFunction)py_image_segment, METH_KEYWORDS, ""},
    {"save", (PyCFunction)py_image_save, METH_KEYWORDS, ""},
    {"openglbitmap", (PyCFunction)py_image_openglbitmap, METH_KEYWORDS, ""},
    {"insert_tiling", (PyCFunction)py_image_insert_tiling, METH_KEYWORDS, ""},
    {"transform", (PyCFunction)py_image_linear_transform, METH_KEYWORDS, ""},
    {"resize", (PyCFunction)py_image_resize, METH_KEYWORDS, ""},
    {"signature", (PyCFunction)py_image_signature, METH_KEYWORDS, ""},
    {0,0,0,0}
};

static PyObject* f_image_load(PyObject* parent, PyObject* args, PyObject* kwargs)
{
    static char *kwlist[] = {"filename", NULL};
    char*filename = 0;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "s", kwlist, &filename))
	return NULL;
    
    ImageObject*self = PyObject_New(ImageObject, &ImageClass);
    self->image = image_load(filename);
    self->strrepr = 0;
    if(!self->image) {
        return PY_ERROR("Couldn't load file %s", filename);
    }
    return (PyObject*)self;
}
static PyObject* f_image_new(PyObject* parent, PyObject* args, PyObject* kwargs)
{
    static char *kwlist[] = {"width", "height", NULL};
    int width=0, height=0;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "ii", kwlist, &width, &height))
	return NULL;
    
    ImageObject*self = PyObject_New(ImageObject, &ImageClass);
    self->image = image_new(width, height);
    self->strrepr = 0;
    return (PyObject*)self;
}
static PyObject* f_image_new_random(PyObject* parent, PyObject* args, PyObject* kwargs)
{
    static char *kwlist[] = {"width", "height", NULL};
    int width=0, height=0;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "ii", kwlist, &width, &height))
	return NULL;
    
    ImageObject*self = PyObject_New(ImageObject, &ImageClass);
    self->image = image_new(width, height);

    int t;
    int size = width*height;
    for(t=0;t<width*height;t++) {
        int r = lrand48();
        self->image->data[t].r = r;
        self->image->data[t].g = r>>8;
        self->image->data[t].b = r>>16;
        self->image->data[t].a = 255;
    }
    self->strrepr = 0;
    return (PyObject*)self;
}
static PyObject* py_image_crop(PyObject* _self, PyObject* args, PyObject* kwargs)
{
    ImageObject*self = (ImageObject*)_self;
    static char *kwlist[] = {"width", "height", NULL};
    int width, height;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "ii", kwlist, &width, &height))
	return NULL;
    image_crop(self->image, width, height);
    return PY_NONE;
}
static PyObject* py_image_paste(PyObject* _self, PyObject* args, PyObject* kwargs)
{
    ImageObject*self = (ImageObject*)_self;
    static char *kwlist[] = {"image", "x", "y", NULL};
    ImageObject*other=0;
    int x,y;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!ii", kwlist, &ImageClass, &other, &x, &y))
	return NULL;
    image_paste(self->image, x, y, other->image);
    return PY_NONE;
}
static PyObject* py_image_linear_transform(PyObject* _self, PyObject* args, PyObject* kwargs)
{
    ImageObject*self = (ImageObject*)_self;
    static char *kwlist[] = {"x0", "x1", NULL};
    ImageObject*other=0;
    float x0, x1;
    if(!PyArg_ParseTupleAndKeywords(args, kwargs, "ff", kwlist, &x0, &x1))
	return NULL;
    int t,size=self->image->width*self->image->height;
    for(t=0;t<size;t++) {
        RGBA*p = &self->image->data[t];
        p->r = clamp00ff(p->r*x1 + x0);
        p->g = clamp00ff(p->g*x1 + x0);
        p->b = clamp00ff(p->b*x1 + x0);
        p->a = clamp00ff(p->a*x1 + x0);
    }
    return PY_NONE;
}
static PyObject* py_image_insert_tiling(PyObject* _self, PyObject* args, PyObject* kwargs)
{
    ImageObject*self = (ImageObject*)_self;
    static char *kwlist[] = {"map", "replacement", "marker", NULL};
    ImageObject*other=0;
    BytearrayObject*map=0;
    int marker;
    if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O!|O!i", kwlist, &ImageClass, &other, &BytearrayClass, &map, &marker))
	return NULL;
    if(map) {
        image_insert_tiling(self->image, map->bytearray, other->image, marker);
    } else {
        image_insert_tiling(self->image, 0, other->image, 0);
    }
    return PY_NONE;
}
static PyObject* py_image_cut(PyObject* _self, PyObject* args, PyObject* kwargs)
{
    ImageObject*self = (ImageObject*)_self;
    int x1=0,y1=0,x2=0,y2=0;
    if (!PyArg_ParseTuple(args, "(iiii)", &x1, &y1, &x2, &y2))
	return NULL;
    image_t*img = image_cut(self->image, x1, y1, x2-x1, y2-y1);
    if(!img)
        return PY_ERROR("invalid cut region");
    ImageObject*other = PyObject_New(ImageObject, &ImageClass);
    other->image = img;
    other->strrepr = 0;
    return (PyObject*)other;
}
static PyObject* py_image_clone(PyObject* _self, PyObject* args, PyObject* kwargs)
{
    ImageObject*self = (ImageObject*)_self;
    if (!PyArg_ParseTuple(args, ""))
	return NULL;
    image_t*img = image_clone(self->image);
    if(!img)
        return PY_ERROR("couldn't clone image");
    ImageObject*other = PyObject_New(ImageObject, &ImageClass);
    other->image = img;
    other->strrepr = 0;
    return (PyObject*)other;
}
static PyObject* py_image_resize(PyObject* _self, PyObject* args, PyObject* kwargs)
{
    ImageObject*self = (ImageObject*)_self;
    int width=0, height=0;
    if (!PyArg_ParseTuple(args, "ii", &width, &height))
	return NULL;
    image_t*img = image_new(width, height);

    double fx = (double)self->image->width / (double)width;
    double fy = (double)self->image->height / (double)height;
    int x,y;
    double py = 0;
    for(y=0;y<height;y++) {
        int yy = (int)py;
        RGBA*row1 = &img->data[y*img->width];
        RGBA*row2 = &self->image->data[yy*self->image->width];
        double px = 0;
        for(x=0;x<width;x++) {
            row1[x] = row2[(int)px];
            px += fx;
        }
        py += fy;
    }
    ImageObject*other = PyObject_New(ImageObject, &ImageClass);
    other->image = img;
    other->strrepr = 0;
    return (PyObject*)other;
}
static PyObject* py_image_compare(PyObject* _self, PyObject* args, PyObject* kwargs)
{
    ImageObject*self = (ImageObject*)_self;
    ImageObject*other = 0;
    if (!PyArg_ParseTuple(args, "O!", &ImageClass, &other))
	return NULL;

    double diff = image_compare(self->image, other->image);
    return (PyObject*)PyFloat_FromDouble(diff);
}
static PyObject* py_image_astexture(PyObject* _self, PyObject* args, PyObject* kwargs)
{
    ImageObject*self = (ImageObject*)_self;
    static char *kwlist[] = {"method", NULL};
    int method=2;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "|i", kwlist, &method))
	return NULL;
    TexturemodelObject*texture = PyObject_New(TexturemodelObject, &TexturemodelClass);
    texture->texturedata = texturedata_fromimage(self->image, method);
    return (PyObject*)texture;
}
static PyObject* py_image_segment(PyObject* _self, PyObject* args, PyObject* kwargs)
{
    ImageObject*self = (ImageObject*)_self;
    static char *kwlist[] = {"centers", "bx", "by", "dist", "flags", NULL};
    int num_centers=0;
    int bx=0;
    int by=0;
    int dist=0;
    int flags=SEGMENT_DEFAULTS;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "iiii|i", kwlist, &num_centers, &bx, &by, &dist, &flags))
	return NULL;
    bytearray_t* centers = 0;
    bytearray_t* b = segment(self->image, num_centers, bx, by, dist, /*&centers*/0, flags);
    BytearrayObject*array = PyObject_New(BytearrayObject, &BytearrayClass);
    array->max = -1;
    array->bytearray = b;
    return (PyObject*)array;
}
static PyObject* py_image_openglbitmap(PyObject* _self, PyObject* args, PyObject* kwargs)
{
    ImageObject*self = (ImageObject*)_self;
    static char *kwlist[] = {NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "", kwlist))
	return NULL;
    return PyString_FromStringAndSize((char*)self->image->data, self->image->width*self->image->height*4);
}
static PyObject* py_image_signature(PyObject* _self, PyObject* args, PyObject* kwargs)
{
    ImageObject*self = (ImageObject*)_self;
    static char *kwlist[] = {NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "", kwlist))
	return NULL;

    if(self->image->width!=256 || self->image->height!=256) {
        return PY_ERROR("image must be 256x256");
    }

    PyObject* dict = PyDict_New();

    int t;
    int base = 0;
    for(t=0;t<=4;t++) {
        int width = 256 >> t;
        int height = 256 >> t;
        int x,y,xx,yy;
        for(xx=0;xx<256;xx+=width)
        for(yy=0;yy<256;yy+=height) {
            int x,y;
            int r=0,g=0,b=0;
            for(y=0;y<height;y++) {
                RGBA*c = &self->image->data[(yy+y)*self->image->width+xx];
                for(x=0;x<width;x++) {
                    r+=c[x].r;
                    g+=c[x].g;
                    b+=c[x].b;
                }
            }
            RGBA rgb;
            rgb.r = r / (width*height);
            rgb.g = g / (width*height);
            rgb.b = b / (width*height);
            HLS hls = rgb2hls(rgb);

            int f = base + (int)(hls.h * 17.999) + 18*(int)(hls.s*2.999) + 18*3*(int)(hls.l*2.999);
            PyObject*o =  PyInt_FromLong(f);
            PyDict_SetItem(dict, o, PY_NONE);
            base += 162;
        }
    }
    return dict;
}
static PyObject* py_image_save(PyObject* _self, PyObject* args, PyObject* kwargs)
{
    ImageObject*self = (ImageObject*)_self;
    static char *kwlist[] = {"filename", NULL};
    char*filename = 0;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "s", kwlist, &filename))
	return NULL;
    image_save(self->image, filename);
    return PY_NONE;
}
static void image_dealloc(PyObject* _self) {
    ImageObject* self = (ImageObject*)_self;
    image_delete(self->image);self->image=0;
    PyObject_Del(self);
}
static PyObject* image_getattr(PyObject * _self, char* a)
{
    ImageObject*self = (ImageObject*)_self;
    if(!strcmp(a, "width")) {
        return PyInt_FromLong(self->image->width);
    } else if(!strcmp(a, "height")) {
        return PyInt_FromLong(self->image->height);
    } else if(!strcmp(a, "rgb")) {
        int l = self->image->width*self->image->height;
        unsigned char*data = malloc(self->image->width*self->image->height*3);
        int s,t;
        for(t=0,s=0;t<l;s+=3,t++) {
            data[s+0] = self->image->data[t].r;
            data[s+1] = self->image->data[t].g;
            data[s+2] = self->image->data[t].b;
        }
        return PyString_FromStringAndSize((char*)data,self->image->width*self->image->height*3);
    }
    return Py_FindMethod(image_methods, _self, a);
}
static int image_setattr(PyObject * _self, char* a, PyObject * o) 
{
    ImageObject*self = (ImageObject*)_self;
    return -1;
}
static int image_print(PyObject * _self, FILE *fi, int flags)
{
    ImageObject*self = (ImageObject*)_self;
    fprintf(fi, "%08x(%d)", (int)_self, _self?_self->ob_refcnt:0);
    return 0;
}

//---------------------------------------------------------------------

staticforward PyObject* py_bytearray_save(PyObject* _self, PyObject* args, PyObject* kwargs);
staticforward PyObject* py_bytearray_save_map(PyObject* _self, PyObject* args, PyObject* kwargs);
staticforward PyObject* py_bytearray_findmaxwindow(PyObject* _self, PyObject* args, PyObject* kwargs);
staticforward PyObject* py_bytearray_fillhline(PyObject* _self, PyObject* args, PyObject* kwargs);
staticforward PyObject* py_bytearray_combine_maps(PyObject* _self, PyObject* args, PyObject* kwargs);
staticforward PyObject* py_bytearray_cut(PyObject* _self, PyObject* args, PyObject* kwargs);

static PyMethodDef bytearray_methods[] =
{
    /* Bytearray functions */
    {"save", (PyCFunction)py_bytearray_save, METH_KEYWORDS, ""},
    {"save_map", (PyCFunction)py_bytearray_save_map, METH_KEYWORDS, ""},
    {"find_max_window", (PyCFunction)py_bytearray_findmaxwindow, METH_KEYWORDS, ""},
    {"fillhline", (PyCFunction)py_bytearray_fillhline, METH_KEYWORDS, ""},
    {"combine_maps", (PyCFunction)py_bytearray_combine_maps, METH_KEYWORDS, ""},
    {"cut", (PyCFunction)py_bytearray_cut, METH_KEYWORDS, ""},
    {0,0,0,0}
};

static PyObject* f_bytearray_new(PyObject* parent, PyObject* args, PyObject* kwargs)
{
    static char *kwlist[] = {"width", "height", NULL};
    int width=0, height=0;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "ii", kwlist, &width, &height))
	return NULL;
    
    BytearrayObject*self = PyObject_New(BytearrayObject, &BytearrayClass);
    self->bytearray = bytearray_new(width, height);
    self->max = -1;
    return (PyObject*)self;
}
static PyObject* f_bytearray_load_map(PyObject* parent, PyObject* args, PyObject* kwargs)
{
    static char *kwlist[] = {"filename", NULL};
    char*filename = 0;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "s", kwlist, &filename))
	return NULL;
    
    BytearrayObject*self = PyObject_New(BytearrayObject, &BytearrayClass);
    self->bytearray = bytearray_load_map(filename);
    self->max = -1;
    return (PyObject*)self;
}
static PyObject* py_bytearray_save(PyObject* _self, PyObject* args, PyObject* kwargs)
{
    BytearrayObject*self = (BytearrayObject*)_self;
    static char *kwlist[] = {"filename", NULL};
    char*filename = 0;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "s", kwlist, &filename))
	return NULL;
    bytearray_save(self->bytearray, filename);
    return PY_NONE;
}

static PyObject*listfromwindow(window_t w)
{
    PyObject*tuple = PyTuple_New(4);
    PyTuple_SetItem(tuple, 0, PyInt_FromLong(w.x1));
    PyTuple_SetItem(tuple, 1, PyInt_FromLong(w.y1));
    PyTuple_SetItem(tuple, 2, PyInt_FromLong(w.x2));
    PyTuple_SetItem(tuple, 3, PyInt_FromLong(w.y2));
    return tuple;
}

static PyObject* py_bytearray_findmaxwindow(PyObject* _self, PyObject* args, PyObject* kwargs)
{
    BytearrayObject*self = (BytearrayObject*)_self;
    static char *kwlist[] = {"marker", NULL};
    int marker = 0;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "i", kwlist, &marker))
	return NULL;
    window_t w = find_max_window(self->bytearray, marker);
    return listfromwindow(w);
}
static PyObject* py_bytearray_fillhline(PyObject* _self, PyObject* args, PyObject* kwargs)
{
    BytearrayObject*self = (BytearrayObject*)_self;
    static char *kwlist[] = {"v", "y", "x1", "x2", NULL};
    int v=0,y=0,x1=0,x2=0;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "iiii", kwlist, &v, &y, &x1, &x2))
	return NULL;

    if(y >= self->bytearray->height || y<0 || x1>x2 || x1<0 || x2>self->bytearray->width)
        return PY_ERROR("Bad hline position in image of size %dx%d", self->bytearray->width, self->bytearray->height);
    if(v >= 255)
        return PY_ERROR("Map value overflow: %d", v);

    int x;
    unsigned char*line = &self->bytearray->data[y*self->bytearray->width];
    for(x=x1;x<x2;x++) {
        line[x] = v;
    }
    return PY_NONE;
}
static PyObject* py_bytearray_cut(PyObject* _self, PyObject* args, PyObject* kwargs)
{
    BytearrayObject*self = (BytearrayObject*)_self;
    int x1=0,y1=0,x2=0,y2=0;
    if (!PyArg_ParseTuple(args, "(iiii)", &x1, &y1, &x2, &y2))
	return NULL;
    bytearray_t*img = bytearray_cut(self->bytearray, x1, y1, x2-x1, y2-y1);
    if(!img)
        return PY_ERROR("invalid cut region");
    BytearrayObject*other = PyObject_New(BytearrayObject, &BytearrayClass);
    other->bytearray = img;
    return (PyObject*)other;
}
static PyObject* py_bytearray_combine_maps(PyObject* _self, PyObject* args, PyObject* kwargs)
{
    BytearrayObject*self = (BytearrayObject*)_self;
    static char *kwlist[] = {"other", NULL};
    BytearrayObject*other = 0;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!", kwlist, &BytearrayClass, &other))
	return NULL;

    bytearray_t*b = bytearray_combine_maps(self->bytearray, other->bytearray);
    if(!b)
	return PY_ERROR("Couldn't combine maps");
    
    BytearrayObject*n = PyObject_New(BytearrayObject, &BytearrayClass);
    n->bytearray = b;
    n->max = -1;
    return (PyObject*)n;
}
static PyObject* py_bytearray_save_map(PyObject* _self, PyObject* args, PyObject* kwargs)
{
    BytearrayObject*self = (BytearrayObject*)_self;
    static char *kwlist[] = {"filename", NULL};
    char*filename = 0;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "s", kwlist, &filename))
	return NULL;
    bytearray_save_map(self->bytearray, filename);
    return PY_NONE;
}
static void bytearray_dealloc(PyObject* _self) {
    BytearrayObject* self = (BytearrayObject*)_self;
    bytearray_delete(self->bytearray);self->bytearray = 0;
    PyObject_Del(self);
}
static PyObject* bytearray_getattr(PyObject * _self, char* a)
{
    BytearrayObject*self = (BytearrayObject*)_self;
    if(!strcmp(a, "width")) {
        return PyInt_FromLong(self->bytearray->width);
    } else if(!strcmp(a, "height")) {
        return PyInt_FromLong(self->bytearray->height);
    } else if(!strcmp(a, "rgb")) {
        int l = self->bytearray->width*self->bytearray->height;
        unsigned char*data = malloc(self->bytearray->width*self->bytearray->height*3);
        int s,t;
        for(t=0,s=0;t<l;s+=3,t++) {
            data[s+0] = (self->bytearray->data[t]<<0)&0xe0;
            data[s+1] = (self->bytearray->data[t]<<3)&0xe0;
            data[s+2] = (self->bytearray->data[t]<<6)&0xc0;
        }
        return PyString_FromStringAndSize((char*)data,self->bytearray->width*self->bytearray->height*3);
    } else if(!strcmp(a, "max")) {
        if(self->max < 0) {
            bytearray_t*b = self->bytearray;
            int size = b->width*b->height;
            int t;
            for(t=0;t<size;t++) {
                if(b->data[t]!=255 && b->data[t] > self->max)
                    self->max = b->data[t];
            }
        }
        return PyInt_FromLong(self->max);
    }
    return Py_FindMethod(bytearray_methods, _self, a);
}
static int bytearray_setattr(PyObject * _self, char* a, PyObject * o) 
{
    BytearrayObject*self = (BytearrayObject*)_self;
    return -1;
}
static int bytearray_print(PyObject * _self, FILE *fi, int flags)
{
    BytearrayObject*self = (BytearrayObject*)_self;
    fprintf(fi, "%08x(%d)", (int)_self, _self?_self->ob_refcnt:0);
    return 0;
}

//---------------------------------------------------------------------

staticforward PyObject* py_texturemodel_distance(PyObject* _self, PyObject* args, PyObject* kwargs);
staticforward PyObject* py_texturemodel_synthesize(PyObject* _self, PyObject* args, PyObject* kwargs);
staticforward PyObject* py_texturemodel_setparam(PyObject* _self, PyObject* args, PyObject* kwargs);

static PyMethodDef texturemodel_methods[] =
{
    /* Texturemodel functions */
    {"distance", (PyCFunction)py_texturemodel_distance, METH_KEYWORDS, ""},
    {"synthesize", (PyCFunction)py_texturemodel_synthesize, METH_KEYWORDS, ""},
    {"setparam", (PyCFunction)py_texturemodel_setparam, METH_KEYWORDS, ""},
    {0,0,0,0}
};

static PyObject* py_texturemodel_distance(PyObject* _self, PyObject* args, PyObject* kwargs)
{
    TexturemodelObject*self = (TexturemodelObject*)_self;
    static char *kwlist[] = {"other", NULL};
    TexturemodelObject*other = 0;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!", kwlist, &TexturemodelClass, &other))
	return NULL;
    double dist = texturedata_dist(self->texturedata, other->texturedata);
    return PyFloat_FromDouble(dist);
}
static PyObject* py_texturemodel_setparam(PyObject* _self, PyObject* args, PyObject* kwargs)
{
    TexturemodelObject*self = (TexturemodelObject*)_self;
    static char *kwlist[] = {"key", "value", NULL};
    char*key, *value;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "ss", kwlist, &key, &value))
	return NULL;
    texturedata_set_parameter(key, value);
    return PY_NONE;
}
static PyObject* py_texturemodel_synthesize(PyObject* _self, PyObject* args, PyObject* kwargs)
{
    TexturemodelObject*self = (TexturemodelObject*)_self;
    static char *kwlist[] = {"width", "height", "method", NULL};
    int width=0, height=0, method=0;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "ii|i", kwlist, &width, &height, &method))
	return NULL;
    if(method)
        self->texturedata->type = method;
    image_t*alpha = image_new(width, height);
    image_t*img = texturedata_synthesize_image(self->texturedata, alpha);
    image_delete(alpha);
    
    ImageObject*image = PyObject_New(ImageObject, &ImageClass);
    image->image = img;
    image->strrepr = 0;
    return (PyObject*)image;
}
static void texturemodel_dealloc(PyObject* _self) {
    TexturemodelObject* self = (TexturemodelObject*)_self;
    texturedata_delete(self->texturedata); self->texturedata = 0;
    PyObject_Del(self);
}
static PyObject* texturemodel_getattr(PyObject * _self, char* a)
{
    TexturemodelObject*self = (TexturemodelObject*)_self;
    return Py_FindMethod(texturemodel_methods, _self, a);
}
static int texturemodel_setattr(PyObject * _self, char* a, PyObject * o) 
{
    TexturemodelObject*self = (TexturemodelObject*)_self;
    return -1;
}
static int texturemodel_print(PyObject * _self, FILE *fi, int flags)
{
    TexturemodelObject*self = (TexturemodelObject*)_self;
    fprintf(fi, "%08x(%d)", (int)_self, _self?_self->ob_refcnt:0);
    return 0;
}

//---------------------------------------------------------------------

static PyTypeObject ImageClass =
{
    PyObject_HEAD_INIT(NULL)
    0,
    tp_name: "Image",
    tp_basicsize: sizeof(ImageObject),
    tp_itemsize: 0,
    tp_dealloc: image_dealloc,
    tp_print: image_print,
    tp_getattr: image_getattr,
    tp_setattr: image_setattr,
};
static PyTypeObject BytearrayClass =
{
    PyObject_HEAD_INIT(NULL)
    0,
    tp_name: "Bytearray",
    tp_basicsize: sizeof(BytearrayObject),
    tp_itemsize: 0,
    tp_dealloc: bytearray_dealloc,
    tp_print: bytearray_print,
    tp_getattr: bytearray_getattr,
    tp_setattr: bytearray_setattr,
};
static PyTypeObject TexturemodelClass =
{
    PyObject_HEAD_INIT(NULL)
    0,
    tp_name: "Texturemodel",
    tp_basicsize: sizeof(TexturemodelObject),
    tp_itemsize: 0,
    tp_dealloc: texturemodel_dealloc,
    tp_print: texturemodel_print,
    tp_getattr: texturemodel_getattr,
    tp_setattr: texturemodel_setattr,
};

//=====================================================================


static PyMethodDef gfximage_methods[] =
{
    /* module functions */
    {"image_load", (PyCFunction)f_image_load, METH_KEYWORDS, ""},
    {"image_new", (PyCFunction)f_image_new, METH_KEYWORDS, ""},
    {"image_new_random", (PyCFunction)f_image_new_random, METH_KEYWORDS, ""},
    {"bytearray_new", (PyCFunction)f_bytearray_new, METH_KEYWORDS, ""},
    {"bytearray_load_map", (PyCFunction)f_bytearray_load_map, METH_KEYWORDS, ""},
    {0, 0, 0, 0}
};

void initgfximage(void)
{
    ImageClass.ob_type = &PyType_Type;
    BytearrayClass.ob_type = &PyType_Type;
    TexturemodelClass.ob_type = &PyType_Type;

    PyObject*module = Py_InitModule("gfximage", gfximage_methods);
    PyModule_AddIntConstant(module, "SEGMENT_UV", SEGMENT_UV);
    PyModule_AddIntConstant(module, "SEGMENT_SPATIAL", SEGMENT_SPATIAL);
    PyModule_AddIntConstant(module, "SEGMENT_MARGINAL", SEGMENT_MARGINAL);
    PyModule_AddIntConstant(module, "SEGMENT_DEFAULTS", SEGMENT_DEFAULTS);

}
