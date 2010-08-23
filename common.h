#ifndef __common_h__
#define __common_h__

typedef         unsigned long   U32;
typedef         signed long     S32;
typedef         unsigned short  U16;
typedef         signed short    S16;
typedef         unsigned char   U8;
typedef         signed char     S8;

typedef struct _RGBA
{ U8    a;
  U8    r;
  U8    g;
  U8    b;
} RGBA;

typedef struct _YUV
{
  U8	y,u,v;
} YUV;

typedef struct _HLS
{
    double h,l,s;
} HLS;

typedef struct _image
{
    RGBA*data;
    int width;
    int height;
} image_t;

typedef struct _cxform
{
    float rr,rg,rb,ra, tr;
    float gr,gg,gb,ga, tg;
    float br,bg,bb,ba, tb;
    float ar,ag,ab,aa, ta;
} cxform_t;

typedef struct _matrix
{
    double*data;
    int width;
    int height;
} matrix_t;

typedef struct _complex
{
    double real;
    double imag;
} complex_t;

typedef struct _complex_matrix
{
    complex_t*data;
    int width;
    int height;
} complex_matrix_t;

typedef struct _matrixset
{
    matrix_t**m;
    int num;
} matrixset_t;

typedef struct _complex_matrixset
{
    complex_matrix_t**m;
    int num;
} complex_matrixset_t;

typedef struct _bytearray
{
    int width, height;
    unsigned char*data;
} bytearray_t;

typedef struct _bytearrayset
{
    bytearray_t**m;
    int num;
} bytearrayset_t;

typedef struct _window
{
    int x1,y1,x2,y2;
    int width, height;
} window_t;

#endif
