/*
 * The Python Imaging Library.
 *
 * the imaging library bindings
 *
 * history:
 * 1995-09-24 fl   Created
 * 1996-03-24 fl   Ready for first public release (release 0.0)
 * 1996-03-25 fl   Added fromstring (for Jack's "img" library)
 * 1996-03-28 fl   Added channel operations
 * 1996-03-31 fl   Added point operation
 * 1996-04-08 fl   Added new/new_block/new_array factories
 * 1996-04-13 fl   Added decoders
 * 1996-05-04 fl   Added palette hack
 * 1996-05-12 fl   Compile cleanly as C++
 * 1996-05-19 fl   Added matrix conversions, gradient fills
 * 1996-05-27 fl   Added display_mode
 * 1996-07-22 fl   Added getbbox, offset
 * 1996-07-23 fl   Added sequence semantics
 * 1996-08-13 fl   Added logical operators, point mode
 * 1996-08-16 fl   Modified paste interface
 * 1996-09-06 fl   Added putdata methods, use abstract interface
 * 1996-11-01 fl   Added xbm encoder
 * 1996-11-04 fl   Added experimental path stuff, draw_lines, etc
 * 1996-12-10 fl   Added zip decoder, crc32 interface
 * 1996-12-14 fl   Added modulo arithmetics
 * 1996-12-29 fl   Added zip encoder
 * 1997-01-03 fl   Added fli and msp decoders
 * 1997-01-04 fl   Added experimental sun_rle and tga_rle decoders
 * 1997-01-05 fl   Added gif encoder, getpalette hack
 * 1997-02-23 fl   Added histogram mask
 * 1997-05-12 fl   Minor tweaks to match the IFUNC95 interface
 * 1997-05-21 fl   Added noise generator, spread effect
 * 1997-06-05 fl   Added mandelbrot generator
 * 1997-08-02 fl   Modified putpalette to coerce image mode if necessary
 * 1998-01-11 fl   Added INT32 support
 * 1998-01-22 fl   Fixed draw_points to draw the last point too
 * 1998-06-28 fl   Added getpixel, getink, draw_ink
 * 1998-07-12 fl   Added getextrema
 * 1998-07-17 fl   Added point conversion to arbitrary formats
 * 1998-09-21 fl   Added support for resampling filters
 * 1998-09-22 fl   Added support for quad transform
 * 1998-12-29 fl   Added support for arcs, chords, and pieslices
 * 1999-01-10 fl   Added some experimental arrow graphics stuff
 * 1999-02-06 fl   Added draw_bitmap, font acceleration stuff
 * 2001-04-17 fl   Fixed some egcs compiler nits
 * 2001-09-17 fl   Added screen grab primitives (win32)
 * 2002-03-09 fl   Added stretch primitive
 * 2002-03-10 fl   Fixed filter handling in rotate
 * 2002-06-06 fl   Added I, F, and RGB support to putdata
 * 2002-06-08 fl   Added rankfilter
 * 2002-06-09 fl   Added support for user-defined filter kernels
 * 2002-11-19 fl   Added clipboard grab primitives (win32)
 * 2002-12-11 fl   Added draw context
 * 2003-04-26 fl   Tweaks for Python 2.3 beta 1
 * 2003-05-21 fl   Added createwindow primitive (win32)
 * 2003-09-13 fl   Added thread section hooks
 * 2003-09-15 fl   Added expand helper
 * 2003-09-26 fl   Added experimental LA support
 * 2004-02-21 fl   Handle zero-size images in quantize
 * 2004-06-05 fl   Added ptr attribute (used to access Imaging objects)
 * 2004-06-05 fl   Don't crash when fetching pixels from zero-wide images
 * 2004-09-17 fl   Added getcolors
 * 2004-10-04 fl   Added modefilter
 * 2005-10-02 fl   Added access proxy
 * 2006-06-18 fl   Always draw last point in polyline
 *
 * Copyright (c) 1997-2006 by Secret Labs AB
 * Copyright (c) 1995-2006 by Fredrik Lundh
 *
 * See the README file for information on usage and redistribution.
 */

#define PY_SSIZE_T_CLEAN
#include "Python.h"
#include "hpy.h"
#include "hpy_utils.h"

#ifdef HAVE_LIBJPEG
#include "jconfig.h"
#endif

#ifdef HAVE_LIBZ
#include "zlib.h"
#endif

#ifdef HAVE_LIBTIFF
#ifndef _TIFFIO_
#include <tiffio.h>
#endif
#endif

#include "libImaging/Imaging.h"

#define _USE_MATH_DEFINES
#include <math.h>

/* Configuration stuff. Feel free to undef things you don't need. */
#define WITH_IMAGECHOPS  /* ImageChops support */
#define WITH_IMAGEDRAW   /* ImageDraw support */
#define WITH_MAPPING     /* use memory mapping to read some file formats */
#define WITH_IMAGEPATH   /* ImagePath stuff */
#define WITH_ARROW       /* arrow graphics stuff (experimental) */
#define WITH_EFFECTS     /* special effects */
#define WITH_QUANTIZE    /* quantization support */
#define WITH_RANKFILTER  /* rank filter */
#define WITH_MODEFILTER  /* mode filter */
#define WITH_THREADING   /* "friendly" threading support */
#define WITH_UNSHARPMASK /* Kevin Cazabon's unsharpmask module */

#undef VERBOSE

#define B16(p, i) ((((int)p[(i)]) << 8) + p[(i) + 1])
#define L16(p, i) ((((int)p[(i) + 1]) << 8) + p[(i)])
#define S16(v) ((v) < 32768 ? (v) : ((v)-65536))

/* -------------------------------------------------------------------- */
/* OBJECT ADMINISTRATION                        */
/* -------------------------------------------------------------------- */

typedef struct {
    PyObject_HEAD Imaging image;
    ImagingAccess access;
} ImagingObject;

static PyTypeObject *Imaging_Type;
static HPyGlobal hg_Imaging_Type;

HPyType_LEGACY_HELPERS(ImagingObject);

#ifdef WITH_IMAGEDRAW

typedef struct {
    /* to write a character, cut out sxy from glyph data, place
       at current position plus dxy, and advance by (dx, dy) */
    int dx, dy;
    int dx0, dy0, dx1, dy1;
    int sx0, sy0, sx1, sy1;
} Glyph;

typedef struct {
    PyObject_HEAD ImagingObject *ref;
    Imaging bitmap;
    int ysize;
    int baseline;
    Glyph glyphs[256];
} ImagingFontObject;

static PyTypeObject ImagingFont_Type;

typedef struct {
    PyObject_HEAD ImagingObject *image;
    UINT8 ink[4];
    int blend;
} ImagingDrawObject;

static PyTypeObject ImagingDraw_Type;

#endif

typedef struct {
    ImagingObject *image;
    int readonly;
} PixelAccessObject;
HPyType_HELPERS(PixelAccessObject)

static HPyGlobal hg_PixelAccess_Type;

PyObject *
PyImagingNew(Imaging imOut) {
    ImagingObject *imagep;

    if (!imOut) {
        return NULL;
    }

    imagep = PyObject_New(ImagingObject, Imaging_Type);
    if (imagep == NULL) {
        ImagingDelete(imOut);
        return NULL;
    }

#ifdef VERBOSE
    printf("imaging %p allocated\n", imagep);
#endif

    imagep->image = imOut;
    imagep->access = ImagingAccessNew(imOut);

    return (PyObject *)imagep;
}

HPy HPyImagingNew(HPyContext *ctx, Imaging imOut) {
    ImagingObject *imagep;

    if (!imOut) {
        return HPy_NULL;
    }

    HPy hImagingType = HPyGlobal_Load(ctx, hg_Imaging_Type);
    HPy hOut = HPy_New(ctx, hImagingType, &imagep);
    HPy_Close(ctx, hImagingType);

#ifdef VERBOSE
    printf("imaging %p allocated\n", imagep);
#endif

    imagep->image = imOut;
    imagep->access = ImagingAccessNew(imOut);

    return hOut;
}

HPyDef_SLOT(Imaging_destroy, HPy_tp_destroy)
static void Imaging_destroy_impl(void *obj) {
#ifdef VERBOSE
    printf("imaging %p deleted\n", imagep);
#endif
    ImagingObject *imagep = (ImagingObject *)obj;
    if (imagep->access) {
        ImagingAccessDelete(imagep->image, imagep->access);
    }
    ImagingDelete(imagep->image);
}

#define PyImaging_Check(op) (Py_TYPE(op) == Imaging_Type)

static inline int
HPyImaging_Check(HPyContext *ctx, HPy obj) {
    HPy type = HPy_Type(ctx, obj);
    assert(!HPy_IsNull(type));
    HPy h_Imaging_Type = HPyGlobal_Load(ctx, hg_Imaging_Type);
    int res = HPy_Is(ctx, type, h_Imaging_Type);
    HPy_Close(ctx, type);
    HPy_Close(ctx, h_Imaging_Type);
    return res;
}

Imaging
PyImaging_AsImaging(PyObject *op) {
    if (!PyImaging_Check(op)) {
        PyErr_BadInternalCall();
        return NULL;
    }

    return ((ImagingObject *)op)->image;
}

Imaging
HPyImaging_AsImaging(HPyContext *ctx, HPy op) {
    if (!HPyImaging_Check(ctx, op)) {
        HPyErr_BadInternalCall(ctx);
        return NULL;
    }

    return ImagingObject_AsStruct(ctx, op)->image;
}

/* -------------------------------------------------------------------- */
/* THREAD HANDLING                                                      */
/* -------------------------------------------------------------------- */

void
ImagingSectionEnter(ImagingSectionCookie *cookie) {
#ifdef WITH_THREADING
    *cookie = (PyThreadState *)PyEval_SaveThread();
#endif
}

void
ImagingSectionLeave(ImagingSectionCookie *cookie) {
#ifdef WITH_THREADING
    PyEval_RestoreThread((PyThreadState *)*cookie);
#endif
}

/* -------------------------------------------------------------------- */
/* BUFFER HANDLING                                                      */
/* -------------------------------------------------------------------- */
/* Python compatibility API */

int
PyImaging_CheckBuffer(PyObject *buffer) {
    return PyObject_CheckBuffer(buffer);
}

int
PyImaging_GetBuffer(PyObject *buffer, Py_buffer *view) {
    /* must call check_buffer first! */
    return PyObject_GetBuffer(buffer, view, PyBUF_SIMPLE);
}

/* -------------------------------------------------------------------- */
/* EXCEPTION REROUTING                                                  */
/* -------------------------------------------------------------------- */

/* error messages */
static const char *must_be_sequence = "argument must be a sequence";
static const char *must_be_two_coordinates =
    "coordinate list must contain exactly 2 coordinates";
static const char *wrong_mode = "unrecognized image mode";
static const char *wrong_raw_mode = "unrecognized raw mode";
static const char *outside_image = "image index out of range";
static const char *outside_palette = "palette index out of range";
static const char *wrong_palette_size = "invalid palette size";
static const char *no_palette = "image has no palette";
static const char *readonly = "image is readonly";
/* static const char* no_content = "image has no content"; */

void *
ImagingError_OSError(void) {
    PyErr_SetString(PyExc_OSError, "error when accessing file");
    return NULL;
}

void *
ImagingError_MemoryError(void) {
    return PyErr_NoMemory();
}

void *
ImagingError_Mismatch(void) {
    PyErr_SetString(PyExc_ValueError, "images do not match");
    return NULL;
}

void *
ImagingError_ModeError(void) {
    PyErr_SetString(PyExc_ValueError, "image has wrong mode");
    return NULL;
}

void *
ImagingError_ValueError(const char *message) {
    PyErr_SetString(
        PyExc_ValueError, (message) ? (char *)message : "unrecognized argument value");
    return NULL;
}

void
ImagingError_Clear(void) {
    PyErr_Clear();
}

/* -------------------------------------------------------------------- */
/* HELPERS                                */
/* -------------------------------------------------------------------- */

static int
getbands(const char *mode) {
    Imaging im;
    int bands;

    /* FIXME: add primitive to libImaging to avoid extra allocation */
    im = ImagingNew(mode, 0, 0);
    if (!im) {
        return -1;
    }

    bands = im->bands;

    ImagingDelete(im);

    return bands;
}

#define TYPE_UINT8 (0x100 | sizeof(UINT8))
#define TYPE_INT32 (0x200 | sizeof(INT32))
#define TYPE_FLOAT16 (0x500 | sizeof(FLOAT16))
#define TYPE_FLOAT32 (0x300 | sizeof(FLOAT32))
#define TYPE_DOUBLE (0x400 | sizeof(double))

static void *
getlist(PyObject *arg, Py_ssize_t *length, const char *wrong_length, int type) {
    /* - allocates and returns a c array of the items in the
          python sequence arg.
       - the size of the returned array is in length
       - all of the arg items must be numeric items of the type
          specified in type
       - sequence length is checked against the length parameter IF
          an error parameter is passed in wrong_length
       - caller is responsible for freeing the memory
    */

    Py_ssize_t i, n;
    int itemp;
    double dtemp;
    FLOAT32 ftemp;
    UINT8 *list;
    PyObject *seq;
    PyObject *op;

    if (!PySequence_Check(arg)) {
        PyErr_SetString(PyExc_TypeError, must_be_sequence);
        return NULL;
    }

    n = PySequence_Size(arg);
    if (length && wrong_length && n != *length) {
        PyErr_SetString(PyExc_ValueError, wrong_length);
        return NULL;
    }

    /* malloc check ok, type & ff is just a sizeof(something)
       calloc checks for overflow */
    list = calloc(n, type & 0xff);
    if (!list) {
        return ImagingError_MemoryError();
    }

    seq = PySequence_Fast(arg, must_be_sequence);
    if (!seq) {
        free(list);
        return NULL;
    }

    for (i = 0; i < n; i++) {
        op = PySequence_Fast_GET_ITEM(seq, i);
        // DRY, branch prediction is going to work _really_ well
        // on this switch. And 3 fewer loops to copy/paste.
        switch (type) {
            case TYPE_UINT8:
                itemp = PyLong_AsLong(op);
                list[i] = CLIP8(itemp);
                break;
            case TYPE_INT32:
                itemp = PyLong_AsLong(op);
                memcpy(list + i * sizeof(INT32), &itemp, sizeof(itemp));
                break;
            case TYPE_FLOAT32:
                ftemp = (FLOAT32)PyFloat_AsDouble(op);
                memcpy(list + i * sizeof(ftemp), &ftemp, sizeof(ftemp));
                break;
            case TYPE_DOUBLE:
                dtemp = PyFloat_AsDouble(op);
                memcpy(list + i * sizeof(dtemp), &dtemp, sizeof(dtemp));
                break;
        }
    }

    Py_DECREF(seq);

    if (PyErr_Occurred()) {
        free(list);
        return NULL;
    }

    if (length) {
        *length = n;
    }

    return list;
}

FLOAT32
float16tofloat32(const FLOAT16 in) {
    UINT32 t1;
    UINT32 t2;
    UINT32 t3;
    FLOAT32 out[1] = {0};

    t1 = in & 0x7fff;  // Non-sign bits
    t2 = in & 0x8000;  // Sign bit
    t3 = in & 0x7c00;  // Exponent

    t1 <<= 13;  // Align mantissa on MSB
    t2 <<= 16;  // Shift sign bit into position

    t1 += 0x38000000;  // Adjust bias

    t1 = (t3 == 0 ? 0 : t1);  // Denormals-as-zero

    t1 |= t2;  // Re-insert sign bit

    memcpy(out, &t1, 4);
    return out[0];
}

static inline HPy
getpixel(HPyContext *ctx, Imaging im, ImagingAccess access, int x, int y) {
    union {
        UINT8 b[4];
        UINT16 h;
        INT32 i;
        FLOAT32 f;
    } pixel;

    if (x < 0) {
        x = im->xsize + x;
    }
    if (y < 0) {
        y = im->ysize + y;
    }

    if (x < 0 || x >= im->xsize || y < 0 || y >= im->ysize) {
        HPyErr_SetString(ctx, ctx->h_IndexError, outside_image);
        return HPy_NULL;
    }

    access->get_pixel(im, x, y, &pixel);

    switch (im->type) {
        case IMAGING_TYPE_UINT8:
            switch (im->bands) {
                case 1:
                    return HPyLong_FromLong(ctx, pixel.b[0]);
                case 2:
                    return HPy_BuildValue(ctx, "II", pixel.b[0], pixel.b[1]);
                case 3:
                    return HPy_BuildValue(ctx, "III", pixel.b[0], pixel.b[1], pixel.b[2]);
                case 4:
                    return HPy_BuildValue(
                        ctx, "IIII", pixel.b[0], pixel.b[1], pixel.b[2], pixel.b[3]);
            }
            break;
        case IMAGING_TYPE_INT32:
            return HPyLong_FromLong(ctx, pixel.i);
        case IMAGING_TYPE_FLOAT32:
            return HPyFloat_FromDouble(ctx, pixel.f);
        case IMAGING_TYPE_SPECIAL:
            if (strncmp(im->mode, "I;16", 4) == 0) {
                return HPyLong_FromLong(ctx, pixel.h);
            }
            break;
    }

    /* unknown type */
    return HPy_Dup(ctx, ctx->h_None);
}

static char *
getink(PyObject *color, Imaging im, char *ink) {
    int g = 0, b = 0, a = 0;
    double f = 0;
    /* Windows 64 bit longs are 32 bits, and 0xFFFFFFFF (white) is a
       python long (not int) that raises an overflow error when trying
       to return it into a 32 bit C long
    */
    PY_LONG_LONG r = 0;
    FLOAT32 ftmp;
    INT32 itmp;

    /* fill ink buffer (four bytes) with something that can
       be cast to either UINT8 or INT32 */

    int rIsInt = 0;
    if (PyTuple_Check(color) && PyTuple_GET_SIZE(color) == 1) {
        color = PyTuple_GetItem(color, 0);
    }
    if (im->type == IMAGING_TYPE_UINT8 || im->type == IMAGING_TYPE_INT32 ||
        im->type == IMAGING_TYPE_SPECIAL) {
        if (PyLong_Check(color)) {
            r = PyLong_AsLongLong(color);
            if (r == -1 && PyErr_Occurred()) {
                return NULL;
            }
            rIsInt = 1;
        } else if (im->type == IMAGING_TYPE_UINT8) {
            if (!PyTuple_Check(color)) {
                PyErr_SetString(PyExc_TypeError, "color must be int or tuple");
                return NULL;
            }
        } else {
            PyErr_SetString(
                PyExc_TypeError, "color must be int or single-element tuple");
            return NULL;
        }
    }

    switch (im->type) {
        case IMAGING_TYPE_UINT8:
            /* unsigned integer */
            if (im->bands == 1) {
                /* unsigned integer, single layer */
                if (rIsInt != 1) {
                    if (PyTuple_GET_SIZE(color) != 1) {
                        PyErr_SetString(PyExc_TypeError, "color must be int or single-element tuple");
                        return NULL;
                    } else if (!PyArg_ParseTuple(color, "L", &r)) {
                        return NULL;
                    }
                }
                ink[0] = (char)CLIP8(r);
                ink[1] = ink[2] = ink[3] = 0;
            } else {
                a = 255;
                if (rIsInt) {
                    /* compatibility: ABGR */
                    a = (UINT8)(r >> 24);
                    b = (UINT8)(r >> 16);
                    g = (UINT8)(r >> 8);
                    r = (UINT8)r;
                } else {
                    int tupleSize = PyTuple_GET_SIZE(color);
                    if (im->bands == 2) {
                        if (tupleSize != 1 && tupleSize != 2) {
                            PyErr_SetString(PyExc_TypeError, "color must be int, or tuple of one or two elements");
                            return NULL;
                        } else if (!PyArg_ParseTuple(color, "L|i", &r, &a)) {
                            return NULL;
                        }
                        g = b = r;
                    } else {
                        if (tupleSize != 3 && tupleSize != 4) {
                            PyErr_SetString(PyExc_TypeError, "color must be int, or tuple of one, three or four elements");
                            return NULL;
                        } else if (!PyArg_ParseTuple(color, "Lii|i", &r, &g, &b, &a)) {
                            return NULL;
                        }
                    }
                }
                ink[0] = (char)CLIP8(r);
                ink[1] = (char)CLIP8(g);
                ink[2] = (char)CLIP8(b);
                ink[3] = (char)CLIP8(a);
            }
            return ink;
        case IMAGING_TYPE_INT32:
            /* signed integer */
            itmp = r;
            memcpy(ink, &itmp, sizeof(itmp));
            return ink;
        case IMAGING_TYPE_FLOAT32:
            /* floating point */
            f = PyFloat_AsDouble(color);
            if (f == -1.0 && PyErr_Occurred()) {
                return NULL;
            }
            ftmp = f;
            memcpy(ink, &ftmp, sizeof(ftmp));
            return ink;
        case IMAGING_TYPE_SPECIAL:
            if (strncmp(im->mode, "I;16", 4) == 0) {
                ink[0] = (UINT8)r;
                ink[1] = (UINT8)(r >> 8);
                ink[2] = ink[3] = 0;
                return ink;
            }
    }

    PyErr_SetString(PyExc_ValueError, wrong_mode);
    return NULL;
}

static char *
hpy_getink(HPyContext *ctx, HPy color, Imaging im, char *ink) {
    PyObject *py_color = HPy_AsPyObject(ctx, color);
    char *res = getink(py_color, im, ink);
    Py_DECREF(py_color);
    return res;
}

/* -------------------------------------------------------------------- */
/* FACTORIES                                */
/* -------------------------------------------------------------------- */


HPyDef_METH(fill, "fill", HPyFunc_VARARGS)
static HPy fill_impl(HPyContext *ctx, HPy self, const HPy *args, size_t nargs) {
    char *mode;
    int xsize, ysize;
    HPy color, h_size;
    char buffer[4];
    Imaging im;

    xsize = ysize = 256;
    color = HPy_NULL;

    if (!HPyArg_Parse(ctx, NULL, args, nargs, "s|OO", &mode, &h_size, &color)) {
        return HPy_NULL;
    }

    xsize = HPy_GetLongItem_i(ctx, h_size, 0);
    ysize = HPy_GetLongItem_i(ctx, h_size, 1);

    im = ImagingNewDirty(mode, xsize, ysize);
    if (!im) {
        return HPy_NULL;
    }

    buffer[0] = buffer[1] = buffer[2] = buffer[3] = 0;
    if (!HPy_IsNull(color)) {
        if (!hpy_getink(ctx, color, im, buffer)) {
            ImagingDelete(im);
            return HPy_NULL;
        }
    }

    (void)ImagingFill(im, buffer);

    return HPyImagingNew(ctx, im);
}

HPyDef_METH(new, "new", HPyFunc_VARARGS)
static HPy new_impl(HPyContext *ctx, HPy self, const HPy *args, size_t nargs) {
    char *mode;
    int xsize, ysize;
    HPy h_size;

    if (!HPyArg_Parse(ctx, NULL, args, nargs, "sO", &mode, &h_size)) {
        return HPy_NULL;
    }

    xsize = HPy_GetLongItem_i(ctx, h_size, 0);
    ysize = HPy_GetLongItem_i(ctx, h_size, 1);

    return HPyImagingNew(ctx, ImagingNew(mode, xsize, ysize));
}

static PyObject *
_new_block(PyObject *self, PyObject *args) {
    char *mode;
    int xsize, ysize;

    if (!PyArg_ParseTuple(args, "s(ii)", &mode, &xsize, &ysize)) {
        return NULL;
    }

    return PyImagingNew(ImagingNewBlock(mode, xsize, ysize));
}

HPyDef_METH(linear_gradient, "linear_gradient", HPyFunc_VARARGS)
static HPy linear_gradient_impl(HPyContext *ctx, HPy self, const HPy *args, size_t nargs) {
    HPy h_mode;

    if (!HPyArg_Parse(ctx, NULL, args, nargs, "O", &h_mode)) {
        return HPy_NULL;
    }

    const char *mode = HPyUnicode_AsUTF8AndSize(ctx, h_mode, NULL);

    HPy hNew = HPyImagingNew(ctx, ImagingFillLinearGradient(mode));

    return hNew;
}

HPyDef_METH(radial_gradient, "radial_gradient", HPyFunc_VARARGS)
static HPy radial_gradient_impl(HPyContext *ctx, HPy self, const HPy *args, size_t nargs) {
    HPy h_mode;

    if (!HPyArg_Parse(ctx, NULL, args, nargs, "O", &h_mode)) {
        return HPy_NULL;
    }

    const char *mode = HPyUnicode_AsUTF8AndSize(ctx, h_mode, NULL);
    return HPyImagingNew(ctx, ImagingFillRadialGradient(mode));
}

HPyDef_METH(alpha_composite, "alpha_composite", HPyFunc_VARARGS)
static HPy alpha_composite_impl(HPyContext *ctx, HPy self, const HPy *args, size_t nargs) {
    HPy h_image1, h_image2;

    if (!HPyArg_Parse(
            ctx, NULL, args, nargs, "OO", &h_image1, &h_image2)) {
        return HPy_NULL;
    }

    ImagingObject *image1 = ImagingObject_AsStruct(ctx, h_image1);
    ImagingObject *image2 = ImagingObject_AsStruct(ctx, h_image2);

    HPy hNew = HPyImagingNew(ctx, ImagingAlphaComposite(image1->image, image2->image));

    return hNew;
}

HPyDef_METH(blend, "blend", HPyFunc_VARARGS)
static HPy blend_impl(HPyContext *ctx, HPy self, const HPy *args, size_t nargs){
    HPy h_image1, h_image2;
    double alpha;

    alpha = 0.5;
    if (!HPyArg_Parse(
            ctx, NULL, args, nargs, "OO|d", &h_image1, &h_image2, &alpha)) {
        return HPy_NULL;
    }

    ImagingObject *image1 = ImagingObject_AsStruct(ctx, h_image1);
    ImagingObject *image2 = ImagingObject_AsStruct(ctx, h_image2);

    HPy hNew = HPyImagingNew(ctx, ImagingBlend(image1->image, image2->image, (float)alpha));

    return hNew;
}

/* -------------------------------------------------------------------- */
/* METHODS                                                              */
/* -------------------------------------------------------------------- */

static INT16 *
_prepare_lut_table(PyObject *table, Py_ssize_t table_size) {
    int i;
    Py_buffer buffer_info;
    INT32 data_type = TYPE_FLOAT32;
    float item = 0;
    void *table_data = NULL;
    int free_table_data = 0;
    INT16 *prepared;

/* NOTE: This value should be the same as in ColorLUT.c */
#define PRECISION_BITS (16 - 8 - 2)

    const char *wrong_size =
        ("The table should have table_channels * "
         "size1D * size2D * size3D float items.");

    if (PyObject_CheckBuffer(table)) {
        if (!PyObject_GetBuffer(table, &buffer_info, PyBUF_CONTIG_RO | PyBUF_FORMAT)) {
            if (buffer_info.ndim == 1 && buffer_info.shape[0] == table_size) {
                if (strlen(buffer_info.format) == 1) {
                    switch (buffer_info.format[0]) {
                        case 'e':
                            data_type = TYPE_FLOAT16;
                            table_data = buffer_info.buf;
                            break;
                        case 'f':
                            data_type = TYPE_FLOAT32;
                            table_data = buffer_info.buf;
                            break;
                        case 'd':
                            data_type = TYPE_DOUBLE;
                            table_data = buffer_info.buf;
                            break;
                    }
                }
            }
            PyBuffer_Release(&buffer_info);
        }
    }

    if (!table_data) {
        free_table_data = 1;
        table_data = getlist(table, &table_size, wrong_size, TYPE_FLOAT32);
        if (!table_data) {
            return NULL;
        }
    }

    /* malloc check ok, max is 2 * 4 * 65**3 = 2197000 */
    prepared = (INT16 *)malloc(sizeof(INT16) * table_size);
    if (!prepared) {
        if (free_table_data) {
            free(table_data);
        }
        return (INT16 *)ImagingError_MemoryError();
    }

    for (i = 0; i < table_size; i++) {
        FLOAT16 htmp;
        double dtmp;
        switch (data_type) {
            case TYPE_FLOAT16:
                memcpy(&htmp, ((char *)table_data) + i * sizeof(htmp), sizeof(htmp));
                item = float16tofloat32(htmp);
                break;
            case TYPE_FLOAT32:
                memcpy(
                    &item, ((char *)table_data) + i * sizeof(FLOAT32), sizeof(FLOAT32));
                break;
            case TYPE_DOUBLE:
                memcpy(&dtmp, ((char *)table_data) + i * sizeof(dtmp), sizeof(dtmp));
                item = (FLOAT32)dtmp;
                break;
        }
        /* Max value for INT16 */
        if (item >= (0x7fff - 0.5) / (255 << PRECISION_BITS)) {
            prepared[i] = 0x7fff;
            continue;
        }
        /* Min value for INT16 */
        if (item <= (-0x8000 + 0.5) / (255 << PRECISION_BITS)) {
            prepared[i] = -0x8000;
            continue;
        }
        if (item < 0) {
            prepared[i] = item * (255 << PRECISION_BITS) - 0.5;
        } else {
            prepared[i] = item * (255 << PRECISION_BITS) + 0.5;
        }
    }

#undef PRECISION_BITS
    if (free_table_data) {
        free(table_data);
    }
    return prepared;
}

static PyObject *
_color_lut_3d(ImagingObject *self, PyObject *args) {
    char *mode;
    int filter;
    int table_channels;
    int size1D, size2D, size3D;
    PyObject *table;

    INT16 *prepared_table;
    Imaging imOut;

    if (!PyArg_ParseTuple(
            args,
            "siiiiiO:color_lut_3d",
            &mode,
            &filter,
            &table_channels,
            &size1D,
            &size2D,
            &size3D,
            &table)) {
        return NULL;
    }

    /* actually, it is trilinear */
    if (filter != IMAGING_TRANSFORM_BILINEAR) {
        PyErr_SetString(PyExc_ValueError, "Only LINEAR filter is supported.");
        return NULL;
    }

    if (1 > table_channels || table_channels > 4) {
        PyErr_SetString(PyExc_ValueError, "table_channels should be from 1 to 4");
        return NULL;
    }

    if (2 > size1D || size1D > 65 || 2 > size2D || size2D > 65 || 2 > size3D ||
        size3D > 65) {
        PyErr_SetString(
            PyExc_ValueError, "Table size in any dimension should be from 2 to 65");
        return NULL;
    }

    prepared_table =
        _prepare_lut_table(table, table_channels * size1D * size2D * size3D);
    if (!prepared_table) {
        return NULL;
    }

    imOut = ImagingNewDirty(mode, self->image->xsize, self->image->ysize);
    if (!imOut) {
        free(prepared_table);
        return NULL;
    }

    if (!ImagingColorLUT3D_linear(
            imOut,
            self->image,
            table_channels,
            size1D,
            size2D,
            size3D,
            prepared_table)) {
        free(prepared_table);
        ImagingDelete(imOut);
        return NULL;
    }

    free(prepared_table);

    return PyImagingNew(imOut);
}

HPyDef_METH(Imaging_convert, "convert", HPyFunc_VARARGS)
static HPy Imaging_convert_impl(HPyContext *ctx, HPy self, const HPy *args, size_t nargs) {
    HPy h_mode, h_paletteimage = HPy_NULL;
    int dither = 0;
    ImagingPalette palette = NULL;
    const char *mode;

    if (!HPyArg_Parse(ctx, NULL, args, nargs, "O|iO", &h_mode, &dither, &h_paletteimage)) {
        return HPy_NULL;
    }
    if (!HPy_IsNull(h_paletteimage)) {
        ImagingObject *paletteimage = ImagingObject_AsStruct(ctx, h_paletteimage);
        if (!PyImaging_Check(paletteimage)) {
            HPyErr_SetString(ctx, ctx->h_ValueError, "palette argument must be image with mode 'P'");
            return HPy_NULL;
        }
        if (paletteimage->image->palette == NULL) {
            HPyErr_SetString(ctx, ctx->h_ValueError, "null palette");
            return HPy_NULL;
        }
    }

    ImagingObject *im_obj = ImagingObject_AsStruct(ctx, self);
    Imaging im = im_obj->image;
    mode = HPyUnicode_AsUTF8AndSize(ctx, h_mode, NULL);

    return HPyImagingNew(ctx, ImagingConvert(im, mode, palette, dither));
}

static PyObject *
_convert2(ImagingObject *self, PyObject *args) {
    ImagingObject *imagep1;
    ImagingObject *imagep2;
    if (!PyArg_ParseTuple(
            args, "O!O!", Imaging_Type, &imagep1, Imaging_Type, &imagep2)) {
        return NULL;
    }

    if (!ImagingConvert2(imagep1->image, imagep2->image)) {
        return NULL;
    }

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
_convert_matrix(ImagingObject *self, PyObject *args) {
    char *mode;
    float m[12];
    if (!PyArg_ParseTuple(args, "s(ffff)", &mode, m + 0, m + 1, m + 2, m + 3)) {
        PyErr_Clear();
        if (!PyArg_ParseTuple(
                args,
                "s(ffffffffffff)",
                &mode,
                m + 0,
                m + 1,
                m + 2,
                m + 3,
                m + 4,
                m + 5,
                m + 6,
                m + 7,
                m + 8,
                m + 9,
                m + 10,
                m + 11)) {
            return NULL;
        }
    }

    return PyImagingNew(ImagingConvertMatrix(self->image, mode, m));
}

static PyObject *
_convert_transparent(ImagingObject *self, PyObject *args) {
    char *mode;
    int r, g, b;
    if (PyArg_ParseTuple(args, "s(iii)", &mode, &r, &g, &b)) {
        return PyImagingNew(ImagingConvertTransparent(self->image, mode, r, g, b));
    }
    PyErr_Clear();
    if (PyArg_ParseTuple(args, "si", &mode, &r)) {
        return PyImagingNew(ImagingConvertTransparent(self->image, mode, r, 0, 0));
    }
    return NULL;
}

HPyDef_METH(Imaging_copy, "copy", HPyFunc_NOARGS)
static HPy Imaging_copy_impl(HPyContext *ctx, HPy self) {
    ImagingObject *im_obj = ImagingObject_AsStruct(ctx, self);
    return HPyImagingNew(ctx, ImagingCopy(im_obj->image));
}


HPyDef_METH(Imaging_crop, "crop", HPyFunc_VARARGS)
static HPy Imaging_crop_impl(HPyContext *ctx, HPy self, const HPy *args, size_t nargs) {
    int x0, y0, x1, y1;
    HPy h_tuple;

    if (!HPyArg_Parse(ctx, NULL, args, nargs, "O", &h_tuple)) {
        return HPy_NULL;
    }

    x0 = HPy_GetLongItem_i(ctx, h_tuple, 0);
    y0 = HPy_GetLongItem_i(ctx, h_tuple, 1);
    x1 = HPy_GetLongItem_i(ctx, h_tuple, 2);
    y1 = HPy_GetLongItem_i(ctx, h_tuple, 3);
    
    ImagingObject *im_obj = ImagingObject_AsStruct(ctx, self);
    Imaging im = im_obj->image;

    return HPyImagingNew(ctx, ImagingCrop(im, x0, y0, x1, y1));
}

HPyDef_METH(Imaging_expand_image, "expand_image", HPyFunc_VARARGS)
static HPy Imaging_expand_image_impl(HPyContext *ctx, HPy self, const HPy *args, size_t nargs) {
    ImagingObject *im_obj = ImagingObject_AsStruct(ctx, self);
    Imaging im = im_obj->image;
    int x, y;
    int mode = 0;
    if (!HPyArg_Parse(ctx, NULL, args, nargs, "ii|i", &x, &y, &mode)) {
        return HPy_NULL;
    }
    return HPy_FromPyObject(ctx, PyImagingNew(ImagingExpand(im, x, y, mode)));
}

HPyDef_METH(Imaging_filter, "filter", HPyFunc_VARARGS)
static HPy Imaging_filter_impl(HPyContext *ctx, HPy self, const HPy *args, size_t nargs) {
    PyObject *imOut;
    Py_ssize_t kernelsize;
    FLOAT32 *kerneldata;

    ImagingObject *im_obj = ImagingObject_AsStruct(ctx, self);
    Imaging im = im_obj->image;

    int xsize, ysize, i;
    float divisor, offset;
    HPy kernel = HPy_NULL, h_size;

    if (!HPyArg_Parse(
            ctx, NULL, args, nargs, "OffO", &h_size, &divisor, &offset, &kernel)) {
        return HPy_NULL;
    }

    PyObject *py_kernel = HPy_AsPyObject(ctx, kernel);

    xsize = HPy_GetLongItem_i(ctx, h_size, 0);
    ysize = HPy_GetLongItem_i(ctx, h_size, 1);

    /* get user-defined kernel */
    kerneldata = getlist(py_kernel, &kernelsize, NULL, TYPE_FLOAT32);
    Py_DECREF(py_kernel);
    if (!kerneldata) {
        return HPy_NULL;
    }
    if (kernelsize != (Py_ssize_t)xsize * (Py_ssize_t)ysize) {
        HPy_Close(ctx, kernel);
        HPyErr_SetString(ctx, ctx->h_ValueError, "bad kernel size");
        return HPy_NULL;
    }

    for (i = 0; i < kernelsize; ++i) {
        kerneldata[i] /= divisor;
    }

    imOut = PyImagingNew(ImagingFilter(im, xsize, ysize, kerneldata, offset));

    free(kerneldata);

    return HPy_FromPyObject(ctx, imOut);
}

#ifdef WITH_UNSHARPMASK
HPyDef_METH(Imaging_gaussian_blur, "gaussian_blur", HPyFunc_VARARGS)
static HPy Imaging_gaussian_blur_impl(HPyContext *ctx, HPy self, const HPy *args, size_t nargs) {
    Imaging imIn;
    Imaging imOut;

    float radius = 0;
    int passes = 3;
    if (!HPyArg_Parse(ctx, NULL, args, nargs, "f|i", &radius, &passes)) {
        return HPy_NULL;
    }

    ImagingObject *im_obj = ImagingObject_AsStruct(ctx, self);
    imIn = im_obj->image;
    imOut = ImagingNewDirty(imIn->mode, imIn->xsize, imIn->ysize);
    if (!imOut) {
        return HPy_NULL;
    }

    if (!ImagingGaussianBlur(imOut, imIn, radius, passes)) {
        ImagingDelete(imOut);
        return HPy_NULL;
    }

    return HPy_FromPyObject(ctx, PyImagingNew(imOut));
}
#endif

HPyDef_METH(Imaging_getpalette, "getpalette", HPyFunc_VARARGS)
static HPy Imaging_getpalette_impl(HPyContext *ctx, HPy self, const HPy *args, size_t nargs) {
    HPy h_palette;
    int palettesize;
    int bits;
    ImagingShuffler pack;

    ImagingObject *im_obj = ImagingObject_AsStruct(ctx, self);
    Imaging im = im_obj->image;

    char *mode = "RGB";
    char *rawmode = "RGB";
    if (!HPyArg_Parse(ctx, NULL, args, nargs, "|ss", &mode, &rawmode)) {
        return HPy_NULL;
    }

    if (!im->palette) {
        HPyErr_SetString(ctx, ctx->h_ValueError, no_palette);
        return HPy_NULL;
    }

    pack = ImagingFindPacker(mode, rawmode, &bits);
    if (!pack) {
        HPyErr_SetString(ctx, ctx->h_ValueError, wrong_raw_mode);
        return HPy_NULL;
    }

    palettesize = im->palette->size;
    UINT8 *buf = (UINT8 *)malloc((palettesize * bits / 8) * sizeof(UINT8));
    if (buf == NULL) {
        return HPy_NULL;
    }

    pack(buf, im->palette->palette, palettesize);

    h_palette = HPyBytes_FromStringAndSize(ctx, (const char *)buf, palettesize * bits / 8);
    free(buf);
    return h_palette;
}

HPyDef_METH(Imaging_getpalettemode, "getpalettemode", HPyFunc_NOARGS)
static HPy Imaging_getpalettemode_impl(HPyContext *ctx, HPy self) {
    Imaging image = HPyImaging_AsImaging(ctx, self);

    if (!image->palette) {
        HPyErr_SetString(ctx, ctx->h_ValueError, no_palette);
        return HPy_NULL;
    }

    return HPyUnicode_FromString(ctx, image->palette->mode);
}

static inline int
_getxy(HPyContext *ctx, HPy xy, int *x, int *y) {
    HPy value;

    if (!HPyTuple_Check(ctx, xy) || HPy_Length(ctx, xy) != 2) {
        goto badarg;
    }

    value = HPy_GetItem_i(ctx, xy, 0);
    if (HPyLong_Check(ctx, value)) {
        *x = HPyLong_AsLong(ctx, value);
    } else if (HPyFloat_Check(ctx, value)) {
        *x = (int)HPyFloat_AsDouble(ctx, value);
    } else {
        HPy args[] = { value };
        HPy int_value = HPy_CallMethod_s(ctx, "__int__", args, 1, HPy_NULL);
        if (!HPy_IsNull(int_value) && HPyLong_Check(ctx, int_value)) {
            *x = HPyLong_AsLong(ctx, int_value);
        } else {
            HPy_Close(ctx, value);
            goto badval;
        }
    }
    HPy_Close(ctx, value);

    value = HPy_GetItem_i(ctx, xy, 1);
    if (HPyLong_Check(ctx, value)) {
        *y = HPyLong_AsLong(ctx, value);
    } else if (HPyFloat_Check(ctx, value)) {
        *y = (int)HPyFloat_AsDouble(ctx, value);
    } else {
        HPy args[] = { value };
        HPy int_value = HPy_CallMethod_s(ctx, "__int__", args, 1, HPy_NULL);
        if (!HPy_IsNull(int_value) && HPyLong_Check(ctx, int_value)) {
            *y = HPyLong_AsLong(ctx, int_value);
        } else {
            HPy_Close(ctx, value);
            goto badval;
        }
    }
    HPy_Close(ctx, value);

    return 0;

badarg:
    HPyErr_SetString(ctx, ctx->h_TypeError, "argument must be sequence of length 2");
    return -1;

badval:
    HPyErr_SetString(ctx, ctx->h_TypeError, "an integer is required");
    return -1;
}

HPyDef_METH(Imaging_getpixel, "getpixel", HPyFunc_VARARGS)
static HPy Imaging_getpixel_impl(HPyContext *ctx, HPy self, const HPy *args, size_t nargs) {
    HPy xy;
    int x, y;

    ImagingObject *im_self = ImagingObject_AsStruct(ctx, self);

    if (nargs != 1) {
        HPyErr_SetString(ctx, ctx->h_TypeError, "argument 1 must be sequence of length 2");
        return HPy_NULL;
    }

    if (!HPyArg_Parse(ctx, NULL, args, nargs, "O", &xy)) {
        return HPy_NULL;
    }

    if (_getxy(ctx, xy, &x, &y)) {
        return HPy_NULL;
    }

    if (im_self->access == NULL) {
        return HPy_Dup(ctx, ctx->h_None);
    }

    return getpixel(ctx, im_self->image, im_self->access, x, y);
}

union hist_extrema {
    UINT8 u[2];
    INT32 i[2];
    FLOAT32 f[2];
};

static union hist_extrema *
parse_histogram_extremap(
    ImagingObject *self, PyObject *extremap, union hist_extrema *ep) {
    int i0, i1;
    double f0, f1;

    if (extremap) {
        switch (self->image->type) {
            case IMAGING_TYPE_UINT8:
                if (!PyArg_ParseTuple(extremap, "ii", &i0, &i1)) {
                    return NULL;
                }
                ep->u[0] = CLIP8(i0);
                ep->u[1] = CLIP8(i1);
                break;
            case IMAGING_TYPE_INT32:
                if (!PyArg_ParseTuple(extremap, "ii", &i0, &i1)) {
                    return NULL;
                }
                ep->i[0] = i0;
                ep->i[1] = i1;
                break;
            case IMAGING_TYPE_FLOAT32:
                if (!PyArg_ParseTuple(extremap, "dd", &f0, &f1)) {
                    return NULL;
                }
                ep->f[0] = (FLOAT32)f0;
                ep->f[1] = (FLOAT32)f1;
                break;
            default:
                return NULL;
        }
    } else {
        return NULL;
    }
    return ep;
}

static PyObject *
_histogram(ImagingObject *self, PyObject *args) {
    ImagingHistogram h;
    PyObject *list;
    int i;
    union hist_extrema extrema;
    union hist_extrema *ep;

    PyObject *extremap = NULL;
    ImagingObject *maskp = NULL;
    if (!PyArg_ParseTuple(args, "|OO!", &extremap, Imaging_Type, &maskp)) {
        return NULL;
    }

    /* Using a var to avoid allocations. */
    ep = parse_histogram_extremap(self, extremap, &extrema);
    h = ImagingGetHistogram(self->image, (maskp) ? maskp->image : NULL, ep);

    if (!h) {
        return NULL;
    }

    /* Build an integer list containing the histogram */
    list = PyList_New(h->bands * 256);
    for (i = 0; i < h->bands * 256; i++) {
        PyObject *item;
        item = PyLong_FromLong(h->histogram[i]);
        if (item == NULL) {
            Py_DECREF(list);
            list = NULL;
            break;
        }
        PyList_SetItem(list, i, item);
    }

    /* Destroy the histogram structure */
    ImagingHistogramDelete(h);

    return list;
}

HPyDef_METH(Imaging_entropy, "entropy", HPyFunc_VARARGS)
static HPy Imaging_entropy_impl(HPyContext *ctx, HPy self, const HPy *args, size_t nargs) {
    ImagingHistogram h;
    int idx, length;
    long sum;
    double entropy, fsum, p;
    union hist_extrema extrema;
    union hist_extrema *ep;

    ImagingObject *im_obj = ImagingObject_AsStruct(ctx, self);
    Imaging im = im_obj->image;

    HPy h_extremap = HPy_NULL, h_maskp = HPy_NULL;

    PyObject *extremap = NULL;
    ImagingObject *maskp = NULL;
    if (!HPyArg_Parse(ctx, NULL, args, nargs, "|OO", &h_extremap, &h_maskp)) {
        return HPy_NULL;
    }

    extremap = HPy_AsPyObject(ctx, h_extremap);
    maskp = ImagingObject_AsStruct(ctx, h_maskp);

    /* Using a local var to avoid allocations. */
    ep = parse_histogram_extremap(im_obj, extremap, &extrema);
    h = ImagingGetHistogram(im, (maskp) ? maskp->image : NULL, ep);

    if (!h) {
        return HPy_NULL;
    }

    /* Calculate the histogram entropy */
    /* First, sum the histogram data */
    length = h->bands * 256;
    sum = 0;
    for (idx = 0; idx < length; idx++) {
        sum += h->histogram[idx];
    }

    /* Next, normalize the histogram data, */
    /* using the histogram sum value */
    fsum = (double)sum;
    entropy = 0.0;
    for (idx = 0; idx < length; idx++) {
        p = (double)h->histogram[idx] / fsum;
        if (p != 0.0) {
            entropy += p * log(p) * M_LOG2E;
        }
    }

    /* Destroy the histogram structure */
    ImagingHistogramDelete(h);

    return HPyFloat_FromDouble(ctx, -entropy);
}

#ifdef WITH_MODEFILTER
static PyObject *
_modefilter(ImagingObject *self, PyObject *args) {
    int size;
    if (!PyArg_ParseTuple(args, "i", &size)) {
        return NULL;
    }

    return PyImagingNew(ImagingModeFilter(self->image, size));
}
#endif

static PyObject *
_offset(ImagingObject *self, PyObject *args) {
    int xoffset, yoffset;
    if (!PyArg_ParseTuple(args, "ii", &xoffset, &yoffset)) {
        return NULL;
    }

    return PyImagingNew(ImagingOffset(self->image, xoffset, yoffset));
}

HPyDef_METH(Imaging_paste, "paste", HPyFunc_VARARGS)
static HPy Imaging_paste_impl(HPyContext *ctx, HPy self, const HPy *args, size_t nargs) {
    int status;
    char ink[4];

    int x0, y0, x1, y1;
    ImagingObject *maskp = NULL, *im_self = ImagingObject_AsStruct(ctx, self);

    HPy h_source=HPy_NULL, h_pixels=HPy_NULL, h_maskp=HPy_NULL;
    if (!HPyArg_Parse(
            ctx, NULL, args, nargs, "OO|O", &h_source, &h_pixels, &h_maskp)) {
        return HPy_NULL;
    }

    maskp = ImagingObject_AsStruct(ctx, h_maskp);

    x0 = HPy_GetLongItem_i(ctx, h_pixels, 0);
    y0 = HPy_GetLongItem_i(ctx, h_pixels, 1);
    x1 = HPy_GetLongItem_i(ctx, h_pixels, 2);
    y1 = HPy_GetLongItem_i(ctx, h_pixels, 3);

    if (HPyImaging_Check(ctx, h_source)) {
        status = ImagingPaste(
            im_self->image,
            HPyImaging_AsImaging(ctx, h_source),
            (maskp) ? maskp->image : NULL,
            x0,
            y0,
            x1,
            y1);

    } else {
        if (!hpy_getink(ctx, h_source, im_self->image, ink)) {
            return HPy_NULL;
        }
        status = ImagingFill2(
            im_self->image, ink, (maskp) ? maskp->image : NULL, x0, y0, x1, y1);
    }

    if (status < 0) {
        return HPy_NULL;
    }

    return HPy_Dup(ctx, ctx->h_None);
}

static PyObject *
_point(ImagingObject *self, PyObject *args) {
    static const char *wrong_number = "wrong number of lut entries";

    Py_ssize_t n;
    int i, bands;
    Imaging im;

    PyObject *list;
    char *mode;
    if (!PyArg_ParseTuple(args, "Oz", &list, &mode)) {
        return NULL;
    }

    if (mode && !strcmp(mode, "F")) {
        FLOAT32 *data;

        /* map from 8-bit data to floating point */
        n = 256;
        data = getlist(list, &n, wrong_number, TYPE_FLOAT32);
        if (!data) {
            return NULL;
        }
        im = ImagingPoint(self->image, mode, (void *)data);
        free(data);

    } else if (!strcmp(self->image->mode, "I") && mode && !strcmp(mode, "L")) {
        UINT8 *data;

        /* map from 16-bit subset of 32-bit data to 8-bit */
        /* FIXME: support arbitrary number of entries (requires API change) */
        n = 65536;
        data = getlist(list, &n, wrong_number, TYPE_UINT8);
        if (!data) {
            return NULL;
        }
        im = ImagingPoint(self->image, mode, (void *)data);
        free(data);

    } else {
        INT32 *data;
        UINT8 lut[1024];

        if (mode) {
            bands = getbands(mode);
            if (bands < 0) {
                return NULL;
            }
        } else {
            bands = self->image->bands;
        }

        /* map to integer data */
        n = 256 * bands;
        data = getlist(list, &n, wrong_number, TYPE_INT32);
        if (!data) {
            return NULL;
        }

        if (mode && !strcmp(mode, "I")) {
            im = ImagingPoint(self->image, mode, (void *)data);
        } else if (mode && bands > 1) {
            for (i = 0; i < 256; i++) {
                lut[i * 4] = CLIP8(data[i]);
                lut[i * 4 + 1] = CLIP8(data[i + 256]);
                lut[i * 4 + 2] = CLIP8(data[i + 512]);
                if (n > 768) {
                    lut[i * 4 + 3] = CLIP8(data[i + 768]);
                }
            }
            im = ImagingPoint(self->image, mode, (void *)lut);
        } else {
            /* map individual bands */
            for (i = 0; i < n; i++) {
                lut[i] = CLIP8(data[i]);
            }
            im = ImagingPoint(self->image, mode, (void *)lut);
        }
        free(data);
    }

    return PyImagingNew(im);
}

static PyObject *
_point_transform(ImagingObject *self, PyObject *args) {
    double scale = 1.0;
    double offset = 0.0;
    if (!PyArg_ParseTuple(args, "|dd", &scale, &offset)) {
        return NULL;
    }

    return PyImagingNew(ImagingPointTransform(self->image, scale, offset));
}

static PyObject *
_putdata(ImagingObject *self, PyObject *args) {
    Imaging image;
    // i & n are # pixels, require py_ssize_t. x can be as large as n. y, just because.
    Py_ssize_t n, i, x, y;

    PyObject *data;
    PyObject *seq = NULL;
    PyObject *op;
    double scale = 1.0;
    double offset = 0.0;

    if (!PyArg_ParseTuple(args, "O|dd", &data, &scale, &offset)) {
        return NULL;
    }

    if (!PySequence_Check(data)) {
        PyErr_SetString(PyExc_TypeError, must_be_sequence);
        return NULL;
    }

    image = self->image;

    n = PyObject_Length(data);
    if (n > (Py_ssize_t)image->xsize * (Py_ssize_t)image->ysize) {
        PyErr_SetString(PyExc_TypeError, "too many data entries");
        return NULL;
    }

#define set_value_to_item(seq, i) \
op = PySequence_Fast_GET_ITEM(seq, i); \
if (PySequence_Check(op)) { \
    PyErr_SetString(PyExc_TypeError, "sequence must be flattened"); \
    return NULL; \
} else { \
    value = PyFloat_AsDouble(op); \
}
    if (image->image8) {
        if (PyBytes_Check(data)) {
            unsigned char *p;
            p = (unsigned char *)PyBytes_AS_STRING(data);
            if (scale == 1.0 && offset == 0.0) {
                /* Plain string data */
                for (i = y = 0; i < n; i += image->xsize, y++) {
                    x = n - i;
                    if (x > (int)image->xsize) {
                        x = image->xsize;
                    }
                    memcpy(image->image8[y], p + i, x);
                }
            } else {
                /* Scaled and clipped string data */
                for (i = x = y = 0; i < n; i++) {
                    image->image8[y][x] = CLIP8((int)(p[i] * scale + offset));
                    if (++x >= (int)image->xsize) {
                        x = 0, y++;
                    }
                }
            }
        } else {
            seq = PySequence_Fast(data, must_be_sequence);
            if (!seq) {
                PyErr_SetString(PyExc_TypeError, must_be_sequence);
                return NULL;
            }
            double value;
            if (scale == 1.0 && offset == 0.0) {
                /* Clipped data */
                for (i = x = y = 0; i < n; i++) {
                    set_value_to_item(seq, i);
                    image->image8[y][x] = (UINT8)CLIP8(value);
                    if (++x >= (int)image->xsize) {
                        x = 0, y++;
                    }
                }

            } else {
                /* Scaled and clipped data */
                for (i = x = y = 0; i < n; i++) {
                    set_value_to_item(seq, i);
                    image->image8[y][x] = CLIP8(value * scale + offset);
                    if (++x >= (int)image->xsize) {
                        x = 0, y++;
                    }
                }
            }
            PyErr_Clear(); /* Avoid weird exceptions */
        }
    } else {
        /* 32-bit images */
        seq = PySequence_Fast(data, must_be_sequence);
        if (!seq) {
            PyErr_SetString(PyExc_TypeError, must_be_sequence);
            return NULL;
        }
        switch (image->type) {
            case IMAGING_TYPE_INT32:
                for (i = x = y = 0; i < n; i++) {
                    double value;
                    set_value_to_item(seq, i);
                    IMAGING_PIXEL_INT32(image, x, y) =
                        (INT32)(value * scale + offset);
                    if (++x >= (int)image->xsize) {
                        x = 0, y++;
                    }
                }
                PyErr_Clear(); /* Avoid weird exceptions */
                break;
            case IMAGING_TYPE_FLOAT32:
                for (i = x = y = 0; i < n; i++) {
                    double value;
                    set_value_to_item(seq, i);
                    IMAGING_PIXEL_FLOAT32(image, x, y) =
                        (FLOAT32)(value * scale + offset);
                    if (++x >= (int)image->xsize) {
                        x = 0, y++;
                    }
                }
                PyErr_Clear(); /* Avoid weird exceptions */
                break;
            default:
                for (i = x = y = 0; i < n; i++) {
                    union {
                        char ink[4];
                        INT32 inkint;
                    } u;

                    u.inkint = 0;

                    op = PySequence_Fast_GET_ITEM(seq, i);
                    if (!op || !getink(op, image, u.ink)) {
                        Py_DECREF(seq);
                        return NULL;
                    }
                    /* FIXME: what about scale and offset? */
                    image->image32[y][x] = u.inkint;
                    if (++x >= (int)image->xsize) {
                        x = 0, y++;
                    }
                }
                PyErr_Clear(); /* Avoid weird exceptions */
                break;
        }
    }

    Py_XDECREF(seq);

    Py_INCREF(Py_None);
    return Py_None;
}

#ifdef WITH_QUANTIZE

HPyDef_METH(Imaging_quantize, "quantize", HPyFunc_VARARGS)
static HPy Imaging_quantize_impl(HPyContext *ctx, HPy self, const HPy *args, size_t nargs) {
    int colours = 256;
    int method = 0;
    int kmeans = 0;
    if (!HPyArg_Parse(ctx, NULL, args, nargs, "|iii", &colours, &method, &kmeans)) {
        return HPy_NULL;
    }

    Imaging im = HPyImaging_AsImaging(ctx, self);

    if (!im->xsize || !im->ysize) {
        /* no content; return an empty image */
        return HPyImagingNew(ctx, ImagingNew("P", im->xsize, im->ysize));
    }

    return HPyImagingNew(ctx, ImagingQuantize(im, colours, method, kmeans));
}
#endif

static PyObject *
_putpalette(ImagingObject *self, PyObject *args) {
    ImagingShuffler unpack;
    int bits;

    char *rawmode, *palette_mode;
    UINT8 *palette;
    Py_ssize_t palettesize;
    if (!PyArg_ParseTuple(args, "sy#", &rawmode, &palette, &palettesize)) {
        return NULL;
    }

    if (strcmp(self->image->mode, "L") && strcmp(self->image->mode, "LA") &&
        strcmp(self->image->mode, "P") && strcmp(self->image->mode, "PA")) {
        PyErr_SetString(PyExc_ValueError, wrong_mode);
        return NULL;
    }

    palette_mode = strncmp("RGBA", rawmode, 4) == 0 ? "RGBA" : "RGB";
    unpack = ImagingFindUnpacker(palette_mode, rawmode, &bits);
    if (!unpack) {
        PyErr_SetString(PyExc_ValueError, wrong_raw_mode);
        return NULL;
    }

    if (palettesize * 8 / bits > 256) {
        PyErr_SetString(PyExc_ValueError, wrong_palette_size);
        return NULL;
    }

    ImagingPaletteDelete(self->image->palette);

    strcpy(self->image->mode, strlen(self->image->mode) == 2 ? "PA" : "P");

    self->image->palette = ImagingPaletteNew(palette_mode);

    self->image->palette->size = palettesize * 8 / bits;
    unpack(self->image->palette->palette, palette, self->image->palette->size);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
_putpalettealpha(ImagingObject *self, PyObject *args) {
    int index;
    int alpha = 0;
    if (!PyArg_ParseTuple(args, "i|i", &index, &alpha)) {
        return NULL;
    }

    if (!self->image->palette) {
        PyErr_SetString(PyExc_ValueError, no_palette);
        return NULL;
    }

    if (index < 0 || index >= 256) {
        PyErr_SetString(PyExc_ValueError, outside_palette);
        return NULL;
    }

    strcpy(self->image->palette->mode, "RGBA");
    self->image->palette->palette[index * 4 + 3] = (UINT8)alpha;

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
_putpalettealphas(ImagingObject *self, PyObject *args) {
    int i;
    UINT8 *values;
    Py_ssize_t length;
    if (!PyArg_ParseTuple(args, "y#", &values, &length)) {
        return NULL;
    }

    if (!self->image->palette) {
        PyErr_SetString(PyExc_ValueError, no_palette);
        return NULL;
    }

    if (length > 256) {
        PyErr_SetString(PyExc_ValueError, outside_palette);
        return NULL;
    }

    strcpy(self->image->palette->mode, "RGBA");
    for (i = 0; i < length; i++) {
        self->image->palette->palette[i * 4 + 3] = (UINT8)values[i];
    }

    Py_INCREF(Py_None);
    return Py_None;
}

HPyDef_METH(Imaging_putpixel, "putpixel", HPyFunc_VARARGS)
static HPy Imaging_putpixel_impl(HPyContext *ctx, HPy self, const HPy *args, size_t nargs) {
    Imaging im;
    char ink[4];

    int x, y;
    HPy h_coords, h_color=HPy_NULL;

    if (!HPyArg_Parse(ctx, NULL, args, nargs, "OO", &h_coords, &h_color)) {
        return HPy_NULL;
    }

    x = HPy_GetLongItem_i(ctx, h_coords, 0);
    y = HPy_GetLongItem_i(ctx, h_coords, 1);

    ImagingObject *im_self = ImagingObject_AsStruct(ctx, self);
    im = im_self->image;

    if (x < 0) {
        x = im->xsize + x;
    }
    if (y < 0) {
        y = im->ysize + y;
    }

    if (x < 0 || x >= im->xsize || y < 0 || y >= im->ysize) {
        HPyErr_SetString(ctx, ctx->h_IndexError, outside_image);
        return HPy_NULL;
    }

    if (!hpy_getink(ctx, h_color, im, ink)) {
        return HPy_NULL;
    }

    if (im_self->access) {
        im_self->access->put_pixel(im, x, y, ink);
    }

    return HPy_Dup(ctx, ctx->h_None);
}

#ifdef WITH_RANKFILTER
static PyObject *
_rankfilter(ImagingObject *self, PyObject *args) {
    int size, rank;
    if (!PyArg_ParseTuple(args, "ii", &size, &rank)) {
        return NULL;
    }

    return PyImagingNew(ImagingRankFilter(self->image, size, rank));
}
#endif

HPyDef_METH(Imaging_resize, "resize", HPyFunc_VARARGS)
static HPy Imaging_resize_impl(HPyContext *ctx, HPy self, const HPy *args, size_t nargs) {
    Imaging imIn;
    Imaging imOut;

    int xsize, ysize;
    int filter = IMAGING_TRANSFORM_NEAREST;
    float box[4] = {0, 0, 0, 0};

    ImagingObject *im_self = ImagingObject_AsStruct(ctx, self);
    imIn = im_self->image;
    box[2] = imIn->xsize;
    box[3] = imIn->ysize;

    HPy h_size, h_box;

    if (!HPyArg_Parse(ctx, NULL, args, nargs, "O|iO", &h_size, &filter, &h_box)) {
        return HPy_NULL;
    } 

    xsize = HPy_GetLongItem_i(ctx, h_size, 0);
    ysize = HPy_GetLongItem_i(ctx, h_size, 1);

    box[0] = HPy_GetLongItem_i(ctx, h_box, 0);
    box[1] = HPy_GetLongItem_i(ctx, h_box, 1);
    box[2] = HPy_GetLongItem_i(ctx, h_box, 2);
    box[3] = HPy_GetLongItem_i(ctx, h_box, 3);

    if (xsize < 1 || ysize < 1) {
        HPyErr_SetString(ctx, ctx->h_ValueError, "height and width must be > 0");
        return HPy_NULL;
    }

    if (box[0] < 0 || box[1] < 0) {
        HPyErr_SetString(ctx, ctx->h_ValueError, "box offset can't be negative");
        return HPy_NULL;
    }

    if (box[2] > imIn->xsize || box[3] > imIn->ysize) {
        HPyErr_SetString(ctx, ctx->h_ValueError, "box can't exceed original image size");
        return HPy_NULL;
    }

    if (box[2] - box[0] < 0 || box[3] - box[1] < 0) {
        HPyErr_SetString(ctx, ctx->h_ValueError, "box can't be empty");
        return HPy_NULL;
    }

    // If box's coordinates are int and box size matches requested size
    if (box[0] - (int)box[0] == 0 && box[2] - box[0] == xsize &&
        box[1] - (int)box[1] == 0 && box[3] - box[1] == ysize) {
        imOut = ImagingCrop(imIn, box[0], box[1], box[2], box[3]);
    } else if (filter == IMAGING_TRANSFORM_NEAREST) {
        double a[6];

        memset(a, 0, sizeof a);
        a[0] = (double)(box[2] - box[0]) / xsize;
        a[4] = (double)(box[3] - box[1]) / ysize;
        a[2] = box[0];
        a[5] = box[1];

        imOut = ImagingNewDirty(imIn->mode, xsize, ysize);

        imOut = ImagingTransform(
            imOut, imIn, IMAGING_TRANSFORM_AFFINE, 0, 0, xsize, ysize, a, filter, 1);
    } else {
        imOut = ImagingResample(imIn, xsize, ysize, filter, box);
    }

    return HPyImagingNew(ctx, imOut);
}

HPyDef_METH(Imaging_reduce, "reduce", HPyFunc_VARARGS)
static HPy Imaging_reduce_impl(HPyContext *ctx, HPy self, const HPy *args, size_t nargs) {
    Imaging imIn;
    Imaging imOut;

    int xscale, yscale;
    int box[4] = {0, 0, 0, 0};

    ImagingObject *im_obj = ImagingObject_AsStruct(ctx, self);
    imIn = im_obj->image;
    box[2] = imIn->xsize;
    box[3] = imIn->ysize;

    HPy h_scale, h_box;

    if (!HPyArg_Parse(ctx, NULL, args, nargs, "O|O", &h_scale, &h_box)) {
        return HPy_NULL;
    }

    xscale = HPy_GetLongItem_i(ctx, h_scale, 0);
    yscale = HPy_GetLongItem_i(ctx, h_scale, 1);

    box[0] = HPy_GetLongItem_i(ctx, h_box, 0);
    box[1] = HPy_GetLongItem_i(ctx, h_box, 1);
    box[2] = HPy_GetLongItem_i(ctx, h_box, 2);
    box[3] = HPy_GetLongItem_i(ctx, h_box, 3);

    if (xscale < 1 || yscale < 1) {
        HPyErr_SetString(ctx, ctx->h_ValueError, "scale must be > 0");
        return HPy_NULL;
    }

    if (box[0] < 0 || box[1] < 0) {
        HPyErr_SetString(ctx, ctx->h_ValueError, "box offset can't be negative");
        return HPy_NULL;
    }

    if (box[2] > imIn->xsize || box[3] > imIn->ysize) {
        HPyErr_SetString(ctx, ctx->h_ValueError, "box can't can't exceed original image size");
        return HPy_NULL;
    }

    if (box[2] <= box[0] || box[3] <= box[1]) {
        HPyErr_SetString(ctx, ctx->h_ValueError, "box can't be empty");
        return HPy_NULL;
    }

    if (xscale == 1 && yscale == 1) {
        imOut = ImagingCrop(imIn, box[0], box[1], box[2], box[3]);
    } else {
        // Change box format: (left, top, width, height)
        box[2] -= box[0];
        box[3] -= box[1];
        imOut = ImagingReduce(imIn, xscale, yscale, box);
    }

    return HPyImagingNew(ctx, imOut);
}

#define IS_RGB(mode) \
    (!strcmp(mode, "RGB") || !strcmp(mode, "RGBA") || !strcmp(mode, "RGBX"))

static PyObject *
im_setmode(ImagingObject *self, PyObject *args) {
    /* attempt to modify the mode of an image in place */

    Imaging im;

    char *mode;
    Py_ssize_t modelen;
    if (!PyArg_ParseTuple(args, "s#:setmode", &mode, &modelen)) {
        return NULL;
    }

    im = self->image;

    /* move all logic in here to the libImaging primitive */

    if (!strcmp(im->mode, mode)) {
        ; /* same mode; always succeeds */
    } else if (IS_RGB(im->mode) && IS_RGB(mode)) {
        /* color to color */
        strcpy(im->mode, mode);
        im->bands = modelen;
        if (!strcmp(mode, "RGBA")) {
            (void)ImagingFillBand(im, 3, 255);
        }
    } else {
        /* trying doing an in-place conversion */
        if (!ImagingConvertInPlace(im, mode)) {
            return NULL;
        }
    }

    if (self->access) {
        ImagingAccessDelete(im, self->access);
    }
    self->access = ImagingAccessNew(im);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
_transform2(ImagingObject *self, PyObject *args) {
    static const char *wrong_number = "wrong number of matrix entries";

    Imaging imOut;
    Py_ssize_t n;
    double *a;

    ImagingObject *imagep;
    int x0, y0, x1, y1;
    int method;
    PyObject *data;
    int filter = IMAGING_TRANSFORM_NEAREST;
    int fill = 1;
    if (!PyArg_ParseTuple(
            args,
            "(iiii)O!iO|ii",
            &x0,
            &y0,
            &x1,
            &y1,
            Imaging_Type,
            &imagep,
            &method,
            &data,
            &filter,
            &fill)) {
        return NULL;
    }

    switch (method) {
        case IMAGING_TRANSFORM_AFFINE:
            n = 6;
            break;
        case IMAGING_TRANSFORM_PERSPECTIVE:
            n = 8;
            break;
        case IMAGING_TRANSFORM_QUAD:
            n = 8;
            break;
        default:
            n = -1; /* force error */
    }

    a = getlist(data, &n, wrong_number, TYPE_DOUBLE);
    if (!a) {
        return NULL;
    }

    imOut = ImagingTransform(
        self->image, imagep->image, method, x0, y0, x1, y1, a, filter, fill);

    free(a);

    if (!imOut) {
        return NULL;
    }

    Py_INCREF(Py_None);
    return Py_None;
}

HPyDef_METH(Imaging_transpose, "transpose", HPyFunc_VARARGS)
static HPy Imaging_transpose_impl(HPyContext *ctx, HPy self, const HPy *args, size_t nargs) {
    Imaging imIn;
    Imaging imOut;

    int op;
    if (!HPyArg_Parse(ctx, NULL, args, nargs, "i", &op)) {
        return HPy_NULL;
    }

    imIn = ImagingObject_AsStruct(ctx, self)->image;

    switch (op) {
        case 0: /* flip left right */
        case 1: /* flip top bottom */
        case 3: /* rotate 180 */
            imOut = ImagingNewDirty(imIn->mode, imIn->xsize, imIn->ysize);
            break;
        case 2: /* rotate 90 */
        case 4: /* rotate 270 */
        case 5: /* transpose */
        case 6: /* transverse */
            imOut = ImagingNewDirty(imIn->mode, imIn->ysize, imIn->xsize);
            break;
        default:
            HPyErr_SetString(ctx, ctx->h_ValueError, "No such transpose operation");
            return HPy_NULL;
    }

    if (imOut) {
        switch (op) {
            case 0:
                (void)ImagingFlipLeftRight(imOut, imIn);
                break;
            case 1:
                (void)ImagingFlipTopBottom(imOut, imIn);
                break;
            case 2:
                (void)ImagingRotate90(imOut, imIn);
                break;
            case 3:
                (void)ImagingRotate180(imOut, imIn);
                break;
            case 4:
                (void)ImagingRotate270(imOut, imIn);
                break;
            case 5:
                (void)ImagingTranspose(imOut, imIn);
                break;
            case 6:
                (void)ImagingTransverse(imOut, imIn);
                break;
        }
    }

    return HPyImagingNew(ctx, imOut);
}

#ifdef WITH_UNSHARPMASK
static PyObject *
_unsharp_mask(ImagingObject *self, PyObject *args) {
    Imaging imIn;
    Imaging imOut;

    float radius;
    int percent, threshold;
    if (!PyArg_ParseTuple(args, "fii", &radius, &percent, &threshold)) {
        return NULL;
    }

    imIn = self->image;
    imOut = ImagingNewDirty(imIn->mode, imIn->xsize, imIn->ysize);
    if (!imOut) {
        return NULL;
    }

    if (!ImagingUnsharpMask(imOut, imIn, radius, percent, threshold)) {
        return NULL;
    }

    return PyImagingNew(imOut);
}
#endif

static PyObject *
_box_blur(ImagingObject *self, PyObject *args) {
    Imaging imIn;
    Imaging imOut;

    float radius;
    int n = 1;
    if (!PyArg_ParseTuple(args, "f|i", &radius, &n)) {
        return NULL;
    }

    imIn = self->image;
    imOut = ImagingNewDirty(imIn->mode, imIn->xsize, imIn->ysize);
    if (!imOut) {
        return NULL;
    }

    if (!ImagingBoxBlur(imOut, imIn, radius, n)) {
        ImagingDelete(imOut);
        return NULL;
    }

    return PyImagingNew(imOut);
}

/* -------------------------------------------------------------------- */

static PyObject *
_isblock(ImagingObject *self) {
    return PyBool_FromLong(self->image->block != NULL);
}

HPyDef_METH(Imaging_getbbox, "getbbox", HPyFunc_NOARGS)
static HPy Imaging_getbbox_impl(HPyContext *ctx, HPy self) {
    int bbox[4];

    ImagingObject *im_self = ImagingObject_AsStruct(ctx, self);
    Imaging im = im_self->image;

    if (!ImagingGetBBox(im, bbox)) {
        return HPy_Dup(ctx, ctx->h_None);
    }

    return HPy_BuildValue(ctx, "iiii",  bbox[0], bbox[1], bbox[2], bbox[3]);
}

HPyDef_METH(Imaging_getcolors, "getcolors", HPyFunc_VARARGS)
static HPy Imaging_getcolors_impl(HPyContext *ctx, HPy self, const HPy *args, size_t nargs) {
    ImagingColorItem *items;
    int i, colors;
    HPy h_out;

    ImagingObject *im_self = ImagingObject_AsStruct(ctx, self);

    int maxcolors = 256;
    if (!HPyArg_Parse(ctx, NULL, args, nargs, "i:getcolors", &maxcolors)) {
        return HPy_NULL;
    }

    items = ImagingGetColors(im_self->image, maxcolors, &colors);
    if (!items) {
        return HPy_NULL;
    }

    if (colors > maxcolors) {
        h_out = HPy_Dup(ctx, ctx->h_None);
    } else {
        HPyListBuilder builder = HPyListBuilder_New(ctx, colors);
        for (i = 0; i < colors; i++) {
            ImagingColorItem *v = &items[i];
            HPy h_item = HPy_BuildValue(
                ctx, "iN", v->count, getpixel(ctx, im_self->image, im_self->access, v->x, v->y));
            if (HPy_IsNull(h_item)) {
                HPyListBuilder_Cancel(ctx, builder);
            }
            HPyListBuilder_Set(ctx, builder, i, h_item);
            HPy_Close(ctx, h_item);
        }
        h_out = HPyListBuilder_Build(ctx, builder);
    }

    free(items);

    return h_out;
}

HPyDef_METH(Imaging_getextrema, "getextrema", HPyFunc_NOARGS)
static HPy Imaging_getextrema_impl(HPyContext *ctx, HPy self) {
    union {
        UINT8 u[2];
        INT32 i[2];
        FLOAT32 f[2];
        UINT16 s[2];
    } extrema;
    int status;
    unsigned s0, s1;

    ImagingObject *im_self = (ImagingObject *) HPy_AsPyObject(ctx, self);

    status = ImagingGetExtrema(im_self->image, &extrema);
    if (status < 0) {
        return HPy_NULL;
    }

    if (status) {
        switch (im_self->image->type) {
            case IMAGING_TYPE_UINT8:
                // return HPy_BuildValue(ctx, "BB", extrema.u[0], extrema.u[1]);
                s0 = extrema.s[0];
                s1 = extrema.s[1];
                return HPy_BuildValue(ctx, "II", s0, s1);
            case IMAGING_TYPE_INT32:
                return HPy_BuildValue(ctx, "ii", extrema.i[0], extrema.i[1]);
            case IMAGING_TYPE_FLOAT32:
                return HPy_BuildValue(ctx, "dd", extrema.f[0], extrema.f[1]);
            case IMAGING_TYPE_SPECIAL:
                if (strcmp(im_self->image->mode, "I;16") == 0) {
                    // return HPy_BuildValue(ctx, "HH", extrema.s[0], extrema.s[1]);
                    s0 = extrema.s[0];
                    s1 = extrema.s[1];
                    return HPy_BuildValue(ctx, "II", s0, s1);
                }
        }
    }

    return HPy_Dup(ctx, ctx->h_None);
}

static PyObject *
_getprojection(ImagingObject *self) {
    unsigned char *xprofile;
    unsigned char *yprofile;
    PyObject *result;

    /* malloc check ok */
    xprofile = malloc(self->image->xsize);
    yprofile = malloc(self->image->ysize);

    if (xprofile == NULL || yprofile == NULL) {
        free(xprofile);
        free(yprofile);
        return ImagingError_MemoryError();
    }

    ImagingGetProjection(
        self->image, (unsigned char *)xprofile, (unsigned char *)yprofile);

    result = Py_BuildValue(
        "y#y#",
        xprofile,
        (Py_ssize_t)self->image->xsize,
        yprofile,
        (Py_ssize_t)self->image->ysize);

    free(xprofile);
    free(yprofile);

    return result;
}

/* -------------------------------------------------------------------- */

HPyDef_METH(Imaging_getband, "getband", HPyFunc_VARARGS)
static HPy Imaging_getband_impl(HPyContext *ctx, HPy self, const HPy *args, size_t nargs) {
    ImagingObject *im_obj = ImagingObject_AsStruct(ctx, self);
    int band;

    if (!HPyArg_Parse(ctx, NULL, args, nargs, "i", &band)) {
        return HPy_NULL;
    }

    return HPyImagingNew(ctx, ImagingGetBand(im_obj->image, band));
}

HPyDef_METH(Imaging_fillband, "fillband", HPyFunc_VARARGS)
static HPy Imaging_fillband_impl(HPyContext *ctx, HPy self, const HPy *args, size_t nargs) {
    ImagingObject *im_obj = ImagingObject_AsStruct(ctx, self);
    int band;
    int color;

    if (!HPyArg_Parse(ctx, NULL, args, nargs, "ii", &band, &color)) {
        return HPy_NULL;
    }

    if (!ImagingFillBand(im_obj->image, band, color)) {
        return HPy_NULL;
    }

    return ctx->h_None;
}

HPyDef_METH(Imaging_putband, "putband", HPyFunc_VARARGS)
static HPy Imaging_putband_impl(HPyContext *ctx, HPy self, const HPy *args, size_t nargs) {
    ImagingObject *im_obj = ImagingObject_AsStruct(ctx, self);

    ImagingObject *imagep;
    int band;
    if (!HPyArg_Parse(ctx, NULL, args, nargs, "Oi", &imagep, &band)) {
        return HPy_NULL;
    }

    if (!ImagingPutBand(im_obj->image, imagep->image, band)) {
        return HPy_NULL;
    }

    return ctx->h_None;
}

HPyDef_METH(merge, "merge", HPyFunc_VARARGS)
static HPy merge_impl(HPyContext *ctx, HPy self, const HPy *args, size_t nargs)
{
    HPy h_mode, h_band0, h_band1, h_band2, h_band3;
    const char *mode;
    Imaging bands[4] = {NULL, NULL, NULL, NULL};

    if (!HPyArg_Parse(ctx, NULL, args, nargs, "sO|OOO", &h_mode, &h_band0, &h_band1, &h_band2, &h_band3)) {
        return HPy_NULL;
    }

    if (!HPy_IsNull(h_band0)) {
        bands[0] = ImagingObject_AsStruct(ctx, h_band0)->image;
    }
    if (!HPy_IsNull(h_band1)) {
        bands[1] = ImagingObject_AsStruct(ctx, h_band1)->image;
    }
    if (!HPy_IsNull(h_band2)) {
        bands[2] = ImagingObject_AsStruct(ctx, h_band2)->image;
    }
    if (!HPy_IsNull(h_band3)) {
        bands[3] = ImagingObject_AsStruct(ctx, h_band3)->image;
    }

    mode = HPyUnicode_AsUTF8AndSize(ctx, h_mode, NULL);
    return HPyImagingNew(ctx, ImagingMerge(mode, bands));
}

HPyDef_METH(Imaging_split, "split", HPyFunc_NOARGS)
static HPy Imaging_split_impl(HPyContext *ctx, HPy self) {
    HPy_ssize_t i;
    HPy h_New;
    Imaging bands[4] = {NULL, NULL, NULL, NULL};

    ImagingObject *im_self = ImagingObject_AsStruct(ctx, self);
    Imaging im = im_self->image;

    if (!ImagingSplit(im, bands)) {
        return HPy_NULL;
    }

    HPyTupleBuilder builder = HPyTupleBuilder_New(ctx, im->bands);
    for (i = 0; i < im->bands; i++) {
        h_New = HPyImagingNew(ctx, bands[i]);
        if (HPy_IsNull(h_New)) {
            HPyTupleBuilder_Cancel(ctx, builder);
            return HPy_NULL;
        }
        HPyTupleBuilder_Set(ctx, builder, i, h_New);
        HPy_Close(ctx, h_New);
    }
    return HPyTupleBuilder_Build(ctx, builder);
}

/* -------------------------------------------------------------------- */

#ifdef WITH_IMAGECHOPS

static PyObject *
_chop_invert(ImagingObject *self) {
    return PyImagingNew(ImagingNegative(self->image));
}

static PyObject *
_chop_lighter(ImagingObject *self, PyObject *args) {
    ImagingObject *imagep;

    if (!PyArg_ParseTuple(args, "O!", Imaging_Type, &imagep)) {
        return NULL;
    }

    return PyImagingNew(ImagingChopLighter(self->image, imagep->image));
}

static PyObject *
_chop_darker(ImagingObject *self, PyObject *args) {
    ImagingObject *imagep;

    if (!PyArg_ParseTuple(args, "O!", Imaging_Type, &imagep)) {
        return NULL;
    }

    return PyImagingNew(ImagingChopDarker(self->image, imagep->image));
}

static PyObject *
_chop_difference(ImagingObject *self, PyObject *args) {
    ImagingObject *imagep;

    if (!PyArg_ParseTuple(args, "O!", Imaging_Type, &imagep)) {
        return NULL;
    }

    return PyImagingNew(ImagingChopDifference(self->image, imagep->image));
}

static PyObject *
_chop_multiply(ImagingObject *self, PyObject *args) {
    ImagingObject *imagep;

    if (!PyArg_ParseTuple(args, "O!", Imaging_Type, &imagep)) {
        return NULL;
    }

    return PyImagingNew(ImagingChopMultiply(self->image, imagep->image));
}

static PyObject *
_chop_screen(ImagingObject *self, PyObject *args) {
    ImagingObject *imagep;

    if (!PyArg_ParseTuple(args, "O!", Imaging_Type, &imagep)) {
        return NULL;
    }

    return PyImagingNew(ImagingChopScreen(self->image, imagep->image));
}

static PyObject *
_chop_add(ImagingObject *self, PyObject *args) {
    ImagingObject *imagep;
    float scale;
    int offset;

    scale = 1.0;
    offset = 0;

    if (!PyArg_ParseTuple(args, "O!|fi", Imaging_Type, &imagep, &scale, &offset)) {
        return NULL;
    }

    return PyImagingNew(ImagingChopAdd(self->image, imagep->image, scale, offset));
}

static PyObject *
_chop_subtract(ImagingObject *self, PyObject *args) {
    ImagingObject *imagep;
    float scale;
    int offset;

    scale = 1.0;
    offset = 0;

    if (!PyArg_ParseTuple(args, "O!|fi", Imaging_Type, &imagep, &scale, &offset)) {
        return NULL;
    }

    return PyImagingNew(ImagingChopSubtract(self->image, imagep->image, scale, offset));
}

static PyObject *
_chop_and(ImagingObject *self, PyObject *args) {
    ImagingObject *imagep;

    if (!PyArg_ParseTuple(args, "O!", Imaging_Type, &imagep)) {
        return NULL;
    }

    return PyImagingNew(ImagingChopAnd(self->image, imagep->image));
}

static PyObject *
_chop_or(ImagingObject *self, PyObject *args) {
    ImagingObject *imagep;

    if (!PyArg_ParseTuple(args, "O!", Imaging_Type, &imagep)) {
        return NULL;
    }

    return PyImagingNew(ImagingChopOr(self->image, imagep->image));
}

static PyObject *
_chop_xor(ImagingObject *self, PyObject *args) {
    ImagingObject *imagep;

    if (!PyArg_ParseTuple(args, "O!", Imaging_Type, &imagep)) {
        return NULL;
    }

    return PyImagingNew(ImagingChopXor(self->image, imagep->image));
}

static PyObject *
_chop_add_modulo(ImagingObject *self, PyObject *args) {
    ImagingObject *imagep;

    if (!PyArg_ParseTuple(args, "O!", Imaging_Type, &imagep)) {
        return NULL;
    }

    return PyImagingNew(ImagingChopAddModulo(self->image, imagep->image));
}

static PyObject *
_chop_subtract_modulo(ImagingObject *self, PyObject *args) {
    ImagingObject *imagep;

    if (!PyArg_ParseTuple(args, "O!", Imaging_Type, &imagep)) {
        return NULL;
    }

    return PyImagingNew(ImagingChopSubtractModulo(self->image, imagep->image));
}

static PyObject *
_chop_soft_light(ImagingObject *self, PyObject *args) {
    ImagingObject *imagep;

    if (!PyArg_ParseTuple(args, "O!", Imaging_Type, &imagep)) {
        return NULL;
    }

    return PyImagingNew(ImagingChopSoftLight(self->image, imagep->image));
}

static PyObject *
_chop_hard_light(ImagingObject *self, PyObject *args) {
    ImagingObject *imagep;

    if (!PyArg_ParseTuple(args, "O!", Imaging_Type, &imagep)) {
        return NULL;
    }

    return PyImagingNew(ImagingChopHardLight(self->image, imagep->image));
}

static PyObject *
_chop_overlay(ImagingObject *self, PyObject *args) {
    ImagingObject *imagep;

    if (!PyArg_ParseTuple(args, "O!", Imaging_Type, &imagep)) {
        return NULL;
    }

    return PyImagingNew(ImagingOverlay(self->image, imagep->image));
}
#endif

/* -------------------------------------------------------------------- */

#ifdef WITH_IMAGEDRAW

static PyObject *
_font_new(PyObject *self_, PyObject *args) {
    ImagingFontObject *self;
    int i, y0, y1;
    static const char *wrong_length = "descriptor table has wrong size";

    ImagingObject *imagep;
    unsigned char *glyphdata;
    Py_ssize_t glyphdata_length;
    if (!PyArg_ParseTuple(
            args, "O!y#", Imaging_Type, &imagep, &glyphdata, &glyphdata_length)) {
        return NULL;
    }

    if (glyphdata_length != 256 * 20) {
        PyErr_SetString(PyExc_ValueError, wrong_length);
        return NULL;
    }

    self = PyObject_New(ImagingFontObject, &ImagingFont_Type);
    if (self == NULL) {
        return NULL;
    }

    /* glyph bitmap */
    self->bitmap = imagep->image;

    y0 = y1 = 0;

    /* glyph glyphs */
    for (i = 0; i < 256; i++) {
        self->glyphs[i].dx = S16(B16(glyphdata, 0));
        self->glyphs[i].dy = S16(B16(glyphdata, 2));
        self->glyphs[i].dx0 = S16(B16(glyphdata, 4));
        self->glyphs[i].dy0 = S16(B16(glyphdata, 6));
        self->glyphs[i].dx1 = S16(B16(glyphdata, 8));
        self->glyphs[i].dy1 = S16(B16(glyphdata, 10));
        self->glyphs[i].sx0 = S16(B16(glyphdata, 12));
        self->glyphs[i].sy0 = S16(B16(glyphdata, 14));
        self->glyphs[i].sx1 = S16(B16(glyphdata, 16));
        self->glyphs[i].sy1 = S16(B16(glyphdata, 18));
        if (self->glyphs[i].dy0 < y0) {
            y0 = self->glyphs[i].dy0;
        }
        if (self->glyphs[i].dy1 > y1) {
            y1 = self->glyphs[i].dy1;
        }
        glyphdata += 20;
    }

    self->baseline = -y0;
    self->ysize = y1 - y0;

    /* keep a reference to the bitmap object */
    Py_INCREF(imagep);
    self->ref = imagep;

    return (PyObject *)self;
}

static void
_font_dealloc(ImagingFontObject *self) {
    Py_XDECREF(self->ref);
    PyObject_Del(self);
}

static inline int
textwidth(ImagingFontObject *self, const unsigned char *text) {
    int xsize;

    for (xsize = 0; *text; text++) {
        xsize += self->glyphs[*text].dx;
    }

    return xsize;
}

void
_font_text_asBytes(PyObject *encoded_string, unsigned char **text) {
    /* Allocates *text, returns a 'new reference'. Caller is required to free */

    PyObject *bytes = NULL;
    Py_ssize_t len = 0;
    char *buffer;

    *text = NULL;

    if (PyUnicode_CheckExact(encoded_string)) {
        bytes = PyUnicode_AsLatin1String(encoded_string);
        if (!bytes) {
            return;
        }
        PyBytes_AsStringAndSize(bytes, &buffer, &len);
    } else if (PyBytes_Check(encoded_string)) {
        PyBytes_AsStringAndSize(encoded_string, &buffer, &len);
    }

    *text = calloc(len + 1, 1);
    if (*text) {
        memcpy(*text, buffer, len);
    } else {
        ImagingError_MemoryError();
    }
    if (bytes) {
        Py_DECREF(bytes);
    }

    return;
}

static PyObject *
_font_getmask(ImagingFontObject *self, PyObject *args) {
    Imaging im;
    Imaging bitmap;
    int x, b;
    int i = 0;
    int status;
    Glyph *glyph;

    PyObject *encoded_string;

    unsigned char *text;
    char *mode = "";

    if (!PyArg_ParseTuple(args, "O|s:getmask", &encoded_string, &mode)) {
        return NULL;
    }

    _font_text_asBytes(encoded_string, &text);
    if (!text) {
        return NULL;
    }

    im = ImagingNew(self->bitmap->mode, textwidth(self, text), self->ysize);
    if (!im) {
        free(text);
        return ImagingError_MemoryError();
    }

    b = 0;
    (void)ImagingFill(im, &b);

    b = self->baseline;
    for (x = 0; text[i]; i++) {
        glyph = &self->glyphs[text[i]];
        bitmap =
            ImagingCrop(self->bitmap, glyph->sx0, glyph->sy0, glyph->sx1, glyph->sy1);
        if (!bitmap) {
            goto failed;
        }
        status = ImagingPaste(
            im,
            bitmap,
            NULL,
            glyph->dx0 + x,
            glyph->dy0 + b,
            glyph->dx1 + x,
            glyph->dy1 + b);
        ImagingDelete(bitmap);
        if (status < 0) {
            goto failed;
        }
        x = x + glyph->dx;
        b = b + glyph->dy;
    }
    free(text);
    return PyImagingNew(im);

failed:
    free(text);
    ImagingDelete(im);
    Py_RETURN_NONE;
}

static PyObject *
_font_getsize(ImagingFontObject *self, PyObject *args) {
    unsigned char *text;
    PyObject *encoded_string;
    PyObject *val;

    if (!PyArg_ParseTuple(args, "O:getsize", &encoded_string)) {
        return NULL;
    }

    _font_text_asBytes(encoded_string, &text);
    if (!text) {
        return NULL;
    }

    val = Py_BuildValue("ii", textwidth(self, text), self->ysize);
    free(text);
    return val;
}

static struct PyMethodDef _font_methods[] = {
    {"getmask", (PyCFunction)_font_getmask, METH_VARARGS},
    {"getsize", (PyCFunction)_font_getsize, METH_VARARGS},
    {NULL, NULL} /* sentinel */
};

/* -------------------------------------------------------------------- */

static PyObject *
_draw_new(PyObject *self_, PyObject *args) {
    ImagingDrawObject *self;

    ImagingObject *imagep;
    int blend = 0;
    if (!PyArg_ParseTuple(args, "O!|i", Imaging_Type, &imagep, &blend)) {
        return NULL;
    }

    self = PyObject_New(ImagingDrawObject, &ImagingDraw_Type);
    if (self == NULL) {
        return NULL;
    }

    /* keep a reference to the image object */
    Py_INCREF(imagep);
    self->image = imagep;

    self->ink[0] = self->ink[1] = self->ink[2] = self->ink[3] = 0;

    self->blend = blend;

    return (PyObject *)self;
}

static void
_draw_dealloc(ImagingDrawObject *self) {
    Py_XDECREF(self->image);
    PyObject_Del(self);
}

extern Py_ssize_t
PyPath_Flatten(PyObject *data, double **xy);

static PyObject *
_draw_ink(ImagingDrawObject *self, PyObject *args) {
    INT32 ink = 0;
    PyObject *color;
    if (!PyArg_ParseTuple(args, "O", &color)) {
        return NULL;
    }

    if (!getink(color, self->image->image, (char *)&ink)) {
        return NULL;
    }

    return PyLong_FromLong((int)ink);
}

static PyObject *
_draw_arc(ImagingDrawObject *self, PyObject *args) {
    double *xy;
    Py_ssize_t n;

    PyObject *data;
    int ink;
    int width = 0;
    float start, end;
    if (!PyArg_ParseTuple(args, "Offi|i", &data, &start, &end, &ink, &width)) {
        return NULL;
    }

    n = PyPath_Flatten(data, &xy);
    if (n < 0) {
        return NULL;
    }
    if (n != 2) {
        PyErr_SetString(PyExc_TypeError, must_be_two_coordinates);
        free(xy);
        return NULL;
    }

    n = ImagingDrawArc(
        self->image->image,
        (int)xy[0],
        (int)xy[1],
        (int)xy[2],
        (int)xy[3],
        start,
        end,
        &ink,
        width,
        self->blend);

    free(xy);

    if (n < 0) {
        return NULL;
    }

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
_draw_bitmap(ImagingDrawObject *self, PyObject *args) {
    double *xy;
    Py_ssize_t n;

    PyObject *data;
    ImagingObject *bitmap;
    int ink;
    if (!PyArg_ParseTuple(args, "OO!i", &data, Imaging_Type, &bitmap, &ink)) {
        return NULL;
    }

    n = PyPath_Flatten(data, &xy);
    if (n < 0) {
        return NULL;
    }
    if (n != 1) {
        PyErr_SetString(
            PyExc_TypeError, "coordinate list must contain exactly 1 coordinate");
        free(xy);
        return NULL;
    }

    n = ImagingDrawBitmap(
        self->image->image, (int)xy[0], (int)xy[1], bitmap->image, &ink, self->blend);

    free(xy);

    if (n < 0) {
        return NULL;
    }

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
_draw_chord(ImagingDrawObject *self, PyObject *args) {
    double *xy;
    Py_ssize_t n;

    PyObject *data;
    int ink, fill;
    int width = 0;
    float start, end;
    if (!PyArg_ParseTuple(args, "Offii|i", &data, &start, &end, &ink, &fill, &width)) {
        return NULL;
    }

    n = PyPath_Flatten(data, &xy);
    if (n < 0) {
        return NULL;
    }
    if (n != 2) {
        PyErr_SetString(PyExc_TypeError, must_be_two_coordinates);
        free(xy);
        return NULL;
    }

    n = ImagingDrawChord(
        self->image->image,
        (int)xy[0],
        (int)xy[1],
        (int)xy[2],
        (int)xy[3],
        start,
        end,
        &ink,
        fill,
        width,
        self->blend);

    free(xy);

    if (n < 0) {
        return NULL;
    }

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
_draw_ellipse(ImagingDrawObject *self, PyObject *args) {
    double *xy;
    Py_ssize_t n;

    PyObject *data;
    int ink;
    int fill = 0;
    int width = 0;
    if (!PyArg_ParseTuple(args, "Oi|ii", &data, &ink, &fill, &width)) {
        return NULL;
    }

    n = PyPath_Flatten(data, &xy);
    if (n < 0) {
        return NULL;
    }
    if (n != 2) {
        PyErr_SetString(PyExc_TypeError, must_be_two_coordinates);
        free(xy);
        return NULL;
    }

    n = ImagingDrawEllipse(
        self->image->image,
        (int)xy[0],
        (int)xy[1],
        (int)xy[2],
        (int)xy[3],
        &ink,
        fill,
        width,
        self->blend);

    free(xy);

    if (n < 0) {
        return NULL;
    }

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
_draw_lines(ImagingDrawObject *self, PyObject *args) {
    double *xy;
    Py_ssize_t i, n;

    PyObject *data;
    int ink;
    int width = 0;
    if (!PyArg_ParseTuple(args, "Oi|i", &data, &ink, &width)) {
        return NULL;
    }

    n = PyPath_Flatten(data, &xy);
    if (n < 0) {
        return NULL;
    }

    if (width <= 1) {
        double *p = NULL;
        for (i = 0; i < n - 1; i++) {
            p = &xy[i + i];
            if (ImagingDrawLine(
                    self->image->image,
                    (int)p[0],
                    (int)p[1],
                    (int)p[2],
                    (int)p[3],
                    &ink,
                    self->blend) < 0) {
                free(xy);
                return NULL;
            }
        }
        if (p) { /* draw last point */
            ImagingDrawPoint(
                self->image->image, (int)p[2], (int)p[3], &ink, self->blend);
        }
    } else {
        for (i = 0; i < n - 1; i++) {
            double *p = &xy[i + i];
            if (ImagingDrawWideLine(
                    self->image->image,
                    (int)p[0],
                    (int)p[1],
                    (int)p[2],
                    (int)p[3],
                    &ink,
                    width,
                    self->blend) < 0) {
                free(xy);
                return NULL;
            }
        }
    }

    free(xy);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
_draw_points(ImagingDrawObject *self, PyObject *args) {
    double *xy;
    Py_ssize_t i, n;

    PyObject *data;
    int ink;
    if (!PyArg_ParseTuple(args, "Oi", &data, &ink)) {
        return NULL;
    }

    n = PyPath_Flatten(data, &xy);
    if (n < 0) {
        return NULL;
    }

    for (i = 0; i < n; i++) {
        double *p = &xy[i + i];
        if (ImagingDrawPoint(
                self->image->image, (int)p[0], (int)p[1], &ink, self->blend) < 0) {
            free(xy);
            return NULL;
        }
    }

    free(xy);

    Py_INCREF(Py_None);
    return Py_None;
}

#ifdef WITH_ARROW

/* from outline.c */
extern ImagingOutline
PyOutline_AsOutline(PyObject *outline);

static PyObject *
_draw_outline(ImagingDrawObject *self, PyObject *args) {
    ImagingOutline outline;

    PyObject *outline_;
    int ink;
    int fill = 0;
    if (!PyArg_ParseTuple(args, "Oi|i", &outline_, &ink, &fill)) {
        return NULL;
    }

    outline = PyOutline_AsOutline(outline_);
    if (!outline) {
        PyErr_SetString(PyExc_TypeError, "expected outline object");
        return NULL;
    }

    if (ImagingDrawOutline(self->image->image, outline, &ink, fill, self->blend) < 0) {
        return NULL;
    }

    Py_INCREF(Py_None);
    return Py_None;
}

#endif

static PyObject *
_draw_pieslice(ImagingDrawObject *self, PyObject *args) {
    double *xy;
    Py_ssize_t n;

    PyObject *data;
    int ink, fill;
    int width = 0;
    float start, end;
    if (!PyArg_ParseTuple(args, "Offii|i", &data, &start, &end, &ink, &fill, &width)) {
        return NULL;
    }

    n = PyPath_Flatten(data, &xy);
    if (n < 0) {
        return NULL;
    }
    if (n != 2) {
        PyErr_SetString(PyExc_TypeError, must_be_two_coordinates);
        free(xy);
        return NULL;
    }

    n = ImagingDrawPieslice(
        self->image->image,
        (int)xy[0],
        (int)xy[1],
        (int)xy[2],
        (int)xy[3],
        start,
        end,
        &ink,
        fill,
        width,
        self->blend);

    free(xy);

    if (n < 0) {
        return NULL;
    }

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
_draw_polygon(ImagingDrawObject *self, PyObject *args) {
    double *xy;
    int *ixy;
    Py_ssize_t n, i;

    PyObject *data;
    int ink;
    int fill = 0;
    int width = 0;
    if (!PyArg_ParseTuple(args, "Oi|ii", &data, &ink, &fill, &width)) {
        return NULL;
    }

    n = PyPath_Flatten(data, &xy);
    if (n < 0) {
        return NULL;
    }
    if (n < 2) {
        PyErr_SetString(
            PyExc_TypeError, "coordinate list must contain at least 2 coordinates");
        free(xy);
        return NULL;
    }

    /* Copy list of vertices to array */
    ixy = (int *)calloc(n, 2 * sizeof(int));
    if (ixy == NULL) {
        free(xy);
        return ImagingError_MemoryError();
    }

    for (i = 0; i < n; i++) {
        ixy[i + i] = (int)xy[i + i];
        ixy[i + i + 1] = (int)xy[i + i + 1];
    }

    free(xy);

    if (ImagingDrawPolygon(self->image->image, n, ixy, &ink, fill, width, self->blend) < 0) {
        free(ixy);
        return NULL;
    }

    free(ixy);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
_draw_rectangle(ImagingDrawObject *self, PyObject *args) {
    double *xy;
    Py_ssize_t n;

    PyObject *data;
    int ink;
    int fill = 0;
    int width = 0;
    if (!PyArg_ParseTuple(args, "Oi|ii", &data, &ink, &fill, &width)) {
        return NULL;
    }

    n = PyPath_Flatten(data, &xy);
    if (n < 0) {
        return NULL;
    }
    if (n != 2) {
        PyErr_SetString(PyExc_TypeError, must_be_two_coordinates);
        free(xy);
        return NULL;
    }

    n = ImagingDrawRectangle(
        self->image->image,
        (int)xy[0],
        (int)xy[1],
        (int)xy[2],
        (int)xy[3],
        &ink,
        fill,
        width,
        self->blend);

    free(xy);

    if (n < 0) {
        return NULL;
    }

    Py_INCREF(Py_None);
    return Py_None;
}

static struct PyMethodDef _draw_methods[] = {
#ifdef WITH_IMAGEDRAW
    /* Graphics (ImageDraw) */
    {"draw_lines", (PyCFunction)_draw_lines, METH_VARARGS},
#ifdef WITH_ARROW
    {"draw_outline", (PyCFunction)_draw_outline, METH_VARARGS},
#endif
    {"draw_polygon", (PyCFunction)_draw_polygon, METH_VARARGS},
    {"draw_rectangle", (PyCFunction)_draw_rectangle, METH_VARARGS},
    {"draw_points", (PyCFunction)_draw_points, METH_VARARGS},
    {"draw_arc", (PyCFunction)_draw_arc, METH_VARARGS},
    {"draw_bitmap", (PyCFunction)_draw_bitmap, METH_VARARGS},
    {"draw_chord", (PyCFunction)_draw_chord, METH_VARARGS},
    {"draw_ellipse", (PyCFunction)_draw_ellipse, METH_VARARGS},
    {"draw_pieslice", (PyCFunction)_draw_pieslice, METH_VARARGS},
    {"draw_ink", (PyCFunction)_draw_ink, METH_VARARGS},
#endif
    {NULL, NULL} /* sentinel */
};

#endif

HPyDef_METH(Imaging_pixel_access_new, "pixel_access", HPyFunc_VARARGS)
static HPy Imaging_pixel_access_new_impl(HPyContext *ctx, HPy self, const HPy *args, size_t nargs) {
    PixelAccessObject *self2;

    int readonly = 0;
    if (!HPyArg_Parse(ctx, NULL, args, nargs, "|i", &readonly)) {
        return HPy_NULL;
    }

    HPy h_PixelAccess_Type = HPyGlobal_Load(ctx, hg_PixelAccess_Type);

    HPy res = HPy_New(ctx, h_PixelAccess_Type, &self2);
    HPy_Close(ctx, h_PixelAccess_Type);
    if (HPy_IsNull(res)) {
        return HPy_NULL;
    }

    /* keep a reference to the image object */
    ImagingObject *imagep = ImagingObject_AsStruct(ctx, self);
    Py_INCREF(imagep);
    self2->image = imagep;

    self2->readonly = readonly;

    return res;
}

HPyDef_SLOT(pixel_access_dealloc, HPy_tp_destroy)
static void
pixel_access_dealloc_impl(void *data) {
    PixelAccessObject *self = (PixelAccessObject *)data;
    Py_XDECREF(self->image);
}

HPyDef_SLOT(pixel_access_getitem, HPy_mp_subscript)
static HPy
pixel_access_getitem_impl(HPyContext *ctx, HPy self, HPy xy) {
    int x, y;
    if (_getxy(ctx, xy, &x, &y)) {
        return HPy_NULL;
    }

    PixelAccessObject *data = PixelAccessObject_AsStruct(ctx, self);
    return getpixel(ctx, data->image->image, data->image->access, x, y);
}

HPyDef_SLOT(pixel_access_setitem, HPy_mp_ass_subscript)
static int
pixel_access_setitem_impl(HPyContext *ctx, HPy self, HPy xy, HPy color) {
    PixelAccessObject *data = PixelAccessObject_AsStruct(ctx, self);
    Imaging im = data->image->image;
    char ink[4];
    int x, y;

    if (data->readonly) {
        (void)ImagingError_ValueError(readonly);
        return -1;
    }

    if (_getxy(ctx, xy, &x, &y)) {
        return -1;
    }

    if (x < 0) {
        x = im->xsize + x;
    }
    if (y < 0) {
        y = im->ysize + y;
    }

    if (x < 0 || x >= im->xsize || y < 0 || y >= im->ysize) {
        HPyErr_SetString(ctx, ctx->h_IndexError, outside_image);
        return -1;
    }

    if (HPy_IsNull(color)) { /* FIXME: raise exception? */
        return 0;
    }

    if (!hpy_getink(ctx, color, im, ink)) {
        return -1;
    }

    data->image->access->put_pixel(im, x, y, ink);

    return 0;
}

/* -------------------------------------------------------------------- */
/* EFFECTS (experimental)                            */
/* -------------------------------------------------------------------- */

#ifdef WITH_EFFECTS

HPyDef_METH(effect_mandelbrot, "effect_mandelbrot", HPyFunc_VARARGS)
static HPy effect_mandelbrot_impl(HPyContext *ctx, HPy self, const HPy *args, size_t nargs) {
    int xsize = 512;
    int ysize = 512;
    double extent[4];
    int quality = 100;

    extent[0] = -3;
    extent[1] = -2.5;
    extent[2] = 2;
    extent[3] = 2.5;

    HPy h_size, h_extent;

    if (!HPyArg_Parse(ctx, NULL, args, nargs, "|OOi", &h_size, &h_extent, &quality)) {
        return HPy_NULL;
    }

    xsize = HPy_GetLongItem_i(ctx, h_size, 0);
    ysize = HPy_GetLongItem_i(ctx, h_size, 1);

    extent[0] = HPy_GetLongItem_i(ctx, h_extent, 0);
    extent[1] = HPy_GetLongItem_i(ctx, h_extent, 1);
    extent[2] = HPy_GetLongItem_i(ctx, h_extent, 2);
    extent[3] = HPy_GetLongItem_i(ctx, h_extent, 3);

    HPy hNew = HPyImagingNew(ctx, ImagingEffectMandelbrot(xsize, ysize, extent, quality));

    return hNew;
}

HPyDef_METH(effect_noise, "effect_noise", HPyFunc_VARARGS)
static HPy effect_noise_impl(HPyContext *ctx, HPy self, const HPy *args, size_t nargs) {
    int xsize, ysize;
    float sigma = 128;
    HPy h_size;
    if (!HPyArg_Parse(ctx, NULL, args, nargs, "O|f", &h_size, &sigma)) {
        return HPy_NULL;
    }

    xsize = HPy_GetLongItem_i(ctx, h_size, 0);
    ysize = HPy_GetLongItem_i(ctx, h_size, 1);

    HPy hNew = HPyImagingNew(ctx, ImagingEffectNoise(xsize, ysize, sigma));

    return hNew;
}

HPyDef_METH(Imaging_effect_spread, "effect_spread", HPyFunc_VARARGS)
static HPy Imaging_effect_spread_impl(HPyContext *ctx, HPy self, const HPy *args, size_t nargs) {
    int dist;
    ImagingObject *im_obj = ImagingObject_AsStruct(ctx, self);
    Imaging im = im_obj->image;

    if (!HPyArg_Parse(ctx, NULL, args, nargs, "i", &dist)) {
        return HPy_NULL;
    }

    HPy hNew = HPyImagingNew(ctx, ImagingEffectSpread(im, dist));

    return hNew;
}

#endif

/* -------------------------------------------------------------------- */
/* UTILITIES                                */
/* -------------------------------------------------------------------- */

HPyDef_METH(Imaging_getcodecstatus, "getcodecstatus", HPyFunc_VARARGS)
static HPy Imaging_getcodecstatus_impl(HPyContext *ctx, HPy self, const HPy *args, size_t nargs) {

    int status;
    char *msg;

    if (!HPyArg_Parse(ctx, NULL, args, nargs, "i", &status)) {
        return HPy_NULL;
    }

    switch (status) {
        case IMAGING_CODEC_OVERRUN:
            msg = "buffer overrun";
            break;
        case IMAGING_CODEC_BROKEN:
            msg = "broken data stream";
            break;
        case IMAGING_CODEC_UNKNOWN:
            msg = "unrecognized data stream contents";
            break;
        case IMAGING_CODEC_CONFIG:
            msg = "codec configuration error";
            break;
        case IMAGING_CODEC_MEMORY:
            msg = "out of memory";
            break;
        default:
            return HPy_Dup(ctx, ctx->h_None);
    }

    return HPyUnicode_FromString(ctx, msg);
}

/* -------------------------------------------------------------------- */
/* DEBUGGING HELPERS                            */
/* -------------------------------------------------------------------- */

static PyObject *
_save_ppm(ImagingObject *self, PyObject *args) {
    char *filename;

    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return NULL;
    }

    if (!ImagingSavePPM(self->image, filename)) {
        return NULL;
    }

    Py_INCREF(Py_None);
    return Py_None;
}

/* -------------------------------------------------------------------- */

/* methods */

static struct PyMethodDef methods[] = {

    /* Put commonly used methods first */

    /* Standard processing methods (Image) */
    {"color_lut_3d", (PyCFunction)_color_lut_3d, METH_VARARGS},
    {"convert2", (PyCFunction)_convert2, METH_VARARGS},
    {"convert_matrix", (PyCFunction)_convert_matrix, METH_VARARGS},
    {"convert_transparent", (PyCFunction)_convert_transparent, METH_VARARGS},
    {"histogram", (PyCFunction)_histogram, METH_VARARGS},
#ifdef WITH_MODEFILTER
    {"modefilter", (PyCFunction)_modefilter, METH_VARARGS},
#endif
    {"offset", (PyCFunction)_offset, METH_VARARGS},
    {"point", (PyCFunction)_point, METH_VARARGS},
    {"point_transform", (PyCFunction)_point_transform, METH_VARARGS},
    {"putdata", (PyCFunction)_putdata, METH_VARARGS},
#ifdef WITH_RANKFILTER
    {"rankfilter", (PyCFunction)_rankfilter, METH_VARARGS},
#endif
    {"transform2", (PyCFunction)_transform2, METH_VARARGS},

    {"isblock", (PyCFunction)_isblock, METH_NOARGS},
    {"getprojection", (PyCFunction)_getprojection, METH_NOARGS},

    {"setmode", (PyCFunction)im_setmode, METH_VARARGS},

    {"putpalette", (PyCFunction)_putpalette, METH_VARARGS},
    {"putpalettealpha", (PyCFunction)_putpalettealpha, METH_VARARGS},
    {"putpalettealphas", (PyCFunction)_putpalettealphas, METH_VARARGS},

#ifdef WITH_IMAGECHOPS
    /* Channel operations (ImageChops) */
    {"chop_invert", (PyCFunction)_chop_invert, METH_NOARGS},
    {"chop_lighter", (PyCFunction)_chop_lighter, METH_VARARGS},
    {"chop_darker", (PyCFunction)_chop_darker, METH_VARARGS},
    {"chop_difference", (PyCFunction)_chop_difference, METH_VARARGS},
    {"chop_multiply", (PyCFunction)_chop_multiply, METH_VARARGS},
    {"chop_screen", (PyCFunction)_chop_screen, METH_VARARGS},
    {"chop_add", (PyCFunction)_chop_add, METH_VARARGS},
    {"chop_subtract", (PyCFunction)_chop_subtract, METH_VARARGS},
    {"chop_add_modulo", (PyCFunction)_chop_add_modulo, METH_VARARGS},
    {"chop_subtract_modulo", (PyCFunction)_chop_subtract_modulo, METH_VARARGS},
    {"chop_and", (PyCFunction)_chop_and, METH_VARARGS},
    {"chop_or", (PyCFunction)_chop_or, METH_VARARGS},
    {"chop_xor", (PyCFunction)_chop_xor, METH_VARARGS},
    {"chop_soft_light", (PyCFunction)_chop_soft_light, METH_VARARGS},
    {"chop_hard_light", (PyCFunction)_chop_hard_light, METH_VARARGS},
    {"chop_overlay", (PyCFunction)_chop_overlay, METH_VARARGS},

#endif

#ifdef WITH_UNSHARPMASK
    /* Kevin Cazabon's unsharpmask extension */
    {"unsharp_mask", (PyCFunction)_unsharp_mask, METH_VARARGS},
#endif

    {"box_blur", (PyCFunction)_box_blur, METH_VARARGS},

    /* Misc. */
    {"new_block", (PyCFunction)_new_block, METH_VARARGS},

    {"save_ppm", (PyCFunction)_save_ppm, METH_VARARGS},

    {NULL, NULL} /* sentinel */
};



/* attributes */

HPyDef_GET(Imaging_getattr_mode, "mode")
static HPy Imaging_getattr_mode_get(HPyContext *ctx, HPy self, void *closure) {
    return HPyUnicode_FromString(ctx, ImagingObject_AsStruct(ctx, self)->image->mode);
}

HPyDef_GET(Imaging_getattr_size, "size")
static HPy Imaging_getattr_size_get(HPyContext *ctx, HPy self, void *closure) {
    Imaging image = ImagingObject_AsStruct(ctx, self)->image;
    return HPy_BuildValue(ctx, "ii", image->xsize, image->ysize);
}

HPyDef_GET(Imaging_getattr_bands, "bands")
static HPy Imaging_getattr_bands_get(HPyContext *ctx, HPy self, void *closure) {
    return HPyLong_FromLong(ctx, ((Imaging) ImagingObject_AsStruct(ctx, self))->bands);
}

HPyDef_GET(Imaging_getattr_id, "id")
static HPy Imaging_getattr_id_get(HPyContext *ctx, HPy self, void *closure) {
    return HPyLong_FromSize_t(ctx, (HPy_ssize_t)ImagingObject_AsStruct(ctx, self)->image);
}

HPyDef_GET(Imaging_getattr_ptr, "ptr")
static HPy Imaging_getattr_ptr_get(HPyContext *ctx, HPy self, void *closure) {
    return HPyCapsule_New(ctx, ImagingObject_AsStruct(ctx, self)->image, IMAGING_MAGIC, NULL);
}

HPyDef_GET(Imaging_getattr_unsafe_ptrs, "unsafe_ptrs")
static HPy Imaging_getattr_unsafe_ptrs_get(HPyContext *ctx, HPy self, void *closure) {
    Imaging im = ImagingObject_AsStruct(ctx, self)->image;
    return HPy_BuildValue(ctx,
            "(sn)(sn)(sn)",
            "image8", im->image8,
            "image32", im->image32,
            "image", im->image8);
};

/* basic sequence semantics */

HPyDef_SLOT(Imaging_image_length, HPy_sq_length)
static HPy_ssize_t Imaging_image_length_impl(HPyContext *ctx, HPy self) {
    Imaging im = ImagingObject_AsStruct(ctx, self)->image;
    return (HPy_ssize_t)im->xsize * im->ysize;
}

HPyDef_SLOT(Imaging_image_item, HPy_sq_item)
static HPy Imaging_image_item_impl(HPyContext *ctx, HPy self, HPy_ssize_t i) {
    int x, y;
    ImagingObject *im_obj = ImagingObject_AsStruct(ctx, self);
    Imaging im = im_obj->image;

    if (im->xsize > 0) {
        x = i % im->xsize;
        y = i / im->xsize;
    } else {
        x = y = 0; /* leave it to getpixel to raise an exception */
    }

    return getpixel(ctx, im, im_obj->access, x, y);
}

/* type description */

static PyType_Slot Imaging_Type_slots[] = {
    {Py_tp_methods, methods},
    {0, NULL}
};

static HPyDef *Imaging_type_defines[]={

    &Imaging_destroy,
    &Imaging_image_length,
    &Imaging_image_item,

    &Imaging_getattr_mode,
    &Imaging_getattr_size,
    &Imaging_getattr_bands,
    &Imaging_getattr_id,
    &Imaging_getattr_ptr,
    &Imaging_getattr_unsafe_ptrs,

    /* Put commonly used methods first */
    &Imaging_pixel_access_new,
    &Imaging_getpixel,
    &Imaging_putpixel,

    /* Standard processing methods (Image) */
    &Imaging_convert,
    &Imaging_copy,
    &Imaging_crop,
    &Imaging_entropy,
    &Imaging_expand_image,
    &Imaging_filter,
    &Imaging_getpalettemode,

    &Imaging_paste,

#ifdef WITH_QUANTIZE
    &Imaging_quantize,
#endif

    &Imaging_resize,
    &Imaging_reduce,

    &Imaging_getextrema,

    &Imaging_getbbox,
    &Imaging_getcolors,
    &Imaging_split,
    &Imaging_transpose,

    &Imaging_getband,
    &Imaging_fillband,
    &Imaging_putband,

    &Imaging_getpalette,

    /* Utilities */
    &Imaging_getcodecstatus,

    /* Kevin Cazabon's unsharpmask extension */
    &Imaging_gaussian_blur,


#ifdef WITH_EFFECTS
    /* Special effects */
    &Imaging_effect_spread,
#endif

    NULL
};

HPyType_Spec Imaging_Type_spec = {
    .name = "ImagingCoreOriginal",
    .basicsize = sizeof(ImagingObject),
    .flags = (HPy_TPFLAGS_DEFAULT | HPy_TPFLAGS_BASETYPE),
    .builtin_shape = SHAPE(ImagingObject),
    .legacy_slots = Imaging_Type_slots,
    .defines = Imaging_type_defines,
};

#ifdef WITH_IMAGEDRAW

static PyTypeObject ImagingFont_Type = {
    PyVarObject_HEAD_INIT(NULL, 0) "ImagingFont", /*tp_name*/
    sizeof(ImagingFontObject),                    /*tp_size*/
    0,                                            /*tp_itemsize*/
    /* methods */
    (destructor)_font_dealloc, /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number */
    0,                         /*tp_as_sequence */
    0,                         /*tp_as_mapping */
    0,                         /*tp_hash*/
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT,        /*tp_flags*/
    0,                         /*tp_doc*/
    0,                         /*tp_traverse*/
    0,                         /*tp_clear*/
    0,                         /*tp_richcompare*/
    0,                         /*tp_weaklistoffset*/
    0,                         /*tp_iter*/
    0,                         /*tp_iternext*/
    _font_methods,             /*tp_methods*/
    0,                         /*tp_members*/
    0,                         /*tp_getset*/
};

static PyTypeObject ImagingDraw_Type = {
    PyVarObject_HEAD_INIT(NULL, 0) "ImagingDraw", /*tp_name*/
    sizeof(ImagingDrawObject),                    /*tp_size*/
    0,                                            /*tp_itemsize*/
    /* methods */
    (destructor)_draw_dealloc, /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number */
    0,                         /*tp_as_sequence */
    0,                         /*tp_as_mapping */
    0,                         /*tp_hash*/
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT,        /*tp_flags*/
    0,                         /*tp_doc*/
    0,                         /*tp_traverse*/
    0,                         /*tp_clear*/
    0,                         /*tp_richcompare*/
    0,                         /*tp_weaklistoffset*/
    0,                         /*tp_iter*/
    0,                         /*tp_iternext*/
    _draw_methods,             /*tp_methods*/
    0,                         /*tp_members*/
    0,                         /*tp_getset*/
};

#endif

/* type description */

static HPyDef *PixelAccess_defines[] = {
    &pixel_access_dealloc,
    &pixel_access_getitem,
    &pixel_access_setitem,
    NULL
};

static HPyType_Spec PixelAccess_Type_spec = {
    .name = "PixelAccess",
    .basicsize = sizeof(PixelAccessObject),
    .builtin_shape = SHAPE(PixelAccessObject),
    .defines = PixelAccess_defines,
};

/* -------------------------------------------------------------------- */

static PyObject *
_get_stats(PyObject *self, PyObject *args) {
    PyObject *d;
    ImagingMemoryArena arena = &ImagingDefaultArena;

    if (!PyArg_ParseTuple(args, ":get_stats")) {
        return NULL;
    }

    d = PyDict_New();
    if (!d) {
        return NULL;
    }
    PyDict_SetItemString(d, "new_count", PyLong_FromLong(arena->stats_new_count));
    PyDict_SetItemString(
        d, "allocated_blocks", PyLong_FromLong(arena->stats_allocated_blocks));
    PyDict_SetItemString(
        d, "reused_blocks", PyLong_FromLong(arena->stats_reused_blocks));
    PyDict_SetItemString(
        d, "reallocated_blocks", PyLong_FromLong(arena->stats_reallocated_blocks));
    PyDict_SetItemString(d, "freed_blocks", PyLong_FromLong(arena->stats_freed_blocks));
    PyDict_SetItemString(d, "blocks_cached", PyLong_FromLong(arena->blocks_cached));
    return d;
}

static PyObject *
_reset_stats(PyObject *self, PyObject *args) {
    ImagingMemoryArena arena = &ImagingDefaultArena;

    if (!PyArg_ParseTuple(args, ":reset_stats")) {
        return NULL;
    }

    arena->stats_new_count = 0;
    arena->stats_allocated_blocks = 0;
    arena->stats_reused_blocks = 0;
    arena->stats_reallocated_blocks = 0;
    arena->stats_freed_blocks = 0;

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
_get_alignment(PyObject *self, PyObject *args) {
    if (!PyArg_ParseTuple(args, ":get_alignment")) {
        return NULL;
    }

    return PyLong_FromLong(ImagingDefaultArena.alignment);
}

static PyObject *
_get_block_size(PyObject *self, PyObject *args) {
    if (!PyArg_ParseTuple(args, ":get_block_size")) {
        return NULL;
    }

    return PyLong_FromLong(ImagingDefaultArena.block_size);
}

static PyObject *
_get_blocks_max(PyObject *self, PyObject *args) {
    if (!PyArg_ParseTuple(args, ":get_blocks_max")) {
        return NULL;
    }

    return PyLong_FromLong(ImagingDefaultArena.blocks_max);
}

static PyObject *
_set_alignment(PyObject *self, PyObject *args) {
    int alignment;
    if (!PyArg_ParseTuple(args, "i:set_alignment", &alignment)) {
        return NULL;
    }

    if (alignment < 1 || alignment > 128) {
        PyErr_SetString(PyExc_ValueError, "alignment should be from 1 to 128");
        return NULL;
    }
    /* Is power of two */
    if (alignment & (alignment - 1)) {
        PyErr_SetString(PyExc_ValueError, "alignment should be power of two");
        return NULL;
    }

    ImagingDefaultArena.alignment = alignment;

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
_set_block_size(PyObject *self, PyObject *args) {
    int block_size;
    if (!PyArg_ParseTuple(args, "i:set_block_size", &block_size)) {
        return NULL;
    }

    if (block_size <= 0) {
        PyErr_SetString(PyExc_ValueError, "block_size should be greater than 0");
        return NULL;
    }

    if (block_size & 0xfff) {
        PyErr_SetString(PyExc_ValueError, "block_size should be multiple of 4096");
        return NULL;
    }

    ImagingDefaultArena.block_size = block_size;

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
_set_blocks_max(PyObject *self, PyObject *args) {
    int blocks_max;
    if (!PyArg_ParseTuple(args, "i:set_blocks_max", &blocks_max)) {
        return NULL;
    }

    if (blocks_max < 0) {
        PyErr_SetString(PyExc_ValueError, "blocks_max should be greater than 0");
        return NULL;
    } else if (
        (unsigned long)blocks_max >
        SIZE_MAX / sizeof(ImagingDefaultArena.blocks_pool[0])) {
        PyErr_SetString(PyExc_ValueError, "blocks_max is too large");
        return NULL;
    }

    if (!ImagingMemorySetBlocksMax(&ImagingDefaultArena, blocks_max)) {
        return ImagingError_MemoryError();
    }

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
_clear_cache(PyObject *self, PyObject *args) {
    int i = 0;

    if (!PyArg_ParseTuple(args, "|i:clear_cache", &i)) {
        return NULL;
    }

    ImagingMemoryClearCache(&ImagingDefaultArena, i);

    Py_INCREF(Py_None);
    return Py_None;
}

/* -------------------------------------------------------------------- */

/* FIXME: this is something of a mess.  Should replace this with
   pluggable codecs, but not before PIL 1.2 */

/* Decoders (in decode.c) */
extern PyObject *
PyImaging_BcnDecoderNew(PyObject *self, PyObject *args);
extern PyObject *
PyImaging_BitDecoderNew(PyObject *self, PyObject *args);
extern PyObject *
PyImaging_FliDecoderNew(PyObject *self, PyObject *args);
extern PyObject *
PyImaging_GifDecoderNew(PyObject *self, PyObject *args);
extern PyObject *
PyImaging_HexDecoderNew(PyObject *self, PyObject *args);
extern PyObject *
PyImaging_JpegDecoderNew(PyObject *self, PyObject *args);
extern PyObject *
PyImaging_Jpeg2KDecoderNew(PyObject *self, PyObject *args);
extern PyObject *
PyImaging_LibTiffDecoderNew(PyObject *self, PyObject *args);
extern PyObject *
PyImaging_PackbitsDecoderNew(PyObject *self, PyObject *args);
extern PyObject *
PyImaging_PcdDecoderNew(PyObject *self, PyObject *args);
extern PyObject *
PyImaging_PcxDecoderNew(PyObject *self, PyObject *args);
extern PyObject *
PyImaging_RawDecoderNew(PyObject *self, PyObject *args);
extern PyObject *
PyImaging_SgiRleDecoderNew(PyObject *self, PyObject *args);
extern PyObject *
PyImaging_SunRleDecoderNew(PyObject *self, PyObject *args);
extern PyObject *
PyImaging_TgaRleDecoderNew(PyObject *self, PyObject *args);
extern PyObject *
PyImaging_XbmDecoderNew(PyObject *self, PyObject *args);
extern PyObject *
PyImaging_ZipDecoderNew(PyObject *self, PyObject *args);

/* Encoders (in encode.c) */
extern PyObject *
PyImaging_EpsEncoderNew(PyObject *self, PyObject *args);
extern PyObject *
PyImaging_GifEncoderNew(PyObject *self, PyObject *args);
extern PyObject *
PyImaging_JpegEncoderNew(PyObject *self, PyObject *args);
extern PyObject *
PyImaging_Jpeg2KEncoderNew(PyObject *self, PyObject *args);
extern PyObject *
PyImaging_PcxEncoderNew(PyObject *self, PyObject *args);
extern PyObject *
PyImaging_RawEncoderNew(PyObject *self, PyObject *args);
extern PyObject *
PyImaging_TgaRleEncoderNew(PyObject *self, PyObject *args);
extern PyObject *
PyImaging_XbmEncoderNew(PyObject *self, PyObject *args);
extern PyObject *
PyImaging_ZipEncoderNew(PyObject *self, PyObject *args);
extern PyObject *
PyImaging_LibTiffEncoderNew(PyObject *self, PyObject *args);

/* Display support etc (in display.c) */
#ifdef _WIN32
extern PyObject *
PyImaging_CreateWindowWin32(PyObject *self, PyObject *args);
extern PyObject *
PyImaging_DisplayWin32(PyObject *self, PyObject *args);
extern PyObject *
PyImaging_DisplayModeWin32(PyObject *self, PyObject *args);
extern PyObject *
PyImaging_GrabScreenWin32(PyObject *self, PyObject *args);
extern PyObject *
PyImaging_GrabClipboardWin32(PyObject *self, PyObject *args);
extern PyObject *
PyImaging_ListWindowsWin32(PyObject *self, PyObject *args);
extern PyObject *
PyImaging_EventLoopWin32(PyObject *self, PyObject *args);
extern PyObject *
PyImaging_DrawWmf(PyObject *self, PyObject *args);
#endif
#ifdef HAVE_XCB
extern PyObject *
PyImaging_GrabScreenX11(PyObject *self, PyObject *args);
#endif

/* Experimental path stuff (in path.c) */
extern PyObject *
PyPath_Create(ImagingObject *self, PyObject *args);

/* Experimental outline stuff (in outline.c) */
extern PyObject *
PyOutline_Create(ImagingObject *self, PyObject *args);

extern PyObject *
PyImaging_MapBuffer(PyObject *self, PyObject *args);

static PyMethodDef functions[] = {

    /* Functions */
    {"convert", (PyCFunction)_convert2, METH_VARARGS},

    /* Codecs */
    {"bcn_decoder", (PyCFunction)PyImaging_BcnDecoderNew, METH_VARARGS},
    {"bit_decoder", (PyCFunction)PyImaging_BitDecoderNew, METH_VARARGS},
    {"eps_encoder", (PyCFunction)PyImaging_EpsEncoderNew, METH_VARARGS},
    {"fli_decoder", (PyCFunction)PyImaging_FliDecoderNew, METH_VARARGS},
    {"gif_decoder", (PyCFunction)PyImaging_GifDecoderNew, METH_VARARGS},
    {"gif_encoder", (PyCFunction)PyImaging_GifEncoderNew, METH_VARARGS},
    {"hex_decoder", (PyCFunction)PyImaging_HexDecoderNew, METH_VARARGS},
    {"hex_encoder", (PyCFunction)PyImaging_EpsEncoderNew, METH_VARARGS}, /* EPS=HEX! */
#ifdef HAVE_LIBJPEG
    {"jpeg_decoder", (PyCFunction)PyImaging_JpegDecoderNew, METH_VARARGS},
    {"jpeg_encoder", (PyCFunction)PyImaging_JpegEncoderNew, METH_VARARGS},
#endif
#ifdef HAVE_OPENJPEG
    {"jpeg2k_decoder", (PyCFunction)PyImaging_Jpeg2KDecoderNew, METH_VARARGS},
    {"jpeg2k_encoder", (PyCFunction)PyImaging_Jpeg2KEncoderNew, METH_VARARGS},
#endif
#ifdef HAVE_LIBTIFF
    {"libtiff_decoder", (PyCFunction)PyImaging_LibTiffDecoderNew, METH_VARARGS},
    {"libtiff_encoder", (PyCFunction)PyImaging_LibTiffEncoderNew, METH_VARARGS},
#endif
    {"packbits_decoder", (PyCFunction)PyImaging_PackbitsDecoderNew, METH_VARARGS},
    {"pcd_decoder", (PyCFunction)PyImaging_PcdDecoderNew, METH_VARARGS},
    {"pcx_decoder", (PyCFunction)PyImaging_PcxDecoderNew, METH_VARARGS},
    {"pcx_encoder", (PyCFunction)PyImaging_PcxEncoderNew, METH_VARARGS},
    {"raw_decoder", (PyCFunction)PyImaging_RawDecoderNew, METH_VARARGS},
    {"raw_encoder", (PyCFunction)PyImaging_RawEncoderNew, METH_VARARGS},
    {"sgi_rle_decoder", (PyCFunction)PyImaging_SgiRleDecoderNew, METH_VARARGS},
    {"sun_rle_decoder", (PyCFunction)PyImaging_SunRleDecoderNew, METH_VARARGS},
    {"tga_rle_decoder", (PyCFunction)PyImaging_TgaRleDecoderNew, METH_VARARGS},
    {"tga_rle_encoder", (PyCFunction)PyImaging_TgaRleEncoderNew, METH_VARARGS},
    {"xbm_decoder", (PyCFunction)PyImaging_XbmDecoderNew, METH_VARARGS},
    {"xbm_encoder", (PyCFunction)PyImaging_XbmEncoderNew, METH_VARARGS},
#ifdef HAVE_LIBZ
    {"zip_decoder", (PyCFunction)PyImaging_ZipDecoderNew, METH_VARARGS},
    {"zip_encoder", (PyCFunction)PyImaging_ZipEncoderNew, METH_VARARGS},
#endif

/* Memory mapping */
#ifdef WITH_MAPPING
    {"map_buffer", (PyCFunction)PyImaging_MapBuffer, METH_VARARGS},
#endif

/* Display support */
#ifdef _WIN32
    {"display", (PyCFunction)PyImaging_DisplayWin32, METH_VARARGS},
    {"display_mode", (PyCFunction)PyImaging_DisplayModeWin32, METH_VARARGS},
    {"grabscreen_win32", (PyCFunction)PyImaging_GrabScreenWin32, METH_VARARGS},
    {"grabclipboard_win32", (PyCFunction)PyImaging_GrabClipboardWin32, METH_VARARGS},
    {"createwindow", (PyCFunction)PyImaging_CreateWindowWin32, METH_VARARGS},
    {"eventloop", (PyCFunction)PyImaging_EventLoopWin32, METH_VARARGS},
    {"listwindows", (PyCFunction)PyImaging_ListWindowsWin32, METH_VARARGS},
    {"drawwmf", (PyCFunction)PyImaging_DrawWmf, METH_VARARGS},
#endif
#ifdef HAVE_XCB
    {"grabscreen_x11", (PyCFunction)PyImaging_GrabScreenX11, METH_VARARGS},
#endif

/* Drawing support stuff */
#ifdef WITH_IMAGEDRAW
    {"font", (PyCFunction)_font_new, METH_VARARGS},
    {"draw", (PyCFunction)_draw_new, METH_VARARGS},
#endif

/* Experimental path stuff */
#ifdef WITH_IMAGEPATH
    {"path", (PyCFunction)PyPath_Create, METH_VARARGS},
#endif

/* Experimental arrow graphics stuff */
#ifdef WITH_ARROW
    {"outline", (PyCFunction)PyOutline_Create, METH_VARARGS},
#endif

    /* Resource management */
    {"get_stats", (PyCFunction)_get_stats, METH_VARARGS},
    {"reset_stats", (PyCFunction)_reset_stats, METH_VARARGS},
    {"get_alignment", (PyCFunction)_get_alignment, METH_VARARGS},
    {"get_block_size", (PyCFunction)_get_block_size, METH_VARARGS},
    {"get_blocks_max", (PyCFunction)_get_blocks_max, METH_VARARGS},
    {"set_alignment", (PyCFunction)_set_alignment, METH_VARARGS},
    {"set_block_size", (PyCFunction)_set_block_size, METH_VARARGS},
    {"set_blocks_max", (PyCFunction)_set_blocks_max, METH_VARARGS},
    {"clear_cache", (PyCFunction)_clear_cache, METH_VARARGS},

    {NULL, NULL} /* sentinel */
};

HPyDef_SLOT(_imaging_exec, HPy_mod_exec)
static int _imaging_exec_impl(HPyContext *ctx, HPy h_module) {

    const char *version = (char *)PILLOW_VERSION;

    HPy h_array_type = HPyType_FromSpec(ctx, &Imaging_Type_spec, NULL);
    if (HPy_IsNull(h_array_type)) {
        return 1;
    }

    HPyGlobal_Store(ctx, &hg_Imaging_Type, h_array_type);
    Imaging_Type = (PyTypeObject *) HPy_AsPyObject(ctx, h_array_type);
    HPy_Close(ctx, h_array_type);

#ifdef WITH_IMAGEDRAW
    if (PyType_Ready(&ImagingFont_Type) < 0) {
        return 1;
    }


    if (PyType_Ready(&ImagingDraw_Type) < 0) {
        return 1;
    }
#endif
    HPy h_PixelAcess_Type = HPyType_FromSpec(ctx, &PixelAccess_Type_spec, NULL);
    if (HPy_IsNull(h_PixelAcess_Type)) {
        return 1;
    }
    HPyGlobal_Store(ctx, &hg_PixelAccess_Type, h_PixelAcess_Type);
    HPy_Close(ctx, h_array_type);

    ImagingAccessInit();

#ifdef HAVE_LIBJPEG
    {
        extern const char *ImagingJpegVersion(void);
        HPy_SetAttrStringConstant_s(ctx, h_module,
            "jpeglib_version", ImagingJpegVersion());
    }
#endif

#ifdef HAVE_OPENJPEG
    {
        extern const char *ImagingJpeg2KVersion(void);
        HPy_SetAttrStringConstant_s(ctx, h_module,
            "jp2klib_version", ImagingJpeg2KVersion());
    }
#endif

    HPy have_libjpegturbo;
#ifdef LIBJPEG_TURBO_VERSION
    have_libjpegturbo = ctx->h_True;
#define tostr1(a) #a
#define tostr(a) tostr1(a)
    HPy_SetAttrStringConstant_s(ctx, h_module,
            "libjpeg_turbo_version", tostr(LIBJPEG_TURBO_VERSION));
#undef tostr
#undef tostr1
#else
    have_libjpegturbo = ctx->h_False;
#endif
    HPy_SetAttr_s(ctx, h_module, "HAVE_LIBJPEGTURBO", have_libjpegturbo);

    HPy have_libimagequant;
#ifdef HAVE_LIBIMAGEQUANT
    have_libimagequant = ctx->h_True;
    {
        extern const char *ImagingImageQuantVersion(void);
        HPy_SetAttrStringConstant_s(ctx, h_module, "imagequant_version", ImagingImageQuantVersion());
    }
#else
    have_libimagequant = ctx->h_False;
#endif
    HPy_SetAttr_s(ctx, h_module, "HAVE_LIBIMAGEQUANT", have_libimagequant);

#ifdef HAVE_LIBZ
    /* zip encoding strategies */
    HPy_SetAttrIntConstant_s(ctx, h_module, "DEFAULT_STRATEGY", Z_DEFAULT_STRATEGY);
    HPy_SetAttrIntConstant_s(ctx, h_module, "FILTERED", Z_FILTERED);
    HPy_SetAttrIntConstant_s(ctx, h_module, "HUFFMAN_ONLY", Z_HUFFMAN_ONLY);
    HPy_SetAttrIntConstant_s(ctx, h_module, "RLE", Z_RLE);
    HPy_SetAttrIntConstant_s(ctx, h_module, "FIXED", Z_FIXED);
    {
        extern const char *ImagingZipVersion(void);
        HPy_SetAttrStringConstant_s(ctx, h_module, "zlib_version", ImagingZipVersion());
    }
#endif

#ifdef HAVE_LIBTIFF
    {
        extern const char *ImagingTiffVersion(void);
        HPy_SetAttrStringConstant_s(ctx, h_module, "libtiff_version", ImagingTiffVersion());

        // Test for libtiff 4.0 or later, excluding libtiff 3.9.6 and 3.9.7
        HPy support_custom_tags;
#if TIFFLIB_VERSION >= 20111221 && TIFFLIB_VERSION != 20120218 && \
    TIFFLIB_VERSION != 20120922
        support_custom_tags = ctx->h_True;
#else
        support_custom_tags = ctx->h_False;
#endif
        HPy_SetAttr_s(ctx, h_module, "libtiff_support_custom_tags", support_custom_tags);
    }
#endif

    HPy have_xcb;
#ifdef HAVE_XCB
    have_xcb = ctx->h_True;
#else
    have_xcb = ctx->h_False;
#endif
    HPy_SetAttr_s(ctx, h_module, "HAVE_XCB", have_xcb);

    HPy_SetAttrStringConstant_s(ctx, h_module, "PILLOW_VERSION", version);

    return 0;
}

static HPyDef *module_defines[] = {
    &_imaging_exec,

    /* Object factories */
    &fill,
    &new,
    &alpha_composite,
    &blend,
    &merge,

    /* Special effects (experimental) */
#ifdef WITH_EFFECTS
    &effect_mandelbrot,
    &effect_noise,
    &linear_gradient,
    &radial_gradient,
#endif

    NULL,
};

static HPyGlobal *module_globals[] = {
    &hg_Imaging_Type,
    &hg_PixelAccess_Type,
    NULL
};

static HPyModuleDef module_def = {
    .doc = NULL,
    .size = 0,
    .legacy_methods = functions,
    .defines = module_defines,
    .globals = module_globals,
};

HPy_MODINIT(_imaging, module_def)
