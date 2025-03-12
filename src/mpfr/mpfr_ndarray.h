/* This header file is not a part of API. */


#ifndef __MPFR_NDARRAY_H__
#define __MPFR_NDARRAY_H__


#include <mpfr.h>


/* "mpfr_ndarray" is a structure representing something that resembles
 * N-dimensional array of "mpfr_t" data type */
typedef struct
{
    mpfr_t *data;
    size_t ndim;
    size_t *shape;
    size_t size;
    _Bool  owner;
} mpfr_ndarray;


#endif

