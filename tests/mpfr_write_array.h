/* This header file is not a part of API. */


#ifndef __MPFR_WRITE_ARRAY_H__
#define __MPFR_WRITE_ARRAY_H__


#include <config.h>
#include <mpfr.h>


#ifdef __cplusplus
extern "C"
{
#endif


extern long int mpfr_write_array(mpfr_t *,
                                 size_t,
                                 FILE *);


#ifdef __cplusplus
}
#endif


#endif
