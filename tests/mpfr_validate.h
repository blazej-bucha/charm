/* This header file is not a part of API. */


#ifndef __MPFR_VALIDATE_H__
#define __MPFR_VALIDATE_H__


#include <config.h>
#include <mpfr.h>


#ifdef __cplusplus
extern "C"
{
#endif


extern long int mpfr_validate(char *,
                              mpfr_t *,
                              size_t,
                              mpfr_t);


#ifdef __cplusplus
}
#endif


#endif
