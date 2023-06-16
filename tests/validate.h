/* This header file is not a part of API. */


#ifndef __VALIDATE_H__
#define __VALIDATE_H__


#include <config.h>
#include "../src/prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern long int validate(char *,
                         REAL *,
                         size_t,
                         REAL);


#ifdef __cplusplus
}
#endif


#endif
