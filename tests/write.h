/* This header file is not a part of API. */


#ifndef __WRITE_H__
#define __WRITE_H__


#include <config.h>
#include "../src/prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern long int write(char *,
                      REAL *,
                      size_t);


#ifdef __cplusplus
}
#endif


#endif
