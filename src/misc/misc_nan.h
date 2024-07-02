/* This header file is not a part of API. */


#ifndef __MISC_NAN_H__
#define __MISC_NAN_H__


#include <config.h>
#include <math.h>
#include "../prec.h"


#ifndef NAN
#   define NAN (PREC(0.0) / PREC(0.0))
#endif


#endif
