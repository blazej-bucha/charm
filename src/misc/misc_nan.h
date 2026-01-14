/* This header file is not a part of API. */


#ifndef __MISC_NAN_H__
#define __MISC_NAN_H__


#include <config.h>
#include <math.h>
#include "../prec.h"


#undef CHARM_NAN
#if CHARM_FLOAT
#   define CHARM_NAN        nanf("")
#elif CHARM_QUAD
#   define CHARM_NAN        nanq("")
#else
#   define CHARM_NAN        nan("")
#endif


#endif
