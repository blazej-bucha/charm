/* This header file is not a part of API. */


#ifndef __SHC_READ_GFC_H__
#define __SHC_READ_GFC_H__


#include <config.h>
#include "../prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


/* Size of char arrays to store value of parameters loaded from the "gfc" file
 * */
#undef SHC_READ_GFC_NSTR
#define SHC_READ_GFC_NSTR (128)


/* Size of the char array to store a single line of the "gfc" file */
#undef SHC_READ_GFC_NLINE
#define SHC_READ_GFC_NLINE (2048)


/* Mandatory keywords of "gfc" files */
/* ------------------------------------------------------------------------- */
#undef SHC_READ_GFC_BOH
#define SHC_READ_GFC_BOH "begin_of_head"


#undef SHC_READ_GFC_EOH
#define SHC_READ_GFC_EOH "end_of_head"


#undef SHC_READ_GFC_NMAX
#define SHC_READ_GFC_NMAX "max_degree"


#undef SHC_READ_GFC_EARTH_GM
#define SHC_READ_GFC_EARTH_GM "earth_gravity_constant"


#undef SHC_READ_GFC_GM
#define SHC_READ_GFC_GM "gravity_constant"


#undef SHC_READ_GFC_ALL_GM_KEYWORDS
#define SHC_READ_GFC_ALL_GM_KEYWORDS SHC_READ_GFC_EARTH_GM "\" or \"" \
                                     SHC_READ_GFC_GM


#undef SHC_READ_GFC_R
#define SHC_READ_GFC_R "radius"


/* Optional keyword of "gfc" files */
#undef SHC_READ_GFC_NORM
#define SHC_READ_GFC_NORM "norm"


/* Mandatory keyword of "gfc" files */
#undef SHC_READ_GFC_ERRORS
#define SHC_READ_GFC_ERRORS "errors"
/* ------------------------------------------------------------------------- */


/* Optional keywords of "gfc" files */
/* ------------------------------------------------------------------------- */
#undef SHC_READ_GFC_FORMAT
#define SHC_READ_GFC_FORMAT "format"


#undef SHC_READ_GFC_FORMAT1d0
#define SHC_READ_GFC_FORMAT1d0 "icgem1.0"


#undef SHC_READ_GFC_FORMAT2d0
#define SHC_READ_GFC_FORMAT2d0 "icgem2.0"
/* ------------------------------------------------------------------------- */


/* Valid values of the "SHC_READ_GFC_ERRORS" */
#undef SHC_READ_GFC_ERRORS_NO
#define SHC_READ_GFC_ERRORS_NO "no"
#undef SHC_READ_GFC_ERRORS_CALIBRATED
#define SHC_READ_GFC_ERRORS_CALIBRATED "calibrated"
#undef SHC_READ_GFC_ERRORS_FORMAL
#define SHC_READ_GFC_ERRORS_FORMAL "formal"
#undef SHC_READ_GFC_ERRORS_CALIBRATED_AND_FORMAL
#define SHC_READ_GFC_ERRORS_CALIBRATED_AND_FORMAL "calibrated_and_formal"


/* *If* the "GFC_NORM" keyword is found, its value has to be "GFC_NORM_FULL" */
#undef SHC_READ_GFC_NORM_FULL
#define SHC_READ_GFC_NORM_FULL "fully_normalized"


/* Keyword from the data section of "gfc" files */
#undef SHC_READ_GFC_GFC
#define SHC_READ_GFC_GFC "gfc"
#undef SHC_READ_GFC_GFCT
#define SHC_READ_GFC_GFCT "gfct"
#undef SHC_READ_GFC_TRND
#define SHC_READ_GFC_TRND "trnd"
#undef SHC_READ_GFC_DOT
#define SHC_READ_GFC_DOT "dot"
#undef SHC_READ_GFC_ASIN
#define SHC_READ_GFC_ASIN "asin"
#undef SHC_READ_GFC_ACOS
#define SHC_READ_GFC_ACOS "acos"


#ifdef __cplusplus
}
#endif


#endif
