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


/* Mandatory keyword of "gfc" files */
#undef SHC_READ_GFC_EOH
#define SHC_READ_GFC_EOH "end_of_head"


/* Mandatory keyword of "gfc" files */
#undef SHC_READ_GFC_NMAX
#define SHC_READ_GFC_NMAX "max_degree"


/* Mandatory keyword of "gfc" files */
#undef SHC_READ_GFC_GM
#define SHC_READ_GFC_GM "earth_gravity_constant"


/* Mandatory keyword of "gfc" files */
#undef SHC_READ_GFC_R
#define SHC_READ_GFC_R "radius"


/* Optional keyword of "gfc" files */
#undef SHC_READ_GFC_NORM
#define SHC_READ_GFC_NORM "norm"


/* *If* the "GFC_NORM" keyword is found, is value has to be "GFC_NORM_FULL" */
#undef SHC_READ_GFC_NORM_FULL
#define SHC_READ_GFC_NORM_FULL "fully_normalized"


#ifdef __cplusplus
}
#endif


#endif
