/* This header file is not a part of API. */


#ifndef __SHC_READ_LINE_H__
#define __SHC_READ_LINE_H__


#include <config.h>
#include "../prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern int CHARM(shc_read_line)(FILE *, unsigned long *, unsigned long *,
                                REAL *, REAL *, int, CHARM(err) *);


/* Symbolic constant to read a line from a table of spherical harmonic
 * coefficients */
#undef SHC_READ_LINE_TBL
#define SHC_READ_LINE_TBL 0


/* Symbolic constant to read a line from a table of spherical harmonic
 * coefficients */
#undef SHC_READ_LINE_GFC
#define SHC_READ_LINE_GFC 1


/* Size of char arrays to store value of parameters loaded from the "gfc" file
 * */
#undef SHC_READ_LINE_NSTR
#define SHC_READ_LINE_NSTR (128)


/* Size of the char array to store a single line of the "gfc" file */
#undef SHC_READ_LINE_NLINE
#define SHC_READ_LINE_NLINE (2048)


/* Keyword from the data section of "gfc" files */
#undef SHC_READ_LINE_GFC_KEYWORD
#define SHC_READ_LINE_GFC_KEYWORD "gfc"
#undef SHC_READ_LINE_GFCT_KEYWORD
#define SHC_READ_LINE_GFCT_KEYWORD "gfct"
#undef SHC_READ_LINE_TRND_KEYWORD
#define SHC_READ_LINE_TRND_KEYWORD "trnd"
#undef SHC_READ_LINE_ASIN_KEYWORD
#define SHC_READ_LINE_ASIN_KEYWORD "asin"
#undef SHC_READ_LINE_ACOS_KEYWORD
#define SHC_READ_LINE_ACOS_KEYWORD "acos"


#ifdef __cplusplus
}
#endif


#endif
