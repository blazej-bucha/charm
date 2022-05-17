/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
/* ------------------------------------------------------------------------- */






/* A small function to print a floating point number to a text file.  Not
 * a part of the CHarm's API. */
int CHARM(misc_fprintf_real)(FILE *stream, const char *format, REAL c)
{
#if CHARM_QUAD


    int flt128_n;
    char buf[256];


    flt128_n = quadmath_snprintf(buf, sizeof(buf), format, 70, c);
    if ((size_t)flt128_n >= 256)
        return -1;


    return fprintf(stream, "%s", buf);


#else


    return fprintf(stream, format, c);


#endif
}
