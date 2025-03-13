/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#if HAVE_MPFR
#   include <gmp.h>
#endif
#include "../prec.h"
#include "misc_buildopt.h"
/* ------------------------------------------------------------------------- */






const char *CHARM(misc_buildopt_version_gmp)(int *major,
                                             int *minor,
                                             int *patch)
{
#if HAVE_MPFR && defined(__GNU_MP_VERSION)
    *major = __GNU_MP_VERSION;
    *minor = __GNU_MP_VERSION_MINOR;
    *patch = __GNU_MP_VERSION_PATCHLEVEL;
    return gmp_version;
#else
    *major = LIB_NA_VAL;
    *minor = LIB_NA_VAL;
    *patch = LIB_NA_VAL;
    return LIB_NA_STR;
#endif
}
