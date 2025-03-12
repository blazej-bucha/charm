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
#else
    *major = LIB_NA_VAL;
#endif


#if HAVE_MPFR && defined(__GNU_MP_VERSION_MINOR)
    *minor = __GNU_MP_VERSION_MINOR;
#else
    *minor = LIB_NA_VAL;
#endif


#if HAVE_MPFR && defined(__GNU_MP_VERSION_PATCHLEVEL)
    *patch = __GNU_MP_VERSION_PATCHLEVEL;
#else
    *patch = LIB_NA_VAL;
#endif


    return
#if HAVE_MPFR
    gmp_version;
#else
    LIB_NA_STR;
#endif
}
