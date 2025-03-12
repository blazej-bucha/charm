/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#if HAVE_MPFR
#   include <mpfr.h>
#endif
#include "../prec.h"
#include "misc_buildopt.h"
/* ------------------------------------------------------------------------- */






const char *CHARM(misc_buildopt_version_mpfr)(int *major,
                                              int *minor,
                                              int *patch)
{
#if HAVE_MPFR && defined(MPFR_VERSION_MAJOR)
    *major = MPFR_VERSION_MAJOR;
#else
    *major = LIB_NA_VAL;
#endif


#if HAVE_MPFR && defined(MPFR_VERSION_MINOR)
    *minor = MPFR_VERSION_MINOR;
#else
    *minor = LIB_NA_VAL;
#endif


#if HAVE_MPFR && defined(MPFR_VERSION_PATCHLEVEL)
    *patch = MPFR_VERSION_PATCHLEVEL;
#else
    *patch = LIB_NA_VAL;
#endif


    return
#if HAVE_MPFR
    mpfr_get_version();
#else
    LIB_NA_STR;
#endif
}
