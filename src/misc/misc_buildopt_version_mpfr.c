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
    *minor = MPFR_VERSION_MINOR;
    *patch = MPFR_VERSION_PATCHLEVEL;
    return mpfr_get_version();
#else
    *major = LIB_NA_VAL;
    *minor = LIB_NA_VAL;
    *patch = LIB_NA_VAL;
    return LIB_NA_STR;
#endif
}
