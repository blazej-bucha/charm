/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <fftw3.h>
#include "../prec.h"
#include "misc_buildopt.h"
/* ------------------------------------------------------------------------- */






const char * CHARM(misc_buildopt_version_fftw)(void)
{
    return
           /* Not sure why, but "FFTW(version)" is not available on Windows, at
            * least not in the binaries precompiled by the FFTW team. */
#if HAVE_FFTW_VERSION
           FFTW(version)
#else
           LIB_NA_STR
#endif
    ;
}
