/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <fftw3.h>
#include "../prec.h"
/* ------------------------------------------------------------------------- */






const char * const CHARM(misc_buildopt_version_fftw)(void)
{
    return FFTW(version);
}
