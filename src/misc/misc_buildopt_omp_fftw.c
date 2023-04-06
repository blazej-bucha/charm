/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "misc_buildopt.h"
/* ------------------------------------------------------------------------- */






int CHARM(misc_buildopt_omp_fftw)(void)
{
#if FFTW3_OMP
    return BUILDOPT_OMP_FFTW;
#else
    return 0;
#endif
}
