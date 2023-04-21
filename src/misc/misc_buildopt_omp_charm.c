/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "misc_buildopt.h"
/* ------------------------------------------------------------------------- */






int CHARM(misc_buildopt_omp_charm)(void)
{
#if CHARM_OPENMP
    return BUILDOPT_OMP_CHARM;
#else
    return 0;
#endif
}
