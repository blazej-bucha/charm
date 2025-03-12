/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "misc_buildopt.h"
/* ------------------------------------------------------------------------- */






int CHARM(misc_buildopt_mpfr)(void)
{
#if HAVE_MPFR
    return BUILDOPT_MPI;
#else
    return 0;
#endif
}
