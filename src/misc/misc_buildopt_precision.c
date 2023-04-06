/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "misc_buildopt.h"
/* ------------------------------------------------------------------------- */






int CHARM(misc_buildopt_precision)(void)
{
#if CHARM_FLOAT
    return BUILDOPT_PRECISION_SINGLE;
#elif CHARM_QUAD
    return BUILDOPT_PRECISION_QUAD;
#else
    return BUILDOPT_PRECISION_DOUBLE;
#endif
}
