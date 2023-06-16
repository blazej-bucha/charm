/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "misc_buildopt.h"
/* ------------------------------------------------------------------------- */






int CHARM(misc_buildopt_isfinite)(void)
{
#if HAVE_ISFINITE
    return BUILDOPT_ISFINITE;
#else
    return 0;
#endif
}
