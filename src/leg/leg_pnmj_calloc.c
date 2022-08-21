/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdlib.h>
#include "../prec.h"
#include "leg_pnmj_alloc.h"
#include "../misc/misc_calloc.h"
/* ------------------------------------------------------------------------- */






CHARM(pnmj) *CHARM(leg_pnmj_calloc)(unsigned long nmax, int ordering)
{
    return CHARM(leg_pnmj_alloc)(nmax, ordering, CHARM(misc_calloc));
}
