/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdlib.h>
#include "../prec.h"
#include "leg_pnmj_alloc.h"
/* ------------------------------------------------------------------------- */






CHARM(pnmj) *CHARM(leg_pnmj_malloc)(unsigned long nmax, int ordering)
{
    return CHARM(leg_pnmj_alloc)(nmax, ordering, malloc);
}
