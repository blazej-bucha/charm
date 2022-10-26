/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
#include "shc_alloc.h"
#include "../misc/misc_calloc.h"
/* ------------------------------------------------------------------------- */






CHARM(shc) *CHARM(shc_calloc)(unsigned long nmax, REAL mu, REAL r)
{
    return CHARM(shc_alloc)(nmax, mu, r, CHARM(misc_calloc));
}
