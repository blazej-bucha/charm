/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
#include "shc_alloc.h"
/* ------------------------------------------------------------------------- */






CHARM(shc) *CHARM(shc_malloc)(unsigned long nmax, REAL mu, REAL r)
{
    return CHARM(shc_alloc)(nmax, mu, r, malloc);
}
