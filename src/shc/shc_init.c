/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
#include "shc_init_chunk.h"
/* ------------------------------------------------------------------------- */






CHARM(shc) *CHARM(shc_init)(unsigned long nmax,
                            REAL mu,
                            REAL r,
                            REAL *c,
                            REAL *s)
{
    CHARM(err) *err = CHARM(err_init)();


    unsigned long chunk[2] = {0, nmax};
    CHARM(shc) *shcs = CHARM(shc_init_chunk)(nmax, mu, r, c, s, 1, chunk, err);
    if (!CHARM(err_isempty)(err) || (shcs == NULL))
    {
        CHARM(shc_free)(shcs);
        shcs = NULL;
        goto EXIT;
    }


EXIT:
    CHARM(err_free)(err);
    return shcs;
}

