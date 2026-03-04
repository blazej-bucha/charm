/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
/* ------------------------------------------------------------------------- */






CHARM(shc) *CHARM(shc_copy)(const CHARM(shc) *shcs,
                            unsigned long nmin,
                            unsigned long nmax,
                            unsigned long nmax_shcs_out)
{
    if (nmin > shcs->nmax)
        return NULL;


    if (nmax > shcs->nmax)
        return NULL;


    if (nmin > nmax)
        return NULL;


    if (nmax > nmax_shcs_out)
        return NULL;


    if (shcs->distributed)
        return NULL;


    CHARM(shc) *shcs_out = CHARM(shc_calloc)(nmax_shcs_out, shcs->mu, shcs->r);
    if (shcs_out == NULL)
        return NULL;


    unsigned long nmm;
    for (unsigned long m = 0; m <= nmax; m++)
    {
        for (unsigned long n = CHARM_MAX(m, nmin); n <= nmax; n++)
        {
            nmm = n - m;
            shcs_out->c[m][nmm] = shcs->c[m][nmm];
            shcs_out->s[m][nmm] = shcs->s[m][nmm];
        }
    }


    return shcs_out;
}
