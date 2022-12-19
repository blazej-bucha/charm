/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
/* ------------------------------------------------------------------------- */






CHARM(shc) *CHARM(shc_init)(unsigned long nmax, REAL mu, REAL r,
                            REAL *c, REAL *s)
{
    /* --------------------------------------------------------------------- */
    CHARM(shc) *shcs = (CHARM(shc) *)malloc(sizeof(CHARM(shc)));
    if (shcs == NULL)
        return shcs;


    shcs->c = shcs->s = NULL;
    /* --------------------------------------------------------------------- */


    /* Save the scalar input values to the "CHARM(shc)" structure */
    /* --------------------------------------------------------------------- */
    shcs->nmax = nmax;
    shcs->mu   = mu;

    if (r <= PREC(0.0))
        goto FAILURE;
    shcs->r = r;

    /* Size of the first dimension of the "shcs->c" and "shcs->s" pointer
     * arrays (number of harmonic orders) */
    size_t nmax1 = (size_t)(nmax) + 1;
    /* Total number of spherical harmonic coefficients to be stored in both
     * "shcs->c[0]" and "shcs->s[0]".  Note that the "S_{n,0}" coefficients,
     * which do not exist, will be stored in "shcs->s[0]" as zeros. */
    size_t len = ((nmax1 + 1) * nmax1) / 2;
    shcs->nc   = len;
    shcs->ns   = len;
    /* --------------------------------------------------------------------- */


    /* Allocate the "shcs->c" and "shcs->s" pointer arrays */
    /* --------------------------------------------------------------------- */
    shcs->c = (REAL **)malloc(nmax1 * sizeof(REAL *));
    if (shcs->c == NULL)
        goto FAILURE;


    shcs->s = (REAL **)malloc(nmax1 * sizeof(REAL *));
    if (shcs->s == NULL)
        goto FAILURE;
    /* --------------------------------------------------------------------- */


    /* Set "shcs->c[m]" and "shcs->s[m]" to point to "C_{m,m}" and "S_{m,m}" */
    /* --------------------------------------------------------------------- */
    shcs->c[0]  = c;
    shcs->s[0]  = s;
    shcs->owner = 0;


    unsigned long nmax1mm = 0;
    for (unsigned long m = 0; m <= nmax; m++)
    {
        /* Now set "shcs->c[m]" and "shcs->s[m]" to point to "C_{m,m}" and
         * "S_{m,m}" */
        shcs->c[m] = shcs->c[0] + nmax1mm;
        shcs->s[m] = shcs->s[0] + nmax1mm;


        nmax1mm += nmax1 - m;
    }
    /* --------------------------------------------------------------------- */


    /* --------------------------------------------------------------------- */
EXIT:
    return shcs;
    /* --------------------------------------------------------------------- */


FAILURE:
    /* --------------------------------------------------------------------- */
    free(shcs->c); free(shcs->s);
    free(shcs);
    shcs = NULL;


    goto EXIT;
    /* --------------------------------------------------------------------- */
}
