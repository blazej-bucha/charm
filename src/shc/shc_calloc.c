/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
/* ------------------------------------------------------------------------- */






CHARM(shc) *CHARM(shc_calloc)(unsigned long nmax, REAL mu, REAL r)
{
    CHARM(shc) *shcs = NULL;
    REAL *c          = NULL;
    REAL *s          = NULL;


    /* Prepare arrays with spherical harmonic coefficients filled with zeros */
    /* --------------------------------------------------------------------- */
    size_t nmax1 = (size_t)(nmax) + 1;
    size_t len   = ((nmax1 + 1) * nmax1) / 2;


    c = (REAL *)calloc(len, sizeof(REAL));
    if (c == NULL)
        goto FAILURE;


    s = (REAL *)calloc(len, sizeof(REAL));
    if (s == NULL)
        goto FAILURE;
    /* --------------------------------------------------------------------- */


    /* --------------------------------------------------------------------- */
    shcs = CHARM(shc_init)(nmax, mu, r, c, s);
    if (shcs == NULL)
        goto FAILURE;


    /* "shc_init" sets "owner" to "False".  But in this case, CHarm allocated
     * the arrays of spherical harmonic coefficients, so we have to set the
     * "owner" to "True". */
    shcs->owner = 1;
    /* --------------------------------------------------------------------- */


    /* --------------------------------------------------------------------- */
EXIT:
    return shcs;
    /* --------------------------------------------------------------------- */


    /* --------------------------------------------------------------------- */
FAILURE:
    free(c);
    free(s);
    goto EXIT;
    /* --------------------------------------------------------------------- */
}
