/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
/* ------------------------------------------------------------------------- */






CHARM(shc) *CHARM(shc_init)(unsigned long nmax, REAL mu, REAL r)
{
    /* --------------------------------------------------------------------- */
    /* Allocate memory for the "CHARM(shc)" data type */
    CHARM(shc) *shcs = (CHARM(shc) *)malloc(sizeof(CHARM(shc)));
    if (shcs == NULL)
        return shcs;


    shcs->c = shcs->s = NULL;
    REAL *cnm0_tmp    = NULL;
    REAL *snm0_tmp    = NULL;
    /* --------------------------------------------------------------------- */


    /* Save the input values to the "CHARM(shc)" structure */
    /* --------------------------------------------------------------------- */
    shcs->nmax = nmax;
    shcs->mu   = mu;

    if (r <= PREC(0.0))
        goto FAILURE;
    shcs->r = r;
    /* --------------------------------------------------------------------- */


    /* Allocate the output "shcs->c" and "shcs->s" arrays */
    /* --------------------------------------------------------------------- */
    /* Get the size of the first dimension of the output "shcs->c" and
     * "shcs->s" arrays (number of harmonic orders) */
    size_t nmax1 = (size_t)(nmax) + 1;


    shcs->c = (REAL **)malloc(nmax1 * sizeof(REAL *));
    if (shcs->c == NULL)
        goto FAILURE;


    shcs->s = (REAL **)malloc(nmax1 * sizeof(REAL *));
    if (shcs->s == NULL)
        goto FAILURE;
    /* --------------------------------------------------------------------- */


    /* Allocate two blocks of memory to store the spherical harmonic
     * coefficients and initialize them to zeros */
    /* --------------------------------------------------------------------- */
    /* Get the total number of coefficients to be stored in "shcs->c" and
     * "shcs->s".  Note that the "S_{n,0}" coefficients, which do not exist,
     * will be stored in "shcs->s" as zeros. */
    size_t len = ((nmax1 + 1) * nmax1) / 2;


    shcs->nc = len;
    shcs->ns = len;


    cnm0_tmp = (REAL *)calloc(len, sizeof(REAL));
    if (cnm0_tmp == NULL)
        goto FAILURE;
    snm0_tmp = (REAL *)calloc(len, sizeof(REAL));
    if (snm0_tmp == NULL)
        goto FAILURE;


    shcs->c[0] = cnm0_tmp;
    shcs->s[0] = snm0_tmp;
    /* --------------------------------------------------------------------- */


    /* Set "shcs->c[m]" and "shcs->s[m]" to point to "C_{m,m}" and "S_{m,m}" */
    /* --------------------------------------------------------------------- */
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
    free(cnm0_tmp);  free(snm0_tmp);
    free(shcs->c); free(shcs->s);
    free(shcs);
    shcs = NULL;


    goto EXIT;
    /* --------------------------------------------------------------------- */
}
