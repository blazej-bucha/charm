/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdlib.h>
#include "../prec.h"
#include "leg_pnmj_init.h"
/* ------------------------------------------------------------------------- */






CHARM(pnmj) *CHARM(leg_pnmj_calloc)(unsigned long nmax, int ordering)
{
    CHARM(pnmj) *pnmj = NULL;
    REAL *pnmj_coeffs = NULL;


    /* Prepare the array to store the Fourier coefficients of Legendre
     * functions */
    /* --------------------------------------------------------------------- */
    /* Get the total number of Fourier coefficients that are associated with
     * the maximum degree "nmax" */
    size_t npnmj = CHARM(leg_pnmj_length)(nmax);


    pnmj_coeffs = (REAL *)calloc(npnmj, sizeof(REAL));
    if (pnmj_coeffs == NULL)
        return NULL;


    pnmj = CHARM(leg_pnmj_init)(nmax, ordering, pnmj_coeffs);
    if (pnmj == NULL)
        goto FAILURE;
    /* --------------------------------------------------------------------- */


    /* --------------------------------------------------------------------- */
EXIT:
    return pnmj;
    /* --------------------------------------------------------------------- */


    /* --------------------------------------------------------------------- */
FAILURE:
    free(pnmj_coeffs);
    goto EXIT;
    /* --------------------------------------------------------------------- */
}