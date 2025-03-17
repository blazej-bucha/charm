/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include "../prec.h"
#include <stdio.h>
#include <mpfr.h>
#include "../err/err_set.h"
#include "mpfr_flush_unreleased_memory.h"
#include "mpfr_ndarray.h"
#include "mpfr_ndarray_check.h"
/* ------------------------------------------------------------------------- */






/* Internal function to compute unnormalized Legendre function of argument
 * "cos(psi)" for all degrees and orders up to "nmax" and "mmax", respectively,
 * using "NBITS" for "mpfr_t".  The output is stored in "pnm"; Legendre
 * function of degree "n" and order "m" is to be accessed as
 *
 *      pnm->data[m * pnm->shape[1] + n]
 *
 * */
void CHARM(mpfr_legendre)(mpfr_ndarray *pnm,
                          unsigned long nmax,
                          unsigned long mmax,
                          const mpfr_t psi,
                          mpfr_prec_t NBITS,
                          CHARM(err) *err)
{
    if (CHARM(mpfr_ndarray_check)(pnm, 2, mmax + 1, nmax + 1))
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Wrong shape of the input \"pnm\" mpfr_ndarray.");
        return;
    }


    /* Set all array elements of "pnm->data" to zero */
    for (size_t i = 0; i < pnm->size; i++)
        mpfr_set_ui(pnm->data[i], 0, MPFR_RNDN);


    mpfr_t cospsi, sinpsi, tmp;
    mpfr_inits2(NBITS, cospsi, sinpsi, tmp, (mpfr_ptr)NULL);


    mpfr_set_ui(pnm->data[0], 1, MPFR_RNDN);
    if (nmax == 0)
        goto EXIT;


    if (mmax > nmax)
        goto EXIT;


    mpfr_cos(cospsi, psi, MPFR_RNDN);
    mpfr_sin(sinpsi, psi, MPFR_RNDN);
    mpfr_set(pnm->data[1], cospsi, MPFR_RNDN);


    /* Zonal Legendre functions */
    for (unsigned long n = 2; n <= nmax; n++)
    {
        mpfr_mul_ui(pnm->data[n], cospsi, 2 * n - 1, MPFR_RNDN);
        mpfr_mul(pnm->data[n], pnm->data[n], pnm->data[n - 1], MPFR_RNDN);
        mpfr_mul_ui(tmp, pnm->data[n - 2], n - 1, MPFR_RNDN);
        mpfr_sub(pnm->data[n], pnm->data[n], tmp, MPFR_RNDN);
        mpfr_div_ui(pnm->data[n], pnm->data[n], n, MPFR_RNDN);
    }


    if (mmax == 0)
        goto EXIT;


    /* The "nmax + 2" index below is derived from "m * (nmax + 1) + n" for "n"
     * and "m" being "1". */
    mpfr_set(pnm->data[nmax + 2], sinpsi, MPFR_RNDN);
    if (nmax == 1)
        goto EXIT;


    size_t idx;
    size_t offset = pnm->shape[1];


    /* Sectorial and tesseral Legendre functions */
    for (unsigned long n = 2; n <= nmax; n++)
    {
        for (unsigned long m = 1; m <= CHARM_MIN(n, mmax); m++)
        {
            idx = m * offset + n;
            mpfr_mul_ui(pnm->data[idx], sinpsi, n - m + 1, MPFR_RNDN);
            mpfr_mul(pnm->data[idx], pnm->data[idx],
                     pnm->data[(m - 1) * offset + n], MPFR_RNDN);
            mpfr_add(pnm->data[idx], pnm->data[idx], pnm->data[idx - 1],
                     MPFR_RNDN);
            mpfr_div(pnm->data[idx], pnm->data[idx], cospsi, MPFR_RNDN);
        }
    }


EXIT:
    mpfr_clears(cospsi, sinpsi, tmp, (mpfr_ptr)NULL);
    mpfr_free_cache();
    FLUSH_UNRELEASED_MEMORY;


    return;
}

