/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpfr.h>
#include "../prec.h"
#include "../err/err_set.h"
#include "../err/err_propagate.h"
#include "../mpfr/mpfr_ndarray.h"
#include "../mpfr/mpfr_ndarray_malloc.h"
#include "../mpfr/mpfr_flush_unreleased_memory.h"
#include "../mpfr/mpfr_ndarray_free.h"
#include "../mpfr/mpfr_fact.h"
#include "../mpfr/mpfr_double_fact.h"
#include "../mpfr/mpfr_legendre.h"
#include "../mpfr/mpfr_binomial.h"
#include "../misc/misc_idx_4d.h"
#include "gfm_cap_q_rpows.h"
#include "gfm_cap_q_dr.h"
#include "gfm_cap_q_dkernel.h"
#include "gfm_cap_quu.h"
#include "gfm_check_p.h"
#include "gfm_check_kminkmax.h"
/* ------------------------------------------------------------------------- */






void CHARM(gfm_cap_quu)(const mpfr_t rref,
                   const mpfr_t r,
                   const mpfr_t psi,
                   unsigned long nmax,
                   unsigned pmax,
                   unsigned kmin,
                   unsigned kmax,
                   unsigned imax,
                   unsigned u,
                   int zone,
                   mpfr_prec_t NBITS,
                   mpfr_t *qkpin,
                   CHARM(err) *err)
{
    /* Error checks */
    /* --------------------------------------------------------------------- */
    if ((u != 1) && (u != 2))
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "The input parameter \"u\" must either be \"1\" or "
                       "\"2\".");
        return;
    }


    CHARM(gfm_check_p)(pmax, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return;
    }


    CHARM(gfm_check_kminkmax)(kmin, kmax, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return;
    }
    /* --------------------------------------------------------------------- */






    /* Inits */
    /* --------------------------------------------------------------------- */
    mpfr_ndarray *pnm         = NULL;
    mpfr_ndarray *fact        = NULL;
    mpfr_ndarray *double_fact = NULL;
    mpfr_ndarray *rpows       = NULL;
    mpfr_ndarray *rwq         = NULL;
    mpfr_ndarray *binomial    = NULL;
    mpfr_ndarray *dkdr        = NULL;
    mpfr_ndarray *dkdrdp      = NULL;
    mpfr_t *qkpin00           = NULL;
    /* --------------------------------------------------------------------- */






    /* Legendre functions */
    /* --------------------------------------------------------------------- */
    /* "u" represents the maximum order of the derivative of the trunction
     * coefficients with respect to "psi" */
    pnm = CHARM(mpfr_ndarray_malloc)(NBITS, 2, u + 1, nmax + 1);
    if (pnm == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        goto EXIT;
    }


    CHARM(mpfr_legendre)(pnm, nmax, u, psi, NBITS, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        goto EXIT;
    }
    /* --------------------------------------------------------------------- */






    /* Factorials and double factorials */
    /* --------------------------------------------------------------------- */
    unsigned dmax = pmax + kmax + u + 1;


    fact = CHARM(mpfr_ndarray_malloc)(NBITS, 1, dmax + 2);
    if (fact == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        goto EXIT;
    }
    for (unsigned i = 0; i < dmax + 2; i++)
        CHARM(mpfr_fact)(i, fact->data[i]);


    double_fact = CHARM(mpfr_ndarray_malloc)(NBITS, 1, 2 * dmax + 1);
    if (double_fact == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        goto EXIT;
    }
    for (unsigned i = 0; i <= 2 * dmax; i++)
        CHARM(mpfr_double_fact)(i, double_fact->data[i]);
    /* --------------------------------------------------------------------- */






    /* Powers of "r" */
    /* --------------------------------------------------------------------- */
    unsigned rpows_max = CHARM_MAX(pmax, kmax + u);
    rpows = CHARM(mpfr_ndarray_malloc)(NBITS, 1, rpows_max + 1);
    if (rpows == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        goto EXIT;
    }


    CHARM(gfm_cap_q_rpows)(rpows, r, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        goto EXIT;
    }
    /* --------------------------------------------------------------------- */






    /* The "R"-terms (Eq. 69 of Bucha et al., 2019b) */
    /* --------------------------------------------------------------------- */
    rwq = CHARM(mpfr_ndarray_malloc)(NBITS, 2, kmax + u + 1, pmax);
    if (rwq == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        goto EXIT;
    }
    CHARM(gfm_cap_q_dr)(rwq, rpows, pmax, kmax + u, NBITS, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        goto EXIT;
    }
    /* --------------------------------------------------------------------- */






    /* Binomial coefficients */
    /* --------------------------------------------------------------------- */
    size_t nbin = CHARM_MAX(kmax + 1, imax + 1);
    binomial = CHARM(mpfr_ndarray_malloc)(NBITS, 2, nbin + 1, nbin + 1);
    if (binomial == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        goto EXIT;
    }


    CHARM(mpfr_binomial)(binomial, nbin, NBITS, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        goto EXIT;
    }
    /* --------------------------------------------------------------------- */






    /* Derivatives of the integral kernels "K" with respect to "r" (Eq. 67 of
     * Bucha et al., 2019b) */
    /* --------------------------------------------------------------------- */
    dkdr = CHARM(mpfr_ndarray_malloc)(NBITS, 3, pmax, imax + 1, kmax + 1);
    if (dkdr == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        goto EXIT;
    }


    CHARM(gfm_cap_q_dkernel)(dkdr, r, rref, psi, rwq, binomial, pmax, kmax,
                             imax, 0, fact, double_fact, NBITS, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        goto EXIT;
    }
    /* --------------------------------------------------------------------- */






    /* Derivatives of the integral kernels "K" with respect to "r" (multiple
     * times) and with respect to "psi" (only once) (Eq. 67 of Bucha et al.,
     * 2019b) */
    /* --------------------------------------------------------------------- */
    if (u == 2)
    {
        dkdrdp = CHARM(mpfr_ndarray_malloc)(NBITS, 3, pmax, imax + 1,
                                            kmax + 1);
        if (dkdrdp == NULL)
        {
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                           CHARM_ERR_MALLOC_FAILURE);
            goto EXIT;
        }


        /* Later in this code, we assume the dimensions of "dkdr" match the
         * dimensions of "dkdrdp", so let's do the check.  There should be no
         * problem, as both arrays are allocated within this function, but one
         * never knows... */
        /* ................................................................. */
        if (dkdr->ndim != dkdrdp->ndim)
            goto EXIT;


        for (size_t l = 0; l < dkdr->ndim; l++)
            if (dkdr->shape[l] != dkdrdp->shape[l])
                goto EXIT;
        /* ................................................................. */


        CHARM(gfm_cap_q_dkernel)(dkdrdp, r, rref, psi, rwq, binomial, pmax,
                                 kmax, imax, 1, fact, double_fact, NBITS,
                                 err);
        if (!CHARM(err_isempty)(err))
        {
            CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
            goto EXIT;
        }
    }


    CHARM(mpfr_ndarray_free)(rwq);
    /* --------------------------------------------------------------------- */






    /* Truncation coefficients (Eq. 66 of Bucha et al., 2019b) */
    /* --------------------------------------------------------------------- */
    /* Compute auxiliary truncation coefficients, from which the sought ones
     * will be obtained */
    /* ..................................................................... */
    /* "Q00" */
    size_t qkpin00_size = (kmax + 1) * pmax * (imax + 1) * (nmax + 1);
    qkpin00 = (mpfr_t *)malloc(qkpin00_size * sizeof(mpfr_t));
    if (qkpin00 == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        goto EXIT;
    }
    for (unsigned v = 0; v < qkpin00_size; v++)
        mpfr_init2(qkpin00[v], NBITS);
    CHARM(gfm_cap_q)(rref, r, psi, nmax, pmax, 0, kmax, imax, zone,
                     CHARM_GFM_Q00, NBITS, qkpin00, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        goto EXIT;
    }
    /* ..................................................................... */


    /* Set all "qkpin" coefficients to zero */
    /* ..................................................................... */
    size_t qkpin_size = CHARM(gfm_cap_nq)(nmax, pmax, kmin, kmax, imax, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        goto EXIT;
    }


    for (size_t v = 0; v < qkpin_size; v++)
        mpfr_set_ui(qkpin[v], 0, MPFR_RNDN);
    /* ..................................................................... */


    /* We use the array of Legendre functions to pre-compute some terms related
     * to Legendre functions.  After this code block, "pnm->data" therefore no
     * longer store Legendre functions. */
    /* ..................................................................... */
    int cj = (zone >= 0) ? 1 : -1;
    mpfr_t sinpsi_cj, tmp1;
    mpfr_inits2(NBITS, sinpsi_cj, tmp1, (mpfr_ptr)NULL);
    mpfr_sin(sinpsi_cj, psi, MPFR_RNDN);
    mpfr_mul_si(sinpsi_cj, sinpsi_cj, cj, MPFR_RNDN);
    size_t idx_pnm1 = pnm->shape[1];


    /* To compute the trunction coefficients, we in fact need only Legendre
     * functions of order "m".  Therefore, we can recompute only the elements
     * of "pnm->data" that are related to this order, as only these elements
     * will be used later.
     *
     * The loop over "n" runs from "1" to avoid division by zero. */
    {
    unsigned long n;
#if CHARM_OPENMP
#   pragma omp parallel for default(shared) private(n)
#endif
    for (n = 1; n <= nmax; n++)
    {
        mpfr_mul(pnm->data[idx_pnm1 + n], pnm->data[idx_pnm1 + n], sinpsi_cj,
                 MPFR_RNDN);
        mpfr_div_ui(pnm->data[idx_pnm1 + n],
                    pnm->data[idx_pnm1 + n], n * (n + 1), MPFR_RNDN);
    }
    }


    size_t idx_pnm2 = 2 * pnm->shape[1];
    if (u == 2)
    {
        {
        unsigned long n;
#if CHARM_OPENMP
#   pragma omp parallel for default(shared) private(n)
#endif
        for (n = 2; n <= nmax; n++)
        {
            mpfr_mul(pnm->data[idx_pnm2 + n], pnm->data[idx_pnm2 + n],
                     sinpsi_cj, MPFR_RNDN);
            mpfr_div_ui(pnm->data[idx_pnm2 + n], pnm->data[idx_pnm2 + n],
                        (n + 2) * (n + 1) * n * (n - 1), MPFR_RNDN);
        }
        }
    }
    /* ..................................................................... */


    /* ..................................................................... */
    int sign;
    size_t idx_dk;
    size_t idx_qkpin, idx_qkpin00;


    for (unsigned p = 1; p <= pmax; p++)
    {
        for (unsigned k = kmin; k <= kmax; k++)
        {
            for (unsigned i = 0; i <= imax; i++)
            {
                idx_qkpin = CHARM(misc_idx_4d)(k - kmin, p - 1, i, 0, pmax,
                                               imax + 1, nmax + 1);


                for (unsigned q = 0; q <= k; q++)
                {
                    sign = ((k - q) % 2) ? -1 : 1;
                    mpfr_mul_si(tmp1,
                                binomial->data[k * binomial->shape[1] + q],
                                sign, MPFR_RNDN);
                    mpfr_mul(tmp1, tmp1, fact->data[k - q + u - 1], MPFR_RNDN);
                    mpfr_div(tmp1, tmp1, rpows->data[k - q + u], MPFR_RNDN);


                    idx_qkpin00 = CHARM(misc_idx_4d)(q, p - 1, i, 0, pmax,
                                                     imax + 1, nmax + 1);
                    idx_dk      = CHARM(misc_idx_4d)(0, p - 1, i, q, pmax,
                                                     imax + 1, kmax + 1);


#if CHARM_OPENMP
#   pragma omp parallel default(shared)
#endif
                    {
                        mpfr_t tmp2;
                        mpfr_init2(tmp2, NBITS);
                        mpfr_t tmp3;
                        mpfr_init2(tmp3, NBITS);


                        if (u == 1)
                        {
                            {
                            unsigned long n;
#if CHARM_OPENMP
#   pragma omp for
#endif
                            for (n = 0; n <= nmax; n++)
                            {
                                mpfr_mul(tmp2, pnm->data[idx_pnm1 + n],
                                         dkdr->data[idx_dk], MPFR_RNDN);
                                mpfr_sub(tmp2, tmp2, qkpin00[idx_qkpin00 + n],
                                         MPFR_RNDN);


                                mpfr_mul(tmp2, tmp2, tmp1, MPFR_RNDN);
                                mpfr_add(qkpin[idx_qkpin + n],
                                         qkpin[idx_qkpin + n], tmp2,
                                         MPFR_RNDN);
                            }
                            }
                        }
                        else if (u == 2)
                        {
                            {
                            unsigned long n;
#if CHARM_OPENMP
#   pragma omp for
#endif
                            for (n = 0; n <= nmax; n++)
                            {
                                mpfr_mul(tmp2, pnm->data[idx_pnm1 + n],
                                         dkdr->data[idx_dk], MPFR_RNDN);
                                mpfr_mul(tmp3, pnm->data[idx_pnm2 + n],
                                         dkdrdp->data[idx_dk], MPFR_RNDN);
                                mpfr_sub(tmp2, tmp2, tmp3, MPFR_RNDN);
                                mpfr_sub(tmp2, tmp2, qkpin00[idx_qkpin00 + n],
                                         MPFR_RNDN);
                                mpfr_neg(tmp2, tmp2, MPFR_RNDN);


                                mpfr_mul(tmp2, tmp2, tmp1, MPFR_RNDN);
                                mpfr_add(qkpin[idx_qkpin + n],
                                         qkpin[idx_qkpin + n], tmp2,
                                         MPFR_RNDN);
                            }
                            }
                        }


                        mpfr_clear(tmp2);
                        mpfr_clear(tmp3);
                        mpfr_free_cache();
                    }
                }
            }
        }
    }


    mpfr_clears(sinpsi_cj, tmp1, (mpfr_ptr)NULL);
    /* ..................................................................... */


    /* Set the non-existing coefficients to "0.0" */
    /* ..................................................................... */
    for (unsigned k = kmin; k <= kmax; k++)
    {
        for (unsigned p = 1; p <= pmax; p++)
        {
            for (unsigned i = 0; i <= imax; i++)
            {
                idx_qkpin = CHARM(misc_idx_4d)(k - kmin, p - 1, i, 0, pmax,
                                               imax + 1, nmax + 1);


                for (unsigned n = 0; n < CHARM_MIN(u, nmax + 1); n++)
                    mpfr_set_ui(qkpin[idx_qkpin + n], 0, MPFR_RNDN);
            }
        }
    }
    /* ..................................................................... */
    /* --------------------------------------------------------------------- */






EXIT:
    /* Clean up */
    /* --------------------------------------------------------------------- */
    if (qkpin00 != NULL)
        for (unsigned v = 0; v < qkpin00_size; v++)
            mpfr_clear(qkpin00[v]);
    free(qkpin00);


    CHARM(mpfr_ndarray_free)(pnm);
    CHARM(mpfr_ndarray_free)(fact);
    CHARM(mpfr_ndarray_free)(double_fact);
    CHARM(mpfr_ndarray_free)(rpows);
    CHARM(mpfr_ndarray_free)(binomial);
    CHARM(mpfr_ndarray_free)(dkdr);
    CHARM(mpfr_ndarray_free)(dkdrdp);


    FLUSH_UNRELEASED_MEMORY;
    /* --------------------------------------------------------------------- */






    return;
}

