/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <mpfr.h>
#if CHARM_OPENMP
#   include <omp.h>
#endif
#include "../prec.h"
#include "../err/err_set.h"
#include "../err/err_propagate.h"
#include "../mpfr/mpfr_check_bits.h"
#include "../mpfr/mpfr_flush_unreleased_memory.h"
#include "../misc/misc_idx_4d.h"
#include "gfm_cap_q_check_radius.h"
#include "gfm_cap_q_check_psi.h"
#include "gfm_check_p.h"
#include "gfm_check_kminkmax.h"
#include "gfm_cap_q_check_type.h"
#include "gfm_cap_q_ref.h"
/* ------------------------------------------------------------------------- */






long CHARM(gfm_cap_q_check_prec)(const mpfr_t rref,
                                 const mpfr_t r,
                                 const mpfr_t psi,
                                 unsigned long nmax,
                                 unsigned pmax,
                                 unsigned kmin,
                                 unsigned kmax,
                                 unsigned imax,
                                 unsigned type,
                                 mpfr_prec_t NBITS,
                                 mpfr_prec_t NBITSREF,
                                 CHARM(err) *err)
{
    /* Some initializations */
    /* --------------------------------------------------------------------- */
    long ret = -1;


    mpfr_t *qkpin_in  = NULL;
    mpfr_t *qkpin_out = NULL;
    mpfr_t *skpin     = NULL;
    long *ret_team    = NULL;
    /* --------------------------------------------------------------------- */






    /* Checks */
    /* --------------------------------------------------------------------- */
    CHARM(mpfr_check_bits)(NBITS, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        goto EXIT;
    }


    CHARM(mpfr_check_bits)(NBITSREF, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        ret = -2;
        goto EXIT;
    }


    if (NBITSREF < NBITS)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"NBITSREF\" cannot be smaller than \"NBITS\".");
        ret = -3;
        goto EXIT;
    }


    CHARM(gfm_cap_q_check_radius)(rref, NBITS, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        ret = -4;
        goto EXIT;
    }


    CHARM(gfm_cap_q_check_radius)(r, NBITS, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        ret = -5;
        goto EXIT;
    }


    CHARM(gfm_cap_q_check_psi)(psi, NBITS, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        ret = -6;
        goto EXIT;
    }


    CHARM(gfm_check_p)(pmax, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        ret = -7;
        goto EXIT;
    }


    CHARM(gfm_check_kminkmax)(kmin, kmax, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        ret = -8;
        goto EXIT;
    }


    CHARM(gfm_cap_q_check_type)(type, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        ret = -9;
        goto EXIT;
    }
    /* --------------------------------------------------------------------- */






    /* --------------------------------------------------------------------- */
    size_t size = CHARM(gfm_cap_nq)(nmax, pmax, kmin, kmax, imax, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        ret = -10;
        goto EXIT;
    }


    qkpin_in = (mpfr_t *)malloc(size * sizeof(mpfr_t));
    if (qkpin_in == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        ret = -11;
        goto EXIT;
    }
    for (size_t v = 0; v < size; v++)
        mpfr_init2(qkpin_in[v], NBITS);


    qkpin_out = (mpfr_t *)malloc(size * sizeof(mpfr_t));
    if (qkpin_out == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        ret = -12;
        goto EXIT;
    }
    for (size_t v = 0; v < size; v++)
        mpfr_init2(qkpin_out[v], NBITS);


    skpin = (mpfr_t *)malloc(size * sizeof(mpfr_t));
    if (skpin == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        ret = -13;
        goto EXIT;
    }
    for (size_t v = 0; v < size; v++)
        mpfr_init2(skpin[v], NBITS);


    /* Compute the near- and far-zone truncation coefficients */
    CHARM(gfm_cap_q)(rref, r, psi, nmax, pmax, kmin, kmax, imax, 1, type,
                     NBITS, qkpin_in, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        ret = -14;
        goto EXIT;
    }
    CHARM(gfm_cap_q)(rref, r, psi, nmax, pmax, kmin, kmax, imax, -1, type,
                     NBITS, qkpin_out, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        ret = -15;
        goto EXIT;
    }


    /* Compute the reference near- plus far-zone coefficients */
    CHARM(gfm_cap_q_ref)(rref, r, nmax, pmax, kmin, kmax, imax, type, NBITSREF,
                         skpin, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        ret = -16;
        goto EXIT;
    }


    /* This is the default value returned if all digits agree to the last bit
     * */
    ret = LONG_MAX;


    /* Initialize "ret" but for all threads in an OpenMP team, "ret_team" */
    int nthreads;
#if CHARM_OPENMP
    nthreads = omp_get_max_threads();
#else
    nthreads = 1;
#endif
    ret_team = (long *)calloc(nthreads, sizeof(long));
    if (ret_team == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        ret = -17;
        goto EXIT;
    }
    for (int myid = 0; myid < nthreads; myid++)
        ret_team[myid] = LONG_MAX;


    size_t idx;
    for (unsigned k = kmin; k <= kmax; k++)
    {
        for (unsigned p = 1; p <= pmax; p++)
        {
            for (unsigned i = 0; i <= imax; i++)
            {
                idx = CHARM(misc_idx_4d)(k - kmin, p - 1, i, 0, pmax, imax + 1,
                                         nmax + 1);


#if CHARM_OPENMP
#pragma omp parallel default(shared)
#endif
                {
                    int myid;
#if CHARM_OPENMP
                    myid = omp_get_thread_num();
#else
                    myid = 0;
#endif
                    long delta;
                    mpfr_t tmp;
                    mpfr_init2(tmp, NBITS);


                    {
                    unsigned long n;
#if CHARM_OPENMP
#pragma omp for
#endif
                    for (n = 0; n <= nmax; n++)
                    {
                        mpfr_add(tmp, qkpin_in[idx + n], qkpin_out[idx + n],
                                 MPFR_RNDN);
                        if (mpfr_zero_p(skpin[idx + n]))
                        {
                            /* In this case, the sum of near- and far-zone
                             * coefficients should be zero to the last bit.
                             * But since it is unlikely to get an exact zero
                             * when summing the near- and far-zone
                             * coefficients, we need a different way to assess
                             * the accuracy.  We use here a simple rule saying
                             * that if the reference value is zero (to the last
                             * bit), then the accuracy is derived from the
                             * exponent of the value that is being validated.
                             * For instance, if the reference value is
                             *
                             *      "xref = 0.0"
                             *
                             * and the value being validated is
                             *
                             *      "x = 1.232e-15",
                             *
                             * the accuracy would be 15 digits.
                             *
                             * This means we need to compute the absolute value
                             * of "tmp" and then its 10-base logarithm.  This
                             * is done below outside this branch.
                             *
                             * This accuracy assessment is likely to be
                             * unsuitable in many situations, but here, it
                             * allows us to detect values that are wrong in
                             * a relative sense. */
                            ;
                        }
                        else
                        {
                            mpfr_sub(tmp, tmp, skpin[idx + n], MPFR_RNDN);
                            mpfr_div(tmp, tmp, skpin[idx + n], MPFR_RNDN);
                        }


                        mpfr_abs(tmp, tmp, MPFR_RNDN);
                        mpfr_log10(tmp, tmp, MPFR_RNDN);
                        mpfr_neg(tmp, tmp, MPFR_RNDN);


                        if (mpfr_inf_p(tmp))
                            /* This is an exact match to the last bit, so let's
                             * return "LONG_MAX" */
                            delta = LONG_MAX;
                        else
                            delta = lround(mpfr_get_d(tmp, MPFR_RNDN));


                        if (delta < ret_team[myid])
                            ret_team[myid] = delta;
                    }
                    }


                    mpfr_clear(tmp);
                    mpfr_free_cache();
                }
            }
        }
    }


    /* Find the minimum value in "ret_team" */
    ret = ret_team[0];
    for (int myid = 1; myid < nthreads; myid++)
        if (ret_team[myid] < ret)
            ret = ret_team[myid];


    if (ret < 0)
        ret = 0;
    /* --------------------------------------------------------------------- */






EXIT:
    /* Clean up */
    /* --------------------------------------------------------------------- */
    if (qkpin_in != NULL)
        for (size_t v = 0; v < size; v++)
            mpfr_clear(qkpin_in[v]);
    free(qkpin_in);


    if (qkpin_out != NULL)
        for (size_t v = 0; v < size; v++)
            mpfr_clear(qkpin_out[v]);
    free(qkpin_out);


    if (skpin != NULL)
        for (size_t v = 0; v < size; v++)
            mpfr_clear(skpin[v]);
    free(skpin);
    free(ret_team);


    mpfr_free_cache();
    FLUSH_UNRELEASED_MEMORY;
    /* --------------------------------------------------------------------- */






    return ret;
}

