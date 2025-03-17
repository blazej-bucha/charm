/* Header files */
/* ------------------------------------------------------------------------- */
#include <stdio.h>
#include <stdlib.h>
#include <mpfr.h>
#include "../err/err_set.h"
#include "../err/err_propagate.h"
#include "../mpfr/mpfr_ndarray.h"
#include "../mpfr/mpfr_ndarray_free.h"
#include "../mpfr/mpfr_ndarray_malloc.h"
#include "../mpfr/mpfr_binomial.h"
#include "../mpfr/mpfr_check_bits.h"
#include "../mpfr/mpfr_fact.h"
#include "../mpfr/mpfr_double_fact.h"
#include "../mpfr/mpfr_flush_unreleased_memory.h"
#include "../misc/misc_idx_4d.h"
#include "gfm_cap_q_rpows.h"
#include "gfm_cap_q_dr.h"
#include "gfm_cap_q_aps.h"
#include "gfm_cap_q_norm_dist.h"
#include "gfm_cap_q_dnorm_dist.h"
#include "gfm_cap_q_check_radius.h"
#include "gfm_cap_q_check_psi.h"
#include "gfm_check_p.h"
#include "gfm_check_kminkmax.h"
#include "gfm_cap_qu0.h"
/* ------------------------------------------------------------------------- */






void CHARM(gfm_cap_qu0)(const mpfr_t rref,
                        const mpfr_t r,
                        const mpfr_t psi,
                        unsigned long nmax,
                        unsigned pmax,
                        unsigned kmin,
                        unsigned kmax,
                        unsigned imax,
                        int zone,
                        mpfr_prec_t NBITS,
                        mpfr_t *qkpin,
                        CHARM(err) *err)
{
    /* Checks */
    /* --------------------------------------------------------------------- */
    CHARM(mpfr_check_bits)(NBITS, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return;
    }


    CHARM(gfm_cap_q_check_radius)(rref, NBITS, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return;
    }


    CHARM(gfm_cap_q_check_radius)(r, NBITS, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return;
    }


    CHARM(gfm_cap_q_check_psi)(psi, NBITS, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
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
    mpfr_t t, t2, cospsi, ul, uu, u0, gu, gl, tmp1, tmp2, tmp3, ttmp;
    mpfr_t tmi1, tmi2, tp2, gammat, it, iim1t, tpi1, tpi2, t2p1, tt2;
    mpfr_t aps, factp;
    mpfr_inits2(NBITS, t, t2, cospsi,
                ul, uu, u0,
                gu, gl,
                tmp1, tmp2, tmp3, ttmp, tmi1, tmi2, tp2,
                gammat, it, iim1t,
                tpi1, tpi2, t2p1, tt2,
                aps, factp,
                (mpfr_ptr)NULL);


    mpfr_ndarray *fact        = NULL;
    mpfr_ndarray *double_fact = NULL;
    mpfr_ndarray *gammau      = NULL;
    mpfr_ndarray *gammal      = NULL;
    mpfr_ndarray *m           = NULL;
    mpfr_ndarray *binomial    = NULL;
    mpfr_ndarray *bin         = NULL;
    mpfr_ndarray *beta0       = NULL;
    mpfr_ndarray *beta1uu     = NULL;
    mpfr_ndarray *beta1ul     = NULL;
    mpfr_ndarray *ain         = NULL;
    mpfr_ndarray *alpha0uu    = NULL;
    mpfr_ndarray *alpha0ul    = NULL;
    mpfr_ndarray *alpha1uu    = NULL;
    mpfr_ndarray *alpha1ul    = NULL;
    mpfr_ndarray *jin         = NULL;
    mpfr_ndarray *ti          = NULL;
    mpfr_ndarray *bell        = NULL;
    mpfr_ndarray *gin         = NULL;
    mpfr_ndarray *rpows       = NULL;
    mpfr_ndarray *rps         = NULL;
    /* --------------------------------------------------------------------- */






    /* Prepare arrays with factorials and double factorials */
    /* --------------------------------------------------------------------- */
    unsigned dmax = pmax + kmax + 1;


    fact = CHARM(mpfr_ndarray_malloc)(NBITS, 1, dmax + 2);
    if (fact == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        goto EXIT;
    }
    {
    unsigned h;
#if HAVE_OPENMP
#pragma omp parallel for default(shared) private(h)
#endif
    for (h = 0; h < fact->size; h++)
        CHARM(mpfr_fact)(h, fact->data[h]);
    }


    double_fact = CHARM(mpfr_ndarray_malloc)(NBITS, 1, 2 * dmax + 1);
    if (double_fact == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        goto EXIT;
    }
    {
    unsigned h;
#if HAVE_OPENMP
#pragma omp parallel for default(shared) private(h)
#endif
    for (h = 0; h < double_fact->size; h++)
        CHARM(mpfr_double_fact)(h, double_fact->data[h]);
    }
    /* --------------------------------------------------------------------- */






    /* Compute the derivatives of "1 / g(t, uu)" and of "1 / g(t, ul)" */
    /* --------------------------------------------------------------------- */
    mpfr_div(t, rref, r, MPFR_RNDN);
    mpfr_mul(t2, t, t, MPFR_RNDN);


    mpfr_cos(cospsi, psi, MPFR_RNDN);
    if (zone >= 0)  /* Near-zone */
    {
        mpfr_set(ul, cospsi, MPFR_RNDN);
        mpfr_set_ui(uu, 1, MPFR_RNDN);
    }
    else
    {
        mpfr_set_si(ul, -1, MPFR_RNDN);
        mpfr_set(uu, cospsi, MPFR_RNDN);
    }


    /* "1 / g(t, uu)" */
    gammau = CHARM(mpfr_ndarray_malloc)(NBITS, 1, dmax + 1);
    if (gammau == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        goto EXIT;
    }


    /* "1 / g(t, ul)" */
    gammal = CHARM(mpfr_ndarray_malloc)(NBITS, 1, dmax + 1);
    if (gammal == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        goto EXIT;
    }


    CHARM(gfm_cap_q_norm_dist)(t, t2, uu, gu);
    CHARM(gfm_cap_q_norm_dist)(t, t2, ul, gl);


    /* Compute the derivatives */
    CHARM(gfm_cap_q_dnorm_dist)(gammau, gu, t, uu, fact, double_fact, dmax,
                                NBITS, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        goto EXIT;
    }
    CHARM(gfm_cap_q_dnorm_dist)(gammal, gl, t, ul, fact, double_fact, dmax,
                                NBITS, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        goto EXIT;
    }
    /* --------------------------------------------------------------------- */






    /* Integrals of Legendre polynomials (Eqs. 69 -- 71 of Bucha et al., 2019)
     * */
    /* --------------------------------------------------------------------- */
    mpfr_ndarray *gamma0;
    if (zone >= 0)  /* Near-zone */
    {
        mpfr_set(u0, ul, MPFR_RNDN);
        gamma0 = gammal;
    }
    else
    {
        mpfr_set(u0, uu, MPFR_RNDN);
        gamma0 = gammau;
    }


    m = CHARM(mpfr_ndarray_malloc)(NBITS, 1, nmax + 1);
    if (m == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        goto EXIT;
    }


    if (zone >= 0)  /* Near-zone */
    {
        mpfr_ui_sub(m->data[0], 1, ul, MPFR_RNDN);


        if (nmax > 0)
        {
            mpfr_mul(m->data[1], ul, ul, MPFR_RNDN);
            mpfr_ui_sub(m->data[1], 1, m->data[1], MPFR_RNDN);
            mpfr_div_ui(m->data[1], m->data[1], 2, MPFR_RNDN);
        }
    }
    else
    {
        mpfr_add_ui(m->data[0], uu, 1, MPFR_RNDN);


        if (nmax > 0)
        {
            mpfr_mul(m->data[1], uu, uu, MPFR_RNDN);
            mpfr_sub_ui(m->data[1], m->data[1], 1, MPFR_RNDN);
            mpfr_div_ui(m->data[1], m->data[1], 2, MPFR_RNDN);
        }
    }


    for (unsigned long n = 2; n <= nmax; n++)
    {
        mpfr_mul(tmp1, m->data[n - 1], u0, MPFR_RNDN);
        mpfr_mul_ui(tmp1, tmp1, 2 * n - 1, MPFR_RNDN);
        mpfr_mul_ui(tmp2, m->data[n - 2], n - 2, MPFR_RNDN);

        mpfr_sub(m->data[n], tmp1, tmp2, MPFR_RNDN);
        mpfr_div_ui(m->data[n], m->data[n], n + 1, MPFR_RNDN);
    }
    /* --------------------------------------------------------------------- */






    /* Compute binomial coefficients */
    /* --------------------------------------------------------------------- */
    unsigned nbinomial = CHARM_MAX(dmax + 1, imax + 1);
    binomial = CHARM(mpfr_ndarray_malloc)(NBITS, 2, nbinomial + 1,
                                          nbinomial + 1);
    if (binomial == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        goto EXIT;
    }
    CHARM(mpfr_binomial)(binomial, nbinomial, NBITS, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        goto EXIT;
    }
    /* --------------------------------------------------------------------- */






    /* Compute the B-terms (Eqs. 64 -- 66 of Bucha et al., 2019) */
    /* --------------------------------------------------------------------- */
    bin = CHARM(mpfr_ndarray_malloc)(NBITS, 2, dmax + 1, nmax + 1);
    if (bin == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        goto EXIT;
    }
    for (unsigned h = 0; h < bin->size; h++)
        mpfr_set_ui(bin->data[h], 0, MPFR_RNDN);


    beta0 = CHARM(mpfr_ndarray_malloc)(NBITS, 1, dmax + 1);
    if (beta0 == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        goto EXIT;
    }


    beta1uu = CHARM(mpfr_ndarray_malloc)(NBITS, 1, dmax + 1);
    if (beta1uu == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        goto EXIT;
    }


    beta1ul = CHARM(mpfr_ndarray_malloc)(NBITS, 1, dmax + 1);
    if (beta1ul == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        goto EXIT;
    }


    /* "(1 + t2) / t" */
    mpfr_add_ui(ttmp, t2, 1, MPFR_RNDN);
    mpfr_div(ttmp, ttmp, t, MPFR_RNDN);


    int signi  =  1;
    int signi1 = -1;
    for (unsigned h = 0; h <= dmax; h++)
    {
        mpfr_pow_si(tmi1, t, -((long)(h + 1)), MPFR_RNDN);
        mpfr_pow_si(tmi2, t, -((long)(h + 2)), MPFR_RNDN);
        mpfr_ui_div(tp2, 1, t2, MPFR_RNDN);
        mpfr_mul_si(beta0->data[h], fact->data[h], signi, MPFR_RNDN);
        mpfr_mul(beta0->data[h], beta0->data[h], tmi1, MPFR_RNDN);


        if (h == 0)
        {
            mpfr_div(beta1uu->data[h], uu, t, MPFR_RNDN);
            mpfr_sub(beta1uu->data[h], tp2, beta1uu->data[h], MPFR_RNDN);
            mpfr_add_ui(beta1uu->data[h], beta1uu->data[h], 1, MPFR_RNDN);


            mpfr_div(beta1ul->data[h], ul, t, MPFR_RNDN);
            mpfr_sub(beta1ul->data[h], tp2, beta1ul->data[h], MPFR_RNDN);
            mpfr_add_ui(beta1ul->data[h], beta1ul->data[h], 1, MPFR_RNDN);
        }
        else
        {
            mpfr_mul_si(tmp1, fact->data[h + 1], signi, MPFR_RNDN);
            mpfr_mul(tmp1, tmp1, tmi2, MPFR_RNDN);
            mpfr_mul_si(tmp2, fact->data[h], signi1, MPFR_RNDN);
            mpfr_mul(tmp2, tmp2, tmi1, MPFR_RNDN);
            mpfr_mul(tmp2, tmp2, uu, MPFR_RNDN);
            mpfr_add(beta1uu->data[h], tmp1, tmp2, MPFR_RNDN);


            mpfr_mul_si(tmp2, fact->data[h], signi1, MPFR_RNDN);
            mpfr_mul(tmp2, tmp2, tmi1, MPFR_RNDN);
            mpfr_mul(tmp2, tmp2, ul, MPFR_RNDN);
            mpfr_add(beta1ul->data[h], tmp1, tmp2, MPFR_RNDN);
        }


        signi  = -signi;
        signi1 = -signi1;
    }


    size_t idx;
    for (unsigned h = 0; h <= dmax; h++)
    {
        idx = h * bin->shape[1];
        mpfr_div(gammat, gamma0->data[h], t, MPFR_RNDN);
        mpfr_ui_div(it, h, t, MPFR_RNDN);


        if (h == 0)
        {
            mpfr_mul(tmp1, t, gu, MPFR_RNDN);
            mpfr_ui_div(tmp1, 1, tmp1, MPFR_RNDN);
            mpfr_mul(tmp2, t, gl, MPFR_RNDN);
            mpfr_ui_div(tmp2, 1, tmp2, MPFR_RNDN);


            mpfr_sub(bin->data[0], tmp1, tmp2, MPFR_RNDN);


            if (nmax > 0)
            {
                mpfr_mul(tmp1, t2, gu, MPFR_RNDN);
                mpfr_ui_div(tmp1, 1, tmp1, MPFR_RNDN);
                mpfr_mul(tmp3, t, uu, MPFR_RNDN);
                mpfr_ui_sub(tmp3, 1, tmp3, MPFR_RNDN);
                mpfr_add(tmp3, tmp3, t2, MPFR_RNDN);
                mpfr_mul(tmp1, tmp1, tmp3, MPFR_RNDN);


                mpfr_mul(tmp2, t2, gl, MPFR_RNDN);
                mpfr_ui_div(tmp2, 1, tmp2, MPFR_RNDN);
                mpfr_mul(tmp3, t, ul, MPFR_RNDN);
                mpfr_ui_sub(tmp3, 1, tmp3, MPFR_RNDN);
                mpfr_add(tmp3, tmp3, t2, MPFR_RNDN);
                mpfr_mul(tmp2, tmp2, tmp3, MPFR_RNDN);


                mpfr_sub(bin->data[1], tmp1, tmp2, MPFR_RNDN);
            }
        }
        else
        {
            for (unsigned j = 0; j <= h; j++)
            {
                mpfr_mul(tmp1, beta0->data[h - j], gammau->data[j], MPFR_RNDN);
                mpfr_mul(tmp2, beta0->data[h - j], gammal->data[j], MPFR_RNDN);
                mpfr_sub(tmp1, tmp1, tmp2, MPFR_RNDN);
                mpfr_mul(tmp1, binomial->data[h * binomial->shape[1] + j],
                         tmp1, MPFR_RNDN);
                mpfr_add(bin->data[idx], bin->data[idx], tmp1, MPFR_RNDN);


                if (nmax > 0)
                {
                    mpfr_mul(tmp1, beta1uu->data[h - j], gammau->data[j],
                             MPFR_RNDN);
                    mpfr_mul(tmp2, beta1ul->data[h - j], gammal->data[j],
                             MPFR_RNDN);
                    mpfr_sub(tmp1, tmp1, tmp2, MPFR_RNDN);
                    mpfr_mul(tmp1, binomial->data[h * binomial->shape[1] + j],
                             tmp1, MPFR_RNDN);
                    mpfr_add(bin->data[idx + 1], bin->data[idx + 1], tmp1,
                             MPFR_RNDN);
                }
            }
        }


        if (h >= 1)
            mpfr_ui_div(iim1t, h * (h - 1), t, MPFR_RNDN);


        for (unsigned long n = 2; n <= nmax; n++)
        {
            mpfr_mul(bin->data[idx + n], ttmp, bin->data[idx + n - 1],
                     MPFR_RNDN);
            mpfr_sub(bin->data[idx + n], bin->data[idx + n],
                     bin->data[idx + n - 2], MPFR_RNDN);
            mpfr_mul(tmp1, m->data[n - 1], gammat, MPFR_RNDN);
            mpfr_sub(bin->data[idx + n], bin->data[idx + n], tmp1, MPFR_RNDN);


            if (h >= 1)
            {
                mpfr_add(tmp1, bin->data[(h - 1) * (nmax + 1) + n],
                         bin->data[(h - 1) * (nmax + 1) + n - 2], MPFR_RNDN);
                mpfr_mul(tmp1, tmp1, it, MPFR_RNDN);
                mpfr_sub(bin->data[idx + n], bin->data[idx + n], tmp1,
                         MPFR_RNDN);
                mpfr_mul_ui(tmp1, bin->data[(h - 1) * (nmax + 1) + n - 1],
                            2 * h, MPFR_RNDN);
                mpfr_add(bin->data[idx + n], bin->data[idx + n], tmp1,
                         MPFR_RNDN);
            }


            if (h >= 2)
            {
                mpfr_mul(tmp1, iim1t, bin->data[(h - 2) * (nmax + 1) + n - 1],
                         MPFR_RNDN);
                mpfr_add(bin->data[idx + n], bin->data[idx + n], tmp1,
                         MPFR_RNDN);
            }

        }
    }
    /* --------------------------------------------------------------------- */






    /* Compute the A-terms (Eqs. 56 -- 58 of Bucha et al., 2019) */
    /* --------------------------------------------------------------------- */
    ain = CHARM(mpfr_ndarray_malloc)(NBITS, 2, dmax + 1, nmax + 1);
    if (ain == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        goto EXIT;
    }
    for (unsigned h = 0; h < ain->size; h++)
        mpfr_set_ui(ain->data[h], 0, MPFR_RNDN);


    alpha0uu = CHARM(mpfr_ndarray_malloc)(NBITS, 1, dmax + 1);
    if (alpha0uu == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        goto EXIT;
    }


    alpha0ul = CHARM(mpfr_ndarray_malloc)(NBITS, 1, dmax + 1);
    if (alpha0ul == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        goto EXIT;
    }


    alpha1uu = CHARM(mpfr_ndarray_malloc)(NBITS, 1, dmax + 1);
    if (alpha1uu == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        goto EXIT;
    }


    alpha1ul = CHARM(mpfr_ndarray_malloc)(NBITS, 1, dmax + 1);
    if (alpha1ul == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        goto EXIT;
    }


    mpfr_add_ui(t2p1, t2, 1, MPFR_RNDN);
    mpfr_mul_ui(tt2, t, 2, MPFR_RNDN);


    signi  = -1;
    signi1 = 1;
    for (unsigned h = 0; h <= dmax; h++)
    {
        mpfr_pow_ui(tpi1, t, h + 1, MPFR_RNDN);
        mpfr_pow_ui(tpi2, t, h + 2, MPFR_RNDN);
        signi  = -signi;
        signi1 = -signi1;


        if (h == 0)
        {
            mpfr_ui_div(tmp1, 1, t, MPFR_RNDN);
            mpfr_add(tmp1, tmp1, t, MPFR_RNDN);
            mpfr_set(alpha0uu->data[h], tmp1, MPFR_RNDN);
            mpfr_set(alpha0ul->data[h], tmp1, MPFR_RNDN);
            mpfr_mul_ui(tmp1, uu, 2, MPFR_RNDN);
            mpfr_sub(alpha0uu->data[h], alpha0uu->data[h], tmp1, MPFR_RNDN);
            mpfr_mul_ui(tmp1, ul, 2, MPFR_RNDN);
            mpfr_sub(alpha0ul->data[h], alpha0ul->data[h], tmp1, MPFR_RNDN);


            mpfr_set(tmp1, t2, MPFR_RNDN);
            mpfr_ui_div(tmp2, 1, t2, MPFR_RNDN);
            mpfr_add(tmp1, tmp1, tmp2, MPFR_RNDN);
            mpfr_add_ui(tmp1, tmp1, 2, MPFR_RNDN);
            mpfr_set(alpha1uu->data[h], tmp1, MPFR_RNDN);
            mpfr_set(alpha1ul->data[h], tmp1, MPFR_RNDN);
            mpfr_mul(tmp1, t, uu, MPFR_RNDN);
            mpfr_sub(alpha1uu->data[h], alpha1uu->data[h], tmp1, MPFR_RNDN);
            mpfr_div(tmp1, uu, t, MPFR_RNDN);
            mpfr_sub(alpha1uu->data[h], alpha1uu->data[h], tmp1, MPFR_RNDN);
            mpfr_mul(tmp1, uu, uu, MPFR_RNDN);
            mpfr_mul_ui(tmp1, tmp1, 2, MPFR_RNDN);
            mpfr_sub(alpha1uu->data[h], alpha1uu->data[h], tmp1, MPFR_RNDN);
            mpfr_mul(tmp1, t, ul, MPFR_RNDN);
            mpfr_sub(alpha1ul->data[h], alpha1ul->data[h], tmp1, MPFR_RNDN);
            mpfr_div(tmp1, ul, t, MPFR_RNDN);
            mpfr_sub(alpha1ul->data[h], alpha1ul->data[h], tmp1, MPFR_RNDN);
            mpfr_mul(tmp1, ul, ul, MPFR_RNDN);
            mpfr_mul_ui(tmp1, tmp1, 2, MPFR_RNDN);
            mpfr_sub(alpha1ul->data[h], alpha1ul->data[h], tmp1, MPFR_RNDN);
        }
        else
        {
            mpfr_mul_si(alpha0uu->data[h], fact->data[h], signi, MPFR_RNDN);
            mpfr_div(alpha0uu->data[h], alpha0uu->data[h], tpi1, MPFR_RNDN);
            mpfr_set(alpha0ul->data[h], alpha0uu->data[h], MPFR_RNDN);


            mpfr_mul_si(tmp1, fact->data[h + 1], signi, MPFR_RNDN);
            mpfr_div(tmp1, tmp1, tpi2, MPFR_RNDN);
            mpfr_set(alpha1uu->data[h], tmp1, MPFR_RNDN);
            mpfr_set(alpha1ul->data[h], tmp1, MPFR_RNDN);


            mpfr_mul_si(tmp1, fact->data[h], signi1, MPFR_RNDN);
            mpfr_mul(tmp1, tmp1, uu, MPFR_RNDN);
            mpfr_div(tmp1, tmp1, tpi1, MPFR_RNDN);
            mpfr_add(alpha1uu->data[h], alpha1uu->data[h], tmp1, MPFR_RNDN);


            mpfr_mul_si(tmp1, fact->data[h], signi1, MPFR_RNDN);
            mpfr_mul(tmp1, tmp1, ul, MPFR_RNDN);
            mpfr_div(tmp1, tmp1, tpi1, MPFR_RNDN);
            mpfr_add(alpha1ul->data[h], alpha1ul->data[h], tmp1, MPFR_RNDN);


            if (h == 1)
            {
                mpfr_add_ui(alpha0uu->data[h], alpha0uu->data[h], 1,
                            MPFR_RNDN);
                mpfr_add_ui(alpha0ul->data[h], alpha0ul->data[h], 1,
                            MPFR_RNDN);


                mpfr_mul_ui(tmp1, t, 2, MPFR_RNDN);
                mpfr_add(alpha1uu->data[h], alpha1uu->data[h], tmp1,
                         MPFR_RNDN);
                mpfr_add(alpha1ul->data[h], alpha1ul->data[h], tmp1,
                         MPFR_RNDN);
                mpfr_sub(alpha1uu->data[h], alpha1uu->data[h], uu, MPFR_RNDN);
                mpfr_sub(alpha1ul->data[h], alpha1ul->data[h], ul, MPFR_RNDN);
            }
            else if (h == 2)
            {
                mpfr_add_ui(alpha1uu->data[h], alpha1uu->data[h], 2,
                            MPFR_RNDN);
                mpfr_add_ui(alpha1ul->data[h], alpha1ul->data[h], 2,
                            MPFR_RNDN);
            }
        }
    }


    for (unsigned h = 0; h <= dmax; h++)
    {
        idx = h * (nmax + 1);


        if (h == 0)
        {
            mpfr_div(ain->data[idx], gl, t, MPFR_RNDN);
            mpfr_div(tmp1, gu, t, MPFR_RNDN);
            mpfr_sub(ain->data[idx], ain->data[idx], tmp1, MPFR_RNDN);


            if (nmax > 0)
            {
                mpfr_mul_ui(tmp1, t2, 3, MPFR_RNDN);
                mpfr_div(tmp1, gl, tmp1, MPFR_RNDN);
                mpfr_mul(tmp2, t, ul, MPFR_RNDN);
                mpfr_add_ui(tmp2, tmp2, 1, MPFR_RNDN);
                mpfr_add(tmp2, tmp2, t2, MPFR_RNDN);
                mpfr_mul(ain->data[idx + 1], tmp1, tmp2, MPFR_RNDN);
                mpfr_mul_ui(tmp1, t2, 3, MPFR_RNDN);
                mpfr_div(tmp1, gu, tmp1, MPFR_RNDN);
                mpfr_mul(tmp2, t, uu, MPFR_RNDN);
                mpfr_add_ui(tmp2, tmp2, 1, MPFR_RNDN);
                mpfr_add(tmp2, tmp2, t2, MPFR_RNDN);
                mpfr_mul(tmp1, tmp1, tmp2, MPFR_RNDN);
                mpfr_sub(ain->data[idx + 1], ain->data[idx + 1], tmp1,
                         MPFR_RNDN);
            }
        }
        else
        {
            for (unsigned j = 0; j <= h; j++)
            {
                mpfr_mul(tmp1, alpha0uu->data[h - j], gammau->data[j],
                         MPFR_RNDN);
                mpfr_mul(tmp2, alpha0ul->data[h - j], gammal->data[j],
                         MPFR_RNDN);
                mpfr_sub(tmp1, tmp1, tmp2, MPFR_RNDN);
                mpfr_mul(tmp1, tmp1,
                         binomial->data[h * binomial->shape[1] + j],
                         MPFR_RNDN);
                mpfr_sub(ain->data[idx], ain->data[idx], tmp1, MPFR_RNDN);


                if (nmax > 0)
                {
                    mpfr_mul(tmp1, alpha1uu->data[h - j], gammau->data[j],
                             MPFR_RNDN);
                    mpfr_mul(tmp2, alpha1ul->data[h - j], gammal->data[j],
                             MPFR_RNDN);
                    mpfr_sub(tmp1, tmp1, tmp2, MPFR_RNDN);
                    mpfr_mul(tmp1, tmp1,
                             binomial->data[h * binomial->shape[1] + j],
                             MPFR_RNDN);
                    mpfr_sub(ain->data[idx + 1], ain->data[idx + 1], tmp1,
                             MPFR_RNDN);
                }
            }


            if (nmax > 0)
                mpfr_div_ui(ain->data[idx + 1], ain->data[idx + 1], 3,
                            MPFR_RNDN);
        }


#if HAVE_OPENMP
#pragma omp parallel default(shared)
#endif
        {
            mpfr_t tmp4, tmp5;
            mpfr_inits2(NBITS, tmp4, tmp5, (mpfr_ptr)NULL);


            {
            unsigned long n;
#if HAVE_OPENMP
#pragma omp for
#endif
            for (n = 2; n <= nmax; n++)
            {
                mpfr_mul(tmp4, t2p1, bin->data[idx + n], MPFR_RNDN);
                mpfr_mul(tmp5, tt2, bin->data[idx + n - 1], MPFR_RNDN);
                mpfr_sub(tmp4, tmp5, tmp4, MPFR_RNDN);
                mpfr_mul_ui(tmp5, m->data[n], 2 * (n + 1), MPFR_RNDN);
                mpfr_mul(tmp5, tmp5, gamma0->data[h], MPFR_RNDN);
                mpfr_add(ain->data[idx + n], tmp4, tmp5, MPFR_RNDN);


                if (h >= 1)
                {
                    mpfr_mul_ui(tmp4, bin->data[(h - 1) * (nmax + 1)  + n], h,
                                MPFR_RNDN);
                    mpfr_mul(tmp4, tmp4, tt2, MPFR_RNDN);
                    mpfr_sub(ain->data[idx + n], ain->data[idx + n], tmp4,
                             MPFR_RNDN);
                    mpfr_mul_ui(tmp4, bin->data[(h - 1) * (nmax + 1) + n - 1],
                                2 * h, MPFR_RNDN);
                    mpfr_add(ain->data[idx + n], ain->data[idx + n], tmp4,
                             MPFR_RNDN);
                }


                if (h >= 2)
                {
                    mpfr_mul_ui(tmp4, bin->data[(h - 2) * (nmax + 1) + n],
                                h * (h - 1), MPFR_RNDN);
                    mpfr_sub(ain->data[idx + n], ain->data[idx + n], tmp4,
                             MPFR_RNDN);
                }


                mpfr_div_ui(ain->data[idx + n], ain->data[idx + n], 2 * n + 1,
                            MPFR_RNDN);
            }
            }


            mpfr_clears(tmp4, tmp5, (mpfr_ptr)NULL);
            mpfr_free_cache();
        }
    }
    /* --------------------------------------------------------------------- */







    /* Compute the J-terms (Eq. 54 of Bucha et al., 2019) */
    /* --------------------------------------------------------------------- */
    jin = CHARM(mpfr_ndarray_malloc)(NBITS, 2, dmax + 1, nmax + 1);
    if (jin == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        goto EXIT;
    }


    /* h == 0 */
    {
    unsigned long n;
#if HAVE_OPENMP
#pragma omp parallel for default(shared) private(n)
#endif
    for (n = 0; n <= nmax; n++)
        mpfr_mul(jin->data[n], t, ain->data[n], MPFR_RNDN);
    }


    /* h > 0 */
    for (unsigned h = 1; h <= dmax; h++)
    {
        idx = h * (nmax + 1);


#if HAVE_OPENMP
#pragma omp parallel default(shared)
#endif
        {
            mpfr_t tmp4;
            mpfr_init2(tmp4, NBITS);


            {
            unsigned long n;
#if HAVE_OPENMP
#pragma omp for
#endif
            for (n = 0; n <= nmax; n++)
            {
                mpfr_mul(jin->data[idx + n], t, ain->data[idx + n], MPFR_RNDN);
                mpfr_mul_ui(tmp4, ain->data[(h - 1) * (nmax + 1) + n], h,
                            MPFR_RNDN);
                mpfr_add(jin->data[idx + n], jin->data[idx + n], tmp4,
                         MPFR_RNDN);
            }
            }


            mpfr_clear(tmp4);
            mpfr_free_cache();
        }
    }
    /* --------------------------------------------------------------------- */






    /* Compute the T-terms (Eq. 50 of Bucha et al., 2019) */
    /* --------------------------------------------------------------------- */
    ti = CHARM(mpfr_ndarray_malloc)(NBITS, 1, dmax + 1);
    if (ti == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        goto EXIT;
    }


    signi = -1;
    mpfr_div(tmp1, rref, r, MPFR_RNDN);
    for (unsigned h = 0; h <= dmax; h++)
    {
        signi = -signi;


        mpfr_mul(ti->data[h], fact->data[h], tmp1, MPFR_RNDN);
        mpfr_mul_si(ti->data[h], ti->data[h], signi, MPFR_RNDN);


        mpfr_div(tmp1, tmp1, r, MPFR_RNDN);
    }
    /* --------------------------------------------------------------------- */






    /* Compute partial Bell Polynomials (Eq. 51 of Bucha et al., 2019) */
    /* --------------------------------------------------------------------- */
    bell = CHARM(mpfr_ndarray_malloc)(NBITS, 2, dmax + 1, dmax + 1);
    if (bell == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        goto EXIT;
    }
    for (unsigned h = 0; h < bell->size; h++)
        mpfr_set_ui(bell->data[h], 0, MPFR_RNDN);


    mpfr_set_ui(bell->data[0], 1, MPFR_RNDN);
    for (unsigned h = 1; h <= dmax; h++)
    {
        idx = h * (dmax + 1);


        for (unsigned j = 1; j <= h; j++)
        {
            for (unsigned l = 0; l <= h - j; l++)
            {
                mpfr_mul(tmp1, ti->data[l + 1],
                         bell->data[(h - 1 - l) * (dmax + 1) + j - 1],
                         MPFR_RNDN);
                mpfr_mul(tmp1, tmp1,
                         binomial->data[(h - 1) * binomial->shape[1] + l],
                         MPFR_RNDN);
                mpfr_add(bell->data[idx + j], bell->data[idx + j], tmp1,
                         MPFR_RNDN);
            }
        }
    }
    /* --------------------------------------------------------------------- */






    /* Compute the G-terms (Eq. 48 of Bucha et al, 2019) */
    /* --------------------------------------------------------------------- */
    gin = CHARM(mpfr_ndarray_malloc)(NBITS, 2, dmax + 1, nmax + 1);
    if (gin == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        goto EXIT;
    }
    for (unsigned h = 0; h < gin->size; h++)
        mpfr_set_ui(gin->data[h], 0, MPFR_RNDN);


    /* h == 0 */
    {
    unsigned long n;
#if HAVE_OPENMP
#pragma omp parallel for default(shared) private(n)
#endif
    for (n = 0; n <= nmax; n++)
        mpfr_mul(gin->data[n], t, ain->data[n], MPFR_RNDN);
    }


    /* h > 0 */
    for (unsigned h = 1; h <= dmax; h++)
    {
        for (unsigned j = 1; j <= h; j++)
        {
            idx = j * (nmax + 1);


#if HAVE_OPENMP
#pragma omp parallel default(shared)
#endif
            {
                mpfr_t tmp4;
                mpfr_init2(tmp4, NBITS);


                {
                unsigned long n;
#if HAVE_OPENMP
#pragma omp for
#endif
                for (n = 0; n <= nmax; n++)
                {
                    mpfr_mul(tmp4, jin->data[idx + n],
                             bell->data[h * (dmax + 1) + j], MPFR_RNDN);
                    mpfr_add(gin->data[h * (nmax + 1) + n],
                             gin->data[h * (nmax + 1) + n], tmp4, MPFR_RNDN);
                }
                }


                mpfr_clear(tmp4);
                mpfr_free_cache();
            }
        }
    }
    /* --------------------------------------------------------------------- */






    /* Compute the trunction coefficients (Eq. 74 of Bucha et al., 2019) */
    /* --------------------------------------------------------------------- */
    /* First, set all truncation coefficients to zero (other parts of the code
     * rely on this) */
    size_t qkpin_size = CHARM(gfm_cap_nq)(nmax, pmax, kmin, kmax, imax, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return;
    }


    for (unsigned v = 0; v < qkpin_size; v++)
        mpfr_set_ui(qkpin[v], 0, MPFR_RNDN);


    rpows = CHARM(mpfr_ndarray_malloc)(NBITS, 1, pmax + 1);
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


    /* Compute the R-terms (Eq. 75 of Bucha et al., 2019) */
    /* ..................................................................... */
    rps = CHARM(mpfr_ndarray_malloc)(NBITS, 2, kmax + 1, pmax);
    if (rps == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        goto EXIT;
    }


    CHARM(gfm_cap_q_dr)(rps, rpows, pmax, kmax, NBITS, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        goto EXIT;
    }
    /* ..................................................................... */


    size_t idx_qkpi0, idx_gin, idx_gip1n, idx_ginps;


    for (unsigned k = kmin; k <= kmax; k++)
    {
        idx_gin   = k * (nmax + 1);
        idx_gip1n = (k + 1) * (nmax + 1);


        for (unsigned p = 1; p <= pmax; p++)
        {
            idx_qkpi0 = CHARM(misc_idx_4d)(k - kmin, p - 1, 0, 0, pmax,
                                           imax + 1, nmax + 1);
            mpfr_set(factp, fact->data[p], MPFR_RNDN);


            if (p == 1)
            {
#if HAVE_OPENMP
#pragma omp parallel default(shared)
#endif
                {
                    {
                    unsigned long n;
#if HAVE_OPENMP
#pragma omp for
#endif
                    for (n = 0; n <= nmax; n++)
                        mpfr_set(qkpin[idx_qkpi0 + n], gin->data[idx_gin + n],
                                 MPFR_RNDN);
                    }


                    mpfr_free_cache();
                }
            }
            else if (p == 2)
            {
#if HAVE_OPENMP
#pragma omp parallel default(shared)
#endif
                {
                    mpfr_t tmp4;
                    mpfr_init2(tmp4, NBITS);


                    {
                    unsigned long n;
#if HAVE_OPENMP
#pragma omp for
#endif
                    for (n = 0; n <= nmax; n++)
                    {
                        mpfr_mul_si(tmp4, gin->data[idx_gin + n],
                                    -((long)k - 1), MPFR_RNDN);
                        mpfr_mul(qkpin[idx_qkpi0 + n],
                                 gin->data[idx_gip1n + n], r, MPFR_RNDN);
                        mpfr_sub(qkpin[idx_qkpi0 + n], tmp4,
                                 qkpin[idx_qkpi0 + n], MPFR_RNDN);
                        mpfr_div(qkpin[idx_qkpi0 + n], qkpin[idx_qkpi0 + n],
                                 factp, MPFR_RNDN);
                    }
                    }


                    mpfr_clear(tmp4);
                    mpfr_free_cache();
                }
            }
            else
            {
                for (unsigned s = 1; s <= p - 2; s++)
                {
                    CHARM(gfm_cap_q_aps)(p, s, fact, aps, err);
                    if (!CHARM(err_isempty)(err))
                    {
                        CHARM(err_propagate)(err, __FILE__, __LINE__,
                                             __func__);
                        goto EXIT;
                    }


                    for (unsigned v = 0; v <= k; v++)
                    {
                        if ((k - v) > (p - s))
                            continue;


                        mpfr_mul(tmp1, aps,
                                 binomial->data[k * binomial->shape[1] + v],
                                 MPFR_RNDN);
                        mpfr_div(tmp1, tmp1, factp, MPFR_RNDN);


                        idx_ginps = (p - s + v) * (nmax + 1);
                        mpfr_mul(tmp2, tmp1,
                                 rps->data[(k - v) * pmax + p - 1 - s],
                                 MPFR_RNDN);


#if HAVE_OPENMP
#pragma omp parallel default(shared)
#endif
                        {
                            mpfr_t tmp4;
                            mpfr_init2(tmp4, NBITS);


                            {
                            unsigned long n;
#if HAVE_OPENMP
#pragma omp for
#endif
                            for (n = 0; n <= nmax; n++)
                            {
                                mpfr_mul(tmp4, tmp2, gin->data[idx_ginps + n],
                                         MPFR_RNDN);
                                mpfr_add(qkpin[idx_qkpi0 + n],
                                         qkpin[idx_qkpi0 + n], tmp4,
                                         MPFR_RNDN);
                            }
                            }


                            mpfr_clear(tmp4);
                            mpfr_free_cache();
                        }
                    }
                }
            }
        }
    }
    /* --------------------------------------------------------------------- */






    /* Compute the truncation coefficients for "i > 0" */
    /* --------------------------------------------------------------------- */
    size_t idx_qks00;
    unsigned smin;


    for (unsigned k = kmin; k <= kmax; k++)
    {
        for (unsigned p = 1; p <= pmax; p++)
        {
            for (unsigned i = 1; i <= imax; i++)
            {
                idx_qkpi0 = CHARM(misc_idx_4d)(k - kmin, p - 1, i, 0, pmax,
                                               imax + 1, nmax + 1);


                /* "CHARM_MAX(1, p - i)", but in a special way, as "p - i" may
                 * be negative and "p" and "i" are unsigned ints. */
                smin = (p < i) ? 1: CHARM_MAX(1, p - i);
                for (unsigned s = smin; s <= p; s++)
                {
                    mpfr_set(tmp1,
                             binomial->data[i * binomial->shape[1] + (p - s)],
                             MPFR_RNDN);
                    mpfr_mul_ui(tmp1, tmp1, s, MPFR_RNDN);
                    mpfr_div_ui(tmp1, tmp1, p, MPFR_RNDN);


                    idx_qks00 = CHARM(misc_idx_4d)(k - kmin, s - 1, 0, 0, pmax,
                                                   imax + 1, nmax + 1);
#if HAVE_OPENMP
#pragma omp parallel default(shared)
#endif
                    {
                        mpfr_t tmp4;
                        mpfr_init2(tmp4, NBITS);
                        {
                        unsigned long n;
#if HAVE_OPENMP
#pragma omp for
#endif
                        for (n = 0; n <= nmax; n++)
                        {
                            mpfr_mul(tmp4, tmp1, qkpin[idx_qks00 + n],
                                     MPFR_RNDN);
                            mpfr_add(qkpin[idx_qkpi0 + n],
                                     qkpin[idx_qkpi0 + n], tmp4, MPFR_RNDN);
                        }
                        }


                        mpfr_clear(tmp4);
                        mpfr_free_cache();
                    }
                }
            }
        }
    }
    /* --------------------------------------------------------------------- */






EXIT:
    /* Clean up */
    /* --------------------------------------------------------------------- */
    CHARM(mpfr_ndarray_free)(fact);
    CHARM(mpfr_ndarray_free)(double_fact);
    CHARM(mpfr_ndarray_free)(gammau);
    CHARM(mpfr_ndarray_free)(gammal);
    CHARM(mpfr_ndarray_free)(m);
    CHARM(mpfr_ndarray_free)(binomial);
    CHARM(mpfr_ndarray_free)(bin);
    CHARM(mpfr_ndarray_free)(beta0);
    CHARM(mpfr_ndarray_free)(beta1uu);
    CHARM(mpfr_ndarray_free)(beta1ul);
    CHARM(mpfr_ndarray_free)(ain);
    CHARM(mpfr_ndarray_free)(alpha0uu);
    CHARM(mpfr_ndarray_free)(alpha0ul);
    CHARM(mpfr_ndarray_free)(alpha1uu);
    CHARM(mpfr_ndarray_free)(alpha1ul);
    CHARM(mpfr_ndarray_free)(jin);
    CHARM(mpfr_ndarray_free)(ti);
    CHARM(mpfr_ndarray_free)(bell);
    CHARM(mpfr_ndarray_free)(gin);
    CHARM(mpfr_ndarray_free)(rpows);
    CHARM(mpfr_ndarray_free)(rps);


    mpfr_clears(t, t2, cospsi, ul, uu, u0, gu, gl, tmp1, tmp2, tmp3, ttmp,
                tmi1, tmi2, tp2, gammat, it, iim1t, tpi1, tpi2, t2p1, tt2,
                aps, factp, (mpfr_ptr)NULL);


    mpfr_free_cache();
    FLUSH_UNRELEASED_MEMORY;
    /* --------------------------------------------------------------------- */


    return;
}

