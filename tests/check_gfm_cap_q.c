/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../src/prec.h"
#include "parameters.h"
#include "error_messages.h"
#include "../src/mpfr/mpfr_flush_unreleased_memory.h"
#ifdef GENREF
#   include "mpfr_array2file.h"
#else
#   include "mpfr_validate.h"
#endif
#include "check_gfm_cap_q.h"
/* ------------------------------------------------------------------------- */





long int check_gfm_cap_q(void)
{
    long int e = 0;


    /* Create CHarm's error structure */
    /* --------------------------------------------------------------------- */
    CHARM(err) *err = CHARM(err_init)();
    if (err == NULL)
    {
        fprintf(stderr, ERR_MSG_ERR);
        exit(CHARM_FAILURE);
    }
    /* --------------------------------------------------------------------- */


    /* Allocate memory for trunction coefficients */
    /* --------------------------------------------------------------------- */
    /* Get the number of trunction coefficients */
    size_t q_size = CHARM(gfm_cap_nq)(GFM_NMAX_POTENTIAL,
                                      GFM_MOON_PMAX,
                                      GFM_CAP_KMIN,
                                      GFM_CAP_KMAX,
                                      GFM_Q_IMAX,
                                      err);
    CHARM(err_handler)(err, 1);


    mpfr_t *q = (mpfr_t *)malloc(q_size * sizeof(mpfr_t));
    if (q == NULL)
    {
        fprintf(stderr, CHARM_ERR_MALLOC_FAILURE);
        exit(CHARM_FAILURE);
    }
    for (size_t j = 0; j < q_size; j++)
        mpfr_init2(q[j], GFM_CAP_NBITS);
    /* --------------------------------------------------------------------- */


    /* Input parameters to "CHARM(gfm_cap_q)" */
    /* --------------------------------------------------------------------- */
    /* Integration zones to be tested */
#define NZONE 2
    int zone[NZONE] = {CHARM_GFM_NEAR_ZONE,
                       CHARM_GFM_FAR_ZONE};


    /* Types of truncation coefficients to be tested */
#define NQTYPE 6
    unsigned q_type[NQTYPE] = {CHARM_GFM_Q00,
                               CHARM_GFM_Q10,
                               CHARM_GFM_Q11,
                               CHARM_GFM_Q20,
                               CHARM_GFM_Q21,
                               CHARM_GFM_Q22};


    mpfr_t rref, r, psi;
    mpfr_inits2(GFM_CAP_NBITS, rref, r, psi, (mpfr_ptr)NULL);


    mpfr_set_d(rref, MOON_RREF, MPFR_RNDN);
    mpfr_set_d(r, MOON_RREF + GFM_CAP_HEIGHT, MPFR_RNDN);
    mpfr_set_d(psi, GFM_CAP_PSI0, MPFR_RNDN);
    /* --------------------------------------------------------------------- */


    /* --------------------------------------------------------------------- */
    for (size_t z = 0; z < NZONE; z++)
    {
    for (size_t t = 0; t < NQTYPE; t++)
    {
        CHARM(gfm_cap_q)(rref,
                         r,
                         psi,
                         GFM_NMAX_POTENTIAL,
                         GFM_MOON_PMAX,
                         GFM_CAP_KMIN,
                         GFM_CAP_KMAX,
                         GFM_Q_IMAX,
                         zone[z],
                         q_type[t],
                         GFM_CAP_NBITS,
                         q,
                         err);
        CHARM(err_handler)(err, 1);


        char file[NSTR_LONG];
        sprintf(file, "%s/gfm_q_t%u_z%d%s", FOLDER, q_type[t], zone[z], FTYPE);
#if GENREF
        e += mpfr_array2file(file, q, q_size);
#else
        mpfr_t eps;
        mpfr_init2(eps, GFM_CAP_NBITS);
        mpfr_set_d(eps, GFM_Q_THRESHOLD, MPFR_RNDN);

        e += mpfr_validate(file, q, q_size, eps);

        mpfr_clear(eps);
#endif
    }
    }
    /* --------------------------------------------------------------------- */


    /* --------------------------------------------------------------------- */
    mpfr_clears(rref, r, psi, (mpfr_ptr)NULL);
    for (size_t j = 0; j < q_size; j++)
        mpfr_clear(q[j]);
    free(q);
    mpfr_free_cache();
    FLUSH_UNRELEASED_MEMORY;
    CHARM(err_free)(err);
    /* --------------------------------------------------------------------- */


    return e;
}
