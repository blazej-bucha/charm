/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../src/prec.h"
#include "parameters.h"
#include "error_messages.h"
#include "../src/mpfr/mpfr_flush_unreleased_memory.h"
#include "check_gfm_cap_q_check_prec.h"
/* ------------------------------------------------------------------------- */





long int check_gfm_cap_q_check_prec(void)
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


    /* Input parameters to "CHARM(gfm_cap_q)" */
    /* --------------------------------------------------------------------- */
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
    /* For the current value of "GFM_CAP_NBITS", "CHARM(gfm_cap_q_check_prec)"
     * should never return value smaller than "ref". */
    const long ref = 50L;


    for (unsigned long n = 0; n <= GFM_NMAX_POTENTIAL; n++)
    {
    for (unsigned p = 1; p <= GFM_MOON_PMAX; p++)
    {
    for (unsigned k = GFM_CAP_KMIN; k <= GFM_CAP_KMAX; k++)
    {
    for (unsigned i = 0; i <= GFM_Q_IMAX; i++)
    {
    for (size_t t = 0; t < NQTYPE; t++)
    {
        long digits = CHARM(gfm_cap_q_check_prec)(rref,
                                                  r,
                                                  psi,
                                                  n,
                                                  p,
                                                  0,
                                                  k,
                                                  i,
                                                  q_type[t],
                                                  GFM_CAP_NBITS,
                                                  2 * GFM_CAP_NBITS,
                                                  err);
        CHARM(err_handler)(err, 1);
        if (digits < ref)
        {
            printf("\n");
            printf("            result:    %ld\n", digits);
            printf("            reference: %ld\n", ref);
            e += 1;
        }
    }
    }
    }
    }
    }
    /* --------------------------------------------------------------------- */


    /* --------------------------------------------------------------------- */
    mpfr_clears(rref, r, psi, (mpfr_ptr)NULL);
    mpfr_free_cache();
    FLUSH_UNRELEASED_MEMORY;
    CHARM(err_free)(err);
    /* --------------------------------------------------------------------- */


    return e;
}
