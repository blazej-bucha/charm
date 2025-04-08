/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../src/prec.h"
#include "cmp_vals.h"
#include "parameters.h"
#include "error_messages.h"
#include "check_integ_pn1m1pn2m2.h"
/* ------------------------------------------------------------------------- */






/* Checks "integ_pn1m1pn2m2" using the orthogonality of Legendre functions.
 *
 * Because of the orthogonality test, the checks are done only for Legendre
 * functions of equal orders.  This test therefore does not cover all possible
 * variations of "n1", "m1", "n2" and "m2". */
long int check_integ_pn1m1pn2m2(void)
{
    /* --------------------------------------------------------------------- */
    CHARM(err) *err = CHARM(err_init)();
    if (err == NULL)
    {
        fprintf(stderr, "%s", ERR_MSG_ERR);
        exit(CHARM_FAILURE);
    }


    CHARM(pnmj) *pnmj;


    long int e = 0;
    REAL ip, ip_ref;
    REAL delta_0m, delta_n1n2;
    /* --------------------------------------------------------------------- */


    /* Loop over the two ordering schemes of the Fourier coefficients */
    /* --------------------------------------------------------------------- */
    unsigned long nmax = 8;


    for (int o = 0; o < 2; o++)
    {
        /* Compute the Fourier coefficients */
        /* ................................................................. */
        pnmj = CHARM(leg_pnmj_calloc)(nmax,
                                      (o == 0) ? CHARM_LEG_PMNJ :
                                                 CHARM_LEG_PMJN);
        if (pnmj == NULL)
        {
            fprintf(stderr, "%s", ERR_MSG_PNMJ);
            exit(CHARM_FAILURE);
        }


        CHARM(leg_pnmj_coeffs)(pnmj, nmax, err);
        CHARM(err_handler)(err, 1);
        /* ................................................................. */


        for (unsigned long m = 0; m <= nmax; m++)
        {
            delta_0m = (m == 0) ? PREC(1.0) : PREC(0.0);


            for (unsigned long n1 = m; n1 <= nmax; n1++)
            {
                for (unsigned long n2 = m; n2 <= nmax; n2++)
                {
                    delta_n1n2 = (n1 == n2) ? PREC(1.0) : PREC(0.0);


                    /* Compute the reference value using the orthogonality
                     * property of Legendre functions */
                    /* ..................................................... */
                    ip_ref = PREC(2.0) * (PREC(2.0) - delta_0m) * delta_n1n2;
                    /* ..................................................... */


                    /* Integral of Legendre functions over the interval of
                     * co-latitudes "[0, pi]" */
                    /* ..................................................... */
                    ip = CHARM(integ_pn1m1pn2m2)(PREC(0.0), PI, n1, m, n2, m,
                                                 pnmj, err);
                    CHARM(err_handler)(err, 1);


                    e += cmp_vals_real(ip, ip_ref,
                                       PREC(10.0) * CHARM(glob_threshold));
                    /* ..................................................... */


                    /* Integrals of Legendre functions over the intervals of
                     * co-latitudes "[0, CLT0]" and "[CLT0, PI]", the sum of
                     * which must be equal to the integral over "[0, PI]",
                     * hence to the reference values. */
                    /* ..................................................... */
                    ip = CHARM(integ_pn1m1pn2m2)(PREC(0.0), CLT0, n1, m, n2, m,
                                                 pnmj, err);
                    CHARM(err_handler)(err, 1);


                    ip += CHARM(integ_pn1m1pn2m2)(CLT0, PI, n1, m, n2, m, pnmj,
                                                  err);
                    CHARM(err_handler)(err, 1);


                    e += cmp_vals_real(ip, ip_ref,
                                       PREC(10.0) * CHARM(glob_threshold));
                    /* ..................................................... */
                }
            }
        }


        CHARM(leg_pnmj_free)(pnmj);
    }
    /* --------------------------------------------------------------------- */


    /* --------------------------------------------------------------------- */
    CHARM(err_free)(err);


    return e;
    /* --------------------------------------------------------------------- */
}
