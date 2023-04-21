/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../src/prec.h"
#include "cmp_vals.h"
#include "parameters.h"
/* ------------------------------------------------------------------------- */






/* Checks "integ_yi1n1m1yi2n2m2" using the orthogonality of spherical
 * harmonics. */
long int check_integ_yi1n1m1yi2n2m2(void)
{
    /* --------------------------------------------------------------------- */
    CHARM(err) *err = CHARM(err_init)();
    if (err == NULL)
    {
        printf("Failed to initialize the \"err\" structure.\n");
        exit(CHARM_FAILURE);
    }


    CHARM(pnmj) *pnmj;


    long int e = 0;
    REAL iy, iy_ref;
    REAL delta_i1i2, delta_m1m2, delta_n1n2;
    /* --------------------------------------------------------------------- */


    /* Loop over the two ordering schemes of the Fourier coefficients */
    /* --------------------------------------------------------------------- */
    unsigned long nmax = 8;


    for (int o = 0; o < 2; o++)
    {
        /* Compute the Fourier coefficients */
        /* ................................................................. */
        pnmj = CHARM(leg_pnmj_calloc)(nmax,
                                      (o == 0) ? CHARM_LEG_PNMJ_ORDER_MNJ :
                                                 CHARM_LEG_PNMJ_ORDER_MJN);
        if (pnmj == NULL)
        {
            fprintf(stderr, "Failed to initialize the \"pnmj\" "
                            "structure.\n");
            exit(CHARM_FAILURE);
        }


        CHARM(leg_pnmj_coeffs)(pnmj, nmax, err);
        CHARM(err_handler)(err, 1);
        /* ................................................................. */


        for (int i1 = 0; i1 < 2; i1++)
        {
        for (int i2 = 0; i2 < 2; i2++)
        {
            delta_i1i2 = (i1 == i2) ? PREC(1.0) : PREC(0.0);


            for (unsigned long m1 = 0; m1 <= nmax; m1++)
            {
            for (unsigned long m2 = 0; m2 <= nmax; m2++)
            {
                if (m1 == m2)
                {
                    if ((i1 * i2 == 1) && (m1 == 0))
                        /* The sine harmonics of order zero do not exist, so
                         * make sure the reference value will be zero. */
                        delta_m1m2 = PREC(0.0);
                    else
                        delta_m1m2 = PREC(1.0);
                }
                else
                    delta_m1m2 = PREC(0.0);


                for (unsigned long n1 = m1; n1 <= nmax; n1++)
                {
                for (unsigned long n2 = m2; n2 <= nmax; n2++)
                {
                    delta_n1n2 = (n1 == n2) ? PREC(1.0) : PREC(0.0);


                    /* Compute the reference value using the orthogonality
                     * property of Legendre functions */
                    /* ................................................. */
                    iy_ref = PREC(4.0) * PI * delta_i1i2 * delta_n1n2 *
                             delta_m1m2;
                    /* ................................................. */


                    /* Integral of Legendre functions over the interval of
                     * co-latitudes "[0, pi]" */
                    /* ................................................. */
                    iy = CHARM(integ_yi1n1m1yi2n2m2)(PREC(0.0), PI,
                                                     PREC(0.0), PREC(2.0) *
                                                                PI,
                                                     i1, n1, m1,
                                                     i2, n2, m2,
                                                     pnmj, err);
                    CHARM(err_handler)(err, 1);


                    e += cmp_vals(iy, iy_ref,
                                  PREC(10.0) * CHARM(glob_threshold));
                    /* ................................................. */


                    /* Split the unit sphere into four spherical rectangles,
                     * compute the integral over each rectangle separately,
                     * then sum the four integrals and compare the result with
                     * "iy_ref". */
                    /* ................................................. */
                    iy = CHARM(integ_yi1n1m1yi2n2m2)(PREC(0.0), CLT0,
                                                     PREC(0.0), LON0,
                                                     i1, n1, m1,
                                                     i2, n2, m2,
                                                     pnmj, err);
                    CHARM(err_handler)(err, 1);


                    iy += CHARM(integ_yi1n1m1yi2n2m2)(PREC(0.0), CLT0,
                                                      LON0, PREC(2.0) * PI,
                                                      i1, n1, m1,
                                                      i2, n2, m2,
                                                      pnmj, err);
                    CHARM(err_handler)(err, 1);


                    iy += CHARM(integ_yi1n1m1yi2n2m2)(CLT0, PI,
                                                      PREC(0.0), LON0,
                                                      i1, n1, m1,
                                                      i2, n2, m2,
                                                      pnmj, err);
                    CHARM(err_handler)(err, 1);


                    iy += CHARM(integ_yi1n1m1yi2n2m2)(CLT0, PI,
                                                      LON0, PREC(2.0) * PI,
                                                      i1, n1, m1,
                                                      i2, n2, m2,
                                                      pnmj, err);
                    CHARM(err_handler)(err, 1);


                    e += cmp_vals(iy, iy_ref,
                                  PREC(10.0) * CHARM(glob_threshold));
                    /* ..................................................... */
                }
                }
            }
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
