/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdlib.h>
#include <math.h>
#include "../prec.h"
#include "../xnum/xnum_xnorm.h"
#include "../xnum/xnum_xlsum2.h"
#include "../err/err_set.h"
#include "leg_pnmj_gpeven.h"
/* ------------------------------------------------------------------------- */






/**
 * Computes the Fourier coefficients \f$P_{nmj}\f$ when degree ``n`` is even
 * and ``m <= (n - 2)``.
 *
 * This function is due to Fukushima (2018), Table 6.
 *
 *
 * @param[in] jmax Maximum value of ``j``
 * @param[in] n Harmonic degree
 * @param[in] m Harmonic order
 * @param[in] xp1 Significand of the X-number representation of \f$P_{n, m + 1,
 * j}\f$
 * @param[in] ip1 Exponent of the X-number representation of \f$P_{n, m + 1,
 * j}\f$
 * @param[in] xp2 Significand of the X-number representation of \f$P_{n, m + 2,
 * j}\f$
 * @param[in] ip2 Exponent of the X-number representation of \f$P_{n, m + 2,
 * j}\f$
 *
 *
 * @param[out] xp0 Significand of the X-number representation of \f$P_{nmj}\f$
 * @param[out] ip0 Exponent of the X-number representation of \f$P_{nmj}\f$
 * @param[out] err Error reported by the function (if any).
 *
 *
 * **References**: Fukushima, T. (2018) Fast computation of sine/cosine series
 * coefficients of associated Legendre functions of arbitrary high degree and
 * order.
 *
 * */
void CHARM(leg_pnmj_gpeven)(unsigned long jmax, unsigned long n,
                            unsigned long m, const REAL *xp2,
                            const REAL *xp1, REAL *xp0, const int *ip2,
                            const int *ip1, int *ip0,
                            CHARM(err) *err)
{
    if (n % 2)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"n\" has to be even.");


        return;
    }


    if ((n - 2) < m)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"n - 2\" cannot be smaller than \"m\".");


        return;
    }


    if (m > n)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Harmonic order \"m\" cannot be larger "
                       "than harmonic degree \"n\".");


        return;
    }


    unsigned long m1 = m + 1;
    unsigned long m2 = m + 2;
    unsigned long modd = m - (m / 2) * 2;


    REAL u;
    if (m == 0)
        u = SQRT(PREC(0.5) / ((REAL)n * (REAL)(n + 1)));
    else
        u = SQRT(PREC(1.0) / ((REAL)(n - m) * (REAL)(n + m1)));


    REAL alpha2 = PREC(4.0) * u;
    REAL beta = SQRT((REAL)(n - m1) * (REAL)(n + m2)) * u;


    /* (75) of Fukushima (2018) */
    xp0[0] = beta * xp2[0];
    ip0[0] = ip2[0];
    CHARM(xnum_xnorm)(&xp0[0], &ip0[0]);


    if (modd == 0)
        /* Eq. (76) -- Note that there appears to be a misprint for the
         * condition on "j" right before this equation. Correctly, it should
         * read (in C syntax) "j != 0". */
        for (unsigned long j = 1; j <= jmax; j++)
            CHARM(xnum_xlsum2)((REAL)(j) * alpha2, xp1[j], beta, xp2[j],
                               &xp0[j], ip1[j], ip2[j], &ip0[j]);
    else
        /* Eq. (77) -- Note that there appears to be a misprint for the
         * condition on "j" right before this equation. Correctly, it should
         * read (in C syntax) "j != 0". */
        for (unsigned long j = 1; j <= jmax; j++)
            CHARM(xnum_xlsum2)(-(REAL)(j) * alpha2, xp1[j], beta, xp2[j],
                               &xp0[j], ip1[j], ip2[j], &ip0[j]);


    return;
}
