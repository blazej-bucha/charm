/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdlib.h>
#include <math.h>
#include "../prec.h"
#include "../xnum/xnum_xlsum2.h"
#include "../xnum/xnum_xnorm.h"
#include "../err/err_set.h"
/* ------------------------------------------------------------------------- */






/**
 * Returns the Fourier coefficients \f$P_{nnj}\f$ and \f$P_{n, n - 1, j}\f$
 * when degree ``n`` is odd and ``n >= 5``.
 *
 * This function is due to Fukushima (2018), Table 7.
 *
 *
 * @param[in] n Harmonic degree.
 *
 * @param[in] xpold Significand of the X-number representation of \f$P_{n-2,
 * n-2, j}\f$.
 *
 * @param[in] ipold Exponent of the X-number representation of \f$P_{n-2, n-2,
 * j}\f$.
 *
 * @param[out] xp Significand of the X-number representation of \f$P_{n, n,
 * j}\f$.
 *
 * @param[out] xp1 Significand of the X-number representation of \f$P_{n, n -
 * 1, j}\f$.
 *
 * @param[out] ip Exponent of the X-number representation of \f$P_{n, n, j}\f$.
 *
 * @param[out] ip1 Exponent of the X-number representation of \f$P_{n, n - 1,
 * j}\f$
 *
 * @param[out] err Error reported by the function (if any).
 *
 * **References**: Fukushima, T. (2018) Fast computation of sine/cosine series
 * coefficients of associated Legendre functions of arbitrary high degree and
 * order.
 *
 * */
void CHARM(leg_pnmj_dpodd)(unsigned long n, const REAL *xpold, REAL *xp,
                           REAL *xp1, const int *ipold, int *ip, int *ip1,
                           CHARM(err) *err)
{
    if ((n % 2) == 0)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"n\" has to be odd.");


        return;
    }


    if (n < 5)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"n\" cannot be smaller than \"6\".");


        return;
    }


    /* Note that "n >= 5", so unsigned longs are safe */
    unsigned long jx = (n - 1) / 2;
    unsigned long jxm2 = jx - 2;
    unsigned long jxm1 = jx - 1;
    unsigned long n2 = 2 * n;


    /* Eq. (19) of Fukushima (2018) */
    REAL gamma = PREC(0.125) * SQRT((REAL)(n2 + 1) * (REAL)(n2 - 1)
                                    / ((REAL)n * (REAL)(n - 1)));

    REAL gamma2 = PREC(2.0) * gamma;


    /* Eq. (70) */
    CHARM(xnum_xlsum2)(PREC(3.0) * gamma, xpold[0], -gamma, xpold[1], &xp[0],
                       ipold[0], ipold[1], &ip[0]);


    /* Eq. (71) */
    unsigned long jm1, jp1;
    REAL xtemp;
    int itemp;
    for (unsigned long j = 1; j <= jxm2; j++)
    {
        jm1 = j - 1;
        jp1 = j + 1;


        CHARM(xnum_xlsum2)(-gamma, xpold[jm1], gamma2, xpold[j], &xtemp,
                           ipold[jm1], ipold[j], &itemp);
        CHARM(xnum_xlsum2)(PREC(1.0), xtemp, -gamma, xpold[jp1], &xp[j],
                           itemp, ipold[jp1], &ip[j]);
    }


    /* Eq. (72) */
    CHARM(xnum_xlsum2)(-gamma, xpold[jxm2], gamma2, xpold[jxm1], &xp[jxm1],
                       ipold[jxm2], ipold[jxm1], &ip[jxm1]);


    /* Eq. (73) */
    xp[jx] = -gamma * xpold[jxm1];
    ip[jx] = ipold[jxm1];

    CHARM(xnum_xnorm)(&xp[jx], &ip[jx]);

    REAL alpha = SQRT(PREC(2.0) / (REAL)n);


    /* Eq. (74) */
    for (unsigned long j = 0; j <= jx; j++)
    {
        xp1[j] = (REAL)(2 * j + 1) * alpha * xp[j];
        ip1[j] = ip[j];

        CHARM(xnum_xnorm)(&xp1[j], &ip1[j]);
    }


    return;
}
