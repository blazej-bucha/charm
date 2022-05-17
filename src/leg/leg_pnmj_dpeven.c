/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdlib.h>
#ifdef _MSC_VER
#   define _USE_MATH_DEFINES
#endif
#include <math.h>
#include "../prec.h"
#include "../err/err_set.h"
/* ------------------------------------------------------------------------- */






/**
 * Computes the Fourier coefficients \f$P_{nnj}\f$ and \f$P_{n, n - 1, j}\f$
 * when degree ``n`` is even and ``n >= 6``.
 *
 * This function is due to Fukushima (2018), Table 4.
 *
 *
 * @param[in] n Harmonic degree
 * @param[in] xpold Significand of the X-number representation of \f$P_{n-2,
 * n-2, j}\f$
 * @param[in] ipold Exponent of the X-number representation of \f$P_{n-2, n-2,
 * j}\f$
 *
 * @param[out] xp Significand of the X-number representation of \f$P_{n, n,
 * j}\f$
 * @param[out] xp1 Significand of the X-number representation of \f$P_{n, n -
 * 1, j}\f$
 * @param[out] ip Exponent of the X-number representation of \f$P_{n, n, j}\f$
 * @param[out] ip1 Exponent of the X-number representation of \f$P_{n, n - 1,
 * j}\f$
 * @param[out] err Error reported by the function (if any).
 *
 *
 * **References**: Fukushima, T. (2018) Fast computation of sine/cosine series
 * coefficients of associated Legendre functions of arbitrary high degree and
 * order.
 *
 * */
void CHARM(leg_pnmj_dpeven)(unsigned long n, const REAL *xpold, REAL *xp,
                            REAL *xp1, const int *ipold, int *ip, int *ip1,
                            CHARM(err) *err)
{
    if (n % 2)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"n\" has to be even.");


        return;
    }


    if (n < 6)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"n\" cannot be smaller than \"6\".");


        return;
    }


    /* Note that "n >= 6", so unsigned longs are safe */
    unsigned long jx   = n / 2;
    unsigned long jxm2 = jx - 2;
    unsigned long jxm1 = jx - 1;
    unsigned long n2   = 2 * n;


    /* Eq. (19) of Fukushima (2018) */
    REAL gamma = ADDP(0.125) * SQRT((REAL)(n2 + 1) * (REAL)(n2 - 1) /
                                   ((REAL)n * (REAL)(n - 1)));

    REAL gamma2 = ADDP(2.0) * gamma;


    /* Eq. (64) */
    CHARM(xnum_xlsum2)(gamma2, xpold[0], -gamma, xpold[1], &xp[0], ipold[0],
                       ipold[1], &ip[0]);


    /* Eq. (65) */
    REAL xtemp;
    int itemp;
    CHARM(xnum_xlsum2)(-gamma2, xpold[0], gamma2, xpold[1], &xtemp,
                       ipold[0], ipold[1], &itemp);
    CHARM(xnum_xlsum2)(ADDP(1.0), xtemp, -gamma, xpold[2], &xp[1], itemp,
                       ipold[2], &ip[1]);


    /* Eq. (66) */
    unsigned long jm1, jp1;
    for (unsigned long j = 2; j <= jxm2; j++)
    {
        jm1 = j - 1;
        jp1 = j + 1;

        CHARM(xnum_xlsum2)(-gamma, xpold[jm1], gamma2, xpold[j], &xtemp,
                           ipold[jm1], ipold[j], &itemp);
        CHARM(xnum_xlsum2)(ADDP(1.0), xtemp, -gamma, xpold[jp1], &xp[j],
                           itemp, ipold[jp1], &ip[j]);
    }


    /* Eq. (67) */
    CHARM(xnum_xlsum2)(-gamma, xpold[jxm2], gamma2, xpold[jxm1], &xp[jxm1],
                       ipold[jxm2], ipold[jxm1], &ip[jxm1]);


    /* Eq. (68) */
    xp[jx] = -gamma * xpold[jxm1];
    ip[jx] = ipold[jxm1];
    CHARM(xnum_xnorm)(&xp[jx], &ip[jx]);


    REAL alpha2 = ADDP(2.0) * SQRT(ADDP(2.0) / (REAL)n);
    xp1[0] = ADDP(0.0);
    ip1[0] = 0;


    /* Eq. (69) */
    for (unsigned long j = 1; j <= jx; j++)
    {
        xp1[j] = -(REAL)j * alpha2 * xp[j];
        ip1[j] = ip[j];

        CHARM(xnum_xnorm)(&xp1[j], &ip1[j]);
    }


    return;
}
