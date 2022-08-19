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
 * Computes Fourier coefficients of a Legendre function of degree ``n`` equal
 * to or less than ``4``.
 *
 * The function is due to Fukushima (2018), Table 3.  In order to allow for an
 * easy transition to quadruple precision, we replaced the explicit values from
 * Table 3 of Fukushima (2018) by the analytical expressions from Table 1 of
 * the same paper.
 *
 * @param[in] n Harmonic degree equal to or less than ``4``.
 * @param[in] m Harmonic order.
 * @param[out] p Fourier coefficients of a Legendre function of degree ``n``
 * and order ``m`` for ``j = 0``, ``1``, ``...``, ``floor(n / 2) + 1``.
 *
 *
 * @param[out] err Error reported by the function (if any).
 *
 *
 * **References**: Fukushima, T. (2018) Fast computation of sine/cosine series
 * coefficients of associated Legendre functions of arbitrary high degree and
 * order.
 *
 * */
void CHARM(leg_pnmj_leq4)(unsigned long n, unsigned long m, REAL *p,
                          CHARM(err) *err)
{
    if (n > 4)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"n\" cannot be larger than \"4\".");


        return;
    }


    if (m > n)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Harmonic order \"m\" cannot be larger "
                       "than harmonic degree \"n\".");


        return;
    }


    if (n == 0)
        p[0] = PREC(1.0);
    else if (n == 1)
        p[0] = SQRT(PREC(3.0));
    else if (n == 2)
    {
        if (m == 0)
        {
            REAL tmp = SQRT(PREC(5.0)) / PREC(4.0);
            p[0] = tmp;
            p[1] = PREC(3.0) * tmp;
        }
        else if (m == 1)
        {
            p[0] = PREC(0.0);
            p[1] = SQRT(PREC(15.0)) / PREC(2.0);
        }
        else if (m == 2)
        {
            REAL tmp = SQRT(PREC(15.0)) / PREC(4.0);
            p[0] =  tmp;
            p[1] = -tmp;
        }
    }
    else if (n == 3)
    {
        if (m == 0)
        {
            REAL tmp = SQRT(PREC(7.0)) / PREC(8.0);
            p[0] = PREC(3.0) * tmp;
            p[1] = PREC(5.0) * tmp;
        }
        else if (m == 1)
        {
            REAL tmp = SQRT(PREC(42.0)) / PREC(16.0);
            p[0] =       tmp;
            p[1] = PREC(5.0) * tmp;
        }
        else if (m == 2)
        {
            REAL tmp = SQRT(PREC(105.0)) / PREC(8.0);
            p[0] =  tmp;
            p[1] = -tmp;
        }
        else if (m == 3)
        {
            REAL tmp = SQRT(PREC(70.0)) / PREC(16.0);
            p[0] = PREC(3.0) * tmp;
            p[1] =      -tmp;
        }
    }
    else if (n == 4)
    {
        if (m == 0)
        {
            p[0] =  PREC(27.0) / PREC(64.0);
            p[1] =  PREC(15.0) / PREC(16.0);
            p[2] = PREC(105.0) / PREC(64.0);
        }
        else if (m == 1)
        {
            REAL tmp = SQRT(PREC(10.0));
            p[0] =  PREC(0.0);
            p[1] =  PREC(3.0) * tmp / PREC(16.0);
            p[2] = PREC(21.0) * tmp / PREC(32.0);
        }
        else if (m == 2)
        {
            REAL tmp = SQRT(PREC(5.0));
            p[0] =   PREC(9.0) * tmp / PREC(32.0);
            p[1] =   PREC(3.0) * tmp / PREC(8.0);
            p[2] = PREC(-21.0) * tmp / PREC(32.0);
        }
        else if (m == 3)
        {
            REAL tmp = PREC(3.0) * SQRT(PREC(70.0));
            p[0] =  PREC(0.0);
            p[1] =  tmp / PREC(16.0);
            p[2] = -tmp / PREC(32.0);
        }
        else if (m == 4)
        {
            REAL tmp = SQRT(PREC(35.0));
            p[0] =  PREC(9.0) * tmp / PREC(64.0);
            p[1] = PREC(-3.0) * tmp / PREC(16.0);
            p[2] =  PREC(3.0) * tmp / PREC(64.0);
        }
    }


    return;
}
