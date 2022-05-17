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
        p[0] = ADDP(1.0);
    else if (n == 1)
        p[0] = SQRT(ADDP(3.0));
    else if (n == 2)
    {
        if (m == 0)
        {
            REAL tmp = SQRT(ADDP(5.0)) / ADDP(4.0);
            p[0] = tmp;
            p[1] = ADDP(3.0) * tmp;
        }
        else if (m == 1)
        {
            p[0] = ADDP(0.0);
            p[1] = SQRT(ADDP(15.0)) / ADDP(2.0);
        }
        else if (m == 2)
        {
            REAL tmp = SQRT(ADDP(15.0)) / ADDP(4.0);
            p[0] =  tmp;
            p[1] = -tmp;
        }
    }
    else if (n == 3)
    {
        if (m == 0)
        {
            REAL tmp = SQRT(ADDP(7.0)) / ADDP(8.0);
            p[0] = ADDP(3.0) * tmp;
            p[1] = ADDP(5.0) * tmp;
        }
        else if (m == 1)
        {
            REAL tmp = SQRT(ADDP(42.0)) / ADDP(16.0);
            p[0] =       tmp;
            p[1] = ADDP(5.0) * tmp;
        }
        else if (m == 2)
        {
            REAL tmp = SQRT(ADDP(105.0)) / ADDP(8.0);
            p[0] =  tmp;
            p[1] = -tmp;
        }
        else if (m == 3)
        {
            REAL tmp = SQRT(ADDP(70.0)) / ADDP(16.0);
            p[0] = ADDP(3.0) * tmp;
            p[1] =      -tmp;
        }
    }
    else if (n == 4)
    {
        if (m == 0)
        {
            p[0] =  ADDP(27.0) / ADDP(64.0);
            p[1] =  ADDP(15.0) / ADDP(16.0);
            p[2] = ADDP(105.0) / ADDP(64.0);
        }
        else if (m == 1)
        {
            REAL tmp = SQRT(ADDP(10.0));
            p[0] =  ADDP(0.0);
            p[1] =  ADDP(3.0) * tmp / ADDP(16.0);
            p[2] = ADDP(21.0) * tmp / ADDP(32.0);
        }
        else if (m == 2)
        {
            REAL tmp = SQRT(ADDP(5.0));
            p[0] =   ADDP(9.0) * tmp / ADDP(32.0);
            p[1] =   ADDP(3.0) * tmp / ADDP(8.0);
            p[2] = ADDP(-21.0) * tmp / ADDP(32.0);
        }
        else if (m == 3)
        {
            REAL tmp = ADDP(3.0) * SQRT(ADDP(70.0));
            p[0] =  ADDP(0.0);
            p[1] =  tmp / ADDP(16.0);
            p[2] = -tmp / ADDP(32.0);
        }
        else if (m == 4)
        {
            REAL tmp = SQRT(ADDP(35.0));
            p[0] =  ADDP(9.0) * tmp / ADDP(64.0);
            p[1] = ADDP(-3.0) * tmp / ADDP(16.0);
            p[2] =  ADDP(3.0) * tmp / ADDP(64.0);
        }
    }


    return;
}
