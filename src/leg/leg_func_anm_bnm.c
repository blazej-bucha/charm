/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdlib.h>
#include "../prec.h"
/* ------------------------------------------------------------------------- */






void CHARM(leg_func_anm_bnm)(unsigned long nmax, unsigned long m,
                             const REAL *r, const REAL *ri,
                             REAL *anm, REAL *bnm)
{
    if (m > nmax)
    {
        /* Technically, this should be an error.  But since
         * "CHARM(leg_func_anm_bnm)" is usually used inside computationally
         * expensive parallel for loops, proper error treatment would most
         * likely require to introduce the "#pragma omp barrier" directive.
         * This, however, might slow down the computations.  Therefore, only
         * a return statement is applied.  This is not 100 % safe (in a very
         * stupid case, one might encounter, for instance, an infinite loop),
         * but should be reasonably safe in most normal circumstances. */
        return;
    }
    /* --------------------------------------------------------------------- */



    /* Let's treat the zero/non-defined coefficients */
    /* --------------------------------------------------------------------- */
    /* "amm", "bmm" */
    anm[0] = bnm[0] = PREC(0.0);


    if (m == nmax)
        return;


    /* "a_{m + 1, m}", "b_{m + 1, m}" */
    unsigned long n = m + 1;
    anm[n] = r[2 * m + 3];
    bnm[n] = PREC(0.0);
    /* --------------------------------------------------------------------- */


    /* "a_{m + 2, m}", ..., "a_{nmax, m}", "b_{m + 2, m}", ..., "b_{nmax, m}"
     * */
    /* --------------------------------------------------------------------- */
    REAL w;
    for (n = (m + 2); n <= nmax; n++)
    {
        w      = r[2 * n + 1] * ri[n - m] * ri[n + m];
        anm[n] = r[2 * n - 1] * w;
        bnm[n] = r[n - m - 1] * r[n + m - 1] * ri[2 * n - 3] * w;
    }
    /* --------------------------------------------------------------------- */


    return;
}
