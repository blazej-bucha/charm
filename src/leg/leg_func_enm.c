/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdlib.h>
#include "../prec.h"
#include "leg_func_enm.h"
/* ------------------------------------------------------------------------- */






void CHARM(leg_func_enm)(unsigned long nmax,
                         unsigned long m,
                         const REAL *r,
                         const REAL *ri,
                         REAL *enm)
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
    /* "emm" */
    enm[m] = PREC(0.0);


    if (m == nmax)
        return;
    /* --------------------------------------------------------------------- */


    /* "e_{m + 1, m}", ..., "e_{nmax, m}" */
    /* --------------------------------------------------------------------- */
    for (unsigned long n = m + 1; n <= nmax; n++)
        enm[n] = r[2 * n + 1] * r[n - m] * r[n + m] * ri[2 * n - 1];
    /* --------------------------------------------------------------------- */


    return;
}
