/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include "../prec.h"
/* ------------------------------------------------------------------------- */






/* Symbolic constants related to this file only (not a part of API) */
/* ------------------------------------------------------------------------- */
/* Maximum number of iterations to seek for a particular latitude of the
 * Gauss-Legendre grid. */
#undef GL_MAX_ITER
#define GL_MAX_ITER 1000
/* ------------------------------------------------------------------------- */






CHARM(point) *CHARM(crd_point_gl)(unsigned long nmax, REAL r)
{
    /* Some simple error checks */
    /* --------------------------------------------------------------------- */
    if (r <= PREC(0.0))
        return NULL;
    /* --------------------------------------------------------------------- */


    /* Initialize a "CHARM(point)" structure and "w" based on the "nmax"
     * value. */
    /* --------------------------------------------------------------------- */
    unsigned long L = nmax + 1;


    CHARM(point) *glg = CHARM(crd_point_calloc)(CHARM_CRD_POINT_GRID_GL,
                                                L, 2 * L);
    if (glg == NULL)
        return NULL;
    /* --------------------------------------------------------------------- */


    /* Some pre-computations */
    /* --------------------------------------------------------------------- */
    REAL L_fp = (REAL)L;


    /* Note that "m" is rounded to the greatest integer less than or equal
     * to "(L + 1) / 2". */
    unsigned long m = (L + 1) / 2;


    REAL c = L_fp + PREC(0.5);
    /* --------------------------------------------------------------------- */


    /* Latitudes */
    /* --------------------------------------------------------------------- */
    unsigned int ERROR_glob = 0;


#if CHARM_PARALLEL
#pragma omp parallel default(none) shared(glg, L, L_fp, m, c, ERROR_glob)
#endif
    {
    REAL z, z1, p1, p2, p3, pp;
    unsigned long it;
    unsigned int ERROR = 0;


#if CHARM_PARALLEL
#pragma omp for
#endif
    for (unsigned long i = 0; i < m; i++)
    {
        z  = COS(PI * ((REAL)(i + 1) - PREC(0.25)) / c);
        it = 0;


        do
        {
            p1 = PREC(1.0);
            p2 = PREC(0.0);


            for (unsigned long j = 1; j <= L; j++)
            {
                p3 = p2;
                p2 = p1;
                p1 = ((REAL)(2 * j - 1) * z * p2 - (REAL)(j - 1) * p3)
                     / (REAL)j;
            }


            pp = L_fp * (z * p1 - p2) / (z * z - PREC(1.0));
            z1 = z;
            z -= p1 / pp;


            it++;
        } while ((FABS(z - z1) > EPS) && (it < GL_MAX_ITER));


        if (it >= GL_MAX_ITER)
        {
            /* Here, we make a note that the iterations did not converge by
             * setting "ERROR" to "1".  After the "for" loop, we
             * check whether or not the iterations converged across the
             * threads.  If the convergence was not achieved for at least one
             * thread and one latitude, we return a NULL pointer. */
            ERROR = 1;
        }


#if HAVE_ISFINITE
        if (!isfinite(z))
            ERROR = 1;
#endif


        glg->lat[L - 1 - i] = -ASIN(z);
        glg->lat[i]         = -glg->lat[L - 1 - i];


        glg->w[L - 1 - i] = PREC(2.0) / ((PREC(1.0) - z * z) * pp * pp);
        glg->w[i]         = glg->w[L - 1 - i];
    }


#if CHARM_PARALLEL
#pragma omp critical
#endif
    {
    if (ERROR)
        ERROR_glob += ERROR;
    }
    }


    if (ERROR_glob > 0)
    {
        CHARM(crd_point_free)(glg);
        return NULL;
    }
    /* --------------------------------------------------------------------- */


    /* Longitudes */
    /* --------------------------------------------------------------------- */
    c = PI / L_fp;


#if CHARM_PARALLEL
#pragma omp parallel for default(none) shared(glg, L, c)
#endif
    for (unsigned long j = 0; j < (2 * L); j++)
        glg->lon[j] = c * (REAL)(j);
    /* --------------------------------------------------------------------- */


    /* Spherical radii */
    /* --------------------------------------------------------------------- */
#if CHARM_PARALLEL
#pragma omp parallel for default(none) shared(glg, L, r)
#endif
    for (unsigned long i = 0; i < L; i++)
        glg->r[i] = r;
    /* --------------------------------------------------------------------- */


    return glg;
}
