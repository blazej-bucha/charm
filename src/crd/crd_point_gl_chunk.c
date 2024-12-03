/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <math.h>
#include "../prec.h"
#include "../misc/misc_is_nearly_equal.h"
#include "../misc/misc_check_radius.h"
#include "../err/err_set.h"
#include "../err/err_propagate.h"
#include "crd_point_quad_l.h"
#include "crd_point_quad_nlat_north.h"
#include "crd_point_gl_chunk.h"
/* ------------------------------------------------------------------------- */






/* Symbolic constants related to this file only (not a part of API) */
/* ------------------------------------------------------------------------- */
/* Maximum number of iterations to seek for a particular latitude of the
 * Gauss-Legendre grid. */
#undef GL_MAX_ITER
#define GL_MAX_ITER 1000
/* ------------------------------------------------------------------------- */






CHARM(point) *CHARM(crd_point_gl_chunk)(unsigned long nmax,
                                        REAL r,
                                        size_t local_nlat,
                                        size_t local_0_start,
                                        CHARM(err) *err)
{
    /* Some simple error checks */
    /* --------------------------------------------------------------------- */
    CHARM(misc_check_radius)(r, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return NULL;
    }
    /* --------------------------------------------------------------------- */


    /* Initialize a "CHARM(point)" structure */
    /* --------------------------------------------------------------------- */
    /* Shape of the full Gauss--Legendre grid */
    size_t nlat, nlon;
    CHARM(crd_point_gl_shape)(nmax, &nlat, &nlon);


    size_t local_nlat_north = CHARM(crd_point_quad_nlat_north)(local_nlat,
                                                    local_0_start, nlat,
                                                    CHARM_CRD_POINT_GRID_GL,
                                                    nmax, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return NULL;
    }


    CHARM(point) *glg = CHARM(crd_point_calloc)(CHARM_CRD_POINT_GRID_GL,
                                                local_nlat, nlon);
    if (glg == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        return NULL;
    }


#if HAVE_MPI
    if (glg->local_nlat == 0)
        /* This is a valid case, so skip the computation of latitudes,
         * integration weights and spherical radii and continue with
         * longitudes.  */
        goto LONGITUDES;
#endif
    /* --------------------------------------------------------------------- */


    /* Some pre-computations */
    /* --------------------------------------------------------------------- */
    const unsigned long L = CHARM(crd_point_quad_l)(nmax);
    const REAL L_fp = (REAL)L;
    const REAL c1 = L_fp + PREC(0.5);
    /* --------------------------------------------------------------------- */


    /* Latitudes */
    /* --------------------------------------------------------------------- */
    unsigned int ERROR_glob = 0;
    unsigned long imin = local_0_start;
    unsigned long imax = local_0_start + local_nlat_north;
    REAL thold = CHARM(glob_threshold);


#if HAVE_OPENMP
#pragma omp parallel default(none) shared(glg, L, L_fp, c1, ERROR_glob) \
    shared(local_nlat, imin, imax, thold)
#endif
    {
    REAL z, z1, p1, p2, p3, pp;
    unsigned long it;
    unsigned int ERROR = 0;
    unsigned long north, south;


    unsigned long i;
#if HAVE_OPENMP
#pragma omp for
#endif
    for (i = imin; i < imax; i++)
    {
        z  = COS(PI * ((REAL)(i + 1) - PREC(0.25)) / c1);
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
        }
        while ((FABS(z - z1) > EPS) && (it < GL_MAX_ITER));


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


        north = i - imin;
        south = local_nlat - 1 - north;


        glg->lat[south] = -ASIN(z);
        glg->w[south] = PREC(2.0) / ((PREC(1.0) - z * z) * pp * pp);


        /* Do not apply the symmetry property at the equator in order not to
         * change the sign of the near-zero latitude. */
        if (CHARM(misc_is_nearly_equal)(glg->lat[south], PREC(0.0), thold))
            continue;


        glg->lat[north] = -glg->lat[south];
        glg->w[north] = glg->w[south];
    }


#if HAVE_OPENMP
#pragma omp critical
#endif
    {
    if (ERROR)
        ERROR_glob += ERROR;
    }
    }


    if (ERROR_glob > 0)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Error in computing latitudes of the "
                       "Gauss--Legendre grid.  Either the maximum "
                       "number of iterations was exceeded or at least "
                       "one latitude was detected to be not a real number.");
        CHARM(crd_point_free)(glg);
        return NULL;
    }
    /* --------------------------------------------------------------------- */


    /* Spherical radii */
    /* --------------------------------------------------------------------- */
    for (unsigned long i = 0; i < local_nlat; i++)
        glg->r[i] = r;
    /* --------------------------------------------------------------------- */


    /* Longitudes */
    /* --------------------------------------------------------------------- */
#if HAVE_MPI
LONGITUDES:
    ;  /* To avoid declaration after a label */
#endif


    const REAL c = PI / (REAL)(CHARM(crd_point_quad_l)(nmax));


    for (unsigned long i = 0; i < glg->nlon; i++)
        glg->lon[i] = c * (REAL)(i);
    /* --------------------------------------------------------------------- */


    return glg;
}
