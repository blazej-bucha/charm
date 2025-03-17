/* This source file should never be compiled directly.  Instead, it is always
 * include in other source files (e.g., "gfm_global_density_3d.c").  The is why
 * it relies on some symbolic constants that neither defined in this file, nor
 * they are defined in the included header files. */




/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <math.h>
#include <string.h>
#if !SGFM_GLOBAL
#   include "../mpfr/mpfr_set_real.h"
#   include "../mpfr/mpfr_get_real.h"
#   include <mpfr.h>
#   include "../mpfr/mpfr_flush_unreleased_memory.h"
#   include "../mpfr/mpfr_check_bits.h"
#endif
#include "../prec.h"
#include "../shc/shc_reset_coeffs.h"
#include "../shc/shc_check_distribution.h"
#include "../misc/misc_is_nearly_equal.h"
#include "../misc/misc_idx_4d.h"
#include "../err/err_set.h"
#include "../err/err_propagate.h"
#include "../err/err_check_distribution.h"
#include "gfm_check_p.h"
#include "gfm_check_kminkmax.h"
/* ------------------------------------------------------------------------- */






/* Macros and symbolic constants */
/* ------------------------------------------------------------------------- */
#ifndef SGFM_COMPILE
    /* If we got here, this means this source file is being compiled directly.
     * This should never be the case though.  But if it indeed happens that we
     * got here, for instance, when compiling CHarm outside the official
     * installation system, we do not want to compile this source code. */
#   define SGFM_COMPILE 0
#endif






#if SGFM_COMPILE

#define GET_NUMBER_OF_DIGITS(x) ((x) > 0) ? (int)ceil(log10((double)((x)))) : 1


/* If "SGFM_RHO_CONST" is "1", a constant mass density is assumed.  Otherwise,
 * the density is allowed to vary in 3D space.
 *
 * Such a conditional compilation is used to lower the memory requirements in
 * case of constant mass density.  Then, we need only a single "REAL" number to
 * represent the density instead of an array filled with the same number.  For
 * high-resolution applications and high integer powers, this noticeably
 * reduces the amount of working memory being used. */
#ifndef SGFM_RHO_CONST
#   define SGFM_RHO_CONST 0
#endif


/* If "SGFM_GLOBAL" is "1", global spectral gravity forward modelling is
 * assumed.  Otherwise, the cap-modified variant is used.  The latter requires
 * MPFR and GMP, so let's set the default to the global variant, where no such
 * requirements are present. */
#ifndef SGFM_GLOBAL
#   define SGFM_GLOBAL 1
#endif


/* If "PYWRAP" is set to "1", compiled will be a wrapper functions for Python
 * wrapper. */
#ifndef PYWRAP
#   define PYWRAP 0
#endif
/* ------------------------------------------------------------------------- */






static _Bool export_powers(const char *string)
{
    return (string != NULL);
}






#if !SGFM_RHO_CONST
static unsigned long get_max_nmax(unsigned long *a, size_t na)
{
    unsigned long max = a[0];
    unsigned long ai;


    for (size_t i = 1; i < na; i++)
    {
        ai = a[i];
        if (ai > max)
            max = ai;
    }


    return max;
}
#endif






static void export_shcs_powers(CHARM(shc) *shcs,
                               const char *path,
                               const char *shcs_file_format,
                               unsigned p,
                               unsigned pmax,
                               unsigned i,
                               unsigned imax,
                               CHARM(err) *err)
{
    /* Get the number of digits in "pmax" and "imax" to be able to add nice
     * leading zeros to the "p" and "i" values in the output file name */
    int ndp = GET_NUMBER_OF_DIGITS(pmax);
    int ndi = GET_NUMBER_OF_DIGITS(imax);


    /* Prepare the output file name */
    /* --------------------------------------------------------------------- */
    char buf[CHARM_ERR_MAX_FILE];
    int ret;
    _Bool format_bin = 0;
    _Bool format_mtx = 0;
    _Bool format_tbl = 0;
    _Bool format_dov = 0;


    if (!strcmp(shcs_file_format, "bin"))
        format_bin = 1;
    else if (!strcmp(shcs_file_format, "mtx"))
        format_mtx = 1;
    else if (!strcmp(shcs_file_format, "tbl"))
        format_tbl = 1;
    else if (!strcmp(shcs_file_format, "dov"))
        format_dov = 1;


    if (format_bin)
        ret = snprintf(buf, sizeof(buf), "%s-p%0*u-i%0*u.shcs",
                       path, ndp, p, ndi, i);
    else
        ret = snprintf(buf, sizeof(buf), "%s-p%0*u-i%0*u.txt",
                       path, ndp, p, ndi, i);


    if ((ret > 0) && ((size_t)ret > sizeof(buf)))
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Too small string buffer for the output file name.");
        return;

    }
    else if (ret <= 0)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Failed to prepare the output file name string.");
        return;
    }
    /* --------------------------------------------------------------------- */


    /* Write the coefficients */
    /* --------------------------------------------------------------------- */
    if (format_bin)
        CHARM(shc_write_bin)(shcs, shcs->nmax, buf, err);
    else if (format_mtx)
        CHARM(shc_write_mtx)(shcs, shcs->nmax, REAL_PRINT_FORMAT, buf, err);
    else if (format_tbl)
        CHARM(shc_write_tbl)(shcs, shcs->nmax, REAL_PRINT_FORMAT,
                             CHARM_SHC_WRITE_N, buf, err);
    else if (format_dov)
        CHARM(shc_write_dov)(shcs, shcs->nmax, REAL_PRINT_FORMAT,
                             CHARM_SHC_WRITE_N, buf, err);


    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return;
    }
    /* --------------------------------------------------------------------- */


    return;
}






#if SGFM_GLOBAL
static void global_vnmc(unsigned long nmax,
                        unsigned p,
                        unsigned i,
                        REAL c,
                        REAL *vnmc)
{
    REAL snpi;
    REAL npip3_real, ptmp_real;
    REAL p_real  = (REAL)p;
    unsigned long ip3 = (unsigned long)(i + 3);


    {
    unsigned long n;
#ifdef CHARM_OPENMP
#pragma omp parallel for default(none) shared(nmax, p, i, c, vnmc, p_real) \
    shared(ip3) private(n, snpi, npip3_real, ptmp_real)
#endif
    for (n = 0; n <= nmax; n++)
    {
        npip3_real = (REAL)(n + ip3);


        /* Compute the binomial coefficient "(n + i + 3)" choose "p" */
        snpi = npip3_real / p_real;
        for (unsigned ptmp = 1; ptmp < p; ptmp++)
        {
            ptmp_real  = (REAL)ptmp;
            snpi      *= (npip3_real - ptmp_real) / (p_real - ptmp_real);
        }


        snpi    *= PREC(2.0) / ((PREC(2.0) * (REAL)n + PREC(1.0)) *
                                npip3_real);
        vnmc[n]  = c * snpi;
    }
    }


    return;
}
#endif






CHARM_EXTERN void CHARM_CDECL
#if SGFM_GLOBAL
#   if SGFM_RHO_CONST
        CHARM(gfm_global_density_const)
#   else
        CHARM(gfm_global_density_3d)
#   endif
#else
#   if PYWRAP
#       if SGFM_RHO_CONST
            CHARM(gfm_cap_density_const_pywrap)
#       else
            CHARM(gfm_cap_density_3d_pywrap)
#       endif
#   else
#       if SGFM_RHO_CONST
            CHARM(gfm_cap_density_const)
#       else
            CHARM(gfm_cap_density_3d)
#       endif
#   endif
#endif
                                (const CHARM(shc) *shape_shcs,
                                 unsigned long shape_nmax,
                                 REAL shape_ref_radius,
#if SGFM_RHO_CONST
                                 REAL density,
#else
                                 CHARM(shc) **density_shcs,
                                 unsigned long *density_nmax,
                                 unsigned density_order,
#endif
                                 REAL grav_const,
                                 REAL mass,
                                 unsigned shape_power_min,
                                 unsigned shape_power_max,
#if !SGFM_GLOBAL
                                 unsigned radial_derivative_min,
                                 unsigned radial_derivative_max,
                                 REAL evaluation_spherical_radius,
                                 REAL integration_radius,
                                 unsigned u,
                                 unsigned v,
                                 int zone,
#   if PYWRAP
                                 int nbits_int,
#   else
                                 mpfr_prec_t nbits,
#   endif
#endif
                                 unsigned long potential_shcs_nmax,
                                 const char *shape_density_shcs_path,
                                 const char *potential_shcs_path,
                                 const char *shcs_file_format,
#if SGFM_GLOBAL
                                 CHARM(shc) *potential_shcs,
#else
                                 CHARM(shc) **potential_shcs,
#endif
                                 CHARM(err) *err)
{
    /* Substitutions */
    /* --------------------------------------------------------------------- */
    /* Minimum and maximum topography powers */
    unsigned pmin = shape_power_min;
    unsigned pmax = shape_power_max;


#if !SGFM_GLOBAL
    /* Minimum and maximum radial derivatives */
    unsigned kmin = radial_derivative_min;
    unsigned kmax = radial_derivative_max;
#endif


    /* Total number of density derivatives */
    unsigned imax;


    /* The maximum degree in "density_nmax" */
    unsigned long density_nmax_max;
    /* --------------------------------------------------------------------- */






    /* Check input arguments */
    /* --------------------------------------------------------------------- */
    CHARM(err_check_distribution)(err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return;
    }


    CHARM(shc_check_distribution)(shape_shcs, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return;
    }


#if !SGFM_RHO_CONST
    for (unsigned i = 0; i <= density_order; i++)
    {
        CHARM(shc_check_distribution)(density_shcs[i], err);
        if (!CHARM(err_isempty)(err))
        {
            CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
            return;
        }
    }
#endif


    if (potential_shcs != NULL)
    {
#if SGFM_GLOBAL
        CHARM(shc_check_distribution)(potential_shcs, err);
        if (!CHARM(err_isempty)(err))
        {
            CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
            return;
        }
#else
        for (unsigned k = kmin; k <= kmax; k++)
        {
            CHARM(shc_check_distribution)(potential_shcs[k - kmin], err);
            if (!CHARM(err_isempty)(err))
            {
                CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
                return;
            }
        }
#endif
    }


    char err_msg[CHARM_ERR_MAX_MSG];


    if (shape_nmax > shape_shcs->nmax)
    {
        sprintf(err_msg, "\"shape_nmax\" is \"%lu\", but it cannot be larger "
                         "than \"shape_shcs->nmax = %lu\".",
                shape_nmax, shape_shcs->nmax);
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       err_msg);
        return;
    }


    if (!CHARM(misc_is_nearly_equal)(shape_shcs->mu, PREC(1.0),
                                     CHARM(glob_threshold)))
    {
        sprintf(err_msg, "\"shape_shcs->mu\" is \"" REAL_PRINT_FORMAT
                         "\", but it has to be \"1.0\".", shape_shcs->mu);
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       err_msg);
        return;
    }


    if (!CHARM(misc_is_nearly_equal)(shape_shcs->r, PREC(1.0),
                                     CHARM(glob_threshold)))
    {
        sprintf(err_msg, "\"shape_shcs->r\" is \"" REAL_PRINT_FORMAT
                         "\", but it has to be \"1.0\".", shape_shcs->r);
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       err_msg);
        return;
    }


    if (shape_ref_radius <= PREC(0.0))
    {
        sprintf(err_msg, "\"shape_ref_radius\" is \"" REAL_PRINT_FORMAT
                         "\", but it has to be positive.", shape_ref_radius);
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       err_msg);
        return;
    }


    if (CHARM(misc_is_nearly_equal)(mass, PREC(0.0), CHARM(glob_threshold)))
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"mass\" cannot be zero.");
        return;
    }


    CHARM(gfm_check_p)(pmin, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return;
    }


    CHARM(gfm_check_p)(pmax, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return;
    }


    if (pmin > pmax)
    {
        sprintf(err_msg, "\"shape_power_min\" is \"%u\", but it cannot be "
                         "larger than \"shape_power_max = %u\".",
                shape_power_min, shape_power_max);
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       err_msg);
        return;
    }


#if !SGFM_GLOBAL
    CHARM(gfm_check_kminkmax)(kmin, kmax, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return;
    }


    if (evaluation_spherical_radius <= PREC(0.0))
    {
        sprintf(err_msg, "\"evaluation_spherical_radius\" is \""
                         REAL_PRINT_FORMAT "\", but it has to be positive.",
                evaluation_spherical_radius);
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       err_msg);
        return;
    }


    if ((integration_radius < PREC(0.0)) || (integration_radius > PI))
    {
        sprintf(err_msg, "\"integration_radius\" is \""
                         REAL_PRINT_FORMAT "\", but it must be from the "
                         "interval \"[0.0, pi]\".", integration_radius);
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       err_msg);
        return;
    }


    if (u > 2)
    {
        sprintf(err_msg, "\"u\" is \"%u\", but it cannot be larger than "
                         "\"2\".", u);
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       err_msg);
        return;
    }


    if (v > u)
    {
        sprintf(err_msg, "\"v\" is \"%u\", but it cannot be larger than "
                         "\"u = %u\".", v, u);
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       err_msg);
        return;
    }


    CHARM(mpfr_check_bits)(
#if PYWRAP
                           nbits_int,
#else
                           nbits,
#endif
                           err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return;
    }


#if PYWRAP
    mpfr_prec_t nbits = nbits_int;
#endif


#endif


#if !SGFM_RHO_CONST
    for (unsigned i = 0; i <= density_order; i++)
    {
        if (density_nmax[i] > density_shcs[i]->nmax)
        {
            sprintf(err_msg, "\"density_nmax[%u]\" is \"%lu\", but it cannot "
                             "be larger than \"density_shcs[%u]->nmax = "
                             "%lu\".",
                    i, density_nmax[i], i, density_shcs[i]->nmax);
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                           err_msg);
            return;
        }


        if (!CHARM(misc_is_nearly_equal)(density_shcs[i]->mu, PREC(1.0),
                                         CHARM(glob_threshold)))
        {
            sprintf(err_msg, "\"density_shcs[%u]->mu\" is \"" REAL_PRINT_FORMAT
                             "\", but it has to be \"1.0\".",
                    i, density_shcs[i]->mu);
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                           err_msg);
            return;
        }


        if (!CHARM(misc_is_nearly_equal)(density_shcs[i]->r, PREC(1.0),
                                         CHARM(glob_threshold)))
        {
            sprintf(err_msg, "\"density_shcs[%u]->r\" is \"" REAL_PRINT_FORMAT
                             "\", but it has to be \"1.0\".",
                    i, density_shcs[i]->r);
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                           err_msg);
            return;
        }
    }
#endif


    if (potential_shcs != NULL)
    {
#if SGFM_GLOBAL
        if (potential_shcs_nmax > potential_shcs->nmax)
        {
            sprintf(err_msg, "\"potential_shcs_nmax\" is \"%lu\", but it "
                             "cannot be larger than "
                             "\"potential_shcs->nmax = %lu\".",
                    potential_shcs_nmax, potential_shcs->nmax);
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                           err_msg);
            return;
        }
#else
        for (unsigned k = kmin; k <= kmax; k++)
        {
            if (potential_shcs_nmax > potential_shcs[k - kmin]->nmax)
            {
                sprintf(err_msg, "\"potential_shcs_nmax\" is \"%lu\", but it "
                                 "cannot be larger than "
                                 "\"potential_shcs[%u]->nmax = %lu\".",
                        potential_shcs_nmax, k - kmin,
                        potential_shcs[k - kmin]->nmax);
                CHARM(err_set)(err, __FILE__, __LINE__, __func__,
                               CHARM_EFUNCARG, err_msg);
                return;
            }
        }
#endif
    }


    /* Are we asked to export spherical harmonic coefficients of shape and
     * density powers? */
    _Bool export_hro_powers = export_powers(shape_density_shcs_path);


    /* And what about exporting the gravitational contributions induced by the
     * individual powers? */
    _Bool export_pot_powers = export_powers(potential_shcs_path);


    /* Check "shcs_file_format", but only if we are exporting spherical
     * harmonic coefficients of topo and density powers to data files */
    if (export_hro_powers || export_pot_powers)
    {
        if ((shcs_file_format == NULL) ||
            ((strcmp(shcs_file_format, "bin") != 0) &&
             (strcmp(shcs_file_format, "mtx") != 0) &&
             (strcmp(shcs_file_format, "tbl") != 0) &&
             (strcmp(shcs_file_format, "dov") != 0)))
        {
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                           "Unsupported value of \"shcs_file_format\".");
            return;
        }
    }
    /* --------------------------------------------------------------------- */


    /* Set all coefficients of "potential_shcs" to zero */
    if (potential_shcs != NULL)
    {
#if SGFM_GLOBAL
        CHARM(shc_reset_coeffs)(potential_shcs);
#else
        for (unsigned k = kmin; k <= kmax; k++)
            CHARM(shc_reset_coeffs)(potential_shcs[k - kmin]);
#endif
    }


#if SGFM_RHO_CONST
    imax             = 0;
    density_nmax_max = 0;
#else
    imax             = density_order;
    density_nmax_max = get_max_nmax(density_nmax, density_order);
#endif


    /* Maximum degree of the output potential */
    unsigned long nmax_pot = potential_shcs_nmax;


    CHARM(point) *grd = NULL;
    REAL *shape       = NULL;
#if SGFM_RHO_CONST
    REAL rho          = density;
#else
    REAL *rho         = NULL;
    /* The "i"-th power of "shape_ref_radius" */
    REAL shape_ref_radius_powi;
#endif
    REAL *vnmc        = NULL;
    CHARM(shc) *shcs  = NULL;
    unsigned long nmaxp, nmaxp_threshold;
    REAL c = (PREC(2.0) * PI * POW(shape_ref_radius, 3)) / mass;
#if SGFM_RHO_CONST
    c *= rho;
#endif
    REAL gm = grav_const * mass;
#if !SGFM_GLOBAL
    mpfr_t *q = NULL;
    REAL shape_ref_radius_powu = POW(shape_ref_radius, u);
    REAL m1powv                = POW(PREC(-1.0), v);
    CHARM(shc) *shcs_tmp = NULL;


    /* Truncation coefficients */
    size_t q_size = CHARM(gfm_cap_nq)(nmax_pot, pmax, kmin, kmax, imax, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return;
    }


    q = (mpfr_t *)malloc(q_size * sizeof(mpfr_t));
    if (q == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        goto EXIT;
    }
    for (size_t j = 0; j < q_size; j++)
        mpfr_init2(q[j], nbits);


    unsigned type = 9999;
    if ((u == 0) && (v == 0))
        type = CHARM_GFM_Q00;
    else if ((u == 1) && (v == 0))
        type = CHARM_GFM_Q10;
    else if ((u == 1) && (v == 1))
        type = CHARM_GFM_Q11;
    else if ((u == 2) && (v == 0))
        type = CHARM_GFM_Q20;
    else if ((u == 2) && (v == 1))
        type = CHARM_GFM_Q21;
    else if ((u == 2) && (v == 2))
        type = CHARM_GFM_Q22;


    mpfr_t rref, r, psi0;
    mpfr_inits2(nbits, rref, r, psi0, (mpfr_ptr)NULL);
    mpfr_set_REAL(rref, shape_ref_radius, MPFR_RNDN);
    mpfr_set_REAL(r, evaluation_spherical_radius, MPFR_RNDN);
    mpfr_set_REAL(psi0, integration_radius, MPFR_RNDN);
    CHARM(gfm_cap_q)(rref, r, psi0, nmax_pot, pmax, kmin, kmax, imax, zone,
                     type, nbits, q, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        goto EXIT;
    }
#endif


    /* Loop over the topography powers.
     *
     * We assume that "shape_nmax" is significantly larger than any maximum
     * degree in "density_nmax".  In that case, it is more efficient to let the
     * loop over the shape-density powers "p" be the outermost one.  This
     * allows to synthesize the shape only once for each power "p" and then
     * reused this array for each density gradient "i".  Nonetheless, if
     * "shape_nmax" is equal to or smaller than max degrees in "density_nmax",
     * everything still works fine, though it is less efficient. */
    for (unsigned p = pmin; p <= pmax; p++)
    {
        /* "nmaxp" is the maximum degree that is required to harmonically
         * analyze the shape-density function without any aliasing.  Since
         * "density_nmax" may vary for each density polynomial order, "nmaxp"
         * should technically be determined inside the loop over "i" as
         *
         *      nmaxp = p * shape_nmax + density_nmax[i];
         *
         * However, this would require to synthesize over and over again the
         * shape, each time up to a different "nmaxp".  To avoid this, we use
         * the maximum of "density_nmax", that is "density_nmax_max", and
         * synthesize the shape only once based on "density_nmax_max".
         *
         * If there are huge differences between the elements of
         * "density_nmax", this may not be the most efficient approach.  Still,
         * since we expect rather similar values inside "density_nmax", this
         * should be efficient. */
        nmaxp           = p * shape_nmax + density_nmax_max;
        nmaxp_threshold = (nmaxp > nmax_pot) ? nmax_pot : nmaxp;


        /* Prepare the Gauss--Legendre grid */
        /* ----------------------------------------------------------------- */
        grd = CHARM(crd_point_gl)(nmaxp, PREC(1.0));
        if (grd == NULL)
        {
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                           "Failed to compute the Gauss--Legendre grid.");
            goto EXIT;
        }
        /* ----------------------------------------------------------------- */


        /* Allocate memory for shape and density grids */
        /* ----------------------------------------------------------------- */
        /* Shape */
        shape = (REAL *)malloc(grd->npoint * sizeof(REAL));
        if (shape == NULL)
        {
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                           CHARM_ERR_MALLOC_FAILURE);
            goto EXIT;
        }
        CHARM(shs_point)(grd, shape_shcs, shape_nmax, shape, err);
        if (!CHARM(err_isempty)(err))
        {
            CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
            goto EXIT;
        }
        /* Remove the reference radius, so that integrated will be only the
         * masses defined by "shape_shcs" that are above this sphere. */
        {
        size_t i;
#if CHARM_OPENMP
#pragma omp parallel for default(shared) private(i)
#endif
        for (i = 0; i < grd->npoint; i++)
            shape[i] -= shape_ref_radius;
        }


#if !SGFM_RHO_CONST
        /* At first, "rho" represents the density polynomial order.  Later, it
         * is used to store the product of the "pth" shape power and of the
         * "i"th polynomial density order. */
        rho = (REAL *)malloc(grd->npoint * sizeof(REAL));
        if (rho == NULL)
        {
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                           CHARM_ERR_MALLOC_FAILURE);
            goto EXIT;
        }


        shape_ref_radius_powi = PREC(1.0);
#endif
        /* ----------------------------------------------------------------- */


        /* Allocate memory for spherical harmonic coefficients of the product
         * of (1) the density, (2) "shape_ref_radius_powi" (if !SGFM_RHO_CONST)
         * and (3) the "p"th shape power */
        /* ----------------------------------------------------------------- */
        shcs = CHARM(shc_calloc)(nmaxp_threshold, PREC(1.0), PREC(1.0));
        /* ----------------------------------------------------------------- */


        /* ----------------------------------------------------------------- */
        vnmc = (REAL *)calloc(nmaxp_threshold + 1, sizeof(REAL));
        if (vnmc == NULL)
        {
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                           CHARM_ERR_MALLOC_FAILURE);
            goto EXIT;
        }
        /* ----------------------------------------------------------------- */


        for (unsigned i = 0; i <= imax; i++)
        {
#if !SGFM_RHO_CONST
            /* Synthesize the "i"th polynomial density order */
            CHARM(shs_point)(grd, density_shcs[i], density_nmax[i], rho, err);
            if (!CHARM(err_isempty)(err))
            {
                CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
                goto EXIT;
            }
#endif


            /* Prepare the product of (1) the "i"th polynomial density order,
             * (2) "shape_ref_radius_powi" (if !SGFM_RHO_CONST) and (3) the
             * "p"th power of the shape */
            /* ------------------------------------------------------------- */
            REAL shape_ij, shape_ij_tmp;
            size_t idx;
            {
            size_t vi;
#if CHARM_OPENMP
#if SGFM_RHO_CONST
#pragma omp parallel for default(none) \
    shared(shape, rho, grd, shape_ref_radius, p) \
    private(vi, idx, shape_ij, shape_ij_tmp)
#else
#pragma omp parallel for default(none) \
    shared(shape, rho, grd, shape_ref_radius, p, shape_ref_radius_powi) \
    private(vi, idx, shape_ij, shape_ij_tmp)
#endif
#endif
            for (vi = 0; vi < grd->nlat; vi++)
            {
                for (size_t vj = 0; vj < grd->nlon; vj++)
                {
                    idx = vi * grd->nlon + vj;
                    shape_ij = shape_ij_tmp = shape[idx] / shape_ref_radius;
                    for (unsigned ptmp = 2; ptmp <= p; ptmp++)
                        shape_ij *= shape_ij_tmp;


#if SGFM_RHO_CONST
                    shape[idx] = shape_ij;
#else
                    /* After the following command, "rho[idx]" represents the
                     * product of (1) the density, (2) "shape_ref_radius_powi"
                     * and (3) the "p"th power of the shape ("Hp_rho") */
                    rho[idx] *= shape_ij * shape_ref_radius_powi;
#endif
                }
            }
            }


#if !SGFM_RHO_CONST
            /* Update "shape_ref_radius_powi" for the next density gradient if
             * any */
            shape_ref_radius_powi *= shape_ref_radius;
#endif
            /* ------------------------------------------------------------- */


            /* Harmonic analysis of Hp_rhoi */
            /* ------------------------------------------------------------- */
            shcs->mu = shcs->r = PREC(1.0);
            CHARM(sha_point)(grd,
#if SGFM_RHO_CONST
                             shape,
#else
                             rho,
#endif
                             nmaxp_threshold, shcs, err);
            if (!CHARM(err_isempty)(err))
            {
                CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
                goto EXIT;
            }


#if !SGFM_GLOBAL
            /* If "v > 0", we have to make a copy of coefficients from "shcs"
             * up to degree "CHARM_MIN(shcs->nmax, v)", because these
             * coefficients in "shcs" will be modified and later the original
             * coefficients have to be restored back.  See below for more
             * details. */
            if (v > 0)
            {
                shcs_tmp = CHARM(shc_calloc)(CHARM_MIN(shcs->nmax,
                                                       (unsigned long)v),
                                             PREC(1.0), PREC(1.0));
                for (unsigned long m = 0; m <= shcs_tmp->nmax; m++)
                {
                    for (unsigned long n = m; n <= shcs_tmp->nmax; n++)
                    {
                        shcs_tmp->c[m][n - m] = shcs->c[m][n - m];
                        shcs_tmp->s[m][n - m] = shcs->s[m][n - m];
                    }
                }
            }
#endif
        /* ----------------------------------------------------------------- */


            if (export_hro_powers)
            {
                export_shcs_powers(shcs, shape_density_shcs_path,
                                   shcs_file_format, p, pmax, i, imax, err);
                if (!CHARM(err_isempty)(err))
                {
                    CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
                    goto EXIT;
                }
            }
            /* ------------------------------------------------------------- */


            /* Compute the gravitational contribution of the "p"th shape power
             * and "i"th polynomial density order */
#if SGFM_GLOBAL
            global_vnmc(shcs->nmax, p, i, c, vnmc);
#else
            for (unsigned k = kmin; k <= kmax; k++)
            {
                size_t idx = CHARM(misc_idx_4d)(k - kmin, p - 1, i, 0, pmax,
                                                imax + 1, nmax_pot + 1);
                for (unsigned long n = 0; n <= shcs->nmax; n++)
                    vnmc[n] = m1powv * c * shape_ref_radius_powu *
                              mpfr_get_REAL(q[idx + n], MPFR_RNDN);
#endif


                /* Final computation of gravitational coefficients due to the
                 * "p"th shape power and "i"th polynomial density order */
                {
                unsigned long m;
#if CHARM_OPENMP
#pragma omp parallel for default(shared) private(m)
#endif
                for (m = 0; m <= shcs->nmax; m++)
                {
                    for (unsigned long n = m; n <= shcs->nmax; n++)
                    {
                        shcs->c[m][n - m] *= vnmc[n];
                        shcs->s[m][n - m] *= vnmc[n];
                    }
                }
                }
                shcs->mu = gm;
                shcs->r  = shape_ref_radius;


                if (export_pot_powers)
                {
                    export_shcs_powers(shcs, potential_shcs_path,
                                       shcs_file_format, p, pmax, i, imax,
                                       err);
                    if (!CHARM(err_isempty)(err))
                    {
                        CHARM(err_propagate)(err, __FILE__, __LINE__,
                                             __func__);
                        goto EXIT;
                    }
                }
                /* --------------------------------------------------------- */


                /* Computation of the output gravitational coefficients that
                 * sum all the partial gravitational contributions from "shcs"
                 * */
                /* --------------------------------------------------------- */
                if (potential_shcs != NULL)
                {
                    {
                    unsigned long m;
#if CHARM_OPENMP
#pragma omp parallel for default(shared) private(m)
#endif
                    for (m = 0; m <= shcs->nmax; m++)
                    {
                        for (unsigned long n = m; n <= shcs->nmax; n++)
                        {
                            potential_shcs
#if !SGFM_GLOBAL
                                          [k - kmin]
#endif
                                          ->c[m][n - m] += shcs->c[m][n - m];
                            potential_shcs
#if !SGFM_GLOBAL
                                          [k - kmin]
#endif
                                          ->s[m][n - m] += shcs->s[m][n - m];
                        }
                    }
                    }
                }
#if !SGFM_GLOBAL

                /* Since we have multiplied "shcs" by "vnmc" above, now we have
                 * to divide "shcs" by "vnmc" to get things back as they were.
                 * This is necessary only in case of cap-modified SGFM and if
                 * "kmin != kmax".  The rational is that for the next radial
                 * derivative, we need a new set of truncation coefficients,
                 * but the original "shcs".  So we either have to store an
                 * additinal copy of "shcs", which requires significant amount
                 * of additional RAM or we now divide "shcs" by "vnmc" and in
                 * the next loop iteration of "k", we will have the same "shcs"
                 * (except perhaps for rounding errors).
                 *
                 * Things get complicated if "v > 0".  In that case, "vnmc" up
                 * to degree "n < 0" are zero.  Therefore, these coefficients
                 * were extracted from "shcs" to "shcs_tmp" before and now they
                 * have to be restored. */
                if (kmin != kmax)
                {
                    {
                    unsigned long m;
#if CHARM_OPENMP
#pragma omp parallel for default(shared) private(m)
#endif
                    for (m = 0; m <= shcs->nmax; m++)
                    {
                        for (unsigned long n = m; n <= shcs->nmax; n++)
                        {
                            shcs->c[m][n - m] /= vnmc[n];
                            shcs->s[m][n - m] /= vnmc[n];
                        }
                    }
                    }


                    /* Now restore the low-degree coefficients */
                    if (v > 0)
                    {
                        for (unsigned long m = 0; m <= shcs_tmp->nmax; m++)
                        {
                            for (unsigned long n = m; n <= shcs_tmp->nmax; n++)
                            {
                                shcs->c[m][n - m] = shcs_tmp->c[m][n - m];
                                shcs->s[m][n - m] = shcs_tmp->s[m][n - m];
                            }
                        }
                    }
                }
            }


            CHARM(shc_free)(shcs_tmp);
#endif
            /* ------------------------------------------------------------- */
        }


        free(shape); shape = NULL;
#if !SGFM_RHO_CONST
        free(rho); rho = NULL;
#endif
        free(vnmc); vnmc = NULL;
        CHARM(crd_point_free)(grd); grd = NULL;
        CHARM(shc_free)(shcs); shcs = NULL;
    }


    /* Finally, set the correct "mu" and "r" values of the output coefficients
     * */
    if (potential_shcs != NULL)
    {
#if SGFM_GLOBAL
        potential_shcs->mu = gm;
        potential_shcs->r  = shape_ref_radius;
#else
        for (unsigned k = kmin; k <= kmax; k++)
        {
            potential_shcs[k - kmin]->mu = gm;
            potential_shcs[k - kmin]->r  = shape_ref_radius;
        }
#endif
    }


EXIT:
    free(shape);
#if !SGFM_RHO_CONST
    free(rho);
#endif
    free(vnmc);
    CHARM(crd_point_free)(grd);
    CHARM(shc_free)(shcs);
#if !SGFM_GLOBAL
    if (q != NULL)
        for (size_t j = 0; j < q_size; j++)
            mpfr_clear(q[j]);
    free(q);

    mpfr_clears(rref, r, psi0, (mpfr_ptr)NULL);
    mpfr_free_cache();
    FLUSH_UNRELEASED_MEMORY;
#endif


    return;
}
#endif

