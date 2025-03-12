/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <mpfr.h>
#include "../prec.h"
#include "../misc/misc_is_nearly_equal.h"
#include "../err/err_set.h"
#include "../err/err_propagate.h"
/* ------------------------------------------------------------------------- */






#if CAP_SGFM_COMPILE
/* ------------------------------------------------------------------------- */
#ifndef PYWRAP
#   define PYWRAP 0
#endif
/* ------------------------------------------------------------------------- */





CHARM_EXTERN void CHARM_CDECL
#if PYWRAP
     CHARM(gfm_cap_density_lateral_pywrap)
#else
     CHARM(gfm_cap_density_lateral)
#endif
                                   (const CHARM(shc) *shape_shcs,
                                    unsigned long shape_nmax,
                                    REAL shape_ref_radius,
                                    CHARM(shc) *density_shcs,
                                    unsigned long density_nmax,
                                    REAL grav_const,
                                    REAL mass,
                                    unsigned shape_power_min,
                                    unsigned shape_power_max,
                                    unsigned radial_derivative_min,
                                    unsigned radial_derivative_max,
                                    REAL evaluation_spherical_radius,
                                    REAL integration_radius,
                                    unsigned u,
                                    unsigned v,
                                    int zone,
#if PYWRAP
                                    int nbits_int,
#else
                                    mpfr_prec_t nbits,
#endif
                                    unsigned long potential_shcs_nmax,
                                    const char *shape_density_shcs_path,
                                    const char *potential_shcs_path,
                                    const char *shcs_file_format,
                                    CHARM(shc) **potential_shcs,
                                    CHARM(err) *err)
{
    /* Check input arguments.  Most of the checks are left to
     * "gfm_cap_density_3d".  Here, we check only "density_nmax" and
     * "density_shcs", so that the error messages make sense with respect to
     * input paramters to this function. */
    /* --------------------------------------------------------------------- */
    if (density_nmax > density_shcs->nmax)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"density_nmax\" cannot be larger than "
                       "\"density_shcs->nmax\".");
        return;
    }


    if (!CHARM(misc_is_nearly_equal)(density_shcs->mu, PREC(1.0),
                                     CHARM(glob_threshold)))
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"density_shcs->mu\" have to be \"1.0\".");
        return;
    }


    if (!CHARM(misc_is_nearly_equal)(density_shcs->r, PREC(1.0),
                                     CHARM(glob_threshold)))
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"density_shcs->r\" have to be \"1.0\".");
        return;
    }
    /* --------------------------------------------------------------------- */


#if PYWRAP
    mpfr_prec_t nbits = nbits_int;
#endif
    CHARM(gfm_cap_density_3d)(shape_shcs,
                              shape_nmax,
                              shape_ref_radius,
                              &density_shcs,
                              &density_nmax,
                              0,
                              grav_const,
                              mass,
                              shape_power_min,
                              shape_power_max,
                              radial_derivative_min,
                              radial_derivative_max,
                              evaluation_spherical_radius,
                              integration_radius,
                              u,
                              v,
                              zone,
                              nbits,
                              potential_shcs_nmax,
                              shape_density_shcs_path,
                              potential_shcs_path,
                              shcs_file_format,
                              potential_shcs,
                              err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return;
    }


    return;
}
#endif
