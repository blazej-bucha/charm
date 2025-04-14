/* IMPORTANT: To compile this code, CHarm must be compiled with the MPFR
 * support, that is, the "--enable-mpfr" installation flag must be used during
 * the installation.  If you didn't enable MPFR, manually delete (or comment
 * out) all lines between the following marks:
 *
 * @@@
 *
 * lines to be deleted or commented out
 *
 * !!!
 *
 * */


#include <stdio.h>
#include <stdlib.h>
/* @@@ */
#include <mpfr.h>
/* !!! */
#include <charm/charm.h>


int main(void)
{
    /* INPUTS */
    /* ===================================================================== */
    /* Shape of the gravitating body */
    /* --------------------------------------------------------------------- */
    /* Upper integration limit in spherical radius (topography) */
    /* ..................................................................... */
    /* Path to spherical harmonic coefficients defining the shape of the Moon
     * */
    char path_shcs_shape[] = "../../data/input/MoonTopo2600p_to10-tbl.txt";


    /* Maximum harmonic degree of the lunar topography defined by coefficients
     * in "path_shcs_shape" */
    unsigned long shape_nmax = 10;
    /* ..................................................................... */


    /* Lower integration limit in spherical radius (sphere) */
    /* ..................................................................... */
    /* Radius of the reference sphere */
    double rref = 1728000.0;
    /* ..................................................................... */
    /* --------------------------------------------------------------------- */


    /* Density of the gravitating body */
    /* --------------------------------------------------------------------- */
    /* Path to spherical harmonic coefficients defining the zero- and
     * first-order polynomial density coefficients */
    char path_shcs_density0[] = "../../data/input/moon-rho0_to10.tbl";
    char path_shcs_density1[] = "../../data/input/moon-rho1_to10.tbl";
    char *path_shcs_density[2];
    path_shcs_density[0] = path_shcs_density0;
    path_shcs_density[1] = path_shcs_density1;


    /* Maximum harmonic degrees of the polynomial density coefficients in
     * "path_shcs_density". */
    unsigned long density_nmax[2] = {10, 10};


    /* Order of the density polynomial */
    unsigned density_order = 1;
    /* --------------------------------------------------------------------- */


    /* Other inputs */
    /* --------------------------------------------------------------------- */
    /* Newton's gravitational constant ("kg^-1 * m^3 * s^-2") */
    double gc = 6.67430e-11;


    /* Mass of the Moon ("kg") */
    double mass = 7.346e22;


    /* Minimum and maximum topography powers.  In theory, "pmax" should be
     * infinitely large, so higher "pmax" means more complete forward modelling
     * (but also longer computation times). */
    unsigned pmin = 1;
    unsigned pmax = 10;


    /* Maximum harmonic degree of the output gravitational potential.  Due to
     * the non-linear relation between topography and its implied gravitational
     * field, the potential series is spectrally unlimited even if "shape_nmax"
     * is finite (the only exception is "shape_nmax = 0").  Therefore, this
     * value should be as high as reasonably possible in order to mitigate
     * forward modelling errors. */
    unsigned long nmax_potential = 100;
    /* --------------------------------------------------------------------- */
    /* ===================================================================== */






    /* ===================================================================== */
    /* --------------------------------------------------------------------- */
    /* Create CHarm's error structure */
    charm_err *err = charm_err_init();
    if (err == NULL)
    {
        fprintf(stderr, "Failed to initialize \"err\".\n");
        exit(CHARM_FAILURE);
    }


    /* Allocate memory for spherical harmonic coefficients that will hold the
     * shape coefficients */
    charm_shc *shape_shcs = charm_shc_malloc(shape_nmax, 1.0, 1.0);
    if (shape_shcs == NULL)
    {
        fprintf(stderr, "Failed to initialize \"shape_shcs\".\n");
        exit(CHARM_FAILURE);
    }


    /* Read the shape coefficients */
    charm_shc_read_tbl(path_shcs_shape, shape_nmax, shape_shcs, err);
    charm_err_handler(err, 1);


    /* Allocate memory for spherical harmonic coefficients that will hold the
     * polynomial density coefficients and read the coefficients from input
     * files */
    charm_shc *density_shcs[2];
    for (unsigned i = 0; i <= density_order; i++)
    {
        density_shcs[i] = charm_shc_malloc(density_nmax[i], 1.0, 1.0);
        if (density_shcs[i] == NULL)
        {
            fprintf(stderr, "Failed to initialize \"density_shcs[i]\".\n");
            exit(CHARM_FAILURE);
        }


        charm_shc_read_tbl(path_shcs_density[i], density_nmax[i],
                           density_shcs[i], err);
        charm_err_handler(err, 1);
    }


    /* Allocate coefficients for the output gravitational potential */
    charm_shc *potential_global_shcs = charm_shc_malloc(nmax_potential,
                                                        1.0, 1.0);
    if (potential_global_shcs == NULL)
    {
        fprintf(stderr, "Failed to initialize \"potential_global_shcs\".\n");
        exit(CHARM_FAILURE);
    }
    /* --------------------------------------------------------------------- */


    /* Global gravity forward modelling using 3D density */
    /* --------------------------------------------------------------------- */
    charm_gfm_global_density_3d(shape_shcs, shape_nmax, rref,
                                density_shcs, density_nmax, density_order,
                                gc, mass,
                                pmin, pmax,
                                nmax_potential,
                                NULL, NULL, NULL,
                                potential_global_shcs,
                                err);
    charm_err_handler(err, 1);


    /* Now we can compute, say, the gravitational potential using
     * "potential_global_shcs".  Let's do this in the Gauss--Legendre grid that
     * passes above all lunar masses. */
    /* ..................................................................... */
    /* Create the Gauss--Legendre grid */
    charm_point *grd = charm_crd_point_gl(potential_global_shcs->nmax,
                                          potential_global_shcs->r + 25000.0);
    if (grd == NULL)
    {
        fprintf(stderr, "Failed to compute the Gauss--Legendre grid "
                        "\"grd\".\n");
        exit(CHARM_FAILURE);
    }


    /* Allocate memory for the gravitational potential to be synthesized */
    double *vgfm = (double *)malloc(grd->npoint * sizeof(double));
    if (vgfm == NULL)
    {
        fprintf(stderr, "Failed to allocate memory for the gravitational "
                        "potential \"vgfm\".\n");
        exit(CHARM_FAILURE);
    }


    /* Do the synthesis of the gravitational potential */
    charm_shs_point(grd, potential_global_shcs, potential_global_shcs->nmax,
                    vgfm, err);
    charm_err_handler(err, 1);
    /* ..................................................................... */
    /* --------------------------------------------------------------------- */


    /* In CHarm, gravity forward modelling using lateral density is fairly
     * similar to using 3D density, so we skip it. */


    /* Global gravity forward modelling using constant density */
    /* --------------------------------------------------------------------- */
    /* Constant density of the lunar crust */
    double density_const = 2550.0;


    charm_gfm_global_density_const(shape_shcs, shape_nmax, rref,
                                   density_const,
                                   gc, mass,
                                   pmin, pmax,
                                   nmax_potential,
                                   NULL, NULL, NULL,
                                   potential_global_shcs,
                                   err);
    charm_err_handler(err, 1);


    /* You could now similarly synthesize the gravitational potential, this
     * time being implied by the constant density model. */
    /* --------------------------------------------------------------------- */


    /* @@@ */
    /* Spatially restricted gravity forward modelling of near-zone masses using
     * 3D density */
    /* --------------------------------------------------------------------- */
    /* Minimum and maximum order of the radial derivatives of the output
     * quantity */
    unsigned kmin = 0;
    unsigned kmax = 3;


    /* Radius of the sphere, on which evaluation points reside */
    double r = rref + 25000.0;


    /* Integration radius in radians */
    double psi0 = 1.0;


    /* Order of the potential derivative ("0" for potential, "1" for quantities
     * related to gravitational vector elements, and "2" for quantities related
     * to gravitational tensor elements) */
    unsigned u = 0;


    /* Order of the potential derivative with respect to the spherical distance
     * */
    unsigned v = 0;


    /* We want to integrate all masses up to distance "psi0" from evalution
     * points.  To integrate masses beyond "psi0", set "zone" to
     * "CHARM_GFM_FAR_ZONE". */
    int zone = CHARM_GFM_NEAR_ZONE;


    /* Number of bits to represent significands of floating points numbers used
     * to evaluate truncation coefficients */
    mpfr_prec_t nbits = 512;


    charm_shc **potential_cap_shcs = malloc((kmax - kmin + 1) *
                                            sizeof(charm_shc *));
    if (potential_cap_shcs == NULL)
    {
        fprintf(stderr, "Failed to initialize \"potential_cap_shcs\".\n");
        exit(CHARM_FAILURE);
    }
    for (unsigned k = kmin; k <= kmax; k++)
    {
        potential_cap_shcs[k - kmin] = charm_shc_malloc(nmax_potential,
                                                        1.0, 1.0);
        if (potential_cap_shcs[k - kmin] == NULL)
        {
            fprintf(stderr, "Failed to initialize "
                            "\"potential_cap_shcs[k - kmin]\".\n");
            exit(CHARM_FAILURE);
        }
    }


    charm_gfm_cap_density_3d(shape_shcs, shape_nmax, rref,
                             density_shcs, density_nmax, density_order,
                             gc, mass,
                             pmin, pmax,
                             kmin, kmax,
                             r,
                             psi0,
                             u, v,
                             zone,
                             nbits,
                             nmax_potential,
                             NULL,
                             NULL,
                             NULL,
                             potential_cap_shcs,
                             err);
    charm_err_handler(err, 1);
    /* --------------------------------------------------------------------- */


    /* Spatially restricted gravity forward modelling of near-zone masses using
     * constant density */
    /* --------------------------------------------------------------------- */
    charm_gfm_cap_density_const(shape_shcs, shape_nmax, rref,
                                density_const,
                                gc, mass,
                                pmin, pmax,
                                kmin, kmax,
                                r,
                                psi0,
                                u, v,
                                zone,
                                nbits,
                                nmax_potential,
                                NULL,
                                NULL,
                                NULL,
                                potential_cap_shcs,
                                err);
    charm_err_handler(err, 1);
    /* --------------------------------------------------------------------- */


    /* Let's compute some Molodensky's truncation coefficients, which are used
     * internally whenever calling the "charm_gfm_cap_density_*" routines. */
    /* --------------------------------------------------------------------- */
    /* Allocate memory for trunction coefficients */
    /* ..................................................................... */
    /* Get the number of trunction coefficients */
    size_t q_size = charm_gfm_cap_nq(nmax_potential, pmax, kmin, kmax,
                                     density_order, err);
    charm_err_handler(err, 1);


    mpfr_t *q = (mpfr_t *)malloc(q_size * sizeof(mpfr_t));
    if (q == NULL)
    {
        fprintf(stderr, CHARM_ERR_MALLOC_FAILURE);
        exit(CHARM_FAILURE);
    }
    for (size_t j = 0; j < q_size; j++)
        mpfr_init2(q[j], nbits);
    /* ..................................................................... */


    /* To compute the truncation coefficients, CHarm requires some floating
     * point numbers to be provided as the "mpfr_t" data type, so let's do the
     * conversion.  We chose "rref", "r" and "psi0" such that they can be
     * represented as "double"s exactly, so using "mpfr_set_d" should be safe.
     * Otherwise, say, if "psi0 = 0.1234567", use "mpfr_set_str" instead. */
    /* ..................................................................... */
    mpfr_t mpfr_rref, mpfr_r, mpfr_psi0;
    mpfr_inits2(nbits, mpfr_rref, mpfr_r, mpfr_psi0, (mpfr_ptr)NULL);
    mpfr_set_d(mpfr_rref, rref, MPFR_RNDN);
    mpfr_set_d(mpfr_r, r, MPFR_RNDN);
    mpfr_set_d(mpfr_psi0, psi0, MPFR_RNDN);
    /* ..................................................................... */


    /* Compute the trucation coefficients */
    /* ..................................................................... */
    charm_gfm_cap_q(mpfr_rref, mpfr_r, mpfr_psi0,
                    nmax_potential, pmax, kmin, kmax, density_order,
                    CHARM_GFM_NEAR_ZONE, CHARM_GFM_Q00, nbits, q, err);
    charm_err_handler(err, 1);
    /* ..................................................................... */


    /* Free the memory associated with the trunction coefficients */
    /* ..................................................................... */
    mpfr_clears(mpfr_rref, mpfr_r, mpfr_psi0, (mpfr_ptr)NULL);
    for (size_t j = 0; j < q_size; j++)
        mpfr_clear(q[j]);
    free(q);


    /* NOTE: At this point, the memory associated with "q" should be released.
     * If you are linking against "glibc" (which is the case on most Linux
     * distributions), the memory may not, however, be released completely.
     * This is not because of some bug in MPFR or CHarm, but instead it seems
     * to be the default behaviour of the GNU's "glibc" library.  It seems that
     * if you allocate a large amount of small memory chunks using malloc, as
     * it happens in our code with
     *
     *      for (size_t j = 0; j < q_size; j++)
     *          mpfr_init2(q[j], nbits);
     *
     * then "glibc" may not release the memory completely when we called
     *
     *      for (size_t j = 0; j < q_size; j++)
     *          mpfr_clear(q[j]);
     *
     * In our example, "q_size" is fairly low though, so this should now be
     * a problem in this case.
     *
     * If you are linking against "glibc", your "q_size" is large and you want
     * to release all the memory, call now
     *
     *      malloc_trim(0);
     *
     * Importantly, "malloc_trim" is a GNU extension, so you can call it only
     * if you are linking agains "glibc".  If you are using, say, "musl", you
     * must not call "malloc_trim".
     *
     * Luckily, CHarm detects on compilation whether we are linking against
     * "glibc".  If so, "malloc_trim" is called before returning from functions
     * dealing internally with truncation coefficients.  Therefore, you do not
     * have to worry that CHarm will leave some unreleased memory.  But if you
     * allocate memory for truncation coefficients outside CHarm, as we did in
     * this example, it is your own responsibility to release or to not release
     * the memory you allocated.
     *
     * For more details, look inside "src/mpfr/mpfr_flush_unreleased_memory.h".
     *
     * */
    /* ..................................................................... */
    /* --------------------------------------------------------------------- */
    /* !!! */


    /* Free the heap memory */
    /* --------------------------------------------------------------------- */
    charm_err_free(err);
    charm_shc_free(shape_shcs);
    for (unsigned i = 0; i <= density_order; i++)
        charm_shc_free(density_shcs[i]);
    charm_shc_free(potential_global_shcs);
    for (unsigned k = kmin; k <= kmax; k++)
        charm_shc_free(potential_cap_shcs[k - kmin]);
    free(potential_cap_shcs);
    charm_crd_point_free(grd);
    free(vgfm);
    /* --------------------------------------------------------------------- */


    /* --------------------------------------------------------------------- */
    printf("\nGreat, all done!\n");
    exit(CHARM_SUCCESS);
    /* --------------------------------------------------------------------- */
    /* ===================================================================== */
}
