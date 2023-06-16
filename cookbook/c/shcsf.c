#include <stdio.h>
#include <stdlib.h>
#include <charm/charmf.h>


int main(void)
{
    /* INPUTS */
    /* ===================================================================== */
    /* Define the path to an input "gfc" text file with spherical harmonic
     * coefficients.  For details on the structure of the "gfc" file, see the
     * description of the "charm_shc_read_gfc" function in the "charm_shc"
     * module. */
    char shcs_in_file[] = "../../data/input/EGM96-degree10.gfc";


    /* Maximum harmonic degree to initialize, read and write spherical harmonic
     * coefficients */
    unsigned long nmax = 10;


    /* Define the path to an output binary file with spherical harmonic
     * coefficients.  For details on the structure of the output binary file,
     * see the description of the "charmf_shc_write_bin" function in the
     * "charmf_shc" module. */
    char shcs_out_file[] = "../../data/output/EGM96-degree10-mtx.shcs";
    /* ===================================================================== */






    /* ===================================================================== */
    /* Read the spherical harmonic coefficients */
    /* --------------------------------------------------------------------- */
    /* Initialize a "charmf_shc" structure up to degree "nmax".  All
     * coefficients will be initialized to zero and both the scaling constant
     * "shcs->mu" and the radius of the reference sphere "shcs->r" will be set
     * to "1.0".  These values will later be overwritten when reading the input
     * file.  In case of failure, returned is "NULL", so always check the
     * returned pointer for the "NULL" value! */
    charmf_shc *shcs = charmf_shc_calloc(nmax, 1.0, 1.0);
    if (shcs == NULL)
    {
        fprintf(stderr, "Failed to initialize the \"charmf_shc\" structure");
        exit(CHARM_FAILURE);
    }


    /* Many CHarm functions take a "charmf_err" structure as their last
     * parameter.  This is to allow the called function to report an error (if
     * encountered) back to the caller, including some useful information, such
     * as the error message, error code, function name and so on.  The error
     * structure is defined in the "charmf_err" module and can be initialized
     * like this: */
    charmf_err *err = charmf_err_init();
    if (err == NULL)
    {
        fprintf(stderr, "Failed to initialize the \"charmf_err\" "
                        "structure.\n");
        exit(CHARM_FAILURE);
    }


    /* Now let's read spherical harmonic coefficients, the scaling constant and
     * the radius of the reference sphere from the input text file. */
    charmf_shc_read_gfc(shcs_in_file, nmax, shcs, err);


    /* At this point, "shcs" should store the loaded spherical harmonic
     * coefficients.  However, before we continue in our program, we should
     * check whether the previous function call was successful or not.  To this
     * end, we shall look into the "err" structure by calling an error handler.
     * The first parameter of the error handler is the error structure itself,
     * "err".  The other parameter says that if there is indeed an error
     * message in "err", the program prints details on the error and
     * subsequently terminates.  If you want to print the error but do not want
     * to terminate your program, replace "1" by "0" (further details in the
     * "charmf_err" module).  After you call the error handler, you can reuse
     * the same "err" structure in your program, since the function resets the
     * error structure to the default empty values.
     *
     * It is not absolutely necessary to call the error handler, but it is
     * really really really recommended to always do so. */
    charmf_err_handler(err, 1);
    /* --------------------------------------------------------------------- */


    /* Now print some more or less randomly chosen spherical harmonic
     * coefficients */
    /* --------------------------------------------------------------------- */
    /* Let's start with zonal coefficients */
    unsigned long n = 9;  /* Harmonic degree */
    unsigned long m = 0;  /* Harmonic order */
    printf("C(%3lu,%3lu) = %0.7e\n", n, m, shcs->c[m][n - m]);
    printf("S(%3lu,%3lu) = %0.7e\n", n, m, shcs->s[m][n - m]);


    /* Now some tesseral coefficients */
    n = 9;
    m = 4;
    printf("C(%3lu,%3lu) = %0.7e\n", n, m, shcs->c[m][n - m]);
    printf("S(%3lu,%3lu) = %0.7e\n", n, m, shcs->s[m][n - m]);


    /* And finally some sectorial coefficients */
    n = 9;
    m = 9;
    printf("C(%3lu,%3lu) = %0.7e\n", n, m, shcs->c[m][n - m]);
    printf("S(%3lu,%3lu) = %0.7e\n", n, m, shcs->s[m][n - m]);
    /* --------------------------------------------------------------------- */


    /* Now let's save the coefficients to a binary file and then read them back
     * to another structure for spherical harmonic cofficients.  Just for the
     * fun...  */
    /* --------------------------------------------------------------------- */
    /* Write the coefficients. */
    charmf_shc_write_bin(shcs, nmax, shcs_out_file, err);


    /* Again, we should call the error handler to see whether the previous call
     * was successful or not. */
    charmf_err_handler(err, 1);


    /* And now read back the coefficients to a new "charmf_shc" structure
     * called "shcs2" */
    charmf_shc *shcs2 = charmf_shc_calloc(nmax, 1.0, 1.0);
    if (shcs2 == NULL)
    {
        fprintf(stderr, "Failed to initialize the \"charmf_shc\" structure.");
        exit(CHARM_FAILURE);
    }
    charmf_shc_read_bin(shcs_out_file, nmax, shcs2, err);
    charmf_err_handler(err, 1);
    /* --------------------------------------------------------------------- */


    /* Compute degree variances from the loaded coefficients */
    /* --------------------------------------------------------------------- */
    /* Compute degree variances of the input signal */
    float *dv = (float *)malloc((nmax + 1) * sizeof(float));
    if (dv == NULL)
    {
        fprintf(stderr, "malloc failure.\n");
        exit(CHARM_FAILURE);
    }


    charmf_shc_dv(shcs, nmax, dv, err);
    charmf_err_handler(err, 1);


    /* Print some degree variances */
    n = 0;
    printf("\n\n");
    printf("Degree variance for harmonic degree %lu = %0.7e\n", n, dv[n]);
    n = 4;
    printf("Degree variance for harmonic degree %lu = %0.7e\n", n, dv[n]);
    n = 10;
    printf("Degree variance for harmonic degree %lu = %0.7e\n", n, dv[n]);
    /* --------------------------------------------------------------------- */


    /* Now check whether "shcs" and "shcs2" contain the same coefficients by
     * computing difference degree variances */
    /* --------------------------------------------------------------------- */
    float *ddv = (float *)malloc((nmax + 1) * sizeof(float));
    if (ddv == NULL)
    {
        fprintf(stderr, "malloc failure.\n");
        exit(CHARM_FAILURE);
    }
    charmf_shc_ddv(shcs, shcs2, nmax, ddv, err);
    charmf_err_handler(err, 1);


    /* Print some difference degree variances */
    n = 0;
    printf("\n\n");
    printf("Difference degree variance for harmonic degree %lu = %0.7e\n", n,
                                                                       ddv[n]);
    n = 4;
    printf("Difference degree variance for harmonic degree %lu = %0.7e\n", n,
                                                                       ddv[n]);
    n = 10;
    printf("Difference degree variance for harmonic degree %lu = %0.7e\n", n,
                                                                       ddv[n]);
    /* --------------------------------------------------------------------- */


    /* --------------------------------------------------------------------- */
    /* Remember: always free the memory associated with the CHarm structures
     * via the special CHarm "charmf_*_free*" functions.  Calling the usual
     * "free" function will not deallocate the memory properly and will lead to
     * memory leaks. */
    charmf_err_free(err);
    charmf_shc_free(shcs);
    charmf_shc_free(shcs2);
    free(dv), free(ddv);
    /* --------------------------------------------------------------------- */


    printf("\nGreat, all done!\n");


    exit(CHARM_SUCCESS);
    /* ===================================================================== */
}
