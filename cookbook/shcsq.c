#include <stdio.h>
#include <stdlib.h>
#include <quadmath.h>
#include <charm/charmq.h>


int main(void)
{
    /* INPUTS */
    /* ===================================================================== */
    /* Define the path to an input "gfc" text file with spherical harmonic
     * coefficients.  For details on the structure of the "gfc" file, see the
     * description of the "charm_shc_read_gfc" function in the "charm_shc"
     * module. */
    char shcs_in_file[] = "../data/input/EGM96-degree10.gfc";


    /* Maximum harmonic degree to initialize, read and write spherical harmonic
     * coefficients */
    unsigned long nmax = 10;


    /* Define the path to an output binary file with spherical harmonic
     * coefficients.  For details on the structure of the output binary file,
     * see the description of the "charmq_shc_write_bin" function in the
     * "charmq_shc" module. */
    char shcs_out_file[] = "../data/output/EGM96-degree10-mtx.shcs";
    /* ===================================================================== */






    /* ===================================================================== */
    /* Read the spherical harmonic coefficients */
    /* --------------------------------------------------------------------- */
    /* Initialize a "charmq_shc" structure up to degree "nmax".  All
     * coefficients will be initialized to zero and both the scaling constant
     * "shcs->mu" and the radius of the reference sphere "shcs->r" will be set
     * to "1.0".  These values will later be overwritten when reading the input
     * file.  In case of failure, returned is "NULL", so always check the
     * returned pointer for the "NULL" value! */
    charmq_shc *shcs = charmq_shc_init(nmax, 1.0, 1.0);
    if (shcs == NULL)
    {
        fprintf(stderr, "Failed to initialize the \"charmq_shc\" structure");
        exit(CHARM_FAILURE);
    }


    /* Many CHarm functions take a "charmq_err" structure as their last
     * parameter.  This is to allow the called function to report an error (if
     * encountered) back to the caller, including some useful information, such
     * as the error message, error code, function name and so on.  The error
     * structure is defined in the "charmq_err" module and can be initialized
     * like this: */
    charmq_err *err = charmq_err_init();
    if (err == NULL)
    {
        fprintf(stderr, "Failed to initialize the \"charmq_err\" "
                        "structure.\n");
        exit(CHARM_FAILURE);
    }


    /* Now let's read spherical harmonic coefficients, the scaling constant and
     * the radius of the reference sphere from the input text file.  We open
     * the input file, check whether the operation was successful, then read
     * the file and, finally, close the stream. */
    FILE *fptr = fopen(shcs_in_file, "r");
    if (fptr == NULL)
    {
        fprintf(stderr, "Failed to open the stream for \"%s\".\n",
                shcs_in_file);
        exit(CHARM_FAILURE);
    }
    charmq_shc_read_gfc(fptr, nmax, shcs, err);
    fclose(fptr);


    /* At this point, "shcs" should store the loaded spherical harmonic
     * coefficients.  However, before we continue in our program, we should
     * check whether the previous function call was successful or not.  To this
     * end, we shall look into the "err" structure by calling an error handler.
     * The first parameter of the error handler is the error structure itself,
     * "err".  The other parameter says that if there is indeed an error
     * message in "err", the program prints details on the error and
     * subsequently terminates.  If you want to print the error but do not want
     * to terminate your program, replace "1" by "0" (further details in the
     * "charmq_err" module).  After you call the error handler, you can reuse
     * the same "err" structure in your program, since the function resets the
     * error structure to the default empty values.
     *
     * It is not absolutely necessary to call the error handler, but it is
     * really really really recommended to always do so. */
    charmq_err_handler(err, 1);
    /* --------------------------------------------------------------------- */


    /* Now print some more or less randomly chosen spherical harmonic
     * coefficients */
    /* --------------------------------------------------------------------- */
    /* Let's start with zonal coefficients */
    unsigned long n = 9;  /* Harmonic degree */
    unsigned long m = 0;  /* Harmonic order */
    printf("C(%3lu,%3lu) = %0.34Qe\n", n, m, shcs->c[m][n - m]);
    printf("S(%3lu,%3lu) = %0.34Qe\n", n, m, shcs->s[m][n - m]);


    /* Now some tesseral coefficients */
    n = 9;
    m = 4;
    printf("C(%3lu,%3lu) = %0.34Qe\n", n, m, shcs->c[m][n - m]);
    printf("S(%3lu,%3lu) = %0.34Qe\n", n, m, shcs->s[m][n - m]);


    /* And finally some sectorial coefficients */
    n = 9;
    m = 9;
    printf("C(%3lu,%3lu) = %0.34Qe\n", n, m, shcs->c[m][n - m]);
    printf("S(%3lu,%3lu) = %0.34Qe\n", n, m, shcs->s[m][n - m]);
    /* --------------------------------------------------------------------- */


    /* Now let's save the coefficients to a binary file and then read them back
     * to another structure for spherical harmonic cofficients.  Just for the
     * fun...  */
    /* --------------------------------------------------------------------- */
    /* Write the coefficients. */
    fptr = fopen(shcs_out_file, "wb");
    if (fptr == NULL)
    {
        fprintf(stderr, "Failed to open the stream for \"%s\".\n",
                shcs_out_file);
        exit(CHARM_FAILURE);
    }
    charmq_shc_write_bin(shcs, nmax, fptr, err);
    fclose(fptr);


    /* Again, we should call the error handler to see whether the previous call
     * was successful or not. */
    charmq_err_handler(err, 1);


    /* And now read back the coefficients to a new "charmq_shc" structure
     * called "shcs2" */
    charmq_shc *shcs2 = charmq_shc_init(nmax, 1.0, 1.0);
    if (shcs2 == NULL)
    {
        fprintf(stderr, "Failed to initialize the \"charmq_shc\" structure.");
        exit(CHARM_FAILURE);
    }
    fptr = fopen(shcs_out_file, "rb");
    if (fptr == NULL)
    {
        fprintf(stderr, "Failed to open the stream for \"%s\".\n",
                shcs_out_file);
        exit(CHARM_FAILURE);
    }
    charmq_shc_read_bin(fptr, nmax, shcs2, err);
    fclose(fptr);
    charmq_err_handler(err, 1);
    /* --------------------------------------------------------------------- */


    /* Compute degree variances from the loaded coefficients */
    /* --------------------------------------------------------------------- */
    /* Compute degree variances of the input signal */
    __float128 *dv = (__float128 *)malloc((nmax + 1) * sizeof(__float128));
    if (dv == NULL)
    {
        fprintf(stderr, "malloc failure.\n");
        exit(CHARM_FAILURE);
    }


    charmq_shc_dv(shcs, nmax, dv, err);
    charmq_err_handler(err, 1);


    /* Print some degree variances */
    n = 0;
    printf("\n\n");
    printf("Degree variance for harmonic degree %lu = %0.34Qe\n", n, dv[n]);
    n = 4;
    printf("Degree variance for harmonic degree %lu = %0.34Qe\n", n, dv[n]);
    n = 10;
    printf("Degree variance for harmonic degree %lu = %0.34Qe\n", n, dv[n]);
    /* --------------------------------------------------------------------- */


    /* Now check whether "shcs" and "shcs2" contain the same coefficients by
     * computing difference degree variances */
    /* --------------------------------------------------------------------- */
    __float128 *ddv = (__float128 *)malloc((nmax + 1) * sizeof(__float128));
    if (ddv == NULL)
    {
        fprintf(stderr, "malloc failure.\n");
        exit(CHARM_FAILURE);
    }
    charmq_shc_ddv(shcs, shcs2, nmax, ddv, err);
    charmq_err_handler(err, 1);


    /* Print some difference degree variances */
    n = 0;
    printf("\n\n");
    printf("Difference degree variance for harmonic degree %lu = %0.34Qe\n",
           n, ddv[n]);
    n = 4;
    printf("Difference degree variance for harmonic degree %lu = %0.34Qe\n",
           n, ddv[n]);
    n = 10;
    printf("Difference degree variance for harmonic degree %lu = %0.34Qe\n",
           n, ddv[n]);
    /* --------------------------------------------------------------------- */


    /* --------------------------------------------------------------------- */
    /* Remember: always free the memory associated with the CHarm structures
     * via the special CHarm "charmq_*_free*" functions.  Calling the usual
     * "free" function will not deallocate the memory properly and will lead to
     * memory leaks. */
    charmq_err_free(err);
    charmq_shc_free(shcs);
    charmq_shc_free(shcs2);
    free(dv), free(ddv);
    /* --------------------------------------------------------------------- */


    printf("\nGreat, all done!\n");


    exit(CHARM_SUCCESS);
    /* ===================================================================== */
}
