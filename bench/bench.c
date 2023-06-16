/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "../src/prec.h"
#include "../src/misc/misc_fprintf_real.h"
/* ------------------------------------------------------------------------- */






/* Symbolic constants */
/* ------------------------------------------------------------------------- */
/* Total number of maximum harmonic degrees, for which the program is
 * executed. */
#undef NMAX
#if CHARM_FLOAT
#   define NMAX 12
#elif CHARM_QUAD
#   define NMAX 13
#else
#   define NMAX 15
#endif


#undef FORMAT
#if CHARM_FLOAT
#   define FORMAT "%0.8e"
#elif CHARM_QUAD
#   define FORMAT "%0.34Qe"
#else
#   define FORMAT "%0.17e"
#endif
/* ------------------------------------------------------------------------- */






int main(void)
{
    /* Maximum harmonic degrees to benchmark CHarm */
#if CHARM_FLOAT
    unsigned long nmax_all[NMAX] = {     10,
                                         25,
                                         50,
                                         75,
                                        100,
                                        250,
                                        500,
                                        750,
                                       1000,
                                       2500,
                                       5000,
                                       7200
                                    };
#elif CHARM_QUAD
    unsigned long nmax_all[NMAX] = {     10,
                                         25,
                                         50,
                                         75,
                                        100,
                                        250,
                                        500,
                                        750,
                                       1000,
                                       2500,
                                       5000,
                                       7500,
                                      10000
                                    };
#else
    unsigned long nmax_all[NMAX] = {     10,
                                         25,
                                         50,
                                         75,
                                        100,
                                        250,
                                        500,
                                        750,
                                       1000,
                                       2500,
                                       5000,
                                       7500,
                                      10000,
                                      25000,
                                      50000
                                    };
#endif


    /* Data folder to save the outputs of the benchmark */
    char path[] = "../data/output";


    /* Pointers to file streams */
    FILE *fid_sha_time;
    FILE *fid_sha_acc;
    FILE *fid_shs_time;


    /* Maximum harmonic degree to be gradually taken from "nmax_all" */
    unsigned long nmax;


    /* Pointer for the signal data to be harmonically analyzed */
    REAL *f;


    /* Some variables to evaluate the accuracy */
    REAL rmse, maxe, diff;


    /* Some variables to measure the execution time */
#if HAVE_CLOCK_GETTIME
    struct timespec t1, t2;
    long sec, nsec;
#endif
    double elapsed;


    /* Initialize the error structure */
    CHARM(err) *err = CHARM(err_init)();
    if (err == NULL)
    {
        fprintf(stderr, "Failed to initialize the error structure.\n");
        exit(CHARM_FAILURE);
    }


    /* Open the file stream to save computation time for SHA */
    char file_sha_time[2048];
#undef SHA_TIME
#if CHARM_FLOAT
#   define SHA_TIME "%s/benchf-sha-time.txt"
#elif CHARM_QUAD
#   define SHA_TIME "%s/benchq-sha-time.txt"
#else
#   define SHA_TIME "%s/bench-sha-time.txt"
#endif
    sprintf(file_sha_time, SHA_TIME, path);
    fid_sha_time = fopen(file_sha_time, "w");


    /* Open the file stream to save computation time for SHS */
    char file_shs_time[2048];
#undef SHS_TIME
#if CHARM_FLOAT
#   define SHS_TIME "%s/benchf-shs-time.txt"
#elif CHARM_QUAD
#   define SHS_TIME "%s/benchq-shs-time.txt"
#else
#   define SHS_TIME "%s/bench-shs-time.txt"
#endif
    sprintf(file_shs_time, SHS_TIME, path);
    fid_shs_time = fopen(file_shs_time, "w");


    /* Open the file stream to save difference degree variances */
    char file_sha_acc[2048];
#undef ACC
#if CHARM_FLOAT
#   define ACC "%s/benchf-acc.txt"
#elif CHARM_QUAD
#   define ACC "%s/benchq-acc.txt"
#else
#   define ACC "%s/bench-acc.txt"
#endif
    sprintf(file_sha_acc, ACC, path);
    fid_sha_acc = fopen(file_sha_acc, "w");


    printf("\n\n");


    /* Loop over all "nmax_all" degrees */
    for (unsigned long nnmax = 0; nnmax < NMAX; nnmax++)
    {
        nmax = nmax_all[nnmax];
        printf("Maximum harmonic degree: %ld\n", nmax);


        printf("    Initializing reference harmonic coefficients...\n");
        CHARM(shc) *shcs_ref = CHARM(shc_calloc)(nmax, PREC(1.0), PREC(1.0));
        if (shcs_ref == NULL)
        {
            fprintf(stderr, "Failed to initialize the shc structure.\n");
            exit(CHARM_FAILURE);
        }
        srand(time(NULL));
        for (unsigned long m = 0; m <= nmax; m++)
        {
            for (unsigned long n = m; n <= nmax; n++)
            {
                /* Generate some more or less random spherical harmonic
                 * coefficients from "-1.0" to "1.0" */
                shcs_ref->c[m][n - m] = PREC(-1.0) + (REAL)rand() /
                                        ((REAL)RAND_MAX / PREC(2.0));
                if (m > 0)
                    shcs_ref->s[m][n - m] = PREC(-1.0) + (REAL)rand() /
                                            ((REAL)RAND_MAX / PREC(2.0));
            }
        }


        printf("    Initializing harmonic coefficients for harmonic "
               "analysis...\n");
        CHARM(shc) *shcs = CHARM(shc_calloc)(nmax, PREC(1.0), PREC(1.0));
        if (shcs == NULL)
        {
            fprintf(stderr, "Failed to initialize the shc structure.\n");
            exit(CHARM_FAILURE);
        }


        printf("    Preparing the Gauss--Legendre grid...\n");
        CHARM(point) *grd = CHARM(crd_point_gl)(nmax, PREC(1.0));
        if (grd == NULL)
        {
            fprintf(stderr, "Failed to compute the Gauss--Legendre grid.\n");
            exit(CHARM_FAILURE);
        }


        printf("    Allocating memory for the signal to be synthesized...\n");
        f = (REAL *)malloc(grd->nlat * grd->nlon * sizeof(REAL));
        if (f == NULL)
        {
            fprintf(stderr, "malloc failure.\n");
            exit(CHARM_FAILURE);
        }


        printf("    Performing spherical harmonic synthesis...\n");
#if HAVE_CLOCK_GETTIME
        clock_gettime(CLOCK_REALTIME, &t1);
#endif
        CHARM(shs_point)(grd, shcs_ref, nmax, f, err);
        CHARM(err_handler)(err, 1);
#if HAVE_CLOCK_GETTIME
        clock_gettime(CLOCK_REALTIME, &t2);
        sec     = t2.tv_sec  - t1.tv_sec;
        nsec    = t2.tv_nsec - t1.tv_nsec;
        elapsed = sec + nsec * 1.0e-9;
#else
        elapsed = 0.0;
#endif


        printf("    Saving wall-clock time to %s...\n", file_shs_time);
        fprintf(fid_shs_time, "%ld %0.17e\n", nmax, elapsed);


        printf("    Performing spherical harmonic analysis...\n");
#if HAVE_CLOCK_GETTIME
        clock_gettime(CLOCK_REALTIME, &t1);
#endif
        CHARM(sha_point)(grd, f, nmax, shcs, err);
        CHARM(err_handler)(err, 1);
#if HAVE_CLOCK_GETTIME
        clock_gettime(CLOCK_REALTIME, &t2);
        sec     = t2.tv_sec  - t1.tv_sec;
        nsec    = t2.tv_nsec - t1.tv_nsec;
        elapsed = sec + nsec * 1.0e-9;
#else
        elapsed = 0.0;
#endif


        printf("    Saving wall-clock time to %s...\n", file_sha_time);
        fprintf(fid_sha_time, "%ld %0.17e\n", nmax, elapsed);


        /* ----------------------------------------------------------------- */
        printf("    Evaluating accuracy...\n");


        maxe = PREC(0.0);
        rmse = PREC(0.0);
        for (unsigned long m = 0; m <= nmax;m ++)
        {
            for (unsigned long n = m; n <= nmax; n++)
            {
                diff = shcs->c[m][n - m] - shcs_ref->c[m][n - m];
                if (FABS(diff) > maxe)
                    maxe = diff;
                rmse += diff * diff;


                diff = shcs->s[m][n - m] - shcs_ref->s[m][n - m];
                if (FABS(diff) > maxe)
                    maxe = diff;
                rmse += diff * diff;
            }
        }


        rmse = SQRT(PREC(2.0) * rmse / ((REAL)((nmax + 1) * (nmax + 2))));
        /* ----------------------------------------------------------------- */


        printf("    Saving accuracy to %s...\n", file_sha_acc);
        fprintf(fid_sha_acc, "%ld ", nmax);
        CHARM(misc_fprintf_real)(fid_sha_acc, FORMAT, maxe);
        fprintf(fid_sha_acc, " ");
        CHARM(misc_fprintf_real)(fid_sha_acc, FORMAT, rmse);
        fprintf(fid_sha_acc, "\n");


        printf("    Freeing the heap memory...\n");
        CHARM(shc_free)(shcs);
        CHARM(shc_free)(shcs_ref);
        CHARM(crd_point_free)(grd);
        free(f);


        printf("\n\n");
    }


    fclose(fid_shs_time);
    fclose(fid_sha_time);
    fclose(fid_sha_acc);
    CHARM(err_free)(err);


    printf("Done.\n");
#if !HAVE_CLOCK_GETTIME
    printf("WARNING: The \"clock_gettime\" function was not found "
           "during the compilation of CHarm, so the benchmark "
           "program couldn't measure the computation time.  The only relevant "
           "output is therefore the computational accuracy.\n");
#endif
    return 0;
}
