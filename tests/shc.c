/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "parameters.h"
#include "../src/prec.h"
#include "cmp_arrays.h"
#include "validate.h"
/* ------------------------------------------------------------------------- */






/* Symbolic constants */
/* ------------------------------------------------------------------------- */
#undef NMAX_TOPO
#define NMAX_TOPO 4


#undef NMAX_POT
#define NMAX_POT 10
/* ------------------------------------------------------------------------- */






/* Function to test the "shc" module */
/* ------------------------------------------------------------------------- */
int shc(unsigned long nmax_topo,
        unsigned long nmax_pot,
        char SHCs_in_topo_mtx_file[],
        char SHCs_in_pot_gfc_file[],
        char SHCs_in_pot_tbl_file[],
        char SHCs_out_topo_bin_file[],
        char SHCs_out_topo_mtx_file[],
        char SHCs_out_pot_tbl_n_file[],
        char SHCs_out_pot_tbl_m_file[])
{
    /* --------------------------------------------------------------------- */
    if (nmax_topo != NMAX_TOPO)
    {
        fprintf(stderr, "\"nmax_topo\" has to be \"%d\".\n", NMAX_TOPO);
        exit(CHARM_FAILURE);
    }


    if (nmax_pot != NMAX_POT)
    {
        fprintf(stderr, "\"nmax_pot\" has to be \"%d\".\n", NMAX_POT);
        exit(CHARM_FAILURE);
    }
    /* --------------------------------------------------------------------- */






    /* Initialize an error structure to be used throughout the function */
    /* --------------------------------------------------------------------- */
    CHARM(err) *err = CHARM(err_init)();
    if (err == NULL)
    {
        fprintf(stderr, "Failed to initialize an \"err\" structure.\n");
        exit(CHARM_FAILURE);
    }
    /* --------------------------------------------------------------------- */






    /* Read spherical harmonic coefficients from "SHCs_in_topo_mtx_file" */
    /* --------------------------------------------------------------------- */
    int errnum = 0;
    printf("    Reading coefficients from a \"mtx\" file...\n");


    CHARM(shc) *shcs_topo_mtx = CHARM(shc_init)(nmax_topo, PREC(1.0),
                                                PREC(1.0));
    if (shcs_topo_mtx == NULL)
    {
        fprintf(stderr, "Failed to initialize a \"shc\" structure.\n");
        exit(CHARM_FAILURE);
    }


    FILE *fptr = fopen(SHCs_in_topo_mtx_file, "r");
    if (fptr == NULL)
    {
        fprintf(stderr, "Failed to open the stream for \"%s\".\n",
                SHCs_in_topo_mtx_file);
        exit(CHARM_FAILURE);
    }
    CHARM(shc_read_mtx)(fptr, nmax_topo, shcs_topo_mtx, err);
    if (!CHARM(err_isempty)(err))
    {
        fprintf(stderr, "        Failed.\n");
        errnum += 1;
    }
    CHARM(err_handler)(err, 0);
    fclose(fptr);
    /* --------------------------------------------------------------------- */






    /* Write spherical harmonic coefficients to a text file
     * "SHCs_out_topo_mtx_file" */
    /* --------------------------------------------------------------------- */
    printf("    Writing coefficients to a \"mtx\" file...\n");


    fptr = fopen(SHCs_out_topo_mtx_file, "w");
    if (fptr == NULL)
    {
        fprintf(stderr, "Failed to open the stream for \"%s\".\n",
                SHCs_out_topo_mtx_file);
        exit(CHARM_FAILURE);
    }


#if CHARM_FLOAT
#   define FORMAT "%0.7e"
#elif CHARM_QUAD
#   define FORMAT "%0.34Qe"
#else
#   define FORMAT "%0.16e"
#endif
    CHARM(shc_write_mtx)(shcs_topo_mtx, nmax_topo, FORMAT, fptr, err);
    if (!CHARM(err_isempty)(err))
    {
        fprintf(stderr, "        Failed.\n");
        errnum += 1;
    }
    CHARM(err_handler)(err, 0);


    fclose(fptr);
    /* --------------------------------------------------------------------- */






    /* Write spherical harmonic coefficients to a binary file
     * "SHCs_out_topo_bin_file" */
    /* --------------------------------------------------------------------- */
    printf("    Writing coefficients to a \"bin\" file...\n");


    fptr = fopen(SHCs_out_topo_bin_file, "wb");
    if (fptr == NULL)
    {
        fprintf(stderr, "Failed to open the stream for \"%s\".\n",
                SHCs_out_topo_bin_file);
        exit(CHARM_FAILURE);
    }


    CHARM(shc_write_bin)(shcs_topo_mtx, nmax_topo, fptr, err);
    if (!CHARM(err_isempty)(err))
    {
        fprintf(stderr, "        Failed.\n");
        errnum += 1;
    }


    CHARM(err_handler)(err, 0);


    fclose(fptr);
    /* --------------------------------------------------------------------- */






    /* Read spherical harmonic coefficients from the saved binary file
     * "SHCs_out_topo_bin_file" */
    /* --------------------------------------------------------------------- */
    printf("    Reading coefficients from a \"bin\" file...\n");


    CHARM(shc) *shcs_topo_bin = CHARM(shc_init)(nmax_topo, PREC(1.0),
                                                PREC(1.0));
    if (shcs_topo_bin == NULL)
    {
        fprintf(stderr, "Failed to initialize a \"shc\" structure");
        exit(CHARM_FAILURE);
    }

    fptr = fopen(SHCs_out_topo_bin_file, "r");
    if (fptr == NULL)
    {
        fprintf(stderr, "Failed to open the stream for \"%s\".\n",
                SHCs_out_topo_bin_file);
        exit(CHARM_FAILURE);
    }
    CHARM(shc_read_bin)(fptr, nmax_topo, shcs_topo_bin, err);
    if (!CHARM(err_isempty)(err))
    {
        fprintf(stderr, "        Failed.\n");
        errnum += 1;
    }
    CHARM(err_handler)(err, 0);
    fclose(fptr);
    /* --------------------------------------------------------------------- */






    /* Read spherical harmonic coefficients from the saved text file
     * "SHCs_out_topo_mtx_file" */
    /* --------------------------------------------------------------------- */
    CHARM(shc) *shcs_topo_mtx_out = CHARM(shc_init)(nmax_topo, PREC(1.0),
                                                    PREC(1.0));
    if (shcs_topo_mtx_out == NULL)
    {
        fprintf(stderr, "Failed to initlize a \"shc\" structure");
        exit(CHARM_FAILURE);
    }

    fptr = fopen(SHCs_out_topo_mtx_file, "r");
    if (fptr == NULL)
    {
        fprintf(stderr, "Failed to open the stream for \"%s\".\n",
                SHCs_out_topo_mtx_file);
        exit(CHARM_FAILURE);
    }
    CHARM(shc_read_mtx)(fptr, nmax_topo, shcs_topo_mtx_out, err);
    CHARM(err_handler)(err, 1);
    fclose(fptr);
    /* --------------------------------------------------------------------- */






    /* Compute degree variances and degree amplitudes */
    /* --------------------------------------------------------------------- */
    char file[NSTR] = "";


    printf("    Degree variances and degree amplitudes...\n");
    for (unsigned long nmax2 = 0; nmax2 <= NMAX; nmax2++)
    {
        /* ................................................................. */
        REAL *dv = (REAL *)malloc((nmax2 + 1) * sizeof(REAL));
        if (dv == NULL)
        {
            fprintf(stderr, "Failed to initialize an array of degree "
                            "variances");
            exit(CHARM_FAILURE);
        }
        CHARM(shc_dv)(shcs_topo_mtx, nmax2, dv, err);
        CHARM(err_handler)(err, 1);


        sprintf(file, "%s/shc_nx%lu_dv%s", FOLDER, (unsigned long)NMAX,
                                           FTYPE);
        errnum += validate(file, dv, nmax2 + 1,
                           PREC(10.0) * CHARM(glob_threshold));


        free(dv);
        /* ................................................................. */


        /* ................................................................. */
        REAL *da = (REAL *)malloc((nmax2 + 1) * sizeof(REAL));
        if (da == NULL)
        {
            fprintf(stderr, "Failed to initialize an array of degree "
                            "variances");
            exit(CHARM_FAILURE);
        }
        CHARM(shc_da)(shcs_topo_mtx, nmax2, da, err);
        CHARM(err_handler)(err, 1);


        sprintf(file, "%s/shc_nx%lu_da%s", FOLDER, (unsigned long)NMAX,
                                           FTYPE);
        errnum += validate(file, dv, nmax2 + 1,
                           PREC(10.0) * CHARM(glob_threshold));


        free(da);
        /* ................................................................. */
    }
    /* --------------------------------------------------------------------- */






    /* Compute difference degree variances and difference degree amplitudes to
     * check reading/writing spherical harmonic coefficients */
    /* --------------------------------------------------------------------- */
    /* Reference values (zeros) */
    REAL *ddv_dda_ref = (REAL *)calloc(nmax_topo + 1, sizeof(REAL));
    if (ddv_dda_ref == NULL)
    {
        fprintf(stderr, "Failed to initialize an array of degree "
                        "amplitudes.\n");
        exit(CHARM_FAILURE);
    }


    /* ..................................................................... */
    REAL *ddv = (REAL *)malloc((nmax_topo + 1) * sizeof(REAL));
    if (ddv == NULL)
    {
        fprintf(stderr, "Failed to initialize an array of degree "
                        "variances.\n");
        exit(CHARM_FAILURE);
    }
    CHARM(shc_ddv)(shcs_topo_mtx, shcs_topo_bin, nmax_topo, ddv, err);
    CHARM(err_handler)(err, 1);


    printf("    Difference degree variances...\n");
    errnum += cmp_arrays(ddv, ddv_dda_ref, nmax_topo + 1,
                         PREC(10.0) * CHARM(glob_threshold));
    /* ..................................................................... */


    /* ..................................................................... */
    REAL *dda = (REAL *)malloc((nmax_topo + 1) * sizeof(REAL));
    if (dda == NULL)
    {
        fprintf(stderr, "Failed to initialize an array of degree "
                        "amplitudes\n");
        exit(CHARM_FAILURE);
    }
    CHARM(shc_dda)(shcs_topo_mtx, shcs_topo_mtx_out, nmax_topo, dda, err);
    CHARM(err_handler)(err, 1);


    printf("    Difference degree amplitudes...\n");
    errnum += cmp_arrays(dda, ddv_dda_ref, nmax_topo + 1,
                         PREC(10.0) * CHARM(glob_threshold));
    /* ..................................................................... */
    /* --------------------------------------------------------------------- */






    /* Read spherical harmonic coefficients from "gfc" and "tbl" files */
    /* --------------------------------------------------------------------- */
    /* ..................................................................... */
    printf("    Reading coefficients from a \"gfc\" file...\n");
    CHARM(shc) *shcs_pot_gfc = CHARM(shc_init)(nmax_pot, PREC(1.0),
                                               PREC(1.0));
    if (shcs_pot_gfc == NULL)
    {
        fprintf(stderr, "Failed to initlize a \"shc\" structure");
        exit(CHARM_FAILURE);
    }

    fptr = fopen(SHCs_in_pot_gfc_file, "r");
    if (fptr == NULL)
    {
        fprintf(stderr, "Failed to open the stream for \"%s\".\n",
                SHCs_in_pot_gfc_file);
        exit(CHARM_FAILURE);
    }
    CHARM(shc_read_gfc)(fptr, nmax_pot, shcs_pot_gfc, err);
    CHARM(err_handler)(err, 1);
    fclose(fptr);
    /* ..................................................................... */


    /* ..................................................................... */
    printf("    Reading coefficients from a \"tbl\" file...\n");
    CHARM(shc) *shcs_pot_tbl = CHARM(shc_init)(nmax_pot, PREC(1.0),
                                               PREC(1.0));
    if (shcs_pot_tbl == NULL)
    {
        fprintf(stderr, "Failed to initlize a \"shc\" structure");
        exit(CHARM_FAILURE);
    }

    fptr = fopen(SHCs_in_pot_tbl_file, "r");
    if (fptr == NULL)
    {
        fprintf(stderr, "Failed to open the stream for \"%s\".\n",
                SHCs_in_pot_tbl_file);
        exit(CHARM_FAILURE);
    }
    CHARM(shc_read_tbl)(fptr, nmax_pot, shcs_pot_tbl, err);
    CHARM(err_handler)(err, 1);
    fclose(fptr);
    /* ..................................................................... */


    /* ..................................................................... */
    printf("    Validating the reading of \"gfc\" and \"tbl\" coefficients "
           "files...\n");
    free(dda);
    dda = (REAL *)malloc((nmax_pot + 1) * sizeof(REAL));
    if (dda == NULL)
    {
        fprintf(stderr, "Failed to initialize an array of degree "
                        "amplitudes\n");
        exit(CHARM_FAILURE);
    }
    CHARM(shc_dda)(shcs_pot_gfc, shcs_pot_tbl, nmax_pot, dda, err);
    CHARM(err_handler)(err, 1);


    /* Reference values (zeros) */
    free(ddv_dda_ref);
    ddv_dda_ref = (REAL *)calloc(nmax_pot + 1, sizeof(REAL));
    if (ddv_dda_ref == NULL)
    {
        fprintf(stderr, "Failed to initialize an array of degree "
                        "amplitudes.\n");
        exit(CHARM_FAILURE);
    }
    errnum += cmp_arrays(dda, ddv_dda_ref, nmax_pot + 1,
                         PREC(10.0) * CHARM(glob_threshold));
    /* ..................................................................... */


    /* Writing to "tbl" with "N" ordering */
    /* ..................................................................... */
    printf("    Writing coefficients to \"tbl\" using the \"N\" ordering "
           "scheme...\n");
    fptr = fopen(SHCs_out_pot_tbl_n_file, "w");
    if (fptr == NULL)
    {
        fprintf(stderr, "Failed to open the stream for \"%s\".\n",
                SHCs_out_pot_tbl_n_file);
        exit(CHARM_FAILURE);
    }
    CHARM(shc_write_tbl)(shcs_pot_tbl, nmax_pot, FORMAT,
                         CHARM_SHC_WRITE_TBL_N, fptr, err);
    if (!CHARM(err_isempty)(err))
    {
        fprintf(stderr, "        Failed.\n");
        errnum += 1;
    }
    CHARM(err_handler)(err, 0);
    fclose(fptr);


    printf("    Validating the writing to \"tbl\" using the \"N\" "
           "ordering scheme...\n");
    CHARM(shc) *shcs_pot_tbl2 = CHARM(shc_init)(nmax_pot, PREC(1.0),
                                                PREC(1.0));
    if (shcs_pot_tbl2 == NULL)
    {
        fprintf(stderr, "Failed to initlize a \"shc\" structure");
        exit(CHARM_FAILURE);
    }
    fptr = fopen(SHCs_out_pot_tbl_n_file, "r");
    if (fptr == NULL)
    {
        fprintf(stderr, "Failed to open the stream for \"%s\".\n",
                SHCs_out_pot_tbl_n_file);
        exit(CHARM_FAILURE);
    }
    CHARM(shc_read_tbl)(fptr, nmax_pot, shcs_pot_tbl2, err);
    CHARM(err_handler)(err, 1);
    fclose(fptr);
    CHARM(shc_dda)(shcs_pot_tbl, shcs_pot_tbl2, nmax_pot, dda, err);
    CHARM(err_handler)(err, 1);
    errnum += cmp_arrays(dda, ddv_dda_ref, nmax_pot + 1,
                         PREC(10.0) * CHARM(glob_threshold));
    /* ..................................................................... */


    /* Writing to "tbl" with "M" ordering */
    /* ..................................................................... */
    printf("    Writing coefficients to \"tbl\" using the \"M\" ordering "
           "scheme...\n");
    fptr = fopen(SHCs_out_pot_tbl_m_file, "w");
    if (fptr == NULL)
    {
        fprintf(stderr, "Failed to open the stream for \"%s\".\n",
                SHCs_out_pot_tbl_m_file);
        exit(CHARM_FAILURE);
    }
    CHARM(shc_write_tbl)(shcs_pot_tbl, nmax_pot, FORMAT,
                         CHARM_SHC_WRITE_TBL_M, fptr, err);
    if (!CHARM(err_isempty)(err))
    {
        fprintf(stderr, "        Failed.\n");
        errnum += 1;
    }
    CHARM(err_handler)(err, 0);
    fclose(fptr);


    printf("    Validating the writing to \"tbl\" using the \"M\" "
           "ordering scheme...\n");
    CHARM(shc) *shcs_pot_tbl3 = CHARM(shc_init)(nmax_pot, PREC(1.0),
                                                PREC(1.0));
    if (shcs_pot_tbl3 == NULL)
    {
        fprintf(stderr, "Failed to initlize a \"shc\" structure");
        exit(CHARM_FAILURE);
    }
    fptr = fopen(SHCs_out_pot_tbl_m_file, "r");
    if (fptr == NULL)
    {
        fprintf(stderr, "Failed to open the stream for \"%s\".\n",
                SHCs_out_pot_tbl_m_file);
        exit(CHARM_FAILURE);
    }
    CHARM(shc_read_tbl)(fptr, nmax_pot, shcs_pot_tbl3, err);
    CHARM(err_handler)(err, 1);
    fclose(fptr);
    CHARM(shc_dda)(shcs_pot_tbl, shcs_pot_tbl3, nmax_pot, dda, err);
    CHARM(err_handler)(err, 1);
    errnum += cmp_arrays(dda, ddv_dda_ref, nmax_pot + 1,
                         PREC(10.0) * CHARM(glob_threshold));
    /* ..................................................................... */
    /* --------------------------------------------------------------------- */






    /* --------------------------------------------------------------------- */
    CHARM(err_free)(err);
    printf("    Freeing \"shc\" structures...\n");
    CHARM(shc_free)(shcs_topo_mtx);
    CHARM(shc_free)(shcs_topo_bin);
    CHARM(shc_free)(shcs_topo_mtx_out);
    CHARM(shc_free)(shcs_pot_gfc);
    CHARM(shc_free)(shcs_pot_tbl);
    CHARM(shc_free)(shcs_pot_tbl2);
    CHARM(shc_free)(shcs_pot_tbl3);
    free(ddv);
    free(dda);
    free(ddv_dda_ref);


    return errnum;
    /* --------------------------------------------------------------------- */
}
/* ------------------------------------------------------------------------- */
