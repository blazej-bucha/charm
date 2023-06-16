/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../src/prec.h"
#include "cmp_arrays.h"
#include "parameters.h"
/* ------------------------------------------------------------------------- */






/* This file checks reading spherical harmonic coefficients.
 *
 * * "shc_read_gfc" -- It is only checked whether or not "shc_read_gfc" threw
 * an error and no checks are done as to whether or not the content of the file
 * was loaded properly.  The latter is checked indirectly when testing reading
 * spherical harmonic coefficients from other text formats.
 *
 * * "shc_read_tbl", "shc_read_mtx" -- These functions are used to load the
 * coefficients, the value of which are compared with respect to coefficients
 * obtained from "shc_read_gfc".
 *
 * * "shc_read_bin" -- Must be called only after "SHCS_OUT_PATH_POT_BIN" was
 * created by "module_shc".  */
long int check_shc_read_type(void (*shc_read_type)(const char *,
                                                   unsigned long,
                                                   CHARM(shc) *,
                                                   CHARM(err) *))
{
    CHARM(err) *err = CHARM(err_init)();
    if (err == NULL)
    {
        fprintf(stderr, "Failed to initialize an \"err\" structure.\n");
        exit(CHARM_FAILURE);
    }


    int long e = 0;


    /* Read the reference "gfc" file */
    CHARM(shc) *shcs_gfc = CHARM(shc_calloc)(SHCS_NMAX_POT, PREC(1.0),
                                             PREC(1.0));
    if (shcs_gfc == NULL)
    {
        fprintf(stderr, "Failed to initlize a \"shc\" structure");
        exit(CHARM_FAILURE);
    }
    CHARM(shc_read_gfc)(SHCS_IN_PATH_POT_GFC, SHCS_NMAX_POT, shcs_gfc, err);
    if (shc_read_type == CHARM(shc_read_gfc))
    {
        /* If testing "shc_read_gfc", check only for potential errors in
         * "err". */
        if (!CHARM(err_isempty)(err))
            e += 1;
        CHARM(err_handler)(err, 0);


        goto EXIT;
    }
    else
        CHARM(err_handler)(err, 1);


    /* Read the coefficients using "tbl", "mtx" and "bin" routines */
    CHARM(shc) *shcs_type = CHARM(shc_calloc)(SHCS_NMAX_POT, PREC(1.0),
                                              PREC(1.0));
    if (shcs_type == NULL)
    {
        fprintf(stderr, "Failed to initlize a \"shc\" structure");
        exit(CHARM_FAILURE);
    }
    if (shc_read_type == CHARM(shc_read_tbl))
        shc_read_type(SHCS_IN_PATH_POT_TBL, SHCS_NMAX_POT, shcs_type, err);
    else if (shc_read_type == CHARM(shc_read_mtx))
        shc_read_type(SHCS_IN_PATH_POT_MTX, SHCS_NMAX_POT, shcs_type, err);
    else if (shc_read_type == CHARM(shc_read_bin))
        shc_read_type(SHCS_OUT_PATH_POT_BIN, SHCS_NMAX_POT, shcs_type, err);
    CHARM(err_handler)(err, 1);


    e += cmp_arrays(shcs_gfc->c[0], shcs_type->c[0], shcs_type->nc,
                    PREC(10.0) * CHARM(glob_threshold));
    e += cmp_arrays(shcs_gfc->s[0], shcs_type->s[0], shcs_type->ns,
                    PREC(10.0) * CHARM(glob_threshold));


    CHARM(shc_free)(shcs_type);


EXIT:
    CHARM(err_free)(err);
    CHARM(shc_free)(shcs_gfc);


    return e;
}
