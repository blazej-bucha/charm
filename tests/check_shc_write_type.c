/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../src/prec.h"
#include "cmp_arrays.h"
#include "parameters.h"
#include "check_shc_read_type.h"
/* ------------------------------------------------------------------------- */






/* Tests writing spherical harmonic coefficients to data files.
 *
 * This function should be called only after "shc_read_gfc" was already
 * checked.
 *
 * As the CHarm's write functions have different APIs, we use integer
 * "shc_write_type" to define the write case we want to test.  "shc_write_type"
 * takes the following values:
 *
 * * "0" -- "shc_write_bin",
 *
 * * "1" -- "shc_write_mtx",
 *
 * * "2" -- "shc_write_tbl" with "CHARM_SHC_WRITE_TBL_N", and
 *
 * * "3" -- "shc_write_tbl" with "CHARM_SHC_WRITE_TBL_M".
 *
 * */
long int check_shc_write_type(int shc_write_type)
{
    CHARM(err) *err = CHARM(err_init)();
    if (err == NULL)
    {
        fprintf(stderr, "Failed to initialize an \"err\" structure.\n");
        exit(CHARM_FAILURE);
    }


    /* Read the reference spherical harmonic coefficients from "gfc" file */
    CHARM(shc) *shcs_gfc = CHARM(shc_calloc)(SHCS_NMAX_POT, PREC(1.0),
                                             PREC(1.0));
    if (shcs_gfc == NULL)
    {
        fprintf(stderr, "Failed to initlize a \"shc\" structure");
        exit(CHARM_FAILURE);
    }
    CHARM(shc_read_mtx)(SHCS_IN_PATH_POT_MTX, SHCS_NMAX_POT, shcs_gfc, err);
    CHARM(err_handler)(err, 1);


    /* Write the reference coefficients to data files */
    if (shc_write_type == 0)
        CHARM(shc_write_bin)(shcs_gfc, SHCS_NMAX_POT, SHCS_OUT_PATH_POT_BIN,
                             err);
    else if (shc_write_type == 1)
        CHARM(shc_write_mtx)(shcs_gfc, SHCS_NMAX_POT, FORMAT,
                             SHCS_OUT_PATH_POT_MTX, err);
    else if (shc_write_type == 2)
        CHARM(shc_write_tbl)(shcs_gfc, SHCS_NMAX_POT, FORMAT,
                             CHARM_SHC_WRITE_TBL_N, SHCS_OUT_PATH_POT_TBL_N,
                             err);
    else if (shc_write_type == 3)
        CHARM(shc_write_tbl)(shcs_gfc, SHCS_NMAX_POT, FORMAT,
                             CHARM_SHC_WRITE_TBL_M, SHCS_OUT_PATH_POT_TBL_M,
                             err);
    CHARM(err_handler)(err, 1);


    /* Read the spherical harmonic coefficients from the data files */
    CHARM(shc) *shcs_type = CHARM(shc_calloc)(SHCS_NMAX_POT, PREC(1.0),
                                               PREC(1.0));
    if (shcs_type == NULL)
    {
        fprintf(stderr, "Failed to initlize a \"shc\" structure");
        exit(CHARM_FAILURE);
    }
    if (shc_write_type == 0)
        CHARM(shc_read_bin)(SHCS_OUT_PATH_POT_BIN, SHCS_NMAX_POT, shcs_type,
                            err);
    else if (shc_write_type == 1)
        CHARM(shc_read_mtx)(SHCS_OUT_PATH_POT_MTX, SHCS_NMAX_POT, shcs_type,
                            err);
    else if (shc_write_type == 2)
        CHARM(shc_read_tbl)(SHCS_OUT_PATH_POT_TBL_N, SHCS_NMAX_POT, shcs_type,
                            err);
    else if (shc_write_type == 3)
        CHARM(shc_read_tbl)(SHCS_OUT_PATH_POT_TBL_M, SHCS_NMAX_POT, shcs_type,
                            err);
    CHARM(err_handler)(err, 1);


    /* Compare the original input coefficients and those obtained by write and
     * read functions */
    int long e = 0;
    e += cmp_arrays(shcs_gfc->c[0], shcs_type->c[0], shcs_type->nc,
                    PREC(10.0) * CHARM(glob_threshold));
    e += cmp_arrays(shcs_gfc->s[0], shcs_type->s[0], shcs_type->ns,
                    PREC(10.0) * CHARM(glob_threshold));


    CHARM(shc_free)(shcs_gfc);
    CHARM(shc_free)(shcs_type);
    CHARM(err_free)(err);


    return e;
}
