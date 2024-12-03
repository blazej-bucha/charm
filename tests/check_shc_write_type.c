/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../src/prec.h"
#include "cmp_vals.h"
#include "cmp_arrays.h"
#include "parameters.h"
#include "error_messages.h"
#include "check_shc_read_type.h"
#include "check_shc_write_type.h"
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
 * * "2" -- "shc_write_tbl" with "CHARM_SHC_WRITE_N",
 *
 * * "3" -- "shc_write_tbl" with "CHARM_SHC_WRITE_M",
 *
 * * "4" -- "shc_write_dov" with "CHARM_SHC_WRITE_N", and
 *
 * * "5" -- "shc_write_dov" with "CHARM_SHC_WRITE_M".
 *
 * */
long int check_shc_write_type(int shc_write_type)
{
    CHARM(err) *err = CHARM(err_init)();
    if (err == NULL)
    {
        fprintf(stderr, ERR_MSG_ERR);
        exit(CHARM_FAILURE);
    }


    /* Read the reference spherical harmonic coefficients from "gfc" file */
    CHARM(shc) *shcs_ref = CHARM(shc_calloc)(SHCS_NMAX_POT, PREC(1.0),
                                             PREC(1.0));
    if (shcs_ref == NULL)
    {
        fprintf(stderr, ERR_MSG_SHC);
        exit(CHARM_FAILURE);
    }
    CHARM(shc_read_gfc)(SHCS_IN_PATH_POT_GFC, SHCS_NMAX_POT, NULL, shcs_ref,
                        err);
    CHARM(err_handler)(err, 1);


    /* Write the reference coefficients to data files */
    if (shc_write_type == 0)
        CHARM(shc_write_bin)(shcs_ref, SHCS_NMAX_POT, SHCS_OUT_PATH_POT_BIN,
                             err);
    else if (shc_write_type == 1)
        CHARM(shc_write_mtx)(shcs_ref, SHCS_NMAX_POT, FORMAT,
                             SHCS_OUT_PATH_POT_MTX, err);
    else if (shc_write_type == 2)
        CHARM(shc_write_tbl)(shcs_ref, SHCS_NMAX_POT, FORMAT,
                             CHARM_SHC_WRITE_N, SHCS_OUT_PATH_POT_TBL_N,
                             err);
    else if (shc_write_type == 3)
        CHARM(shc_write_tbl)(shcs_ref, SHCS_NMAX_POT, FORMAT,
                             CHARM_SHC_WRITE_M, SHCS_OUT_PATH_POT_TBL_M,
                             err);
    else if (shc_write_type == 4)
        CHARM(shc_write_dov)(shcs_ref, SHCS_NMAX_POT, FORMAT,
                             CHARM_SHC_WRITE_N, SHCS_OUT_PATH_POT_DOV_N,
                             err);
    else if (shc_write_type == 5)
        CHARM(shc_write_dov)(shcs_ref, SHCS_NMAX_POT, FORMAT,
                             CHARM_SHC_WRITE_M, SHCS_OUT_PATH_POT_DOV_M,
                             err);
    CHARM(err_handler)(err, 1);


    /* Read the spherical harmonic coefficients from the data files */
    CHARM(shc) *shcs_type = CHARM(shc_calloc)(SHCS_NMAX_POT, PREC(1.0),
                                              PREC(1.0));
    if (shcs_type == NULL)
    {
        fprintf(stderr, ERR_MSG_SHC);
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
    else if (shc_write_type == 4)
        CHARM(shc_read_dov)(SHCS_OUT_PATH_POT_DOV_N, SHCS_NMAX_POT, shcs_type,
                            err);
    else if (shc_write_type == 5)
        CHARM(shc_read_dov)(SHCS_OUT_PATH_POT_DOV_M, SHCS_NMAX_POT, shcs_type,
                            err);
    CHARM(err_handler)(err, 1);


    /* Compare the metadata of the input shcs and those obtained by write and
     * read functions */
    int long e = 0;
    e += (shcs_ref->nmax != shcs_type->nmax) ? 1 : 0;
    e += cmp_vals_real(shcs_ref->mu, shcs_type->mu,
                       PREC(10.0) * CHARM(glob_threshold));
    e += cmp_vals_real(shcs_ref->r, shcs_type->r,
                       PREC(10.0) * CHARM(glob_threshold));



    /* Compare the original input coefficients and those obtained by write and
     * read functions */
    e += cmp_arrays(shcs_ref->c[0], shcs_type->c[0], shcs_type->nc,
                    PREC(10.0) * CHARM(glob_threshold));
    e += cmp_arrays(shcs_ref->s[0], shcs_type->s[0], shcs_type->ns,
                    PREC(10.0) * CHARM(glob_threshold));


    CHARM(shc_free)(shcs_ref);
    CHARM(shc_free)(shcs_type);
    CHARM(err_free)(err);


    return e;
}
