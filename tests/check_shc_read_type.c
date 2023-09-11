/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../src/prec.h"
#include "cmp_vals.h"
#include "cmp_arrays.h"
#include "parameters.h"
/* ------------------------------------------------------------------------- */






/* This file checks reading spherical harmonic coefficients.
 *
 *  * "shc_read_tbl", "shc_read_dov", "shc_read_mtx" -- These functions are
 *  used to load the coefficients, the value of which are compared with respect
 *  to coefficients obtained from "shc_read_gfc".
 *
 * * "shc_read_bin" -- Must be called only after "SHCS_OUT_PATH_POT_BIN" was
 * created by "module_shc".  */
long int check_shc_read_type(unsigned long (*shc_read_type)(const char *,
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
    CHARM(shc) *shcs_ref = CHARM(shc_calloc)(SHCS_NMAX_POT, PREC(1.0),
                                             PREC(1.0));
    if (shcs_ref == NULL)
    {
        fprintf(stderr, "Failed to initlize a \"shc\" structure");
        exit(CHARM_FAILURE);
    }
    CHARM(shc_read_gfc)(SHCS_IN_PATH_POT_GFC, SHCS_NMAX_POT, NULL, shcs_ref,
                        err);
    CHARM(err_handler)(err, 1);


    /* Read the coefficients using "tbl", "mtx", "dov" and "bin" routines */
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
    else if (shc_read_type == CHARM(shc_read_dov))
        shc_read_type(SHCS_IN_PATH_POT_DOV, SHCS_NMAX_POT, shcs_type, err);
    CHARM(err_handler)(err, 1);


    e += (shcs_ref->nmax != shcs_type->nmax) ? 1 : 0;
    e += cmp_vals(shcs_ref->mu, shcs_type->mu,
                  PREC(10.0) * CHARM(glob_threshold));
    e += cmp_vals(shcs_ref->r, shcs_type->r,
                  PREC(10.0) * CHARM(glob_threshold));
    e += cmp_arrays(shcs_ref->c[0], shcs_type->c[0], shcs_type->nc,
                    PREC(10.0) * CHARM(glob_threshold));
    e += cmp_arrays(shcs_ref->s[0], shcs_type->s[0], shcs_type->ns,
                    PREC(10.0) * CHARM(glob_threshold));


    CHARM(shc_free)(shcs_type);


    CHARM(err_free)(err);
    CHARM(shc_free)(shcs_ref);


    return e;
}
