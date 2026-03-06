/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../src/prec.h"
#include "parameters.h"
#include "shc_touch_array_elements.h"
#include "cmp_vals.h"
#include "check_struct.h"
#include "error_messages.h"
#include "check_shc_copy.h"
/* ------------------------------------------------------------------------- */






/* Checks "shc_init".  Assumes that "shc_calloc" and "shc_malloc" have already
 * been tested, as it checks only the features specifically related to
 * "shc_init" */
long int check_shc_copy(void)
{
    unsigned long nmax = SHCS_NMAX_POT;
    REAL mu    = PREC(1.1);
    REAL r     = PREC(2.2);
    long int e = 0;


    /* Reference coefficients */
    /* --------------------------------------------------------------------- */
    CHARM(shc) *shcs_ref = CHARM(shc_calloc)(nmax, mu, r);
    if (shcs_ref == NULL)
    {
        fprintf(stderr, "%s", ERR_MSG_SHC);
        exit(CHARM_FAILURE);
    }


    for (unsigned long m = 0; m <= shcs_ref->nmax; m++)
    {
        for (unsigned long n = m; n <= shcs_ref->nmax; n++)
        {
            shcs_ref->c[m][n - m] = (REAL)(n);
            shcs_ref->s[m][n - m] = (REAL)(n - m);
        }
    }
    /* --------------------------------------------------------------------- */


    /* Check "shc_copy" */
    /* --------------------------------------------------------------------- */
    unsigned long nmax_shcs_out = shcs_ref->nmax;
    unsigned long nmin_all[3] = {0, 2, 5};
    unsigned long nmax_all[3] = {SHCS_NMAX_POT, 2, 7};
    CHARM(shc) *shcs;


    REAL cnmref, snmref, cnm, snm;
    for (size_t i = 0; i < 3; i++)
    {
        shcs = CHARM(shc_copy)(shcs_ref, nmin_all[i], nmax_all[i],
                               nmax_shcs_out);
        if (shcs == NULL)
        {
            fprintf(stderr, "%s", ERR_MSG_SHC);
            exit(CHARM_FAILURE);
        }


        for (unsigned long m = 0; m <= shcs_ref->nmax; m++)
        {
            for (unsigned long n = m; n <= shcs_ref->nmax; n++)
            {
                if ((n < nmin_all[i]) || (n > nmax_all[i]))
                    cnmref = snmref = PREC(0.0);
                else
                {
                    cnmref = shcs_ref->c[m][n - m];
                    snmref = shcs_ref->s[m][n - m];
                }


                cnm = shcs->c[m][n - m];
                snm = shcs->s[m][n - m];


                e += cmp_vals_real(cnm, cnmref, CHARM(glob_threshold));
                e += cmp_vals_real(snm, snmref, CHARM(glob_threshold));
            }
        }


        CHARM(shc_free)(shcs);
    }
    /* --------------------------------------------------------------------- */


    /* Check that the members of "shcs" are properly set */
    /* --------------------------------------------------------------------- */
    char func_call_str[NSTR_LONG];
    char func[NSTR_SHORT];
    snprintf(func, NSTR_SHORT, "shc_copy");


    unsigned long nmin = shcs_ref->nmax + 1;
    nmax = shcs_ref->nmax;
    nmax_shcs_out = shcs_ref->nmax;
    shcs = CHARM(shc_copy)(shcs_ref, nmin, nmax, nmax_shcs_out);
    snprintf(func_call_str, NSTR_LONG, "%s(shcs_ref, %lu, %lu, %lu)",
             func, nmin, nmax, nmax_shcs_out);


    e += check_struct_ptr(shcs, NULL, NEQ, INVALID, func_call_str,
                          "didn't returned NULL pointer");


    shc_touch_array_elements(shcs);
    CHARM(shc_free)(shcs);


    nmin = shcs_ref->nmax;
    nmax = shcs_ref->nmax + 1;
    nmax_shcs_out = shcs_ref->nmax;
    shcs = CHARM(shc_copy)(shcs_ref, nmin, nmax, nmax_shcs_out);
    snprintf(func_call_str, NSTR_LONG, "%s(shcs_ref, %lu, %lu, %lu)",
             func, nmin, nmax, nmax_shcs_out);


    e += check_struct_ptr(shcs, NULL, NEQ, INVALID, func_call_str,
                          "didn't returned NULL pointer");


    shc_touch_array_elements(shcs);
    CHARM(shc_free)(shcs);


    nmin = shcs_ref->nmax;
    nmax = shcs_ref->nmax - 1;
    nmax_shcs_out = shcs_ref->nmax;
    shcs = CHARM(shc_copy)(shcs_ref, nmin, nmax, nmax_shcs_out);
    snprintf(func_call_str, NSTR_LONG, "%s(shcs_ref, %lu, %lu, %lu)",
             func, nmin, nmax, nmax_shcs_out);


    e += check_struct_ptr(shcs, NULL, NEQ, INVALID, func_call_str,
                          "didn't returned NULL pointer");


    shc_touch_array_elements(shcs);
    CHARM(shc_free)(shcs);


    nmin = shcs_ref->nmax;
    nmax = shcs_ref->nmax + 1;
    nmax_shcs_out = shcs_ref->nmax;
    shcs = CHARM(shc_copy)(shcs_ref, nmin, nmax, nmax_shcs_out);
    snprintf(func_call_str, NSTR_LONG, "%s(shcs_ref, %lu, %lu, %lu)",
             func, nmin, nmax, nmax_shcs_out);


    e += check_struct_ptr(shcs, NULL, NEQ, INVALID, func_call_str,
                          "didn't returned NULL pointer");


    shc_touch_array_elements(shcs);
    CHARM(shc_free)(shcs);


    nmin = 0;
    nmax = shcs_ref->nmax;
    nmax_shcs_out = shcs_ref->nmax;
    shcs = CHARM(shc_copy)(shcs_ref, nmin, nmax, nmax_shcs_out);
    snprintf(func_call_str, NSTR_LONG, "%s(shcs_ref, %lu, %lu, %lu)",
             func, nmin, nmax, nmax_shcs_out);


    e += check_struct_ulong(shcs->nmax, shcs_ref->nmax, NEQ, VALID,
                            func_call_str, "returned wrong value of \"nmax\"");


    e += check_struct_REAL(shcs->r, shcs_ref->r, NEQ, VALID,
                           func_call_str, "returned wrong value of \"r\"");


    e += check_struct_REAL(shcs->mu, shcs_ref->mu, NEQ, VALID,
                           func_call_str, "returned wrong value of \"mu\"");


    e += check_struct_ptr(shcs->c, NULL, EQ, VALID, func_call_str,
                          "returned NULL pointer for \"c\"");


    e += check_struct_ptr(shcs->c[0], NULL, EQ, VALID, func_call_str,
                          "returned NULL pointer for \"c[0]\"");


    e += check_struct_ptr(shcs->s, NULL, EQ, VALID, func_call_str,
                          "returned NULL pointer for \"s\"");


    e += check_struct_ptr(shcs->s[0], NULL, EQ, VALID, func_call_str,
                          "returned NULL pointer for \"s[0]\"");


    e += check_struct_size_t(shcs->nc, shcs_ref->nc, NEQ, VALID, func_call_str,
                             "returned wrong value of \"nc\"");


    e += check_struct_size_t(shcs->ns, shcs_ref->ns, NEQ, VALID, func_call_str,
                             "returned wrong value of \"ns\"");


    e += check_struct__Bool(shcs->owner, 1, NEQ, VALID, func_call_str,
                            "returned wrong value of \"owner\"");


    e += check_struct__Bool(shcs->distributed, 0, NEQ, VALID, func_call_str,
                            "returned wrong value of \"distributed\"");


    shc_touch_array_elements(shcs);
    /* --------------------------------------------------------------------- */


    /* --------------------------------------------------------------------- */
    CHARM(shc_free)(shcs_ref);
    CHARM(shc_free)(shcs);
    /* --------------------------------------------------------------------- */


    return e;
}
