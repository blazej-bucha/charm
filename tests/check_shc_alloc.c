/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../src/prec.h"
#include "parameters.h"
#include "shc_touch_array_elements.h"
#include "check_struct.h"
#include "check_shc_alloc.h"
/* ------------------------------------------------------------------------- */






long int check_shc_alloc(CHARM(shc) *(*shc_alloc)(unsigned long,
                                                  REAL,
                                                  REAL))
{
    unsigned long nmax;
    REAL mu    = PREC(1.1);
    REAL r     = PREC(2.2);
    REAL zero  = PREC(0.0);
    long int e = 0;


    char func[NSTR_SHORT];
    if (shc_alloc == CHARM(shc_malloc))
        sprintf(func, "shc_malloc");
    else if (shc_alloc == CHARM(shc_calloc))
        sprintf(func, "shc_calloc");


    char func_call_str[NSTR_LONG];
    CHARM(shc) *shcs;


    /* Check allocation for various maximum harmonic degrees */
    /* --------------------------------------------------------------------- */
    for (nmax = 0; nmax <= SHCS_NMAX_POT; nmax++)
    {
        shcs = shc_alloc(nmax, mu, r);
        sprintf(func_call_str,
                "%s(%lu, " REAL_PRINT_FORMAT ", " REAL_PRINT_FORMAT ")",
                func, nmax, mu, r);


        e += check_struct_ptr(shcs, NULL, EQ, VALID, func_call_str,
                              "returned NULL pointer");


        shc_touch_array_elements(shcs);
        CHARM(shc_free)(shcs);
    }
    /* --------------------------------------------------------------------- */


    /* Check that the NULL pointer is returned if the radius is zero or
     * negative */
    /* --------------------------------------------------------------------- */
    nmax = SHCS_NMAX_POT;


    shcs = shc_alloc(nmax, mu, zero);
    sprintf(func_call_str,
            "%s(%lu, " REAL_PRINT_FORMAT ", " REAL_PRINT_FORMAT ")",
            func, nmax, mu, zero);


    e += check_struct_ptr(shcs, NULL, NEQ, INVALID, func_call_str,
                          "didn't returned NULL pointer");


    shc_touch_array_elements(shcs);
    CHARM(shc_free)(shcs);


    r = -r;


    shcs = shc_alloc(nmax, mu, r);
    sprintf(func_call_str,
            "%s(%lu, " REAL_PRINT_FORMAT ", " REAL_PRINT_FORMAT ")",
            func, nmax, mu, r);


    e += check_struct_ptr(shcs, NULL, NEQ, INVALID, func_call_str,
                          "didn't returned NULL pointer");


    r = -r;


    shc_touch_array_elements(shcs);
    CHARM(shc_free)(shcs);
    /* --------------------------------------------------------------------- */


    /* Check that the members of "shcs" are properly set */
    /* --------------------------------------------------------------------- */
    shcs = shc_alloc(nmax, mu, r);
    sprintf(func_call_str,
            "%s(%lu, " REAL_PRINT_FORMAT ", " REAL_PRINT_FORMAT ")",
            func, nmax, mu, r);


    e += check_struct_ulong(shcs->nmax, nmax, NEQ, VALID, func_call_str,
                            "returned wrong value of \"nmax\"");


    e += check_struct_REAL(shcs->mu, mu, NEQ, VALID, func_call_str,
                           "returned wrong value of \"mu\"");


    e += check_struct_REAL(shcs->r, r, NEQ, VALID, func_call_str,
                           "returned wrong value of \"r\"");


    e += check_struct_ptr(shcs->c, NULL, EQ, VALID, func_call_str,
                          "returned NULL pointer for \"c\"");


    e += check_struct_ptr(shcs->c[0], NULL, EQ, VALID, func_call_str,
                          "returned NULL pointer for \"c[0]\"");


    size_t ncs = ((nmax + 2) * (nmax + 1)) / 2;
    e += check_struct_size_t(shcs->nc, ncs, NEQ, VALID, func_call_str,
                             "returned wrong value of \"nc\"");


    e += check_struct_ptr(shcs->s, NULL, EQ, VALID, func_call_str,
                          "returned NULL pointer for \"s\"");


    e += check_struct_ptr(shcs->s[0], NULL, EQ, VALID, func_call_str,
                          "returned NULL pointer for \"s[0]\"");


    e += check_struct_size_t(shcs->ns, ncs, NEQ, VALID, func_call_str,
                             "returned wrong value of \"ns\"");


    e += check_struct__Bool(shcs->owner, 1, NEQ, VALID, func_call_str,
                            "returned wrong value of \"owner\"");


    e += check_struct__Bool(shcs->distributed, 0, NEQ, VALID, func_call_str,
                            "returned wrong value of \"distributed\"");


    shc_touch_array_elements(shcs);
    CHARM(shc_free)(shcs);
    /* --------------------------------------------------------------------- */


    return e;
}

