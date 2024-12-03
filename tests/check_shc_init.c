/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../src/prec.h"
#include "parameters.h"
#include "shc_touch_array_elements.h"
#include "check_struct.h"
#include "check_shc_init.h"
/* ------------------------------------------------------------------------- */






/* Checks "shc_init".  Assumes that "shc_calloc" and "shc_malloc" have already
 * been tested, as it checks only the features specifically related to
 * "shc_init" */
long int check_shc_init(void)
{
    unsigned long nmax    = SHCS_NMAX_POT;
    REAL mu               = PREC(1.1);
    REAL r                = PREC(2.2);
    size_t ncs            = ((nmax + 2) * (nmax + 1)) / 2;
    long int e            = 0;
    char func[NSTR_SHORT] = "shc_init";
    char func_call_str[NSTR_LONG];


    REAL *c = (REAL *)malloc(ncs * sizeof(REAL));
    if (c == NULL)
    {
        fprintf(stderr, CHARM_ERR_MALLOC_FAILURE"\n");
        exit(CHARM_FAILURE);
    }


    REAL *s = (REAL *)malloc(ncs * sizeof(REAL));
    if (s == NULL)
    {
        fprintf(stderr, CHARM_ERR_MALLOC_FAILURE"\n");
        exit(CHARM_FAILURE);
    }


    CHARM(shc) *shcs = CHARM(shc_init)(NMAX, mu, r, c, s);
    sprintf(func_call_str, "of %s", func);


    e += check_struct_ptr(shcs, NULL, EQ, VALID, func_call_str,
                          "returned NULL pointer");


    e += check_struct_ptr(shcs->c[0], c, NEQ, VALID, func_call_str,
                          "returned wrong value of \"c\"");


    e += check_struct_ptr(shcs->s[0], s, NEQ, VALID, func_call_str,
                          "returned wrong value of \"s\"");


    e += check_struct__Bool(shcs->owner, 0, NEQ, VALID, func_call_str,
                            "returned wrong value of \"owner\"");


    e += check_struct__Bool(shcs->distributed, 0, NEQ, VALID, func_call_str,
                            "returned wrong value of \"distributed\"");


    shc_touch_array_elements(shcs);
    CHARM(shc_free)(shcs);
    free(c);
    free(s);


    return e;
}
