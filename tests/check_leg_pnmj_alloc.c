/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../src/prec.h"
#include "../src/leg/leg_pnmj_length.h"
#include "parameters.h"
#include "check_struct.h"
#include "check_leg_pnmj_alloc.h"
/* ------------------------------------------------------------------------- */






long int check_leg_pnmj_alloc(CHARM(pnmj) *(*pnmj_alloc)(unsigned long,
                                                         int))
{
    unsigned long nmax;
    long int e = 0;


    char func[NSTR_SHORT];
    if (pnmj_alloc == CHARM(leg_pnmj_malloc))
        sprintf(func, "leg_pnmj_malloc");
    else if (pnmj_alloc == CHARM(leg_pnmj_calloc))
        sprintf(func, "leg_pnmj_calloc");


    CHARM(pnmj) *pnmj;
    char func_call_str[NSTR_LONG];


    /* Check allocation for various maximum harmonic degrees */
    /* --------------------------------------------------------------------- */
    {
    int ordering[2] = {CHARM_LEG_PMNJ, CHARM_LEG_PMJN};


    for (nmax = 0; nmax <= SHCS_NMAX_POT; nmax++)
    {
        for (unsigned o = 0; o < 2; o++)
        {
            pnmj = pnmj_alloc(nmax, ordering[o]);
            sprintf(func_call_str, "%s(%lu, %i)", func, nmax, ordering[o]);


            e += check_struct_ptr(pnmj, NULL, EQ, VALID, func_call_str,
                                  "returned NULL pointer");


            CHARM(leg_pnmj_free)(pnmj);
        }
    }
    }
    /* --------------------------------------------------------------------- */


    /* Check that unsupported values of "ordering" cause the allocation
     * function to return NULL pointer */
    /* --------------------------------------------------------------------- */
    int ordering = 9999;
    nmax         = SHCS_NMAX_POT;


    pnmj = pnmj_alloc(nmax, ordering);
    sprintf(func_call_str, "%s(%lu, %i)", func, nmax, ordering);


    e += check_struct_ptr(pnmj, NULL, NEQ, INVALID, func_call_str,
                          "didn't return NULL pointer");


    CHARM(leg_pnmj_free)(pnmj);
    /* --------------------------------------------------------------------- */


    /* Check that the members of "pnmj" are properly set */
    /* --------------------------------------------------------------------- */
    ordering = CHARM_LEG_PMNJ;


    pnmj = pnmj_alloc(nmax, ordering);
    sprintf(func_call_str, "%s(%lu, %i)", func, nmax, ordering);


    e += check_struct_ulong(pnmj->nmax, nmax, NEQ, VALID, func_call_str,
                            "returned wrong value of \"nmax\"");


    e += check_struct_int(pnmj->ordering, ordering, NEQ, VALID, func_call_str,
                          "returned wrong value of \"ordering\"");


    size_t npnmj = CHARM(leg_pnmj_length)(nmax);
    e += check_struct_size_t(pnmj->npnmj, npnmj, NEQ, VALID, func_call_str,
                             "returned wrong value of \"npnmj\"");


    e += check_struct_ptr(pnmj->pnmj, NULL, EQ, VALID, func_call_str,
                          "returned NULL pointer for \"pnmj\"");


    e += check_struct_ptr(pnmj->pnmj[0], NULL, EQ, VALID, func_call_str,
                          "returned NULL pointer for \"pnmj[0]\"");


    e += check_struct_ptr(pnmj->pnmj[0][0], NULL, EQ, VALID, func_call_str,
                          "returned NULL pointer for \"pnmj[0][0]\"");


    CHARM(leg_pnmj_free)(pnmj);
    /* --------------------------------------------------------------------- */


    return e;
}
