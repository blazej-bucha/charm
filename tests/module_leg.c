/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../src/prec.h"
#include "check_leg_pnmj_coeffs.h"
#include "check_leg_pnmj_alloc.h"
#include "check_func.h"
#include "check_outcome.h"
/* ------------------------------------------------------------------------- */






/* Tests the "leg" module */
long int module_leg(void)
{
    long int e    = 0;
    long int esum = 0;


    check_func("leg_pnmj_malloc");
    e = check_leg_pnmj_alloc(CHARM(leg_pnmj_malloc));
    check_outcome(e);
    esum += e;


    check_func("leg_pnmj_calloc");
    e = check_leg_pnmj_alloc(CHARM(leg_pnmj_calloc));
    check_outcome(e);
    esum += e;


    check_func("leg_pnmj_coeffs");
    e = check_leg_pnmj_coeffs();
    check_outcome(e);
    esum += e;


    return esum;
}

