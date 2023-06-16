/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "check_func.h"
#include "check_outcome.h"
#include "check_shs_point.h"
#include "check_shs_cell.h"
#include "check_shs_cell_isurf.h"
/* ------------------------------------------------------------------------- */






/* Tests the "shs" module */
long int module_shs(void)
{
    long int e    = 0;
    long int esum = 0;


    check_func("shs_point");
    e = check_shs_point();
    check_outcome(e);
    esum += e;


    check_func("shs_cell");
    e = check_shs_cell();
    check_outcome(e);
    esum += e;


    check_func("shs_cell_isurf");
    e = check_shs_cell_isurf();
    check_outcome(e);
    esum += e;


    return esum;
}

