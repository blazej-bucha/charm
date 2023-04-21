/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "check_func.h"
#include "check_outcome.h"
#include "check_sha_point.h"
#include "check_sha_cell.h"
/* ------------------------------------------------------------------------- */






/* Tests the "sha" module */
long int module_sha(void)
{
    long int e    = 0;
    long int esum = 0;


    check_func("sha_point");
    e = check_sha_point();
    check_outcome(e);
    esum += e;


    check_func("sha_cell");
    e = check_sha_cell();
    check_outcome(e);
    esum += e;


    return esum;
}

