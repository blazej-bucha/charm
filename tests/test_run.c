/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../src/prec.h"
#include "testing_msg.h"
#include "test_suite_start.h"
#include "test_suite_end.h"
#include "module_shs.h"
#include "module_shc.h"
#include "module_crd.h"
#include "module_sha.h"
#include "module_leg.h"
#include "module_integ.h"
#include "misc.h"
/* ------------------------------------------------------------------------- */






/* Program to test CHarm.
 *
 * If compiled with the symbolic constant "GENREF" defined, instead of a check
 * run, the program generates new reference data instead.  The data files are
 * saved to "FOLDER" from "parameters.h". */
int main(void)
{
    test_suite_start();


    TESTING_MSG("shc");
    long int e = module_shc();


    TESTING_MSG("crd");
    e += module_crd();


    TESTING_MSG("shs");
    e += module_shs();


    TESTING_MSG("sha");
    e += module_sha();


    TESTING_MSG("leg");
    e += module_leg();


    TESTING_MSG("integ");
    e += module_integ();


    printf("Testing miscellaneous functions and macros "
           "(may not be a part of the API)...\n");
    e += misc();


    printf("\n");


    test_suite_end(e);


    return CHARM_SUCCESS;
}

