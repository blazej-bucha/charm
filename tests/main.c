/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../src/prec.h"
#include "module_shs.h"
#include "module_shc.h"
#include "module_crd.h"
#include "module_sha.h"
#include "module_leg.h"
#include "module_integ.h"
#include "misc.h"
/* ------------------------------------------------------------------------- */






/* ------------------------------------------------------------------------- */
#undef TESTING_MSG
#define TESTING_MSG(x) printf("Testing the " x " module...\n");
/* ------------------------------------------------------------------------- */






/* Program to test CHarm.
 *
 * If compiled with the symbolic constant "GENREF" defined, instead of a check
 * run, the program generates new reference data instead.  The data files are
 * saved to "FOLDER" from "parameters.h". */
int main(void)
{
    printf("\n==================================\n");
    printf("\nStart of the test suite\n\n");


    /* Print version number, etc. of the compiled CHarm library */
    printf("Printing info on CHarm...\n");
    printf("..................................\n");
    CHARM(misc_print_version());
    printf("..................................\n\n");


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
    if (e)
        printf("%ld unexpected %s in total\n",
               e, e == 1 ? "result" : "results");
    else
        printf("All tests passed\n");
    printf("\n");
    printf("End of the test suite\n");
    printf("\n");
    printf("==================================\n\n");


    return CHARM_SUCCESS;
}

