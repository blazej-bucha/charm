/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../src/prec.h"
/* ------------------------------------------------------------------------- */






/* ------------------------------------------------------------------------- */
void test_suite_start(void)
{
    printf("\n==================================\n");
    printf("\nStart of the test suite\n\n");


    /* Print version number, etc. of the compiled CHarm library */
    printf("Printing info on CHarm...\n");
    printf("..................................\n");
    CHARM(misc_print_info());
    printf("..................................\n\n");


    return;
}
/* ------------------------------------------------------------------------- */
