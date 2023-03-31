/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../src/prec.h"
#include "check_integ_pn1m1pn2m2.h"
#include "check_integ_yi1n1m1yi2n2m2.h"
/* ------------------------------------------------------------------------- */






/* Tests the "integ" module */
int integ(void)
{
    /* The "nmax" value should be larger than 7, but not too high so that the
     * checks are fast enough. */
    unsigned long nmax = 8;


    int errnum         = 0;



    printf("    Testing \"integ_pn1m1pn2m2\"...\n");
    errnum += check_integ_pn1m1pn2m2(nmax);


    printf("    Testing \"integ_yi1n1m1yi2n2m2\"...\n");
    errnum += check_integ_yi1n1m1yi2n2m2(nmax);


    return errnum;
}

