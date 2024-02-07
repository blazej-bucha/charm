/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "check_integ_pn1m1pn2m2.h"
#include "check_integ_yi1n1m1yi2n2m2.h"
#include "check_func.h"
#include "check_outcome.h"
#include "module_integ.h"
/* ------------------------------------------------------------------------- */






/* Tests the "integ" module */
long int module_integ(void)
{
    long int e    = 0;
    long int esum = 0;


    check_func("integ_pn1m1pn2m2");
    e = check_integ_pn1m1pn2m2();
    check_outcome(e);
    esum += e;


    check_func("integ_yi1n1m1yi2n2m2");
    e = check_integ_yi1n1m1yi2n2m2();
    check_outcome(e);
    esum += e;


    return esum;
}

