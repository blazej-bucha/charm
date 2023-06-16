/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "parameters.h"
#include "../src/prec.h"
#include "check_func.h"
#include "check_outcome.h"
#include "check_shc_read_type.h"
#include "check_shc_write_type.h"
#include "check_shc_dav.h"
#include "check_shc_ddav.h"
#include "check_shc_rescale.h"
#include "check_shc_alloc.h"
#include "check_shc_init.h"
/* ------------------------------------------------------------------------- */






long int module_shc(void)
{
    long int e    = 0;
    long int esum = 0;


    check_func("shc_malloc");
    e = check_shc_alloc(CHARM(shc_malloc));
    check_outcome(e);
    esum += e;


    check_func("shc_calloc");
    e = check_shc_alloc(CHARM(shc_calloc));
    check_outcome(e);
    esum += e;


    check_func("shc_init");
    e = check_shc_init();
    check_outcome(e);
    esum += e;


    check_func("shc_read_gfc");
    e = check_shc_read_type(CHARM(shc_read_gfc));
    check_outcome(e);
    esum += e;


    check_func("shc_read_mtx");
    e = check_shc_read_type(CHARM(shc_read_mtx));
    check_outcome(e);
    esum += e;


    check_func("shc_read_tbl");
    e = check_shc_read_type(CHARM(shc_read_tbl));
    check_outcome(e);
    esum += e;


    /* The following test must run before testing "shc_read_bin" */
    check_func("shc_write_bin");
    e = check_shc_write_type(0);
    check_outcome(e);
    esum += e;


    check_func("shc_read_bin");
    e = check_shc_read_type(CHARM(shc_read_bin));
    check_outcome(e);
    esum += e;


    check_func("shc_write_mtx");
    e = check_shc_write_type(1);
    check_outcome(e);
    esum += e;


    check_func("shc_write_tbl");
    e = check_shc_write_type(2);
    e = check_shc_write_type(3);
    check_outcome(e);
    esum += e;


    check_func("shc_da");
    e = check_shc_dav(CHARM(shc_da));
    check_outcome(e);
    esum += e;


    check_func("shc_dv");
    e = check_shc_dav(CHARM(shc_dv));
    check_outcome(e);
    esum += e;


    check_func("shc_dda");
    e = check_shc_ddav(CHARM(shc_dda));
    check_outcome(e);
    esum += e;


    check_func("shc_ddv");
    e = check_shc_ddav(CHARM(shc_ddv));
    check_outcome(e);
    esum += e;


    check_func("shc_rescale");
    e = check_shc_rescale();
    check_outcome(e);
    esum += e;


    return esum;
}
