/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../src/prec.h"
#include "check_crd_point_quad.h"
#include "check_crd_point_alloc.h"
#include "check_crd_point_init.h"
#include "check_crd_cell_alloc.h"
#include "check_crd_cell_init.h"
#include "check_func.h"
#include "check_outcome.h"
#include "module_crd.h"
/* ------------------------------------------------------------------------- */






/* Tests the "crd" module */
long int module_crd(void)
{
    long int e    = 0;
    long int esum = 0;


    check_func("crd_point_malloc");
    e = check_crd_point_alloc(CHARM(crd_point_malloc));
    check_outcome(e);
    esum += e;


    check_func("crd_point_calloc");
    e = check_crd_point_alloc(CHARM(crd_point_calloc));
    check_outcome(e);
    esum += e;


    check_func("crd_point_init");
    e = check_crd_point_init();
    check_outcome(e);
    esum += e;


    check_func("crd_cell_malloc");
    e = check_crd_cell_alloc(CHARM(crd_cell_malloc));
    check_outcome(e);
    esum += e;


    check_func("crd_cell_calloc");
    e = check_crd_cell_alloc(CHARM(crd_cell_calloc));
    check_outcome(e);
    esum += e;


    check_func("crd_cell_init");
    e = check_crd_cell_init();
    check_outcome(e);
    esum += e;


    check_func("crd_point_gl");
    e = check_crd_point_quad(CHARM(crd_point_gl));
    check_outcome(e);
    esum += e;


    check_func("crd_point_dh1");
    e = check_crd_point_quad(CHARM(crd_point_dh1));
    check_outcome(e);
    esum += e;


    check_func("crd_point_dh2");
    e = check_crd_point_quad(CHARM(crd_point_dh2));
    check_outcome(e);
    esum += e;


    return esum;
}

