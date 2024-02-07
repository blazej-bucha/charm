/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../src/prec.h"
#include "parameters.h"
#include "check_struct.h"
#include "check_crd_cell_init.h"
/* ------------------------------------------------------------------------- */






/* Checks "crd_cell_init".  Assumes that "crd_cell_calloc" and
 * "crd_cell_malloc" have already been tested, as it checks only the features
 * specifically related to "crd_cell_init" */
long int check_crd_cell_init(void)
{
    size_t nlat           = 5;
    size_t nlon           = 10;
    long int e            = 0;
    char func[NSTR_SHORT] = "crd_cell_init";
    char func_call_str[NSTR_LONG];


    REAL *latmin = (REAL *)malloc(nlat * sizeof(REAL));
    if (latmin == NULL)
    {
        fprintf(stderr, "malloc failure.\n");
        exit(CHARM_FAILURE);
    }


    REAL *latmax = (REAL *)malloc(nlat * sizeof(REAL));
    if (latmax == NULL)
    {
        fprintf(stderr, "malloc failure.\n");
        exit(CHARM_FAILURE);
    }


    REAL *lonmin = (REAL *)malloc(nlon * sizeof(REAL));
    if (lonmin == NULL)
    {
        fprintf(stderr, "malloc failure.\n");
        exit(CHARM_FAILURE);
    }


    REAL *lonmax = (REAL *)malloc(nlon * sizeof(REAL));
    if (lonmax == NULL)
    {
        fprintf(stderr, "malloc failure.\n");
        exit(CHARM_FAILURE);
    }


    REAL *r = (REAL *)malloc(nlat * sizeof(REAL));
    if (r == NULL)
    {
        fprintf(stderr, "malloc failure.\n");
        exit(CHARM_FAILURE);
    }


    CHARM(cell) *cell = CHARM(crd_cell_init)(CHARM_CRD_CELL_GRID, nlat, nlon,
                                             latmin, latmax, lonmin, lonmax,
                                             r);
    sprintf(func_call_str, "of %s", func);


    e += check_struct_ptr(cell, NULL, EQ, VALID, func_call_str,
                          "returned a NULL pointer");


    e += check_struct_ptr(cell->latmin, latmin, NEQ, VALID, func_call_str,
                          "returned a wrong value of \"latmin\"");


    e += check_struct_ptr(cell->latmax, latmax, NEQ, VALID, func_call_str,
                          "returned a wrong value of \"latmax\"");


    e += check_struct_ptr(cell->lonmin, lonmin, NEQ, VALID, func_call_str,
                          "returned a wrong value of \"lonmin\"");


    e += check_struct_ptr(cell->lonmax, lonmax, NEQ, VALID, func_call_str,
                          "returned a wrong value of \"lonmax\"");


    e += check_struct_ptr(cell->r, r, NEQ, VALID, func_call_str,
                          "returned a wrong value of \"r\"");


    e += check_struct__Bool(cell->owner, 0, NEQ, VALID, func_call_str,
                            "returned a wrong value of \"owner\"");


    CHARM(crd_cell_free)(cell);
    free(latmin);
    free(latmax);
    free(lonmin);
    free(lonmax);
    free(r);


    return e;
}
