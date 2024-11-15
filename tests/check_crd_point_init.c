/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../src/prec.h"
#include "parameters.h"
#include "check_struct.h"
#include "point_touch_array_elements.h"
#include "check_crd_point_init.h"
/* ------------------------------------------------------------------------- */






/* Checks "crd_point_init".  Assumes that "crd_point_calloc" and
 * "crd_point_malloc" have already been tested, as it checks only the features
 * specifically related to "crd_point_init" */
long int check_crd_point_init(void)
{
    size_t nlat           = 5;
    size_t nlon           = 10;
    long int e            = 0;
    char func[NSTR_SHORT] = "crd_point_init";
    char func_call_str[NSTR_LONG];


    REAL *lat = (REAL *)malloc(nlat * sizeof(REAL));
    if (lat == NULL)
    {
        fprintf(stderr, CHARM_ERR_MALLOC_FAILURE"\n");
        exit(CHARM_FAILURE);
    }


    REAL *lon = (REAL *)malloc(nlon * sizeof(REAL));
    if (lon == NULL)
    {
        fprintf(stderr, CHARM_ERR_MALLOC_FAILURE"\n");
        exit(CHARM_FAILURE);
    }


    REAL *r = (REAL *)malloc(nlat * sizeof(REAL));
    if (r == NULL)
    {
        fprintf(stderr, CHARM_ERR_MALLOC_FAILURE"\n");
        exit(CHARM_FAILURE);
    }


    CHARM(point) *pnt = CHARM(crd_point_init)(CHARM_CRD_POINT_GRID, nlat, nlon,
                                              lat, lon, r);
    sprintf(func_call_str, "of %s", func);


    e += check_struct_ptr(pnt, NULL, EQ, VALID, func_call_str,
                          "returned NULL pointer");


    e += check_struct_ptr(pnt->lat, lat, NEQ, VALID, func_call_str,
                          "returned wrong value of \"lat\"");


    e += check_struct_ptr(pnt->lon, lon, NEQ, VALID, func_call_str,
                          "returned wrong value of \"lon\"");


    e += check_struct_ptr(pnt->r, r, NEQ, VALID, func_call_str,
                          "returned wrong value of \"r\"");


    e += check_struct__Bool(pnt->owner, 0, NEQ, VALID, func_call_str,
                            "returned wrong value of \"owner\"");


    e += check_struct__Bool(pnt->distributed, 0, NEQ, VALID, func_call_str,
                            "returned wrong value of \"distributed\"");


    point_touch_array_elements(pnt);
    CHARM(crd_point_free)(pnt);
    free(lat);
    free(lon);
    free(r);


    return e;
}
