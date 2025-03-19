/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../src/prec.h"
#include "parameters.h"
#include "check_struct.h"
#include "cell_touch_array_elements.h"
#include "check_crd_cell_alloc.h"
/* ------------------------------------------------------------------------- */






/* Total number of supported values of "charm_cell->type" */
#undef CELL_TYPES
#define CELL_TYPES (2)






long int check_crd_cell_alloc(CHARM(cell) *(*crd_cell_alloc)(int,
                                                             size_t,
                                                             size_t))
{
    size_t nlat = 5;
    size_t nlon = 10;
    long int e  = 0;


    char func[NSTR_SHORT];
    if (crd_cell_alloc == CHARM(crd_cell_malloc))
        snprintf(func, NSTR_SHORT, "crd_cell_malloc");
    else if (crd_cell_alloc == CHARM(crd_cell_calloc))
        snprintf(func, NSTR_SHORT, "crd_cell_calloc");


    int type;
    int types[CELL_TYPES] = {CHARM_CRD_CELL_SCATTERED,
                             CHARM_CRD_CELL_GRID};


    char func_call_str[NSTR_LONG];
    CHARM(cell) *cell;


    /* --------------------------------------------------------------------- */
    for (int t = 0; t < CELL_TYPES; t++)
    {
        type = types[t];


        for (size_t i = 0; i < nlat; i++)
        {
            for (size_t j = 0; j < nlon; j++)
            {
                cell = crd_cell_alloc(type, i, j);
                snprintf(func_call_str, NSTR_LONG, "%s(%d, %zu, %zu)",
                         func, type, i, j);


                /* For scattered points, "nlat" must be equal to "nlon" */
                if (type == CHARM_CRD_CELL_SCATTERED)
                {
                    if (i != j)
                    {
                        e += check_struct_ptr(cell, NULL, NEQ, INVALID,
                                              func_call_str,
                                              "didn't return NULL pointer");
                        continue;
                    }
                }


                /* In this case, a valid pointer should be returned */
                e += check_struct_ptr(cell, NULL, EQ, VALID, func_call_str,
                                      "returned NULL pointer");


                cell_touch_array_elements(cell);
                CHARM(crd_cell_free)(cell);
            }
        }
    }
    /* --------------------------------------------------------------------- */


    /* Check that invalid value of "type" causes the allocation function to
     * return NULL pointer */
    /* --------------------------------------------------------------------- */
    type = 9999;
    cell = crd_cell_alloc(type, nlat, nlon);
    snprintf(func_call_str, NSTR_LONG, "%s(%d, %zu, %zu)", func, type, nlat,
             nlon);


    e += check_struct_ptr(cell, NULL, NEQ, INVALID, func_call_str,
                          "didn't return NULL pointer");


    CHARM(crd_cell_free)(cell);
    /* --------------------------------------------------------------------- */


    /* Check that the members of "charm_cell" are properly set */
    /* --------------------------------------------------------------------- */
    type = CHARM_CRD_CELL_GRID;
    cell = crd_cell_alloc(type, nlat, nlon);
    snprintf(func_call_str, NSTR_LONG, "%s(%d, %zu, %zu)", func, type, nlat,
             nlon);


    e += check_struct_int(cell->type, CHARM_CRD_CELL_GRID, NEQ, VALID,
                          func_call_str, "returned wrong value of \"type\"");


    e += check_struct_size_t(cell->nlat, nlat, NEQ, VALID, func_call_str,
                             "returned wrong value of \"nlat\"");


    e += check_struct_size_t(cell->nlon, nlon, NEQ, VALID, func_call_str,
                             "returned wrong value of \"nlon\"");


    e += check_struct_size_t(cell->ncell, nlat * nlon, NEQ, VALID,
                             func_call_str,
                             "returned wrong value of \"ncell\"");


    e += check_struct_ptr(cell->latmin, NULL, EQ, VALID, func_call_str,
                          "returned NULL pointer for \"latmin\"");


    e += check_struct_ptr(cell->latmax, NULL, EQ, VALID, func_call_str,
                          "returned NULL pointer for \"latmax\"");


    e += check_struct_ptr(cell->lonmin, NULL, EQ, VALID, func_call_str,
                          "returned NULL pointer for \"lonmin\"");


    e += check_struct_ptr(cell->lonmax, NULL, EQ, VALID, func_call_str,
                          "returned NULL pointer for \"lonmax\"");


    e += check_struct_ptr(cell->r, NULL, EQ, VALID, func_call_str,
                          "returned NULL pointer for \"r\"");


    e += check_struct__Bool(cell->owner, 1, NEQ, VALID, func_call_str,
                            "returned wrong value of \"owner\"");


    CHARM(crd_cell_free)(cell);
    /* --------------------------------------------------------------------- */


    return e;
}

