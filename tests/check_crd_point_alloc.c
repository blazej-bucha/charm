/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../src/prec.h"
#include "parameters.h"
#include "check_struct.h"
#include "check_crd_point_alloc.h"
/* ------------------------------------------------------------------------- */






/* Total number of supported values of "charm_point->type" */
#undef POINT_TYPES
#define POINT_TYPES (5)






long int check_crd_point_alloc(CHARM(point) *(*crd_point_alloc)(int,
                                                                size_t,
                                                                size_t))
{
    size_t nlat = 5;
    size_t nlon = 10;
    long int e  = 0;


    char func[NSTR_SHORT];
    if (crd_point_alloc == CHARM(crd_point_malloc))
        sprintf(func, "crd_point_malloc");
    else if (crd_point_alloc == CHARM(crd_point_calloc))
        sprintf(func, "crd_point_calloc");


    int type;
    int types[POINT_TYPES] = {CHARM_CRD_POINT_SCATTERED,
                              CHARM_CRD_POINT_GRID,
                              CHARM_CRD_POINT_GRID_GL,
                              CHARM_CRD_POINT_GRID_DH1,
                              CHARM_CRD_POINT_GRID_DH2};


    char func_call_str[NSTR_LONG];
    CHARM(point) *pnt;


    /* --------------------------------------------------------------------- */
    for (int t = 0; t < POINT_TYPES; t++)
    {
        type = types[t];


        for (size_t i = 0; i < nlat; i++)
        {
            for (size_t j = 0; j < nlon; j++)
            {
                pnt = crd_point_alloc(type, i, j);
                sprintf(func_call_str, "%s(%d, %zu, %zu)", func, type, i, j);


                /* For scattered points, "nlat" must be equal to "nlon" */
                if (type == CHARM_CRD_POINT_SCATTERED)
                {
                    if (i != j)
                    {
                        e += check_struct_ptr(pnt, NULL, NEQ, INVALID,
                                              func_call_str,
                                              "didn't return a NULL pointer");
                        continue;
                    }
                }


                if ((i == 0) || (j == 0))
                    /* In this case, "NULL" should be returned */
                    e += check_struct_ptr(pnt, NULL, NEQ, INVALID,
                                          func_call_str,
                                          "didn't return a NULL pointer");
                else
                    /* In this case, a valid pointer should be returned */
                    e += check_struct_ptr(pnt, NULL, EQ, VALID, func_call_str,
                                          "returned a NULL pointer");


                CHARM(crd_point_free)(pnt);
            }
        }
    }
    /* --------------------------------------------------------------------- */


    /* Check that invalid value of "type" causes the allocation function to
     * return a NULL pointer */
    /* --------------------------------------------------------------------- */
    type = 9999;
    pnt  = crd_point_alloc(type, nlat, nlon);
    sprintf(func_call_str, "%s(%d, %zu, %zu)", func, type, nlat, nlon);


    e += check_struct_ptr(pnt, NULL, NEQ, INVALID, func_call_str,
                          "didn't return a NULL pointer");


    CHARM(crd_point_free)(pnt);
    /* --------------------------------------------------------------------- */


    /* Check that the members of "charm_point" are properly set */
    /* --------------------------------------------------------------------- */
    type = CHARM_CRD_POINT_GRID;
    pnt = crd_point_alloc(type, nlat, nlon);
    sprintf(func_call_str, "%s(%d, %zu, %zu)", func, type, nlat, nlon);


    e += check_struct_int(pnt->type, CHARM_CRD_POINT_GRID, NEQ, VALID,
                          func_call_str, "returned a wrong value of \"type\"");


    e += check_struct_size_t(pnt->nlat, nlat, NEQ, VALID, func_call_str,
                             "returned a wrong value of \"nlat\"");


    e += check_struct_size_t(pnt->nlon, nlon, NEQ, VALID, func_call_str,
                             "returned a wrong value of \"nlon\"");


    e += check_struct_size_t(pnt->npoint, nlat * nlon, NEQ, VALID,
                             func_call_str,
                             "returned a wrong value of \"npoint\"");


    e += check_struct_ptr(pnt->lat, NULL, EQ, VALID, func_call_str,
                          "returned a NULL pointer for \"lat\"");


    e += check_struct_ptr(pnt->lon, NULL, EQ, VALID, func_call_str,
                          "returned a NULL pointer for \"lon\"");


    e += check_struct_ptr(pnt->r, NULL, EQ, VALID, func_call_str,
                          "returned a NULL pointer for \"r\"");


    e += check_struct__Bool(pnt->owner, 1, NEQ, VALID, func_call_str,
                            "returned a wrong value of \"owner\"");


    CHARM(crd_point_free)(pnt);


    /* For quadrature point types, "w" should not be NULL; otherwise, it should
     * be NULL */
    for (int t = 0; t < POINT_TYPES; t++)
    {
        type = types[t];


        pnt = crd_point_alloc(type, nlat, nlat);
        sprintf(func_call_str, "%s(%d, %zu, %zu)", func, type, nlat, nlat);


        if ((type == CHARM_CRD_POINT_SCATTERED) ||
            (type == CHARM_CRD_POINT_GRID))
            /* "w" should be NULL */
            e += check_struct_ptr(pnt->w, NULL, NEQ, VALID, func_call_str,
                                  "didn't return a NULL pointer for \"w\"");
        else
            /* "w" should not be NULL */
            e += check_struct_ptr(pnt->w, NULL, EQ, VALID, func_call_str,
                                  "returned a NULL pointer for \"w\"");


        CHARM(crd_point_free)(pnt);
    }
    /* --------------------------------------------------------------------- */


    return e;
}

