/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../src/prec.h"
#include "generate_cell.h"
#include "parameters.h"
#include "error_messages.h"
#ifdef GENREF
#   include "array2file.h"
#else
#   include "validate.h"
#endif
#include "modify_low_degree_coefficients.h"
#include "check_shs_cell_isurf.h"
/* ------------------------------------------------------------------------- */






long int check_shs_cell_isurf(void)
{
    /* Read reference potential coefficients */
    /* --------------------------------------------------------------------- */
    CHARM(shc) *shcs_pot = CHARM(shc_calloc)(SHCS_NMAX_POT, PREC(1.0),
                                             PREC(1.0));
    if (shcs_pot == NULL)
    {
        fprintf(stderr, ERR_MSG_SHC);
        exit(CHARM_FAILURE);
    }


    /* Error structure */
    CHARM(err) *err = CHARM(err_init)();
    if (err == NULL)
    {
        fprintf(stderr, ERR_MSG_ERR);
        exit(CHARM_FAILURE);
    }


    CHARM(shc_read_mtx)(SHCS_IN_PATH_POT_MTX, SHCS_NMAX_POT, shcs_pot, err);
    CHARM(err_handler)(err, 1);


    /* Modify coefficients of degrees "0" and "1" to allow for an accurate
     * validation in all precisions. */
    modify_low_degree_coefficients(shcs_pot);
    /* --------------------------------------------------------------------- */


    /* Read reference topo coefficients */
    /* --------------------------------------------------------------------- */
    CHARM(shc) *shcs_topo = CHARM(shc_calloc)(SHCS_NMAX_TOPO, PREC(1.0),
                                              PREC(1.0));
    if (shcs_topo == NULL)
    {
        fprintf(stderr, ERR_MSG_SHC);
        exit(CHARM_FAILURE);
    }


    CHARM(shc_read_mtx)(SHCS_IN_PATH_TOPO_MTX, SHCS_NMAX_TOPO, shcs_topo, err);
    CHARM(err_handler)(err, 1);
    /* --------------------------------------------------------------------- */


    /* Check whether zero cells in "charm_cell" do not produce any kind of
     * error (see "check_shs_point_all.c" for details) */
    /* --------------------------------------------------------------------- */
    long int e = 0;
    REAL *f    = NULL;


    CHARM(cell) *grd_0cells = NULL;
    grd_0cells = CHARM(crd_cell_malloc)(CHARM_CRD_CELL_GRID, 0, 0);


    CHARM(shs_cell_isurf)(grd_0cells, shcs_pot, NMAX, shcs_topo,
                          NMAX, NMAX2, NMAX2, f, err);
    if (!CHARM(err_isempty)(err))
    {
        printf("\n        WARNING: Synthesis with zero cells didn't pass!\n");
        e += 1;
    }


    CHARM(err_reset)(err);
    CHARM(crd_cell_free)(grd_0cells);
    /* --------------------------------------------------------------------- */


    /* SHS of area-mean values on irregular surfaces */
    /* --------------------------------------------------------------------- */
    size_t nlat[NCUSTOM_GRD_ISURF] = {1, 2, 10};
    size_t nlon[NCUSTOM_GRD_ISURF] = {1, 6, 22};


    CHARM(cell) *grd    = NULL;
    REAL rref           = (REAL)(RREF);
    shcs_topo->c[0][0] += rref; /* Reference the topography to a sphere */
    char file[NSTR_LONG];


    for (unsigned long nmax_p = 0; nmax_p <= NMAX; nmax_p++)
    {
        for (unsigned long nmax_t = 0; nmax_t <= NMAX; nmax_t++)
        {
            for (size_t i = 0; i < NCUSTOM_GRD_ISURF; i++)
            {
                grd = CHARM(crd_cell_calloc)(CHARM_CRD_CELL_GRID,
                                             nlat[i], nlon[i]);
                if (grd == NULL)
                {
                    fprintf(stderr, ERR_MSG_POINT);
                    exit(CHARM_FAILURE);
                }


                CHARM(generate_cell)(grd, PREC(0.0), PI, PREC(2.0) * PI);


                /* Generate output file name */
                snprintf(file, NSTR_LONG, "%s/shs_c_is_nxp%lu_nxt%lu_n%zu%s",
                              FOLDER, nmax_p, nmax_t, i, FTYPE);


                f = (REAL *)malloc(grd->ncell * sizeof(REAL));
                if (f == NULL)
                {
                    fprintf(stderr, CHARM_ERR_MALLOC_FAILURE"\n");
                    exit(CHARM_FAILURE);
                }


                CHARM(shs_cell_isurf)(grd, shcs_pot, nmax_p, shcs_topo,
                                      nmax_t, NMAX2, NMAX2, f, err);
                CHARM(err_handler)(err, 1);


#ifdef GENREF
                e += array2file(file, f, grd->ncell);
#else
                e += validate(file, f, grd->ncell, CHARM(glob_threshold2));
#endif


                CHARM(crd_cell_free)(grd);
                free(f);
            }
        }
    }


    shcs_topo->c[0][0] -= rref;
    /* --------------------------------------------------------------------- */


    /* --------------------------------------------------------------------- */
    CHARM(err_free)(err);
    CHARM(shc_free)(shcs_pot);
    CHARM(shc_free)(shcs_topo);


    return e;
    /* --------------------------------------------------------------------- */
}
