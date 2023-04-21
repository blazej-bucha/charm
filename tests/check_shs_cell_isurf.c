/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../src/prec.h"
#include "generate_cell.h"
#include "parameters.h"
#include "validate.h"
/* ------------------------------------------------------------------------- */






long int check_shs_cell_isurf(void)
{
    /* Read reference potential coefficients */
    /* --------------------------------------------------------------------- */
    CHARM(shc) *shcs_pot = CHARM(shc_calloc)(SHCS_NMAX_POT, PREC(1.0),
                                             PREC(1.0));
    if (shcs_pot == NULL)
    {
        fprintf(stderr, "Failed to initialize a \"shc\" structure.\n");
        exit(CHARM_FAILURE);
    }


    /* Error structure */
    CHARM(err) *err = CHARM(err_init)();
    if (err == NULL)
    {
        fprintf(stderr, "Failed to initialize a \"err\" structure.\n");
        exit(CHARM_FAILURE);
    }


    CHARM(shc_read_mtx)(SHCS_IN_PATH_POT_MTX, SHCS_NMAX_POT, shcs_pot, err);
    CHARM(err_handler)(err, 1);


    /* Modify coefficients of degrees "0" and "1" to allow for an accurate
     * validation in all precisions. */
    shcs_pot->c[0][0 - 0] = (REAL)C00;
    shcs_pot->c[0][1 - 0] = (REAL)C10;
    shcs_pot->c[1][1 - 1] = (REAL)C11;
    shcs_pot->s[1][1 - 1] = (REAL)S11;
    /* --------------------------------------------------------------------- */


    /* Read reference topo coefficients */
    /* --------------------------------------------------------------------- */
    CHARM(shc) *shcs_topo = CHARM(shc_calloc)(SHCS_NMAX_TOPO, PREC(1.0),
                                              PREC(1.0));
    if (shcs_topo == NULL)
    {
        fprintf(stderr, "Failed to initialize a \"shc\" structure.\n");
        exit(CHARM_FAILURE);
    }


    CHARM(shc_read_mtx)(SHCS_IN_PATH_TOPO_MTX, SHCS_NMAX_TOPO, shcs_topo, err);
    CHARM(err_handler)(err, 1);
    /* --------------------------------------------------------------------- */


    /* SHS of area-mean values on irregular surfaces */
    /* --------------------------------------------------------------------- */
    size_t nlat[NCUSTOM_GRD_ISURF] = {1, 2, 10};
    size_t nlon[NCUSTOM_GRD_ISURF] = {1, 6, 22};


    long int e          = 0;
    REAL *f             = NULL;
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
                    fprintf(stderr, "Failed to initialize a " "\"crd\" "
                                    "structure\n");
                    exit(CHARM_FAILURE);
                }


                CHARM(generate_cell)(grd, PREC(0.0), PI, PREC(2.0) * PI);


                /* Generate output file name */
                sprintf(file, "%s/shs_c_is_nxp%lu_nxt%lu_n%zu%s",
                              FOLDER, nmax_p, nmax_t, i, FTYPE);


                f = (REAL *)malloc(grd->nlat * grd->nlon * sizeof(REAL));
                if (f == NULL)
                {
                    fprintf(stderr, "malloc failure.\n");
                    exit(CHARM_FAILURE);
                }


                CHARM(shs_cell_isurf)(grd, shcs_pot, nmax_p, shcs_topo,
                                      nmax_t, NMAX2, NMAX2, f, err);
                CHARM(err_handler)(err, 1);


                e += validate(file, f, grd->nlat * grd->nlon,
#if CHARM_FLOAT
                              CHARM(glob_threshold2)
#else
                              PREC(100.0) * CHARM(glob_threshold)
#endif
                              );


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
