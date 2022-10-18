/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../src/prec.h"
#include "cmp_arrays.h"
#include "parameters.h"
#include "validate.h"
#include "generate_crd.h"
/* ------------------------------------------------------------------------- */






/* Tests the "shs" module */
/* ------------------------------------------------------------------------- */
int shs(unsigned long nmax_topo, char SHCs_topo_file[],
        unsigned long nmax_pot, char SHCs_pot_file[])
{
    /* --------------------------------------------------------------------- */
    if (nmax_topo != 4)
    {
        fprintf(stderr, "\"nmax_topo\" has to be \"4\".\n");
        exit(CHARM_FAILURE);
    }


    if (nmax_pot != 10)
    {
        fprintf(stderr, "\"nmax_pot\" has to be \"10\".\n");
        exit(CHARM_FAILURE);
    }


    int errnum = 0;
    /* --------------------------------------------------------------------- */






    /* Read reference potential coefficients */
    /* --------------------------------------------------------------------- */
    CHARM(shc) *shcs_pot = CHARM(shc_calloc)(nmax_pot, PREC(1.0), PREC(1.0));
    if (shcs_pot == NULL)
    {
        fprintf(stderr, "Failed to initialize a \"shc\" structure.\n");
        exit(1);
    }


    /* Error structure */
    CHARM(err) *err = CHARM(err_init)();
    if (err == NULL)
    {
        fprintf(stderr, "Failed to initialize a \"err\" structure.\n");
        exit(1);
    }


    CHARM(shc_read_mtx)(SHCs_pot_file, nmax_pot, shcs_pot, err);
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
    CHARM(shc) *shcs_topo = CHARM(shc_calloc)(nmax_topo, PREC(1.0), PREC(1.0));
    if (shcs_topo == NULL)
    {
        fprintf(stderr, "Failed to initialize a \"shc\" structure.\n");
        exit(1);
    }


    CHARM(shc_read_mtx)(SHCs_topo_file, nmax_topo, shcs_topo, err);
    CHARM(err_handler)(err, 1);
    /* --------------------------------------------------------------------- */






    /* GL, DH1 and DH2 point grids */
    /* ..................................................................... */
    {
    int grd_types[3] = {CHARM_CRD_POINTS_GRID_GL,
                                  CHARM_CRD_POINTS_GRID_DH1,
                                  CHARM_CRD_POINTS_GRID_DH2};
    for (int g = 0; g < 3; g++)
    {
        int grd_type = grd_types[g];
        if (grd_type == CHARM_CRD_POINTS_GRID_GL)
            printf("    Gauss--Legendre grid...\n");
        else if (grd_type == CHARM_CRD_POINTS_GRID_DH1)
            printf("    Driscoll--Healy grid (DH1)...\n");
        else if (grd_type == CHARM_CRD_POINTS_GRID_DH2)
            printf("    Driscoll--Healy grid (DH2)...\n");
        else
        {
            printf("Incorrect grid type.\n");
            exit(1);
        }


        for (unsigned long nmax = 0; nmax <= NMAX; nmax++)
        {
            for (int deltar = 0; deltar < NDELTAR; deltar++)
            {
                REAL rref = shcs_pot->r + (REAL)(DELTAR) * (REAL)deltar;


                CHARM(crd) *grd = NULL;
                if (grd_type == CHARM_CRD_POINTS_GRID_GL)
                    grd = CHARM(crd_gl)(nmax, rref);
                else if (grd_type == CHARM_CRD_POINTS_GRID_DH1)
                    grd = CHARM(crd_dh1)(nmax, rref);
                else if (grd_type == CHARM_CRD_POINTS_GRID_DH2)
                    grd = CHARM(crd_dh2)(nmax, rref);
                if (grd == NULL)
                {
                    fprintf(stderr, "Failed to initialize a \"crd\" "
                            "structure\n");
                    exit(1);
                }


                /* Generate output file name */
                char file[NSTR] = "";
                char grd_str[NSTR2] = "";
                if (g == 0)
                    strcpy(grd_str, "gl");
                else if (g == 1)
                    strcpy(grd_str, "dh1");
                else if (g == 2)
                    strcpy(grd_str, "dh2");
                sprintf(file, "%s/shs_p_nx%lu_dr%d_fft1_%s%s",
                        FOLDER, nmax, deltar, grd_str, FTYPE);


                REAL *f = (REAL *)malloc(grd->nlat * grd->nlon * sizeof(REAL));
                if (f == NULL)
                {
                    fprintf(stderr, "malloc failure.\n");
                    exit(1);
                }


                CHARM(shs_point)(grd, shcs_pot, nmax, f, err);
                CHARM(err_handler)(err, 1);


                errnum += validate(file, f, grd->nlat * grd->nlon,
                                   PREC(10.0) * CHARM(glob_threshold));


                CHARM(crd_free)(grd);
                free(f);
            }
        }
    }
    }
    /* ..................................................................... */






    /* Custom grid of points and cells */
    /* ..................................................................... */
    {
    size_t nlat[NCUSTOM_GRD] = {1, 1, 3, 10};
    size_t nlon[NCUSTOM_GRD] = {1, 2, 8, 22};


    for (int pc = 0; pc < 2; pc++)
    {
        if (pc == 0)
            printf("    Custom point grids...\n");
        else if (pc == 1)
            printf("    Custom cell grids...\n");


        for (unsigned long nmax = 0; nmax <= NMAX; nmax++)
        {
            for (size_t i = 0; i < NCUSTOM_GRD; i++)
            {
                for (int fft = 0; fft < 2; fft++)
                {
                    if (fft == 1)
                    {
                        if (((nlon[i] - 1) / 2 < nmax) || (nlon[i] < 2))
                            continue;
                    }


                    for (int s = 0; s < 2; s++) /* "s = 0" for non-symm grids 
                                                 */
                    {
                        if ((nlat[i] == 1) && (s == 1))
                            /* If there is only a single latitude, the grid is
                             * automatically non-symmetric */
                            continue;


                        for (int deltar = 0; deltar < NDELTAR; deltar++)
                        {
                            REAL r = shcs_pot->r + (REAL)(DELTAR) *
                                     (REAL)deltar;


                            CHARM(crd) *grd = NULL;
                            if (pc == 0)
                                grd = CHARM(crd_init)(CHARM_CRD_POINTS_GRID,
                                                      nlat[i], nlon[i]);
                            else
                                grd = CHARM(crd_init)(CHARM_CRD_CELLS_GRID,
                                                      nlat[i], nlon[i]);
                            if (grd == NULL)
                            {
                                fprintf(stderr, "Failed to initialize a "
                                                "\"crd\" structure\n");
                                exit(1);
                            }


                            if (fft == 0)
                                CHARM(generate_crd)(grd, r, PI, PI);
                            else
                                CHARM(generate_crd)(grd, r, PI,
                                                    PREC(2.0) * PI);


                            /* To get a non-symmetric grid, we simply add some 
                             * more or less random number to the first latitude 
                             * */
                            REAL break_symm = PREC(0.0);
                            if (s == 0)
                                break_symm = (REAL)(BREAK_SYMM);
                            grd->lat[0] += break_symm;


                            /* Generate output file name */
                            char file[NSTR] = "";
                            sprintf(file, "%s/shs_%s_nx%lu_n%zu_dr%d_fft%d"
                                              "_s%d%s",
                                    FOLDER, (pc == 0) ? "p": "c", nmax, i,
                                    deltar, fft, (s == 0) ? 0 : 1, FTYPE);


                            REAL *f = (REAL *)malloc(grd->nlat * grd->nlon *
                                                     sizeof(REAL));
                            if (f == NULL)
                            {
                                fprintf(stderr, "malloc failure.\n");
                                exit(1);
                            }


                            if (pc == 0)
                            {
                                CHARM(shs_point)(grd, shcs_pot, nmax, f, err);
                                CHARM(err_handler)(err, 1);
                                errnum += validate(file, f,
                                                   grd->nlat * grd->nlon,
                                                   PREC(10.0) *
                                                   CHARM(glob_threshold));
                            }
                            else
                            {
                                CHARM(shs_cell)(grd, shcs_pot, nmax, f, err);
                                CHARM(err_handler)(err, 1);
                                errnum += validate(file, f,
                                                   grd->nlat * grd->nlon,
                                                   CHARM(glob_threshold2));
                            }


                            CHARM(crd_free)(grd);
                            free(f);
                        }
                    }
                }
            }
        }
    }
    }
    /* ..................................................................... */






    /* Scattered points/cells */
    /* ..................................................................... */
    {
    size_t nlat[NSCTR] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    size_t nlon[NSCTR] = {1, 2, 3, 4, 5, 6, 7, 8, 9};


    for (int pc = 0; pc < 2; pc++)
    {
        if (pc == 0)
            printf("    Scattered points...\n");
        else
            printf("    Scattered cells...\n");


        for (unsigned long nmax = 0; nmax <= NMAX; nmax++)
        {
            for (size_t i = 0; i < NSCTR; i++)
            {
                for (int deltar = 0; deltar < NDELTAR; deltar++)
                {
                    REAL r = shcs_pot->r + (REAL)(DELTAR) * (REAL)deltar;


                    CHARM(crd) *sctr = NULL;
                    if (pc == 0)
                        sctr = CHARM(crd_init)(CHARM_CRD_POINTS_SCATTERED,
                                               nlat[i], nlon[i]);
                    else
                        sctr = CHARM(crd_init)(CHARM_CRD_CELLS_SCATTERED,
                                               nlat[i], nlon[i]);
                    if (sctr == NULL)
                    {
                        fprintf(stderr, "Failed to initialize a \"crd\" "
                                "structure\n");
                        exit(1);
                    }


                    CHARM(generate_crd)(sctr, r, PI, PREC(2.0) * PI);


                    /* Generate output file name */
                    char file[NSTR] = "";
                    sprintf(file, "%s/shs_%s_nx%lu_n%zu_dr%d_sctr%s",
                            FOLDER, (pc == 0) ? "p": "c", nmax, i, deltar,
                            FTYPE);


                    REAL *f = (REAL *)malloc(sctr->nlat * sizeof(REAL));
                    if (f == NULL)
                    {
                        fprintf(stderr, "malloc failure.\n");
                        exit(1);
                    }


                    if (pc == 0)
                    {
                        CHARM(shs_point)(sctr, shcs_pot, nmax, f, err);
                        CHARM(err_handler)(err, 1);
                        errnum += validate(file, f, sctr->nlat,
                                           PREC(10.0) * CHARM(glob_threshold));
                    }
                    else
                    {
                        CHARM(shs_cell)(sctr, shcs_pot, nmax, f, err);
                        CHARM(err_handler)(err, 1);
                        errnum += validate(file, f, sctr->nlat,
                                           CHARM(glob_threshold2));
                    }


                    CHARM(crd_free)(sctr);
                    free(f);
                }
            }
        }
    }
    }
    /* ..................................................................... */






    /* SHS of area-mean values on irregular surfaces */
    /* ..................................................................... */
    {
    printf("    Grid cells on an irregular surface...\n");


    size_t nlat[NCUSTOM_GRD_ISURF] = {1, 2, 10};
    size_t nlon[NCUSTOM_GRD_ISURF] = {1, 6, 22};


    REAL rref = (REAL)(RREF);
    shcs_topo->c[0][0] += rref; /* Reference the topography to a sphere */


    for (unsigned long nmax_p = 0; nmax_p <= NMAX; nmax_p++)
    {
        for (unsigned long nmax_t = 0; nmax_t <= NMAX; nmax_t++)
        {
            for (size_t i = 0; i < NCUSTOM_GRD_ISURF; i++)
            {
                CHARM(crd) *grd = CHARM(crd_init)(CHARM_CRD_CELLS_GRID,
                                                  nlat[i], nlon[i]);
                if (grd == NULL)
                {
                    fprintf(stderr, "Failed to initialize a " "\"crd\" "
                                    "structure\n");
                    exit(1);
                }


                CHARM(generate_crd)(grd, PREC(0.0), PI, PREC(2.0) * PI);


                /* Generate output file name */
                char file[NSTR] = "";
                sprintf(file, "%s/shs_c_is_nxp%lu_nxt%lu_n%zu%s",
                              FOLDER, nmax_p, nmax_t, i, FTYPE);


                REAL *f = (REAL *)malloc(grd->nlat * grd->nlon * sizeof(REAL));
                if (f == NULL)
                {
                    fprintf(stderr, "malloc failure.\n");
                    exit(1);
                }


                CHARM(shs_cell_isurf)(grd, shcs_pot, nmax_p, shcs_topo,
                                      nmax_t, NMAX2, NMAX2, f, err);
                CHARM(err_handler)(err, 1);


                errnum += validate(file, f, grd->nlat * grd->nlon,
#if CHARM_FLOAT
                                   CHARM(glob_threshold2)
#else
                                   PREC(100.0) * CHARM(glob_threshold)
#endif
                                   );


                CHARM(crd_free)(grd);
                free(f);
            }
        }
    }


    shcs_topo->c[0][0] -= rref;
    }
    /* ..................................................................... */
    /* --------------------------------------------------------------------- */






    CHARM(err_free)(err);
    CHARM(shc_free)(shcs_pot);
    CHARM(shc_free)(shcs_topo);


    return errnum;
}
/* ------------------------------------------------------------------------- */

