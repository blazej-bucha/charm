/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../src/prec.h"
#include "../src/simd/simd.h"
#include "../src/shs/shs_check_single_derivative.h"
#include "generate_point.h"
#include "parameters.h"
#include "error_messages.h"
#ifdef GENREF
#   include "array2file.h"
#else
#   include "validate.h"
#endif
#include "modify_low_degree_coefficients.h"
#include "check_shs_point_all.h"
/* ------------------------------------------------------------------------- */





#if COMPILE_SHS == 1


#if (GRAD_0 != 0) && (GRAD_0 != 1)
#   error "Wrong value of \"GRAD_0\"."
#endif
#if (GRAD_1 != 0) && (GRAD_1 != 1)
#   error "Wrong value of \"GRAD_1\"."
#endif
#if (GRAD_2 != 0) && (GRAD_2 != 1)
#   error "Wrong value of \"GRAD_2\"."
#endif
#if (GURU != 0) && (GURU != 1)
#   error "Wrong value of \"GRAD_2\"."
#endif
#if (GRAD_0 + GRAD_1 + GRAD_2 + GURU) != 1
#   error "Wrong combination of \"GRAD_0\", \"GRAD_1\", \"GRAD_2\" and \"GURU\"."
#endif


#undef NPAR
#if GRAD_0 == 1
#   define NPAR (1)
#elif GRAD_1 == 1
#   define NPAR (3)
#elif GRAD_2 == 1
#   define NPAR (6)
#elif GURU == 1
#   define NPAR (10)
#endif


#undef GURU_LOOP
#define GURU_LOOP(x)                                                          \
    unsigned gn = 0;                                                          \
    for (dr = 0; dr <= SHS_MAX_DERIVATIVE; dr++)                              \
    {                                                                         \
        for (dlat = 0; dlat <= SHS_MAX_DERIVATIVE; dlat++)                    \
        {                                                                     \
            for (dlon = 0; dlon <= SHS_MAX_DERIVATIVE; dlon++)                \
            {                                                                 \
                if ((dr + dlat + dlon) <= SHS_MAX_DERIVATIVE)                 \
                {                                                             \
                    x;                                                        \
                    gn++;                                                     \
                }                                                             \
            }                                                                 \
        }                                                                     \
    }


static void allocate_signal(REAL **f, size_t npoint)
{
    for (size_t i = 0; i < NPAR; i++)
    {
        f[i] = (REAL *)malloc(npoint * sizeof(REAL));
        if (f[i] == NULL)
        {
            fprintf(stderr, "%s", CHARM_ERR_MALLOC_FAILURE"\n");
            exit(CHARM_FAILURE);
        }
    }


    return;
}


static void free_signal(REAL **f)
{
    for (size_t i = 0; i < NPAR; i++)
        free(f[i]);


    return;
}


#if GURU == 1
/* For a given "j", this function sets "drp", "dlatp" and "dlonp" to the
 * correspondinf "dr", "dlat" and "dlon" values. */
static void stop_at_j_corresponding_to_dr_dlat_dlon(size_t j,
                                                    unsigned *drp,
                                                    unsigned *dlatp,
                                                    unsigned *dlonp)
{
    unsigned dr   = *drp;
    unsigned dlat = *dlatp;
    unsigned dlon = *dlonp;


    /* Let the "GURU_LOOP" run and stop it at that values of "dr", "dlat" and
     * "dlon" that correspond to "j".  The last values of "dlat" and "dlon"
     * then say whether or not the signal should be set to zeros. */
    {
    GURU_LOOP(
              if (j == gn)
                  goto J_EQ_GN;
             );
    }


    /* If we got here, we are in serious troubles, as this implies wrong
     * implementation */
    fprintf(stderr, "%s", "Wrong implementation of "
                          "\"stop_at_j_corresponding_to_dr_dlat_dlon\".");
    exit(CHARM_FAILURE);


J_EQ_GN:
    *drp   = dr;
    *dlatp = dlat;
    *dlonp = dlon;


    return;
}
#endif


/* For low values of "nmax", some elements of "grad1" and "grad2" are zero by
 * definition, particularly for "nmax = 0" and "nmax = 1.  Then, we want the
 * *exported* values to be exactly zero instead of values of orders like "-20"
 * in double precision (or, say, order "-40" in quadruple precision).
 * Otherwise, these values of low orders are not treated as zeros in validation
 * tests as should be the case and false warnings are printed.
 *
 * This function therefore sets "npoints" elements of "f" to zero based on the
 * the symbolic constants "GRAD_0", "GRAD_1" and "GRAD_2", the element of the
 * gradient being treated ("j") and based on "nmax". */
#ifdef GENREF
static void zeros_for_low_nmax(REAL *f,
                               size_t npoint,
                               size_t j,
                               unsigned long nmax)
{
    _Bool zeros = 0;


#if GRAD_0 == 1

    /* No zeroing for "potential" regardless of "nmax" */
    zeros = 0;

#elif GRAD_1 == 1

    if (nmax == 0)
        if ((j == 0) || (j == 1))  /* The "x" and "y" elements of grad1,
                                    * respectively */
            zeros = 1;

#elif GRAD_2 == 1

    if (nmax == 0)
    {
        if ((j == 1) || (j == 2) || (j == 4))  /* The "xy", "xz" and "yz"
                                                * elements of grad2,
                                                * respectively */
            zeros = 1;
    }
    else if (nmax == 1)
    {
        if (j == 1)  /* "xy" */
            zeros = 1;
    }

#elif GURU == 1

    unsigned dr, dlat, dlon;
    stop_at_j_corresponding_to_dr_dlat_dlon(j, &dr, &dlat, &dlon);


    if (nmax == 0)
    {
        if ((dlat > 0) || (dlon > 0))
            zeros = 1;
        else
            zeros = 0;
    }


    /* Keep in mind that that the "(nmax == 1) && (dlat == 1) && (dlon == 1)"
     * case must not be zeroed out. */
    zeros = 0;

#endif


    if (zeros)
    {
        for (size_t i = 0; i < npoint; i++)
            f[i] = PREC(0.0);
    }


    return;
}
#endif


/* For all points with latitudes in absolute value larger than "LAT_THRESHOLD",
 * sets the signal in "f" to zero. */
static void treat_singularity(CHARM(point) *pnt, REAL *f)
{
    REAL LAT_THRESHOLD = PI_2 - PREC(100.0) * FABS(CHARM(glob_threshold));


    if (pnt->type == CHARM_CRD_POINT_SCATTERED)
    {
        for (size_t i = 0; i < pnt->npoint; i++)
            if (FABS(pnt->lat[i]) > LAT_THRESHOLD)
                f[i] = PREC(0.0);
    }
    else
    {
        for (size_t i = 0; i < pnt->nlat; i++)
        {
            if (FABS(pnt->lat[i]) > LAT_THRESHOLD)
            {
                for (size_t j = 0; j < pnt->nlon; j++)
                    f[i * pnt->nlon + j] = PREC(0.0);
            }
        }

    }


    return;
}


/* Sets to zero signal at points near the polar areas to avoid validating
 * inaccurate results that are subject to issues related to the floating point
 * arithmetic. */
static void zeros_for_singularities(CHARM(point) *pnt,
                                    REAL *f,
                                    size_t j)
{
    _Bool sing = 0;


#if GRAD_0 == 1

    /* No singularities for "potentials" */
    sing = 0;

#elif GRAD_1 == 1

    if (j != 2)  /* All elements of grad1 except for "z" */
        sing = 1;

#elif GRAD_2 == 1

    if (j != 5)  /* All elements of grad2 except for "zz" */
        sing = 1;

#elif GURU == 1

    unsigned dr, dlat, dlon;
    stop_at_j_corresponding_to_dr_dlat_dlon(j, &dr, &dlat, &dlon);


    if ((dlat > 0) || (dlon > 0))
        sing = 1;

#endif


    if (sing)
        treat_singularity(pnt, f);


    return;
}


#if GRAD_0 == 1
long int check_shs_point(void)
#elif GRAD_1 == 1
long int check_shs_point_grad1(void)
#elif GRAD_2 == 1
long int check_shs_point_grad2(void)
#elif GURU == 1
long int check_shs_point_guru(void)
#endif
{
    /* Read reference potential coefficients */
    /* --------------------------------------------------------------------- */
    CHARM(shc) *shcs = CHARM(shc_calloc)(SHCS_NMAX_POT, PREC(1.0), PREC(1.0));
    if (shcs == NULL)
    {
        fprintf(stderr, "%s", ERR_MSG_SHC);
        exit(CHARM_FAILURE);
    }


    /* Error structure */
    CHARM(err) *err = CHARM(err_init)();
    if (err == NULL)
    {
        fprintf(stderr, "%s", ERR_MSG_ERR);
        exit(CHARM_FAILURE);
    }


    CHARM(shc_read_mtx)(SHCS_IN_PATH_POT_MTX, SHCS_NMAX_POT, shcs, err);
    CHARM(err_handler)(err, 1);


    /* Modify coefficients of degrees "0" and "1" to allow for an accurate
     * validation in all precisions. */
    modify_low_degree_coefficients(shcs);
    /* --------------------------------------------------------------------- */






    /* Calling any synthesis routine with zero points in "charm_point" is valid
     * in CHarm and must not produce any error or fail, so let's check this. */
    /* --------------------------------------------------------------------- */
    long int e = 0;


    CHARM(point) *sctr_0points = NULL;
    sctr_0points = CHARM(crd_point_malloc)(CHARM_CRD_POINT_SCATTERED, 0, 0);


    REAL **f = (REAL **)malloc(NPAR * sizeof(REAL *));
    if (f == NULL)
    {
        fprintf(stderr, "%s", CHARM_ERR_MALLOC_FAILURE"\n");
        exit(CHARM_FAILURE);
    }
    for (size_t i = 0; i < NPAR; i++)
        f[i] = NULL;


#if GRAD_0 == 1
    CHARM(shs_point)(sctr_0points, shcs, shcs->nmax, f[0], err);
#elif GRAD_1 == 1
    CHARM(shs_point_grad1)(sctr_0points, shcs, shcs->nmax, f, err);
#elif GRAD_2 == 1
    CHARM(shs_point_grad2)(sctr_0points, shcs, shcs->nmax, f, err);
#elif GURU == 1
    unsigned dr, dlat, dlon;
    {
    GURU_LOOP(
              CHARM(shs_point_guru)(sctr_0points, shcs, shcs->nmax,
                                    dr, dlat, dlon, f[gn], err);
              goto SHS_0POINTS_ERROR;
             );
    }
SHS_0POINTS_ERROR:
#endif
    if (!CHARM(err_isempty)(err))
    {
        printf("\n        WARNING: Synthesis with zero points didn't pass!\n");
        e += 1;
    }


    /* If we get here, this means that "err" was empty after calling each
     * synthesis routine and none of the call produced, e.g., segfault. */
    CHARM(err_reset)(err);
    CHARM(crd_point_free)(sctr_0points);
    /* --------------------------------------------------------------------- */






    /* GL, DH1 and DH2 point grids */
    /* ..................................................................... */
    CHARM(point) *grd_pnt = NULL;
    char file[10][NSTR_LONG];  /* "10", because "CHARM(shs_point_guru)" offers
                                * "10" combinations of "dr", "dlat" and "dlon".
                                * */
    char elem[10][NSTR_SHORT];
    char grd_str[NSTR_SHORT];


#if GRAD_0 == 1
    strcpy(elem[0], "");
#elif GRAD_1 == 1
    strcpy(elem[0], "_x");
    strcpy(elem[1], "_y");
    strcpy(elem[2], "_z");
#elif GRAD_2 == 1
    strcpy(elem[0], "_xx");
    strcpy(elem[1], "_xy");
    strcpy(elem[2], "_xz");
    strcpy(elem[3], "_yy");
    strcpy(elem[4], "_yz");
    strcpy(elem[5], "_zz");
#elif GURU == 1
    char guru_ijk[9];
    {
    GURU_LOOP(
              snprintf(guru_ijk, 9, "_guru%u%u%u", dr, dlat, dlon);
              strcpy(elem[gn], guru_ijk);
             );
    }
#endif


    {
    int grd_types[3] = {CHARM_CRD_POINT_GRID_GL,
                        CHARM_CRD_POINT_GRID_DH1,
                        CHARM_CRD_POINT_GRID_DH2};


#if GRAD_0 == 1
    for (int g = 0; g < 3; g++)
#else
    /* No need for such a large number of combination in this case, as many
     * code branches are shared among the functions */
    for (int g = 0; g < 1; g++)
#endif
    {
        int grd_type = grd_types[g];


        for (unsigned long nmax = 0; nmax <= NMAX; nmax++)
        {
            for (int deltar = 0; deltar < NDELTAR; deltar++)
            {
                REAL rref = shcs->r + (REAL)(DELTAR) * (REAL)deltar;


                if (grd_type == CHARM_CRD_POINT_GRID_GL)
                    grd_pnt = CHARM(crd_point_gl)(nmax, rref);
                else if (grd_type == CHARM_CRD_POINT_GRID_DH1)
                    grd_pnt = CHARM(crd_point_dh1)(nmax, rref);
                else if (grd_type == CHARM_CRD_POINT_GRID_DH2)
                    grd_pnt = CHARM(crd_point_dh2)(nmax, rref);
                if (grd_pnt == NULL)
                {
                    fprintf(stderr, "%s", ERR_MSG_POINT);
                    exit(CHARM_FAILURE);
                }


                /* Generate output file name */
                if (g == 0)
                    strcpy(grd_str, "gl");
                else if (g == 1)
                    strcpy(grd_str, "dh1");
                else if (g == 2)
                    strcpy(grd_str, "dh2");
                for (size_t j = 0; j < NPAR; j++)
                    snprintf(file[j], NSTR_LONG,
                             "%s/shs_p_nx%lu_dr%d_fft1_%s%s%s",
                             FOLDER, nmax, deltar, grd_str, elem[j], FTYPE);


                allocate_signal(f, grd_pnt->npoint);


#if GRAD_0 == 1
                CHARM(shs_point)(grd_pnt, shcs, nmax, f[0], err);
#elif GRAD_1 == 1
                CHARM(shs_point_grad1)(grd_pnt, shcs, nmax, f, err);
#elif GRAD_2 == 1
                CHARM(shs_point_grad2)(grd_pnt, shcs, nmax, f, err);
#elif GURU == 1
                GURU_LOOP(
                          CHARM(shs_point_guru)(grd_pnt, shcs, nmax,
                                                dr, dlat, dlon, f[gn], err);
                          CHARM(err_handler)(err, 1);
                         );
#endif
                CHARM(err_handler)(err, 1);


                for (size_t j = 0; j < NPAR; j++)
                {
                    zeros_for_singularities(grd_pnt, f[j], j);
#ifdef GENREF
                    zeros_for_low_nmax(f[j], grd_pnt->npoint, j, nmax);
                    e += array2file(file[j], f[j], grd_pnt->npoint);
#else
                    e += validate(file[j], f[j], grd_pnt->npoint,
#   if GRAD_0 == 1
                                  PREC(10.0)
#   else
                                  PREC(1000.0)
#   endif
                                  * CHARM(glob_threshold));
#endif
                }


                CHARM(crd_point_free)(grd_pnt);
                free_signal(f);
            }
        }
    }
    }
    /* ..................................................................... */






    /* Custom point grids */
    /* ..................................................................... */
    /* If we got here without any errors in "CHARM(shs_point_guru)", the tests
     * that follow add little value, as the some code branch is use as with
     * "CHARM(shs_point)". */
#if GURU != 1
    {
    size_t nlat[NCUSTOM_GRD] = {1, 1, 3, 10};
    size_t nlon[NCUSTOM_GRD] = {1, 2, 8, 22};


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
                        REAL r = shcs->r + (REAL)(DELTAR) * (REAL)deltar;


                        grd_pnt = CHARM(crd_point_calloc)(CHARM_CRD_POINT_GRID,
                                                          nlat[i], nlon[i]);
                        if (grd_pnt == NULL)
                        {
                            fprintf(stderr, "%s", ERR_MSG_POINT);
                            exit(CHARM_FAILURE);
                        }


                        if (fft == 0)
                            CHARM(generate_point)(grd_pnt, r, PI, PI);
                        else
                            CHARM(generate_point)(grd_pnt, r, PI,
                                                  PREC(2.0) * PI);


                        /* To get a non-symmetric grid, we simply add
                         * some more or less random number to the first
                         * latitude */
                        REAL break_symm = PREC(0.0);
                        if (s == 0)
                            break_symm = (REAL)(BREAK_SYMM);
                        grd_pnt->lat[0] -= break_symm;


                        /* Generate output file name */
                        for (size_t j = 0; j < NPAR; j++)
                            snprintf(file[j], NSTR_LONG,
                                              "%s/shs_%s_nx%lu_n%zu_dr%d_fft%d"
                                              "_s%d%s%s",
                                     FOLDER, "p", nmax, i, deltar, fft,
                                     (s == 0) ? 0 : 1, elem[j], FTYPE);


                        allocate_signal(f, grd_pnt->npoint);


#   if GRAD_0 == 1
                        CHARM(shs_point)(grd_pnt, shcs, nmax, f[0], err);
#   elif GRAD_1 == 1
                        CHARM(shs_point_grad1)(grd_pnt, shcs, nmax, f, err);
#   elif GRAD_2 == 1
                        CHARM(shs_point_grad2)(grd_pnt, shcs, nmax, f, err);
#   endif
                        CHARM(err_handler)(err, 1);


                        for (size_t j = 0; j < NPAR; j++)
                        {
                            zeros_for_singularities(grd_pnt, f[j], j);
#   ifdef GENREF
                            zeros_for_low_nmax(f[j], grd_pnt->npoint, j, nmax);
                            e += array2file(file[j], f[j], grd_pnt->npoint);
#   else
                            e += validate(file[j], f[j], grd_pnt->npoint,
#      if GRAD_0 == 1
                                          PREC(100.0) * CHARM(glob_threshold));
#      else
                                          PREC(1000.0) *
                                          CHARM(glob_threshold2));
#      endif
#   endif
                        }


                        CHARM(crd_point_free)(grd_pnt);


                        free_signal(f);
                    }
                }
            }
        }
    }
    }
#endif
    /* ..................................................................... */






    /* Scattered points */
    /* ..................................................................... */
    {
    size_t nlat[NSCTR] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 20, 30, 30 + 1};
    size_t nlon[NSCTR];
    for (size_t i = 0; i < NSCTR; i++)
        nlon[i] = nlat[i];


    CHARM(point) *sctr_pnt = NULL;
#if GURU == 0
    for (unsigned long nmax = 0; nmax <= NMAX; nmax++)
    {
        for (size_t i = 0; i < NSCTR; i++)
        {
            for (int deltar = 0; deltar < NDELTAR; deltar++)
            {
#elif GURU == 1
    /* For the guru function, we do not need many of these tests */
    for (unsigned long nmax = NMAX; nmax <= NMAX; nmax++)
    {
        for (size_t i = NSCTR - 1; i < NSCTR; i++)
        {
            for (int deltar = NDELTAR - 1; deltar < NDELTAR; deltar++)
            {
#endif
                REAL r = shcs->r + (REAL)(DELTAR) * (REAL)deltar;


                sctr_pnt = CHARM(crd_point_malloc)(CHARM_CRD_POINT_SCATTERED,
                                                   nlat[i], nlon[i]);
                if (sctr_pnt == NULL)
                {
                    fprintf(stderr, "%s", ERR_MSG_POINT);
                    exit(CHARM_FAILURE);
                }


                CHARM(generate_point)(sctr_pnt, r, PI, PREC(2.0) * PI);


                /* Generate output file name */
                for (size_t j = 0; j < NPAR; j++)
                    snprintf(file[j], NSTR_LONG,
                             "%s/shs_%s_nx%lu_n%zu_dr%d_sctr%s%s",
                             FOLDER, "p", nmax, i, deltar, elem[j], FTYPE);


                allocate_signal(f, sctr_pnt->npoint);


#if GRAD_0 == 1
                CHARM(shs_point)(sctr_pnt, shcs, nmax, f[0], err);
#elif GRAD_1 == 1
                CHARM(shs_point_grad1)(sctr_pnt, shcs, nmax, f, err);
#elif GRAD_2 == 1
                CHARM(shs_point_grad2)(sctr_pnt, shcs, nmax, f, err);
#elif GURU == 1
                GURU_LOOP(
                          CHARM(shs_point_guru)(sctr_pnt, shcs, nmax,
                                                dr, dlat, dlon, f[gn], err);
                          CHARM(err_handler)(err, 1);
                         );
#endif
                CHARM(err_handler)(err, 1);


                for (size_t j = 0; j < NPAR; j++)
                {
                    zeros_for_singularities(sctr_pnt, f[j], j);
#ifdef GENREF
                    zeros_for_low_nmax(f[j], sctr_pnt->npoint, j, nmax);
                    e += array2file(file[j], f[j], sctr_pnt->npoint);
#else
                    e += validate(file[j], f[j], sctr_pnt->npoint,
#   if GRAD_0 == 1
                                  PREC(10.0)
#   else
                                  PREC(1000.0)
#   endif
                                  * CHARM(glob_threshold));
#endif
                }


                CHARM(crd_point_free)(sctr_pnt);


                free_signal(f);
            }
        }
    }
    }


    CHARM(shc_free)(shcs);
    shcs = NULL;
    /* --------------------------------------------------------------------- */






    /* Dynamical switching and loop unrolling at custom point grids.  Similarly
     * as in the tests of harmonic analysis, we use some fake coefficients.
     * See "check_sha_point.c" for further details. */
    /* ..................................................................... */
    /* If we got here without any errors in "CHARM(shs_point_guru)", the tests
     * that follow add little value, as the some code branch is use as with
     * "CHARM(shs_point)". */
#if GURU != 1
    {
    size_t nlat[NCUSTOM_GRD] = {10, 11, 12, 13};
    size_t nlon[NCUSTOM_GRD] = {22, 24, 26, 28};


    for (unsigned long nmax = NMAX_DS_MIN; nmax <= NMAX_DS_MAX; nmax++)
    {
        shcs = CHARM(shc_calloc)(nmax, PREC(1.0), PREC(1.0));
        if (shcs == NULL)
        {
            fprintf(stderr, "%s", ERR_MSG_SHC);
            exit(CHARM_FAILURE);
        }
        for (size_t m = 0; m <= nmax; m++)
        {
            for (size_t n = m; n <= nmax; n++)
            {
                shcs->c[m][n - m] = PREC(1.0) / (REAL)(n + 1);
                if (m > 0)
                    shcs->s[m][n - m] = PREC(1.0) / (REAL)(n + 1);
            }
        }


        for (size_t i = 0; i < NCUSTOM_GRD; i++)
        {
            for (int fft = 0; fft < 2; fft++)
            {
                if (fft == 1)
                {
                    if (((nlon[i] - 1) / 2 < nmax) || (nlon[i] < 2))
                        continue;
                }


                for (int s = 0; s < 2; s++) /* "s = 0" for non-symm grids */
                {
                    if ((nlat[i] == 1) && (s == 1))
                        /* If there is only a single latitude, the grid is
                         * automatically non-symmetric */
                        continue;


                    for (int deltar = 0; deltar < NDELTAR; deltar++)
                    {
                        REAL r = shcs->r + (REAL)(DELTAR) * (REAL)deltar;


                        grd_pnt = CHARM(crd_point_gl)(nmax, r);
                        if (grd_pnt == NULL)
                        {
                            fprintf(stderr, "%s", ERR_MSG_POINT);
                            exit(CHARM_FAILURE);
                        }
                        /* Given that we are testing features that are related
                         * to latitudes only, we can do a nasty thing and
                         * modify the number of longitudes in "grd_pnt" to "1",
                         * so that the synthesis will be done only for the
                         * first meridian in "grd_pnt".  This makes the tests
                         * faster and avoids having too large reference data
                         * files in "../data/tests". */
#if HAVE_MPI
                        grd_pnt->nlon = 1;
                        grd_pnt->local_nlon = 1;
#else
                        grd_pnt->nlon = 1;
#endif


                        /* To get a non-symmetric grid, we simply add
                         * some more or less random number to the first
                         * latitude */
                        REAL break_symm = PREC(0.0);
                        if (s == 0)
                            break_symm = (REAL)(BREAK_SYMM);
                        grd_pnt->lat[0] -= break_symm;


                        /* Generate output file name */
                        for (size_t j = 0; j < NPAR; j++)
                            snprintf(file[j], NSTR_LONG,
                                              "%s/shs_%s_nx%lu_n%zu_dr%d_fft%d"
                                              "_s%d_dsun%s%s",
                                     FOLDER, "p", nmax, i, deltar, fft,
                                     (s == 0) ? 0 : 1, elem[j], FTYPE);


                        allocate_signal(f, grd_pnt->npoint);


#   if GRAD_0 == 1
                        CHARM(shs_point)(grd_pnt, shcs, nmax, f[0], err);
#   elif GRAD_1 == 1
                        CHARM(shs_point_grad1)(grd_pnt, shcs, nmax, f, err);
#   elif GRAD_2 == 1
                        CHARM(shs_point_grad2)(grd_pnt, shcs, nmax, f, err);
#   endif
                        CHARM(err_handler)(err, 1);


                        for (size_t j = 0; j < NPAR; j++)
                        {
                            zeros_for_singularities(grd_pnt, f[j], j);
#   ifdef GENREF
                            /* Do not use "grd_pnt->npoint" here, because
                             * "grd_pnt->nlon" has been modified. */
                            zeros_for_low_nmax(f[j],
                                               grd_pnt->nlat * grd_pnt->nlon,
                                               j, nmax);
                            e += array2file(file[j], f[j],
                                            grd_pnt->nlat * grd_pnt->nlon);
#   else
                            /* Do not use "grd_pnt->npoint" here, because
                             * "grd_pnt->nlon" has been modified. */
                            e += validate(file[j], f[j],
                                          grd_pnt->nlat * grd_pnt->nlon,
#      if GRAD_1 == 1
                                          PREC(10.0) *
#      elif GRAD_2 == 1
                                          PREC(1000.0) *
#      endif
                                          CHARM(glob_threshold2));
#   endif
                        }


                        CHARM(crd_point_free)(grd_pnt);


                        free_signal(f);
                    }
                }
            }
        }


        CHARM(shc_free)(shcs);
    }
    }
#endif
    /* ..................................................................... */






    /* --------------------------------------------------------------------- */
    CHARM(err_free)(err);
    free(f);


    return e;
    /* --------------------------------------------------------------------- */
}


#endif
