/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../src/prec.h"
#include "parameters.h"
#include "validate.h"
/* ------------------------------------------------------------------------- */






/* Tests the "crd" module */
int crd(void)
{
    int errnum = 0;
    char file_lat[NSTR], file_lon[NSTR], file_r[NSTR];
    char grd_type[NSTR2];
    CHARM(crd) * grd;


    /* Loop over quadrature grid types */
    for (int g = 0; g < 3; g++)
    {
        if (g == 0)
            printf("    Gauss--Legendre grid...\n");
        else if (g == 1)
            printf("    Driscoll--Healy grid (DH1)...\n");
        else if (g == 2)
            printf("    Driscoll--Healy grid (DH2)...\n");


        for (unsigned long nmax = 0; nmax <= NMAX; nmax++)
        {
            for (int deltar = 0; deltar < NDELTAR; deltar++)
            {
                REAL r = PREC(1.0) + DELTAR * (REAL)deltar;


                if (g == 0)
                    grd = CHARM(crd_gl)(nmax, r);
                else if (g == 1)
                    grd = CHARM(crd_dh1)(nmax, r);
                else if (g == 2)
                    grd = CHARM(crd_dh2)(nmax, r);
                if (grd == NULL)
                {
                    fprintf(stderr, "Failed to initialize a crd "
                            "structure\n");
                    exit(1);
                }


                if (g == 0)
                    strcpy(grd_type, "gl");
                else if (g == 1)
                    strcpy(grd_type, "dh1");
                else if (g == 2)
                    strcpy(grd_type, "dh2");
                sprintf(file_lat, "%s/crd_nx%lu_dr%d_%s_lat%s",
                        FOLDER, nmax, deltar, grd_type, FTYPE);
                sprintf(file_lon, "%s/crd_nx%lu_dr%d_%s_lon%s",
                        FOLDER, nmax, deltar, grd_type, FTYPE);
                sprintf(file_r, "%s/crd_nx%lu_dr%d_%s_r%s",
                        FOLDER, nmax, deltar, grd_type, FTYPE);


                errnum += validate(file_lat, grd->lat, grd->nlat,
                                   PREC(10.0) * CHARM(glob_threshold));
                errnum += validate(file_lon, grd->lon, grd->nlon,
                                   PREC(10.0) * CHARM(glob_threshold));
                errnum += validate(file_r,   grd->r,   grd->nlat,
                                   PREC(10.0) * CHARM(glob_threshold));


                CHARM(crd_free)(grd);
            }
        }
    }


    return errnum;
}
