/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <string.h>
#include "../src/prec.h"
#include "parameters.h"
#ifdef GENREF
#   include "write.h"
#else
#   include "validate.h"
#endif
/* ------------------------------------------------------------------------- */






long int check_crd_point_quad(CHARM(point) *(*quad)(unsigned long, REAL))
{
    long int e = 0;
    char file_lat[NSTR_LONG];
    char file_lon[NSTR_LONG];
    char file_r[NSTR_LONG];
    char grd_type[NSTR_SHORT];
    CHARM(point) *grd = NULL;


    for (unsigned long nmax = 0; nmax <= NMAX; nmax++)
    {
        for (int deltar = 0; deltar < NDELTAR; deltar++)
        {
            REAL r = PREC(1.0) + DELTAR * (REAL)deltar;


            grd = quad(nmax, r);
            if (grd == NULL)
            {
                fprintf(stderr, "Failed to initialize a crd structure\n");
                exit(CHARM_FAILURE);
            }


            if (quad == CHARM(crd_point_gl))
                strcpy(grd_type, "gl");
            else if (quad == CHARM(crd_point_dh1))
                strcpy(grd_type, "dh1");
            else if (quad == CHARM(crd_point_dh2))
                strcpy(grd_type, "dh2");
            sprintf(file_lat, "%s/crd_nx%lu_dr%d_%s_lat%s",
                    FOLDER, nmax, deltar, grd_type, FTYPE);
            sprintf(file_lon, "%s/crd_nx%lu_dr%d_%s_lon%s",
                    FOLDER, nmax, deltar, grd_type, FTYPE);
            sprintf(file_r, "%s/crd_nx%lu_dr%d_%s_r%s",
                    FOLDER, nmax, deltar, grd_type, FTYPE);


#ifdef GENREF
            e += write(file_lat, grd->lat, grd->nlat);
            e += write(file_lon, grd->lon, grd->nlon);
            e += write(file_r,   grd->r,   grd->nlat);
#else
            e += validate(file_lat, grd->lat, grd->nlat,
                          PREC(10.0) * CHARM(glob_threshold));
            e += validate(file_lon, grd->lon, grd->nlon,
                          PREC(10.0) * CHARM(glob_threshold));
            e += validate(file_r,   grd->r,   grd->nlat,
                          PREC(10.0) * CHARM(glob_threshold));
#endif


            CHARM(crd_point_free)(grd);
        }
    }


    return e;
}
