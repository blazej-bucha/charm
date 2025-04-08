/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <string.h>
#include "../src/prec.h"
#include "parameters.h"
#include "error_messages.h"
#ifdef GENREF
#   include "array2file.h"
#else
#   include "validate.h"
#endif
#include "cmp_vals.h"
#include "point_touch_array_elements.h"
#include "check_crd_point_quad.h"
/* ------------------------------------------------------------------------- */






long int check_crd_point_quad(CHARM(point) *(*quad)(unsigned long, REAL))
{
    long int e = 0;
    char file_lat[NSTR_LONG];
    char file_lon[NSTR_LONG];
    char file_r[NSTR_LONG];
    char file_w[NSTR_LONG];
    char grd_type[NSTR_SHORT];
    CHARM(point) *grd = NULL;


    for (unsigned long nmax = 0; nmax <= NMAX_QUAD; nmax++)
    {
        for (int deltar = 0; deltar < NDELTAR; deltar++)
        {
            REAL r = PREC(1.0) + DELTAR * (REAL)deltar;


            grd = quad(nmax, r);
            if (grd == NULL)
            {
                fprintf(stderr, "%s", ERR_MSG_POINT);
                exit(CHARM_FAILURE);
            }


            if (quad == CHARM(crd_point_gl))
                strcpy(grd_type, "gl");
            else if (quad == CHARM(crd_point_dh1))
                strcpy(grd_type, "dh1");
            else if (quad == CHARM(crd_point_dh2))
                strcpy(grd_type, "dh2");
            snprintf(file_lat, NSTR_LONG, "%s/crd_nx%lu_dr%d_%s_lat%s",
                     FOLDER, nmax, deltar, grd_type, FTYPE);
            snprintf(file_lon, NSTR_LONG, "%s/crd_nx%lu_dr%d_%s_lon%s",
                     FOLDER, nmax, deltar, grd_type, FTYPE);
            snprintf(file_r, NSTR_LONG, "%s/crd_nx%lu_dr%d_%s_r%s",
                     FOLDER, nmax, deltar, grd_type, FTYPE);
            snprintf(file_w, NSTR_LONG, "%s/crd_nx%lu_dr%d_%s_w%s",
                     FOLDER, nmax, deltar, grd_type, FTYPE);


#ifdef GENREF
            e += array2file(file_lat, grd->lat, grd->nlat);
            e += array2file(file_lon, grd->lon, grd->nlon);
            e += array2file(file_r,   grd->r,   grd->nlat);
            e += array2file(file_w,   grd->w,   grd->nlat);
#else
            e += validate(file_lat, grd->lat, grd->nlat,
                          PREC(10.0) * CHARM(glob_threshold));
            e += validate(file_lon, grd->lon, grd->nlon,
                          PREC(10.0) * CHARM(glob_threshold));
            e += validate(file_r,   grd->r,   grd->nlat,
                          PREC(10.0) * CHARM(glob_threshold));
            e += validate(file_w,   grd->w,   grd->nlat,
                          PREC(10.0) * CHARM(glob_threshold));


            /* The sum of integration weights must equal "2.0" */
            REAL w_sum = PREC(0.0);
            for (size_t i = 0; i < grd->nlat; i++)
                w_sum += grd->w[i];
            e += cmp_vals_real(w_sum, PREC(2.0),
                               PREC(10.0) * CHARM(glob_threshold));
#endif


            point_touch_array_elements(grd);
            CHARM(crd_point_free)(grd);
        }
    }


    return e;
}
