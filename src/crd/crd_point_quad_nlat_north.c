/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdint.h>
#include "../prec.h"
#include "../err/err_set.h"
#include "crd_point_isGLGrid.h"
#include "crd_point_isDHGrid.h"
#include "crd_point_quad_nlat_north.h"
/* ------------------------------------------------------------------------- */







size_t CHARM(crd_point_quad_nlat_north)(size_t local_nlat,
                                        size_t local_0_start,
                                        size_t nlat,
                                        int grd_type,
                                        unsigned long nmax,
                                        CHARM(err) *err)
{
    /* --------------------------------------------------------------------- */
    char err_msg[CHARM_ERR_MAX_MSG];


    if (local_nlat > nlat)
    {
        sprintf(err_msg, "\"local_nlat = %zu\" cannot be larger than the "
                         "number of latitudes \"%zu\".", local_nlat, nlat);
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       err_msg);
        return SIZE_MAX;
    }
    /* --------------------------------------------------------------------- */






    /* Total number of latitudes of the full grid over the northern hemisphere
     * (including the equator) for a given "nmax" */
    /* --------------------------------------------------------------------- */
    size_t nlat_north;
    if (CHARM(crd_point_isGLGrid)(grd_type))
        nlat_north = (nlat + 1) / 2;
    else if (CHARM(crd_point_isDHGrid)(grd_type))
        nlat_north = nlat / 2 + 1;
    else
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Wrong value of \"grd_type\".");
        return SIZE_MAX;
    }
    /* --------------------------------------------------------------------- */






    /* --------------------------------------------------------------------- */
    /* The same as "nlat_north" but for a chunk of the grid based on
     * "local_0_start". */
    size_t local_nlat_north;


    /* Index of the equator */
    size_t equator;


    if ((local_0_start == 0) && (local_nlat == nlat))
    {
        /* We are computing the full GL grid */
        local_nlat_north = nlat_north;
    }
    else
    {
        /* We are computing a portion of the GL grid only. */


        /* The cases when "local_nlat" is odd are treated below */
        local_nlat_north = local_nlat / 2;


        if (CHARM(crd_point_isGLGrid)(grd_type))
        {
            if (nmax % 2)
            {
                /* "nmax" is odd, so the GL grid does not include the equator.
                 * "local_nlat_north" must therefore always be even, because we
                 * have the same latitudes over both hemispheres but with
                 * different signs.  So if "local_nlat" is odd, this is an
                 * error. */
                if (local_nlat % 2)
                {
                    sprintf(err_msg, "For odd \"nmax = %lu\", "
                                     "\"local_nlat = %zu\" must be even.",
                                     nmax, local_nlat);
                    CHARM(err_set)(err, __FILE__, __LINE__, __func__,
                                   CHARM_EFUNCARG, err_msg);
                    return SIZE_MAX;
                }
            }
            else
            {
                /* "nmax" is even, so the GL grid includes the equator.  Based
                 * on "local_0_start" and "local_nlat", we can determine
                 * whether the equator is in this latitudinal chunk.  If this
                 * is the case, then "local_nlat" must be odd because of the
                 * additional zero latitude. */


                equator = nmax / 2;


                if (local_nlat % 2)
                {
                    /* If "local_nlat" is odd, the last latitude of the chunk
                     * must be zero (the equator). */
                    if ((local_0_start + (local_nlat / 2)) != equator)
                    {
                        sprintf(err_msg, "Wrong latitudinal chunk "
                                         " \"local_0_start = %zu\" and "
                                         "\"local_nlat = %zu\".  "
                                         "The last latitude of the chunk with "
                                         "odd \"local_nlat\" must be the "
                                         "equator.  For this grid, the "
                                         "index of the equator is \"%zu\", "
                                         "but \"local_0_start + local_nlat "
                                         "/ 2 = %zu\".",
                                         local_0_start, local_nlat, equator,
                                         local_0_start + local_nlat / 2);
                        CHARM(err_set)(err, __FILE__, __LINE__, __func__,
                                       CHARM_EFUNCARG, err_msg);
                        return SIZE_MAX;
                    }


                    local_nlat_north += 1;
                }
            }
        }
        else if (CHARM(crd_point_isDHGrid)(grd_type))
        {
            _Bool local_nlat_odd = local_nlat % 2;
            equator = nlat / 2;


            if (local_0_start == 0)
            {
                /* "local_nlat" must be odd if "local_0_start" is zero */
                if ((local_nlat > 0) && !local_nlat_odd)
                {
                    sprintf(err_msg, "For \"local_0_start = %zu\", "
                                     "\"local_nlat = %zu\" must be odd.",
                                     local_0_start, local_nlat);
                    CHARM(err_set)(err, __FILE__, __LINE__, __func__,
                                   CHARM_EFUNCARG, err_msg);
                    return SIZE_MAX;
                }
            }


            if (local_nlat_odd && (local_0_start != 0))
            {
                /* If the local number of latitudes is odd and "local_0_start"
                 * is not zero, this means the chunk must contain the equator
                 * */
                if ((local_0_start + (local_nlat / 2)) != equator)
                {
                    sprintf(err_msg, "Wrong latitudinal chunk "
                                     "\"local_0_start = %zu\" and "
                                     "\"local_nlat = %zu\".  "
                                     "The last latitude of the chunk must be "
                                     "the equator.  For this grid, the "
                                     "index of the equator is \"%zu\", "
                                     "but the index of the last latitude "
                                     "\"local_0_start + local_nlat "
                                     "/ 2\" is \"%zu\".",
                                     local_0_start, local_nlat, equator,
                                     local_0_start + local_nlat / 2);
                    CHARM(err_set)(err, __FILE__, __LINE__, __func__,
                                   CHARM_EFUNCARG, err_msg);

                    return SIZE_MAX;
                }
            }


            if (local_nlat_odd)
                local_nlat_north += 1;
        }
    }
    /* --------------------------------------------------------------------- */






    /* Check that we are not reaching non-existing latitudes */
    /* --------------------------------------------------------------------- */
    if ((local_nlat > 0) && (local_0_start + local_nlat_north > nlat_north))
    {
        sprintf(err_msg, "One or more latitudes in the chunk "
                         "\"local_nlat = %zu\" and \"local_0_start = %zu\" "
                         "exceed the total number of latitudes \"%zu\".",
                         local_nlat, local_0_start, nlat);
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       err_msg);
        return SIZE_MAX;
    }
    /* --------------------------------------------------------------------- */






    return local_nlat_north;
}
