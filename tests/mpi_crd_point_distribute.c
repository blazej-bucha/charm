/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdint.h>
#if HAVE_MPI
#   include <mpi.h>
#else
#   error "mpi.h required"
#endif
#include "../src/prec.h"
#include "../src/crd/crd_point_quad_equator.h"
#include "../src/crd/crd_point_quad_get_nmax_from_nlat.h"
#include "../src/crd/crd_point_isGLGrid.h"
#include "../src/crd/crd_point_isDHGrid.h"
#include "partition_interval.h"
#include "mpi_crd_point_distribute.h"
/* ------------------------------------------------------------------------- */






void mpi_crd_point_distribute(size_t nlat,
                              int grd_type,
                              int rank,
                              int size,
                              size_t *local_nlat,
                              size_t *local_0_start)
{
    size_t x1 = 0;
    size_t x2 = (nlat - (size_t)CHARM(crd_point_isGLGrid)(grd_type)) / 2;


    size_t start, stop;
    partition_interval_size_t(x1, x2, size, rank, &start, &stop);


    *local_nlat = 0;
    *local_0_start = 0;
    if ((start != SIZE_MAX) && (stop != SIZE_MAX))
    {
        *local_0_start = start;
        *local_nlat    = 2 * (stop - start + 1);


        const unsigned long nmax_grd =
                      CHARM(crd_point_quad_get_nmax_from_nlat)(grd_type, nlat);
        size_t equator = CHARM(crd_point_quad_equator)(grd_type, nmax_grd);
        if (((nlat + CHARM(crd_point_isDHGrid)(grd_type)) % 2) &&
            ((*local_0_start + *local_nlat / 2 - 1) == equator))
            *local_nlat -= 1;


        if (CHARM(crd_point_isDHGrid)(grd_type) &&
            (*local_0_start == 0) && (*local_nlat > 0))
            *local_nlat -= 1;
    }


    return;
}

