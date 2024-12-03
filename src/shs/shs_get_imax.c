/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#if HAVE_MPI
#   include <mpi.h>
#endif
#include "../prec.h"
#include "../simd/simd.h"
#if HAVE_MPI
#   include "../mpi/mpi_size_t.h"
#endif
#include "shs_get_imax.h"
/* ------------------------------------------------------------------------- */






/* Returns the number of times the loop over latitudes should run in synthesis
 * and analysis. */
size_t CHARM(shs_get_imax)(size_t nlatdo,
                           size_t block,
                           const CHARM(point) *pnt)
{
    size_t imax;


#if HAVE_MPI
    size_t nlatdo_max;
    if (pnt->distributed)
        MPI_Allreduce(&nlatdo, &nlatdo_max, 1, CHARM_MPI_SIZE_T, MPI_MAX,
                      pnt->comm);
    else
        nlatdo_max = nlatdo;


    /* With MPI, the loop over latitudes ("i") must run the same number of
     * times for all MPI processes to avoid potential deadlocks, thus we have
     * to use "nlatdo_max". */
    imax = SIMD_MULTIPLE(nlatdo_max, SIMD_SIZE * block);
#else
    imax = SIMD_MULTIPLE(nlatdo, SIMD_SIZE * block);
#endif


    return imax;
}
