/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#if HAVE_MPI
#   include <mpi.h>
#else
#   error "mpi.h required"
#endif
#include "../prec.h"
#include "../crd/crd_point_isEmpty.h"
#include "../crd/crd_point_issymm.h"
#include "../crd/crd_point_isQuadGrid.h"
#include "../crd/crd_point_isCustGrid.h"
#include "../crd/crd_point_isSctr.h"
#include "../crd/crd_isPoint.h"
#include "mpi_crd_point_issymm.h"
/* ------------------------------------------------------------------------- */






_Bool CHARM(mpi_crd_point_issymm)(const CHARM(point) *pnt)
{
    if (pnt == NULL)
        return 0;


    _Bool symm;


    if (CHARM(crd_point_isSctr)(pnt->type))
        symm = 0;
    else if (CHARM(crd_point_isQuadGrid)(pnt->type))
        symm = 1;
    else if (CHARM(crd_point_isCustGrid)(pnt->type))
    {
        /* Check the symmetry for a given process */
        _Bool local_symm = CHARM(crd_point_issymm)(pnt);


        /* For processes with empty "pnt", set the symmetry to "1", so that we
         * can easily use "MPI_Allreduce" with "MPI_LAND".  In other words,
         * empty "pnt" is not an obstacle for "pnt" to be symmetric across all
         * processes. */
        if (CHARM(crd_point_isEmpty)(pnt))
            local_symm = 1;


        /* Do logical "and" across all processes */
        MPI_Allreduce(&local_symm, &symm, 1, MPI_C_BOOL, MPI_LAND, pnt->comm);
    }
    else
        symm = 0;


    return symm;
}
