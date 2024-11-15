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
/* ------------------------------------------------------------------------- */







_Bool CHARM(mpi_err_isempty)(const CHARM(err) *err)
{
    _Bool isempty;
    if (err->distributed)
    {
        _Bool local_isempty = CHARM(err_isempty)(err);
        MPI_Allreduce(&local_isempty, &isempty, 1, MPI_C_BOOL, MPI_LAND,
                      err->comm);
    }
    else
        isempty = CHARM(err_isempty)(err);


    return isempty;
}
