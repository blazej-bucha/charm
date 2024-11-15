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
#include "../err/err_is_null_ptr.h"
/* ------------------------------------------------------------------------- */







CHARM(err) *CHARM(mpi_err_init)(MPI_Comm comm)
{
    CHARM(err) *err = CHARM(err_init)();
    if (CHARM(err_is_null_ptr)(err, 1, comm))
    {
        CHARM(err_free)(err);
        err = NULL;
        goto EXIT;
    }


    err->distributed = 1;
    err->comm        = comm;


EXIT:
    return err;
}
