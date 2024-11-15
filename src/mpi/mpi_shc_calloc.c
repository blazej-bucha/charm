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
#include "../err/err_propagate.h"
#include "../misc/misc_calloc.h"
#include "mpi_shc_alloc.h"
#include "mpi_err_gather.h"
/* ------------------------------------------------------------------------- */






CHARM(shc) *CHARM(mpi_shc_calloc)(unsigned long nmax,
                                  REAL mu,
                                  REAL r,
                                  size_t nchunk,
                                  const unsigned long *order,
                                  MPI_Comm comm,
                                  CHARM(err) *err)
{
    CHARM(shc) *shcs = CHARM(mpi_shc_alloc)(nmax, mu, r, nchunk, order, comm,
                                            CHARM(misc_calloc), err);
    if (!CHARM(err_isempty)(err))
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);


    if (!CHARM(mpi_err_isempty)(err))
    {
        CHARM(shc_free)(shcs);
        shcs = NULL;
    }


    CHARM(mpi_err_gather)(err);


    return shcs;
}
