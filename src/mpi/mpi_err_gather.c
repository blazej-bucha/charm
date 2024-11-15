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
#include "../err/err_set.h"
#include "mpi_err_gather.h"
/* ------------------------------------------------------------------------- */






/* If at least one MPI process stores an error in "err", this function sets
 * "err" to "ERR_GATHER_MSG" at MPI processes that reported no error.
 *
 * The function should be used before leaving all functions from public MPI
 * that use MPI calls and have "charm_err" as their input parameter. */
void CHARM(mpi_err_gather)(CHARM(err) *err)
{
    /* Do nothing if "err" is not distributed */
    if (!err->distributed)
        return;


    /* Do nothing if all "err" are empty */
    if (CHARM(mpi_err_isempty)(err))
        return;


    /* At least one MPI process holds an error, so write a general error
     * message to all MPI processes that do not report an error. */
    if (CHARM(err_isempty)(err))
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMPIPROCESS,
                       ERR_GATHER_MSG);


    return;
}
