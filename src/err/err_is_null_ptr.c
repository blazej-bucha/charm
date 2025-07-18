/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#if HAVE_MPI
#   include <mpi.h>
#endif
#include "../prec.h"
#include "err_is_null_ptr.h"
/* ------------------------------------------------------------------------- */






/* Returns "1" if "err" is "NULL" on any process.
 *
 * CHarm uses internally the "charm_err" structure to check that memory
 * allocations, structure initializations, etc. were successful.  Obviously, we
 * cannot use, however, "charm_err" to check that "charm_err" itself was
 * initialized properly.  To this end, we have this function.
 *
 * The function is useful mostly in internal functions that call some
 * collective MPI routines after "charm_err" has been initialized.  If there
 * are no collective MPI calls after "charm_err" is initialized, one can use
 * the simply check "err == NULL".
 *
 * If this function is not compiled with the MPI support, it simply returns the
 * result of "err == NULL".
 *
 * Set "distributed" to "1" if the error structure "err" that is to be checked
 * is supposed to be distributed and "0" otherwise. */
_Bool CHARM(err_is_null_ptr)(CHARM(err) *err
#if HAVE_MPI
                             , _Bool distributed,
                             MPI_Comm comm
#endif
                            )
{
    _Bool error = err == NULL;
    _Bool error_any = error;
#if HAVE_MPI
    if (distributed)
        MPI_Allreduce(&error, &error_any, 1, MPI_C_BOOL, MPI_LOR, comm);
#endif

    return error_any;
}
