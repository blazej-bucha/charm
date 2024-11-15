/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <limits.h>
#include "../prec.h"
#include "../err/err_set.h"
#include "mpi_count_aint.h"
#include "mpi_size_t2charm_mpi_count.h"
/* ------------------------------------------------------------------------- */






/* MPI standards prior to "4.0" require counts to be integers.  This function
 * thus casts "x" of the "size_t" datatype to "int".  If "x" exceeds the limit
 * of "int", an error is written to "err" and "-1" is returned. */
CHARM_MPI_COUNT CHARM(mpi_size_t2charm_mpi_count)(size_t x,
                                                  CHARM(err) *err)
{
#if USE__C_MPI_ROUTINES
    return (CHARM_MPI_COUNT)x;
#else
    if (x > INT_MAX)
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Failed to cast a \"size_t\" value to \"int\", because "
                       "it exceeds the \"INT_MAX\" limit.");


    return (x <= INT_MAX) ? (int)x : -1;
#endif
}
