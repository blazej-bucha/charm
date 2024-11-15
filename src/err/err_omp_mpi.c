/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "err_set.h"
#include "err_isempty_all_mpi_processes.h"
#if HAVE_MPI
#   include "../mpi/mpi_err_gather.h"
#endif
#include "err_omp_mpi.h"
/* ------------------------------------------------------------------------- */






/* Detects an error across all OpenMP threads and MPI processes.  Returns
 * non-zero value if an error was found on at least one thread or processes and
 * zero otherwise.
 *
 * "err_glob" must be an OpenMP shared variable initialized to zero.
 *
 * "err_priv" indicates an error on an OpenMP thread, non-zero value in case of
 *            an error and zero otherwise.
 *
 * "err_msg" and "err_code" are error messages to be written to "err".
 *
 * */
int CHARM(err_omp_mpi)(int *err_glob,
                       int *err_priv,
                       const char *err_msg,
                       int err_code,
                       CHARM(err) *err)
{
    /* First, detect an error across OpenMP threads and write it to "err" if
     * any */
    /* --------------------------------------------------------------------- */
#if HAVE_OPENMP
#pragma omp critical
#endif
    *err_glob += *err_priv;


#if HAVE_OPENMP
#pragma omp barrier
#endif
    if (*err_glob && CHARM(err_isempty)(err))
#if HAVE_OPENMP
#pragma omp master
#endif
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, err_code, err_msg);
    /* --------------------------------------------------------------------- */


#if HAVE_MPI
    /* Now do the same but across MPI processes */
    /* --------------------------------------------------------------------- */
    *err_priv = 0;
#   if HAVE_OPENMP
#pragma omp master
#   endif
    {
        *err_priv = !CHARM_ERR_ISEMPTY_ALL_MPI_PROCESSES(err);
        if (err_priv)
            CHARM(mpi_err_gather)(err);
    }


#   if HAVE_OPENMP
#pragma omp critical
#   endif
    *err_glob += *err_priv;
    /* --------------------------------------------------------------------- */
#endif


#if HAVE_OPENMP
#pragma omp barrier
#endif
    return *err_glob;
}
