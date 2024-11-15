/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#if HAVE_MPI
#   include <mpi.h>
#else
#   error "mpi.h required"
#endif
#include "../prec.h"
#include "../misc/misc_is_nearly_equal.h"
#include "../err/err_set.h"
#include "mpi_size_t.h"
#include "mpi_real.h"
#include "mpi_allequal.h"
/* ------------------------------------------------------------------------- */






#undef X1
#define X1(TYPE, xMPI_TYPE)                                                   \
        _Bool ret = 1;                                                        \
                                                                              \
                                                                              \
        int rank, size;                                                       \
        MPI_Comm_rank(comm, &rank);                                           \
        MPI_Comm_size(comm, &size);                                           \
                                                                              \
                                                                              \
        TYPE *x = (TYPE *)malloc(size * sizeof(TYPE));                        \
        if (x == NULL)                                                        \
        {                                                                     \
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,     \
                           CHARM_ERR_MALLOC_FAILURE);                         \
        }                                                                     \
                                                                              \
                                                                              \
        if (!CHARM(mpi_err_isempty)(err))                                     \
        {                                                                     \
            goto EXIT;                                                        \
        }                                                                     \
                                                                              \
                                                                              \
        MPI_Gather(&(x_local), 1, xMPI_TYPE, x, 1, xMPI_TYPE, 0, comm);


#undef X2
#define X2                                                                    \
        if (rank == 0)                                                        \
        {                                                                     \
            for (int i = 1; i < size; i++)                                    \
            {                                                                 \
                if (x[i] != x[0])                                             \
                {                                                             \
                    ret = 0;                                                  \
                    break;                                                    \
                }                                                             \
            }                                                                 \
        }


#undef X3
#define X3                                                                    \
        MPI_Bcast(&ret, 1, MPI_C_BOOL, 0, comm);                              \
                                                                              \
                                                                              \
    EXIT:                                                                     \
        free(x);                                                              \
                                                                              \
                                                                              \
        return ret;


_Bool CHARM(mpi_allequal__Bool)(_Bool x_local,
                                MPI_Comm comm,
                                CHARM(err) *err)
{
    X1(_Bool, MPI_C_BOOL);
    X2;
    X3;
}


_Bool CHARM(mpi_allequal_int)(int x_local,
                              MPI_Comm comm,
                              CHARM(err) *err)
{
    X1(int, MPI_INT);
    X2;
    X3;
}


_Bool CHARM(mpi_allequal_ulong)(unsigned long x_local,
                                MPI_Comm comm,
                                CHARM(err) *err)
{
    X1(unsigned long, MPI_UNSIGNED_LONG);
    X2;
    X3;
}


_Bool CHARM(mpi_allequal_size_t)(size_t x_local,
                                 MPI_Comm comm,
                                 CHARM(err) *err)
{
    X1(size_t, CHARM_MPI_SIZE_T);
    X2;
    X3;
}


_Bool CHARM(mpi_allequal_real)(REAL x_local,
                               MPI_Comm comm,
                               CHARM(err) *err)
{
    X1(REAL, CHARM_MPI_REAL);


    if (rank == 0)
    {
        for (int i = 1; i < size; i++)
        {
            if (!CHARM(misc_is_nearly_equal)(x[i], x[0],
                                             CHARM(glob_threshold)))
            {
                ret = 0;
                break;
            }
        }
    }


    X3;
}
