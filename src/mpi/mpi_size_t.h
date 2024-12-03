/* This header file is not a part of API.
 *
 * The MPI standard does not define any "MPI_Datatype" that would correspond to
 * "size_t".  This header file therefore defines "MPI_SIZE_T" as an alias to
 * that standard unsigned integer "MPI_Datatype", the maximum value of which
 * corresponds the maximum value that can be stored by "size_t" on the host
 * system.
 *
 * */


#ifndef __MPI_SIZE_T_H__
#define __MPI_SIZE_T_H__


#include <config.h>
#include <stdint.h>
#include <limits.h>
#if HAVE_MPI
#   include <mpi.h>
#else
#   error "mpi.h required"
#endif


#undef CHARM_MPI_SIZE_T
#if SIZE_MAX == UCHAR_MAX
#   define CHARM_MPI_SIZE_T MPI_UNSIGNED_CHAR
#elif SIZE_MAX == USHRT_MAX
#   define CHARM_MPI_SIZE_T MPI_UNSIGNED_SHORT
#elif SIZE_MAX == UINT_MAX
#   define CHARM_MPI_SIZE_T MPI_UNSIGNED
#elif SIZE_MAX == ULONG_MAX
#   define CHARM_MPI_SIZE_T MPI_UNSIGNED_LONG
#elif SIZE_MAX == ULLONG_MAX
#   define CHARM_MPI_SIZE_T MPI_UNSIGNED_LONG_LONG
#elif SIZE_MAX == UINT8_MAX
#   define CHARM_MPI_SIZE_T MPI_UINT8_T
#elif SIZE_MAX == UINT16_MAX
#   define CHARM_MPI_SIZE_T MPI_UINT16_T
#elif SIZE_MAX == UINT32_MAX
#   define CHARM_MPI_SIZE_T MPI_UINT32_T
#elif SIZE_MAX == UINT64_MAX
#   define CHARM_MPI_SIZE_T MPI_UINT64_T
#else
#   error "couldn't find any MPI_Datatype that would correspond to size_t"
#endif


#endif

