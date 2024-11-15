/* This header file is not a part of API. */


#ifndef __CHECK_MPI_CRD_POINT_QUAD_H__
#define __CHECK_MPI_CRD_POINT_QUAD_H__


#include <config.h>
#include "../src/prec.h"
#if HAVE_MPI
#   include <mpi.h>
#else
#   error "mpi.h required"
#endif


#ifdef __cplusplus
extern "C"
{
#endif


extern long int check_mpi_crd_point_quad(CHARM(point) *(*)(unsigned long,
                                                           REAL,
                                                           size_t,
                                                           size_t,
                                                           MPI_Comm,
                                                           CHARM(err) *));


#ifdef __cplusplus
}
#endif


#endif
