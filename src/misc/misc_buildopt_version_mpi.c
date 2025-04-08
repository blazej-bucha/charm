/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#if HAVE_MPI
#   include <mpi.h>
#endif
#include "../prec.h"
#include "misc_buildopt.h"
/* ------------------------------------------------------------------------- */






void CHARM(misc_buildopt_version_mpi)(int *major_header,
                                      int *minor_header,
                                      int *major_lib,
                                      int *minor_lib)
{
#if HAVE_MPI
    *major_header = MPI_VERSION;
    *minor_header = MPI_SUBVERSION;
    MPI_Get_version(major_lib, minor_lib);
#else
    *major_header = *minor_header = *major_lib = *minor_lib = LIB_NA_VAL;
#endif


    return;
}
