/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#if HAVE_MPI
#   include <mpi.h>
#endif
#include "parameters.h"
#include "../src/prec.h"
#include "check_func.h"
#include "check_outcome.h"
#include "check_mpi_shc_alloc.h"
#include "check_mpi_shc_init.h"
#include "check_mpi_shc_local_ncs.h"
#include "check_mpi_crd_point_init.h"
#include "check_mpi_crd_point_quad.h"
#include "check_mpi_crd_point_alloc.h"
#include "check_mpi_shs_point.h"
#include "check_mpi_sha_point.h"
#include "module_mpi.h"
/* ------------------------------------------------------------------------- */






long int module_mpi(void)
{
    long int e    = 0;
    long int esum = 0;


    check_func("mpi_shc_malloc");
    e = check_mpi_shc_alloc(CHARM(mpi_shc_malloc));
    check_outcome(e);
    esum += e;


    check_func("mpi_shc_calloc");
    e = check_mpi_shc_alloc(CHARM(mpi_shc_calloc));
    check_outcome(e);
    esum += e;


    check_func("mpi_shc_local_ncs");
    e = check_mpi_shc_local_ncs();
    check_outcome(e);
    esum += e;


    check_func("mpi_shc_init");
    e = check_mpi_shc_init();
    check_outcome(e);
    esum += e;


    check_func("mpi_crd_point_init");
    e = check_mpi_crd_point_init();
    check_outcome(e);
    esum += e;


    check_func("mpi_crd_point_malloc");
    e = check_mpi_crd_point_alloc(CHARM(mpi_crd_point_malloc));
    check_outcome(e);
    esum += e;


    check_func("mpi_crd_point_calloc");
    e = check_mpi_crd_point_alloc(CHARM(mpi_crd_point_calloc));
    check_outcome(e);
    esum += e;


    check_func("mpi_crd_point_gl");
    e = check_mpi_crd_point_quad(CHARM(mpi_crd_point_gl));
    check_outcome(e);
    esum += e;


    check_func("mpi_crd_point_dh1");
    e = check_mpi_crd_point_quad(CHARM(mpi_crd_point_dh1));
    check_outcome(e);
    esum += e;


    check_func("mpi_crd_point_dh2");
    e = check_mpi_crd_point_quad(CHARM(mpi_crd_point_dh2));
    check_outcome(e);
    esum += e;


    check_func("shs_point");
    e = check_mpi_shs_point();
    check_outcome(e);
    esum += e;


    check_func("sha_point");
    e = check_mpi_sha_point();
    check_outcome(e);
    esum += e;


#if HAVE_MPI
    int initialized;
    MPI_Initialized(&initialized);
    if (initialized)
    {
        long int esum_sum;
        MPI_Allreduce(&esum, &esum_sum, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
        esum = esum_sum;
    }
#endif


    return esum;
}
