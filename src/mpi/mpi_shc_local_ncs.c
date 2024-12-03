/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "../shc/shc_check_chunk_orders.h"
#include "../shc/shc_local_ncs.h"
#include "../err/err_propagate.h"
#include "mpi_err_gather.h"
#include "mpi_err_isdistributed.h"
/* ------------------------------------------------------------------------- */






size_t CHARM(mpi_shc_local_ncs)(unsigned long nmax,
                                size_t nchunk,
                                const unsigned long *order,
                                CHARM(err) *err)
{
    if (!CHARM(mpi_err_isdistributed)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        goto BARRIER;
    }


    CHARM(shc_check_chunk_orders)(nmax, nchunk, order, 0, err);
    if (!CHARM(err_isempty)(err))
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);


    CHARM(mpi_err_gather)(err);


BARRIER:
    if (!CHARM(mpi_err_isempty)(err))
        return 0;


    return CHARM(shc_local_ncs)(nmax, nchunk, order);
}
