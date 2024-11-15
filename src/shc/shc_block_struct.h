/* This header file is not a part of API. */


#ifndef __SHC_BLOCK_STRUCT_H__
#define __SHC_BLOCK_STRUCT_H__


#include <config.h>
#if HAVE_MPI
#   include <mpi.h>
#endif
#include "../prec.h"


/* Structure to store a block of spherical harmonic coefficients */
typedef struct
{
    /* C_{mfirst, mfirst}, C_{mfirst + 1, mfirst}, ..., C_{nmax, mfirst}, ...,
     * C_{mfirst + 1, mfirst + 1}, C_{mfirst + 2, mfirst + 1}, ..., C{nmax,
     * mlast} */
    REAL *c;


    /* The same as "c" but for "s" */
    REAL *s;


    /* Rank of the process that holds the current chunk of coefficients */
    int root;


    /* The maximum harmonic degree of the "charm_shc" structure, from which the
     * block was created */
    unsigned long nmax;


    /* Minimum order of the block */
    unsigned long mfirst;


    /* Maximum order of the block */
    unsigned long mlast;


    /* Total number of coefficients associated with "c" */
    size_t nc;


    /* Total number of coefficients associated with "s" */
    size_t ns;


    /* The maximum number of coefficients that can be stored in "c" and "s".
     * "nc" and "ns" must not be larger than "ncs_max" */
    size_t ncs_max;


    /* "1" if owns the memory associated with "c" and "s" and "0" otherwise.
     * "owner" is not used with "have_m_all", which is always owned, hence
     * deallocated, too. */
    _Bool owner;


#if HAVE_MPI
    /* Temporary array used to identify the MPI process, which stores
     * coefficients of some particular order "m" */
    _Bool *have_m_all;


    /* Takes the same value as "shc_charm->distributed", from which "shc_block"
     * is derived */
    _Bool distributed;


    /* MPI Communicator */
    MPI_Comm comm;
#endif
} CHARM(shc_block);


#endif
