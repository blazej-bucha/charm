/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "../simd/simd.h"
#include "shs_check_symm_simd.h"
/* ------------------------------------------------------------------------- */






/* An internal function to determine whether or not the symmetry of Legendre
 * functions has to be applied for a given SIMD vector (or a scalar, depending
 * on the "SIMD_REAL" macro). */
_Bool CHARM(shs_check_symm_simd)(REAL_SIMD v)
{
    if (!MASK_TRUE_ALL((EQ_R(v, SET_ZERO_R))))
        return 1;
    else
        return 0;
}

