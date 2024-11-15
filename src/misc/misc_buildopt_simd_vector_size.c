/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "../simd/simd.h"
/* ------------------------------------------------------------------------- */






int CHARM(misc_buildopt_simd_vector_size)(void)
{
    return SIMD_SIZE;
}
