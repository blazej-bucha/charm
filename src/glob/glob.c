/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include "../prec.h"
#include "../simd/simd.h"
/* ------------------------------------------------------------------------- */






/* Define global variables */
/* ------------------------------------------------------------------------- */
/* Thresholds for floating points comparisons */
REAL CHARM(glob_threshold)  = PREC(100.0)    * EPS;
REAL CHARM(glob_threshold2) = PREC(100000.0) * EPS;


/* Polar optimization (negative "glob_polar_optimization_a2" means no polar
 * optimization by default) */
unsigned long CHARM(glob_polar_optimization_a1) = 100;
REAL CHARM(glob_polar_optimization_a2)          = PREC(-1.0);


#if HAVE_MPI
unsigned long CHARM(glob_shc_block_nmax_multiplier) = 1000;
size_t CHARM(glob_sha_block_lat_multiplier) = SIMD_BLOCK_A;
size_t CHARM(glob_shs_block_lat_multiplier) = SIMD_BLOCK_S;
#endif
/* ------------------------------------------------------------------------- */
