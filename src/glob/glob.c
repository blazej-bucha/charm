/* Header files */
/* ------------------------------------------------------------------------- */
#include "../prec.h"
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
/* ------------------------------------------------------------------------- */
