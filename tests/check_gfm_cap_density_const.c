/* "0" for 3D density */
#undef GFM_DENSITY
#define GFM_DENSITY (2)


/* "1" for global GFM, "0" for cap GFM */
#undef GFM_GLOBAL
#define GFM_GLOBAL (0)


/* The function in "check_gfm_global_all.c" will compile only if
 * "COMPILE_GFM_GLOBAL" is "1". */
#undef COMPILE_GFM
#define COMPILE_GFM (1)


#include "check_gfm_sgfm_all.c"
