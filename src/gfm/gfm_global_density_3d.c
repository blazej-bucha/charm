/* "SGFM_RHO_CONST" is used within the source file "gfm_sgfm.c", which is
 * included below */
#undef SGFM_RHO_CONST
#define SGFM_RHO_CONST 0


/* Global spectral gravity forward modelling */
#undef SGFM_GLOBAL
#define SGFM_GLOBAL 1


/* If "SGFM_COMPILE" is zero, "gfm_sgfm.c" will not compile any of the sgfm
 * variants (global, cap-modified, constant, 3D, etc.).  This is mainly to
 * ensure easy compilation of CHarm outside the official installation
 * mechanism. */
#undef SGFM_COMPILE
#define SGFM_COMPILE 1


#undef GFM_SGFM_COMPILE
#define GFM_SGFM_COMPILE 1
/* Not a typo, c-file is needed here */
#include "gfm_sgfm.c"


