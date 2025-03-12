/* If "CAP_SGFM_COMPILE" is zero, "gfm_sgfm_cap_lateral.c" will not compile any
 * of the sgfm variants (global, cap-modified, constant, 3D, etc.).  This is
 * mainly to ensure easy compilation of CHarm outside the official installation
 * mechanism. */
#undef CAP_SGFM_COMPILE
#define CAP_SGFM_COMPILE 1


/* Compile the Python wrapper */
#undef PYWRAP
#define PYWRAP 1


/* Not a typo, c-file is needed here */
#include "gfm_sgfm_cap_lateral.c"
