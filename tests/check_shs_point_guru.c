#undef GURU
#define GURU (1)


/* The function in "check_shs_point_all.c" will compile only if "COMPILE_SHS"
 * is "1". */
#undef COMPILE_SHS
#define COMPILE_SHS (1)


#include "check_shs_point_all.c"
