#include "shs_point_gradn.h"


#undef COMPILE_GRADS
#define COMPILE_GRADS 1


#undef GRADN
#define GRADN GRAD_0
#include "shs_point_gradn.c"


#undef GRADN
#define GRADN GRAD_1
#include "shs_point_gradn.c"


#undef GRADN
#define GRADN GRAD_2
#include "shs_point_gradn.c"
