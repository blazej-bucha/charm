#undef COMPILE_KERNELS
#define COMPILE_KERNELS 1


#include "shs_point_gradn.h"


#undef DR
#define DR 0
#undef DLAT
#define DLAT 0
#undef DLON
#define DLON 0
#include "shs_point_kernel.c"


#undef DR
#define DR 1
#undef DLAT
#define DLAT 0
#undef DLON
#define DLON 0
#include "shs_point_kernel.c"


#undef DR
#define DR 2
#undef DLAT
#define DLAT 0
#undef DLON
#define DLON 0
#include "shs_point_kernel.c"


#undef DR
#define DR 0
#undef DLAT
#define DLAT 1
#undef DLON
#define DLON 0
#include "shs_point_kernel.c"


#undef DR
#define DR 0
#undef DLAT
#define DLAT 2
#undef DLON
#define DLON 0
#include "shs_point_kernel.c"


#undef DR
#define DR 0
#undef DLAT
#define DLAT 0
#undef DLON
#define DLON 1
#include "shs_point_kernel.c"


#undef DR
#define DR 0
#undef DLAT
#define DLAT 0
#undef DLON
#define DLON 2
#include "shs_point_kernel.c"


#undef DR
#define DR 1
#undef DLAT
#define DLAT 1
#undef DLON
#define DLON 0
#include "shs_point_kernel.c"


#undef DR
#define DR 1
#undef DLAT
#define DLAT 0
#undef DLON
#define DLON 1
#include "shs_point_kernel.c"


#undef DR
#define DR 0
#undef DLAT
#define DLAT 1
#undef DLON
#define DLON 1
#include "shs_point_kernel.c"


#undef DR
#define DR GRAD_1
#undef DLAT
#define DLAT GRAD_1
#undef DLON
#define DLON GRAD_1
#include "shs_point_kernel.c"


#undef DR
#define DR GRAD_2
#undef DLAT
#define DLAT GRAD_2
#undef DLON
#define DLON GRAD_2
#include "shs_point_kernel.c"
