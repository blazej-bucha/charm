/* This header file is not a part of API. */


#ifndef __SHS_POINT_GRADN_H__
#define __SHS_POINT_GRADN_H__


/* ......................................................................... */
/* Zero-order derivative */
#undef GRAD_0
#define GRAD_0 (0)


/* To indicate the synthesis of all first-order derivatives */
#undef GRAD_1
#define GRAD_1 (-1)


/* To indicate the synthesis of all second-order derivatives */
#undef GRAD_2
#define GRAD_2 (-2)


#undef GRAD_LL
#define GRAD_LL (0)


#undef GRAD_LR
#define GRAD_LR (1)


#undef GRAD_LP
#define GRAD_LP (2)


#undef GRAD_RR
#define GRAD_RR (3)


#undef GRAD_RP
#define GRAD_RP (4)


#undef GRAD_PP
#define GRAD_PP (5)


#undef GRAD_P
#define GRAD_P GRAD_LP


#undef GRAD_L
#define GRAD_L GRAD_LL


#undef GRAD_R
#define GRAD_R GRAD_LR
/* ......................................................................... */


#endif
