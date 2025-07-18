/* This header file is not a part of API. */


#ifndef __SHS_POINT_KERNELS_H__
#define __SHS_POINT_KERNELS_H__


#include <config.h>
#include "../prec.h"
#include "../simd/simd.h"
#include "../shc/shc_block_struct.h"
#include "shs_lc_struct.h"


#ifdef __cplusplus
extern "C"
{
#endif


#undef IO_PARS
#define IO_PARS (unsigned long nmax,                                          \
                 unsigned long m,                                             \
                 const CHARM(shc_block) *shcs_block,                          \
                 _Bool is_ratio_one,                                          \
                 const REAL *anm,                                             \
                 const REAL *bnm,                                             \
                 const REAL *enm,                                             \
                 const REAL_SIMD *t,                                          \
                 const REAL_SIMD *u,                                          \
                 const REAL *ps,                                              \
                 const INT *ips,                                              \
                 REAL_SIMD *rpows,                                            \
                 REAL_SIMD *rpows2,                                           \
                 const REAL_SIMD *symm_simd,                                  \
                 unsigned dorder,                                             \
                 CHARM(lc) *lc)


extern void CHARM(shs_point_kernel_dr0_dlat0_dlon0) IO_PARS;
extern void CHARM(shs_point_kernel_dr1_dlat0_dlon0) IO_PARS;
extern void CHARM(shs_point_kernel_dr2_dlat0_dlon0) IO_PARS;
extern void CHARM(shs_point_kernel_dr0_dlat1_dlon0) IO_PARS;
extern void CHARM(shs_point_kernel_dr0_dlat2_dlon0) IO_PARS;
extern void CHARM(shs_point_kernel_dr0_dlat0_dlon1) IO_PARS;
extern void CHARM(shs_point_kernel_dr0_dlat0_dlon2) IO_PARS;
extern void CHARM(shs_point_kernel_dr1_dlat1_dlon0) IO_PARS;
extern void CHARM(shs_point_kernel_dr1_dlat0_dlon1) IO_PARS;
extern void CHARM(shs_point_kernel_dr0_dlat1_dlon1) IO_PARS;
extern void CHARM(shs_point_kernel_grad1) IO_PARS;
extern void CHARM(shs_point_kernel_grad2) IO_PARS;


#ifdef __cplusplus
}
#endif


#endif
