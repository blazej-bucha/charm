/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <math.h>
#include "../prec.h"
#include "../simd/simd.h"
#include "../crd/crd_cell_isGrid.h"
#include "shs_lc_struct.h"
#include "shs_grd_lr.h"
#include "shs_max_npar.h"
/* ------------------------------------------------------------------------- */






/* Macros */
/* ------------------------------------------------------------------------- */
#undef LC_CS_init
#define LC_CS_init(x, par, i, term, offset)                                   \
    par_simd_blk = (par) * simd_blk;                                          \
    par_nfi_1par = (par) * nfi_1par + offset;                                 \
    for (size_t l = 0; l < simd_blk; l++)                                     \
    {                                                                         \
        lss = l * SIMD_SIZE;                                                  \
        CAT2(dm, term, i)[par_simd_blk + l] =                                 \
                                 ADD_R(MUL_R(lc->CAT2(a, x, i)[l], clontmp),  \
                                       MUL_R(lc->CAT2(b, x, i)[l], slontmp)); \
        STORE_R(&CAT(fi, i)[par_nfi_1par + lss],                              \
                ADD_R(LOAD_R(&CAT(fi, i)[par_nfi_1par + lss]),                \
                      CAT2(dm, term, i)[par_simd_blk + l]));                  \
    }


#undef LC_CS
#define LC_CS(x, par, i)                                                      \
    par_simd_blk = (par) * simd_blk;                                          \
    par_nfi_1par = (par) * nfi_1par + jsize_blk;                              \
    for (size_t l = 0; l < simd_blk; l++)                                     \
    {                                                                         \
        lss = l * SIMD_SIZE;                                                  \
        CAT(dm2, i)[par_simd_blk + l] =                                       \
                                SUB_R(MUL_R(cmdlon2,                          \
                                            CAT(dm1, i)[par_simd_blk + l]),   \
                                      CAT(dm0, i)[par_simd_blk + l]);         \
        STORE_R(&CAT(fi, i)[par_nfi_1par + lss],                              \
                ADD_R(LOAD_R(&CAT(fi, i)[par_nfi_1par + lss]),                \
                      CAT(dm2, i)[par_simd_blk + l]));                        \
        CAT(dm0, i)[par_simd_blk + l] = CAT(dm1, i)[par_simd_blk + l];        \
        CAT(dm1, i)[par_simd_blk + l] = CAT(dm2, i)[par_simd_blk + l];        \
    }
/* ------------------------------------------------------------------------- */






/* An internal function to perform the recursion over grid longitudes of the
 * "i"th latitude parallel (the "LR" part of the "PSLR" algorithm). */
void CHARM(shs_grd_lr)(unsigned long m,
                       REAL lon0,
                       REAL deltalon,
                       size_t nlon,
                       int grd_type,
                       int grad,
                       size_t nfi_1par,
                       CHARM(lc) *lc,
                       _Bool symm,
                       REAL *fi,
                       REAL *fi2)
{
    _Bool is_cell_grd = CHARM(crd_cell_isGrid)(grd_type);
    size_t simd_blk = is_cell_grd ? 1 : SIMD_BLOCK_S;
    size_t size_blk = SIMD_SIZE * simd_blk;
    size_t lss, jsize_blk;


    if (is_cell_grd)
    {
        /* Before applying the PSLR algorithm from Balmino et al. (2012), the
         * lumped coefficients need to be multiplied by additional terms, which
         * stem from the integration in the longitudinal direction */


        REAL m2sm2dl;
        if (m == 0)
            m2sm2dl = deltalon;
        else
        {
            REAL m2 = PREC(2.0) / (REAL)m;
            m2sm2dl = m2 * SIN(deltalon / m2);
        }


        /* With cells, "a", "b", "a2" and "b2" implicitly assume that
         * "SIMD_BLOCK_S" is "1". */
        lc->a[0] = MUL_R(lc->a[0], SET1_R(m2sm2dl));
        lc->b[0] = MUL_R(lc->b[0], SET1_R(m2sm2dl)); /* Remember that "b = 0.0"
                                          * for "m == 0.0", so it can safely be
                                          * multiplied by "m2sm2dl
                                          * * = deltalon" */

        if (symm)
        {
            lc->a2[0] = MUL_R(lc->a2[0], SET1_R(m2sm2dl));
            lc->b2[0] = MUL_R(lc->b2[0], SET1_R(m2sm2dl));  /* Remember that
                                                 * "b2 = 0.0" for "m == 0.0",
                                                 * so it can safely be
                                                 * multiplied * by "m2sm2dl
                                                 * = deltalon" */
        }
    }


    /* The PSLR algorithm from Balmino et al. (2012) */
    /* --------------------------------------------------------------------- */
    /* The first longitude point/cell */
    /* ..................................................................... */
    REAL       lontmp = (REAL)m * lon0;
    REAL_SIMD clontmp = SET1_R(COS(lontmp));
    REAL_SIMD slontmp = SET1_R(SIN(lontmp));
    REAL_SIMD dm0[SHS_MAX_NPAR * SIMD_BLOCK_S];
    size_t par_simd_blk, par_nfi_1par;
    LC_CS_init( , 0, , 0, 0);
    if (grad > 0)
    {
        LC_CS_init(r, 1, , 0, 0);
        LC_CS_init(p, 2, , 0, 0);
    }
    if (grad > 1)
    {
        LC_CS_init(rr, 3, , 0, 0);
        LC_CS_init(rp, 4, , 0, 0);
        LC_CS_init(pp, 5, , 0, 0);
    }


    REAL_SIMD dm02[SHS_MAX_NPAR * SIMD_BLOCK_S];
    if (symm)
    {
        LC_CS_init( , 0, 2, 0, 0);
        if (grad > 0)
        {
            LC_CS_init(r, 1, 2, 0, 0);
            LC_CS_init(p, 2, 2, 0, 0);
        }
        if (grad > 1)
        {
            LC_CS_init(rr, 3, 2, 0, 0);
            LC_CS_init(rp, 4, 2, 0, 0);
            LC_CS_init(pp, 5, 2, 0, 0);
        }
    }


    if (nlon == 1)
        return;
    /* ..................................................................... */


    /* The second longitude point/cell */
    /* ..................................................................... */
     lontmp = (REAL)m * (lon0 + deltalon);
    clontmp = SET1_R(COS(lontmp));
    slontmp = SET1_R(SIN(lontmp));
    REAL_SIMD dm1[SHS_MAX_NPAR * SIMD_BLOCK_S];
    LC_CS_init( , 0, , 1, size_blk);
    if (grad > 0)
    {
        LC_CS_init(r, 1, , 1, size_blk);
        LC_CS_init(p, 2, , 1, size_blk);
    }
    if (grad > 1)
    {
        LC_CS_init(rr, 3, , 1, size_blk);
        LC_CS_init(rp, 4, , 1, size_blk);
        LC_CS_init(pp, 5, , 1, size_blk);
    }


    REAL_SIMD dm12[SHS_MAX_NPAR * SIMD_BLOCK_S];
    if (symm)
    {
        LC_CS_init( , 0, 2, 1, size_blk);
        if (grad > 0)
        {
            LC_CS_init(r, 1, 2, 1, size_blk);
            LC_CS_init(p, 2, 2, 1, size_blk);
        }
        if (grad > 1)
        {
            LC_CS_init(rr, 3, 2, 1, size_blk);
            LC_CS_init(rp, 4, 2, 1, size_blk);
            LC_CS_init(pp, 5, 2, 1, size_blk);
        }
    }


    if (nlon == 2)
        return;
    /* ..................................................................... */


    /* The third and all the remaining longitude points/cells */
    /* ..................................................................... */
    REAL_SIMD cmdlon2 = SET1_R(PREC(2.0) * COS((REAL)m * deltalon));
    REAL_SIMD dm2[SHS_MAX_NPAR * SIMD_BLOCK_S];
    for (size_t j = 2; j < nlon; j++)
    {
        jsize_blk = j * size_blk;
        LC_CS( , 0, );
        if (grad > 0)
        {
            LC_CS(r, 1, );
            LC_CS(p, 2, );
        }
        if (grad > 1)
        {
            LC_CS(rr, 3, );
            LC_CS(rp, 4, );
            LC_CS(pp, 5, );
        }
    }


    REAL_SIMD dm22[SHS_MAX_NPAR * SIMD_BLOCK_S];
    if (symm)
    {
        for (size_t j = 2; j < nlon; j++)
        {
            jsize_blk = j * size_blk;
            LC_CS( , 0, 2);
            if (grad > 0)
            {
                LC_CS(r, 1, 2);
                LC_CS(p, 2, 2);
            }
            if (grad > 1)
            {
                LC_CS(rr, 3, 2);
                LC_CS(rp, 4, 2);
                LC_CS(pp, 5, 2);
            }
        }
    }


    return;
    /* ..................................................................... */
    /* --------------------------------------------------------------------- */


}
