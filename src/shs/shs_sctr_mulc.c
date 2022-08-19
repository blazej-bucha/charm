/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "../simd/simd.h"
/* ------------------------------------------------------------------------- */






void CHARM(shs_sctr_mulc)(size_t i, size_t n, REAL_SIMD c, REAL_SIMD tmp,
                          REAL *tmpv, REAL_SIMD fi, REAL *f)
{
    size_t ipv;


    tmp = MUL_R(c, fi);
    if ((i + SIMD_SIZE) <= n)
        STOREU_R(&f[i], tmp);
    else
    {
        STORE_R(&tmpv[0], tmp);
        for (size_t v = 0; v < SIMD_SIZE; v++)
        {
            ipv = i + v;
            if (ipv < n)
                f[ipv] = tmpv[v];
            else
                continue;
        }
    }
}
