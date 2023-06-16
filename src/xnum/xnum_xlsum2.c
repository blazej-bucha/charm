/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include "../prec.h"
#include "xnum_xnorm.h"
/* ------------------------------------------------------------------------- */






void CHARM(xnum_xlsum2)(REAL f, REAL x, REAL g, REAL y, REAL *z,
                        int ix, int iy, int *iz)
{
    int id = ix - iy;


    if (id == 0)
    {
        *z  = f * x + g * y;
        *iz = ix;
    }
    else if (id == 1)
    {
        *z  = f * x + g * (y * BIGI);
        *iz = ix;
    }
    else if (id == -1)
    {
        *z  = g * y + f * (x * BIGI);
        *iz = iy;
    }
    else if (id > 1)
    {
        *z  = f * x;
        *iz = ix;
    }
    else
    {
        *z  = g * y;
        *iz = iy;
    }


    CHARM(xnum_xnorm)(z, iz);


    return;
}
