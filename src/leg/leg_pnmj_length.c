/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdlib.h>
#include "../prec.h"
#include "leg_pnmj_length.h"
/* ------------------------------------------------------------------------- */






size_t CHARM(leg_pnmj_length)(unsigned long nmax)
{
    size_t npnmj = 0;
    unsigned long jmax;
    for (unsigned long n = 0; n <= nmax; n++)
    {
        jmax = (n / 2);
        for (unsigned long m = 0; m <= n; m++)
            for (unsigned long j = 0; j <= jmax; j++)
                npnmj += 1;
    }


    return npnmj;
}
