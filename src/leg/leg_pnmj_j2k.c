/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include "../prec.h"
/* ------------------------------------------------------------------------- */






unsigned long CHARM(leg_pnmj_j2k)(unsigned long n, unsigned long j)
{
    return (n % 2) + 2 * j;
}
