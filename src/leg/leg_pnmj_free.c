/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdlib.h>
#include "../prec.h"
/* ------------------------------------------------------------------------- */






void CHARM(leg_pnmj_free)(CHARM(pnmj) *pnmj)
{
    if (pnmj == NULL)
        return;


    free(pnmj->pnmj[0][0]);
    for (unsigned long m = 0; m <= pnmj->nmax; m++)
        free(pnmj->pnmj[m]);
    free(pnmj->pnmj);
    free(pnmj);


    return;
}
