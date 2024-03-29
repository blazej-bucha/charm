/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include "../prec.h"
#include "leg_pnmj_check_ordering.h"
/* ------------------------------------------------------------------------- */






/** Returns zero if "ordering" has a valid value of the "ordering" member of
 * the "CHARM(pnmj)" structure and non-zero otherwise. */
int CHARM(leg_pnmj_check_ordering)(int ordering)
{
    return ((ordering == CHARM_LEG_PMNJ) ||
            (ordering == CHARM_LEG_PMJN)) ? 0 : 1;
}
