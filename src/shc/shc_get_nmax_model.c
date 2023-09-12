/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
/* ------------------------------------------------------------------------- */






/* Private function to return "CHARM_SHC_NMAX_MODEL".  The function is used
 * only with the Python wrapper and should not be used with CHarm, so is not
 * a part of API. */
unsigned long CHARM(shc_get_nmax_model)(void)
{
    return CHARM_SHC_NMAX_MODEL;
}
