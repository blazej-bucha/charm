/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include "../prec.h"
#include "shc_read_nmax_only.h"
/* ------------------------------------------------------------------------- */






_Bool CHARM(shc_read_nmax_only)(unsigned long nmax, CHARM(shc) *shcs)
{
    return (nmax == CHARM_SHC_NMAX_MODEL) && (shcs == NULL);
}
