/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
/* ------------------------------------------------------------------------- */






void CHARM(shc_free)(CHARM(shc) *shcs)
{
    if (shcs == NULL)
        return;


    if (shcs->owner)
    {
#if HAVE_MPI
        if (shcs->local_nchunk > 0)
        {
            /* Get the spherical harmonic order, to which the blocks of "cnm"
             * and "snm" cofficients are attached */
            size_t m_free = shcs->local_order[0];
            free(shcs->c[m_free]);
            free(shcs->s[m_free]);
        }
#else
        free(shcs->c[0]);
        free(shcs->s[0]);
#endif
    }
#if HAVE_MPI
    /* "local_chunk" is always deallocated regardless of "shcs->owner" and
     * "shcs->distributed" */
    free(shcs->local_order);
#endif
    free(shcs->c);
    free(shcs->s);
    free(shcs);


    return;
}
