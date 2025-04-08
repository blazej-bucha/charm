/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../src/prec.h"
#include "parameters.h"
#include "error_messages.h"
#include "cmp_arrays.h"
#ifdef GENREF
#   include "array2file.h"
#else
#   include "validate.h"
#endif
#include "check_shc_rescale.h"
/* ------------------------------------------------------------------------- */






/* Checks the routine to rescale spherical harmonic coefficients by
 *
 * 1) loading a reference set of coefficients,
 *
 * 2) synthesizing the signal from the reference coefficients,
 *
 * 3) rescaling the coefficients,
 *
 * 4) synthesizing the signal from the rescaled coefficients, and
 *
 * 5) comparing the two signals.  */
long int check_shc_rescale(void)
{
    CHARM(err) *err = CHARM(err_init)();
    if (err == NULL)
    {
        fprintf(stderr, "%s", ERR_MSG_ERR);
        exit(CHARM_FAILURE);
    }


    CHARM(shc) *shcs = CHARM(shc_calloc)(SHCS_NMAX_POT, PREC(1.0), PREC(1.0));
    if (shcs == NULL)
    {
        fprintf(stderr, "%s", ERR_MSG_SHC);
        exit(CHARM_FAILURE);
    }


    CHARM(shc_read_mtx)(SHCS_IN_PATH_POT_MTX, SHCS_NMAX_POT, shcs, err);
    CHARM(err_handler)(err, 1);


    CHARM(shc_rescale)(shcs, shcs->mu * SHCS_RESCALE_MU_FACTOR,
                       shcs->r * SHCS_RESCALE_R_FACTOR, err);
    CHARM(err_handler)(err, 1);


    long int e = 0;
    char filec[NSTR_LONG];
    char files[NSTR_LONG];
    snprintf(filec, NSTR_LONG, "%s/shc_rescale_c%s", FOLDER, FTYPE);
    snprintf(files, NSTR_LONG, "%s/shc_rescale_s%s", FOLDER, FTYPE);
#ifdef GENREF
    e += array2file(filec, shcs->c[0], shcs->nc);
    e += array2file(files, shcs->s[0], shcs->ns);
#else
    e += validate(filec, shcs->c[0], shcs->nc, CHARM(glob_threshold));
    e += validate(files, shcs->s[0], shcs->ns, CHARM(glob_threshold));
#endif


    CHARM(err_free)(err);
    CHARM(shc_free)(shcs);


    return e;
}
