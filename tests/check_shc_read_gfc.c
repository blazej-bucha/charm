/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../src/prec.h"
#include "cmp_vals.h"
#include "validate.h"
#include "array2file.h"
#include "parameters.h"
#include "error_messages.h"
#include "check_shc_read_gfc.h"
/* ------------------------------------------------------------------------- */






long int check_shc_read_gfc(void)
{
    CHARM(err) *err = CHARM(err_init)();
    if (err == NULL)
    {
        fprintf(stderr, ERR_MSG_ERR);
        exit(CHARM_FAILURE);
    }


    int long e = 0;
    char file_c[NSTR_LONG];
    char file_s[NSTR_LONG];


    CHARM(shc) *shcs = CHARM(shc_calloc)(SHCS_NMAX_POT, PREC(1.0), PREC(1.0));
    if (shcs == NULL)
    {
        fprintf(stderr, ERR_MSG_SHC);
        exit(CHARM_FAILURE);
    }


    char *stat = "static";
    char *tvg = "tvg";
    char *stat_tvg;
    for (int i = 0; i < 2; i++)
    {
        if (i == 0)
            CHARM(shc_read_gfc)(SHCS_IN_PATH_POT_GFC, SHCS_NMAX_POT, NULL,
                                shcs, err);
        else if (i == 1)
            CHARM(shc_read_gfc)(SHCS_IN_PATH_POT_GFC_TVG, SHCS_NMAX_POT,
                                TVG_EPOCH, shcs, err);
        CHARM(err_handler)(err, 1);


        if (i == 0)
            stat_tvg = stat;
        else if (i == 1)
            stat_tvg = tvg;


        snprintf(file_c, NSTR_LONG, "%s/shc_nx%lu_c_gfc_%s%s",
                 FOLDER, shcs->nmax, stat_tvg, FTYPE);
        snprintf(file_s, NSTR_LONG, "%s/shc_nx%lu_s_gfc_%s%s",
                 FOLDER, shcs->nmax, stat_tvg, FTYPE);


#ifdef GENREF
        e += array2file(file_c, shcs->c[0], shcs->nc);
        e += array2file(file_s, shcs->s[0], shcs->ns);
#else
        e += validate(file_c, shcs->c[0], shcs->nc,
                      PREC(10.0) * CHARM(glob_threshold));
        e += validate(file_s, shcs->s[0], shcs->ns,
                      PREC(10.0) * CHARM(glob_threshold));
#endif
    }
    /* --------------------------------------------------------------------- */


    /* Check reading the maximum degree of the file only */
    /* --------------------------------------------------------------------- */
    unsigned long nmax_out = CHARM(shc_read_gfc)(SHCS_IN_PATH_POT_GFC,
                                                 CHARM_SHC_NMAX_MODEL,
                                                 NULL, NULL, err);
    CHARM(err_handler)(err, 1);
    e += cmp_vals_ulong(nmax_out, SHCS_NMAX_POT);
    /* --------------------------------------------------------------------- */


    /* --------------------------------------------------------------------- */
    CHARM(shc_free)(shcs);
    CHARM(err_free)(err);


    return e;
    /* --------------------------------------------------------------------- */
}
