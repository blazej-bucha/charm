/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../src/prec.h"
#include "parameters.h"
#include "error_messages.h"
#ifdef GENREF
#   include "array2file.h"
#else
#   include "validate.h"
#endif
#include "check_shc_dav.h"
/* ------------------------------------------------------------------------- */






long int check_shc_dav(void (*shc_dav)(const CHARM(shc) *,
                                       unsigned long,
                                       REAL *,
                                       CHARM(err) *))
{
    CHARM(err) *err = CHARM(err_init)();
    if (err == NULL)
    {
        fprintf(stderr, ERR_MSG_ERR);
        exit(CHARM_FAILURE);
    }


    CHARM(shc) *shcs = CHARM(shc_calloc)(SHCS_NMAX_TOPO, PREC(1.0), PREC(1.0));
    if (shcs == NULL)
    {
        fprintf(stderr, ERR_MSG_SHC);
        exit(CHARM_FAILURE);
    }


    CHARM(shc_read_mtx)(SHCS_IN_PATH_TOPO_MTX, SHCS_NMAX_TOPO, shcs, err);
    CHARM(err_handler)(err, 1);


    long int e = 0;
    REAL *f    = NULL;


    char file[NSTR_LONG];
    char av[NSTR_SHORT];
    if (shc_dav == CHARM(shc_da))
        snprintf(av, NSTR_SHORT, "a");
    else if (shc_dav == CHARM(shc_dv))
        snprintf(av, NSTR_SHORT, "v");
    snprintf(file, NSTR_LONG, "%s/shc_nx%lu_d%s%s", FOLDER,
             (unsigned long)NMAX, av, FTYPE);


    for (unsigned long nmax = 0; nmax <= NMAX; nmax++)
    {
        f = (REAL *)malloc((nmax + 1) * sizeof(REAL));
        if (f == NULL)
        {
            fprintf(stderr, CHARM_ERR_MALLOC_FAILURE"\n");
            exit(CHARM_FAILURE);
        }
        shc_dav(shcs, nmax, f, err);
        CHARM(err_handler)(err, 1);


#ifdef GENREF
        e += array2file(file, f, nmax + 1);
#else
        e += validate(file, f, nmax + 1, PREC(10.0) * CHARM(glob_threshold));
#endif


        free(f);
    }


    CHARM(err_free)(err);
    CHARM(shc_free)(shcs);


    return e;
}
