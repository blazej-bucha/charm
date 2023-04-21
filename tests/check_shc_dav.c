/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../src/prec.h"
#include "parameters.h"
#include "validate.h"
/* ------------------------------------------------------------------------- */






long int check_shc_dav(void (*shc_dav)(const CHARM(shc) *,
                                       unsigned long,
                                       REAL *,
                                       CHARM(err) *))
{
    CHARM(err) *err = CHARM(err_init)();
    if (err == NULL)
    {
        printf("Failed to initialize the \"err\" structure.\n");
        exit(CHARM_FAILURE);
    }


    CHARM(shc) *shcs = CHARM(shc_calloc)(SHCS_NMAX_TOPO, PREC(1.0), PREC(1.0));
    if (shcs == NULL)
    {
        fprintf(stderr, "Failed to initialize a \"shc\" structure.\n");
        exit(CHARM_FAILURE);
    }


    CHARM(shc_read_mtx)(SHCS_IN_PATH_TOPO_MTX, SHCS_NMAX_TOPO, shcs, err);
    CHARM(err_handler)(err, 1);


    long int e = 0;
    REAL *f    = NULL;


    char file[NSTR_LONG];
    char av[NSTR_SHORT];
    if (shc_dav == CHARM(shc_da))
        sprintf(av, "a");
    else if (shc_dav == CHARM(shc_dv))
        sprintf(av, "v");
    sprintf(file, "%s/shc_nx%lu_d%s%s", FOLDER, (unsigned long)NMAX, av,
            FTYPE);


    for (unsigned long nmax = 0; nmax <= NMAX; nmax++)
    {
        f = (REAL *)malloc((nmax + 1) * sizeof(REAL));
        if (f == NULL)
        {
            fprintf(stderr, "Failed to initialize an array of degree "
                            "variances");
            exit(CHARM_FAILURE);
        }
        shc_dav(shcs, nmax, f, err);
        CHARM(err_handler)(err, 1);


        e += validate(file, f, nmax + 1, PREC(10.0) * CHARM(glob_threshold));


        free(f);
    }


    CHARM(err_free)(err);
    CHARM(shc_free)(shcs);


    return e;
}
