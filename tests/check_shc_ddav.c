/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../src/prec.h"
#include "cmp_arrays.h"
#include "parameters.h"
#include "validate.h"
/* ------------------------------------------------------------------------- */






long int check_shc_ddav(void (*shc_ddav)(const CHARM(shc) *,
                                         const CHARM(shc) *,
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
    REAL *fref = (REAL *)calloc(NMAX + 1, sizeof(REAL));
    if (fref == NULL)
    {
        fprintf(stderr, "Failed to initialize an array of degree "
                        "variances");
        exit(CHARM_FAILURE);
    }
    REAL *f    = NULL;


    for (unsigned long nmax = 0; nmax <= NMAX; nmax++)
    {
        f = (REAL *)malloc((nmax + 1) * sizeof(REAL));
        if (f == NULL)
        {
            fprintf(stderr, "Failed to initialize an array of degree "
                            "variances");
            exit(CHARM_FAILURE);
        }
        shc_ddav(shcs, shcs, nmax, f, err);
        CHARM(err_handler)(err, 1);


        e += cmp_arrays(f, fref, nmax + 1, PREC(10.0) * CHARM(glob_threshold));


        free(f);
    }


    CHARM(err_free)(err);
    CHARM(shc_free)(shcs);
    free(fref);


    return e;
}
