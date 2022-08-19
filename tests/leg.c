/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../src/prec.h"
#include "cmp_arrays.h"
#include "validate.h"
#include "parameters.h"
/* ------------------------------------------------------------------------- */






/* Tests the "leg" module */
int leg(void)
{
    /* --------------------------------------------------------------------- */
    /* Maximum degree to test the Fourier coefficients of Legendre functions */
    unsigned long nmax = 10;


    /* Error structure */
    CHARM(err) *err = CHARM(err_init)();
    if (err == NULL)
    {
        printf("Failed to initialize the \"err\" structure.\n");
        exit(1);
    }


    int errnum = 0;
    char file[NSTR];
    CHARM(pnmj) *pnmj;
    /* --------------------------------------------------------------------- */






    /* Loop over the two ordering schemes of the Fourier coefficients */
    for (int o = 0; o < 2; o++)
    {
        if (o == 0)
            printf("    Testing the \"CHARM_LEG_PNMJ_ORDER_MNJ\" ordering "
                   "scheme...\n");
        else
            printf("    Testing the \"CHARM_LEG_PNMJ_ORDER_MJN\" ordering "
                   "scheme...\n");



        for (unsigned long nmax2 = 0; nmax2 <= nmax; nmax2++)
        {
            /* Initialize a structure to store the Fourier coefficients */
            pnmj = CHARM(leg_pnmj_init)(nmax2,
                                        (o == 0) ? CHARM_LEG_PNMJ_ORDER_MNJ :
                                                   CHARM_LEG_PNMJ_ORDER_MJN);
            if (pnmj == NULL)
            {
                fprintf(stderr, "Failed to initialize the \"pnmj\" "
                                "structure.\n");
                exit(CHARM_FAILURE);
            }


            /* Compute the Fourier coefficients */
            CHARM(leg_pnmj_coeffs)(pnmj, nmax2, err);
            CHARM(err_handler)(err, 1);


            if (o == 0)
            {
                for (unsigned long m = 0; m <= nmax2; m++)
                    for (unsigned long n = m; n <= nmax2; n++)
                    {
                        sprintf(file, "%s/leg_pnmj_nx%lu_mnj_m%lu_n%lu%s",
                                FOLDER, nmax, m, n, FTYPE);


                        errnum += validate(file, pnmj->pnmj[m][n - m], (n / 2),
                                           PREC(10.0) * CHARM(glob_threshold));
                    }
            }
            else
            {
                for (unsigned long m = 0; m <= nmax2; m++)
                    for (unsigned long j = 0; j <= (nmax2 / 2); j++)
                    {
                        sprintf(file, "%s/leg_pnmj_nx%lu_mjn_m%lu_j%lu%s",
                                FOLDER, nmax, m, j, FTYPE);


                        errnum += validate(file, pnmj->pnmj[m][j],
                                           nmax2 - CHARM_MAX(m, 2 * j) + 1,
                                           PREC(10.0) * CHARM(glob_threshold));
                    }
            }


            CHARM(leg_pnmj_free)(pnmj);
        }
    }






    CHARM(err_free)(err);


    return errnum;
}

