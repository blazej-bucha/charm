/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../src/prec.h"
#include "parameters.h"
#ifdef GENREF
#   include "array2file.h"
#else
#   include "validate.h"
#endif
#include "check_leg_pnmj_coeffs.h"
/* ------------------------------------------------------------------------- */






long int check_leg_pnmj_coeffs(void)
{
    /* --------------------------------------------------------------------- */
    CHARM(err) *err = CHARM(err_init)();
    if (err == NULL)
    {
        printf("Failed to initialize the \"err\" structure.\n");
        exit(CHARM_FAILURE);
    }


    long int e = 0;
    char file[NSTR_LONG];
    CHARM(pnmj) *pnmj;
    /* --------------------------------------------------------------------- */






    /* Loop over the two ordering schemes of the Fourier coefficients:
     *
     * * "o == 0" for "CHARM_LEG_PMNJ", and
     *
     * * "o == 1" for "CHARM_LEG_PMJN". */
    /* --------------------------------------------------------------------- */
    for (int o = 0; o < 2; o++)
    {
        for (unsigned long nmax2 = 0; nmax2 <= SHCS_NMAX_POT; nmax2++)
        {
            /* Initialize a structure to store the Fourier coefficients */
            pnmj = CHARM(leg_pnmj_calloc)(nmax2,
                                          (o == 0) ? CHARM_LEG_PMNJ :
                                                     CHARM_LEG_PMJN);
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
                                FOLDER, (unsigned long)SHCS_NMAX_POT, m, n,
                                FTYPE);


#ifdef GENREF
                        e += array2file(file, pnmj->pnmj[m][n - m],
                                        (n / 2) + 1);
#else
                        e += validate(file, pnmj->pnmj[m][n - m], (n / 2) + 1,
                                      PREC(10.0) * CHARM(glob_threshold));
#endif
                    }
            }
            else
            {
                for (unsigned long m = 0; m <= nmax2; m++)
                    for (unsigned long j = 0; j <= (nmax2 / 2); j++)
                    {
                        sprintf(file, "%s/leg_pnmj_nx%lu_mjn_m%lu_j%lu%s",
                                FOLDER, (unsigned long)SHCS_NMAX_POT, m, j,
                                FTYPE);


#ifdef GENREF
                        e += array2file(file, pnmj->pnmj[m][j],
                                        nmax2 - CHARM_MAX(m, 2 * j) + 1);
#else
                        e += validate(file, pnmj->pnmj[m][j],
                                      nmax2 - CHARM_MAX(m, 2 * j) + 1,
                                      PREC(10.0) * CHARM(glob_threshold));
#endif
                    }
            }


            CHARM(leg_pnmj_free)(pnmj);
        }
    }
    /* --------------------------------------------------------------------- */






    /* --------------------------------------------------------------------- */
    CHARM(err_free)(err);


    return e;
    /* --------------------------------------------------------------------- */
}
