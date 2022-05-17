#include <stdio.h>
#include <stdlib.h>
#include <charm/charm.h>


int main(void)
{
    /* Maximum harmonic degree to compute the Fourier coefficients of Legendre
     * functions */
    unsigned long nmax = 500;


    /* Initialize a structure to store the Fourier coefficients of Legendre
     * functions */
    charm_pnmj *pnmj = charm_leg_pnmj_init(nmax, CHARM_LEG_PNMJ_ORDER_MNJ);
    if (pnmj == NULL)
    {
        fprintf(stderr, "Failed to initialize the \"charm_pnmj\" "
                        "structure.\n");
        exit(CHARM_FAILURE);
    }


    /* Initialize an "charm_err" structure */
    charm_err *err = charm_err_init();
    if (err == NULL)
    {
        fprintf(stderr, "Failed to initialize the \"charm_err\" structure.\n");
        exit(CHARM_FAILURE);
    }


    /* Compute the Fourier coefficients */
    charm_leg_pnmj_coeffs(pnmj, nmax, err);
    charm_err_handler(err, 1);


    /* Print some Fourier coefficients */
    /* Harmonic degree */
    unsigned long n = 123;
    /* Harmonic order */
    unsigned long m = 23;
    /* Wave-number-related variable */
    unsigned long j = 12;
    /* Wave-number (computed from "j"; for details, see the "charm_leg"
     * module)*/
    unsigned long k = charm_leg_pnmj_j2k(n, j);
    printf("Fourier coefficient for degree %lu, order %lu and "
           "wave-number %lu = %0.16e\n", n, m, k, pnmj->pnmj[m][n - m][j]);


    n = 360; m = 358; j = 101; k = charm_leg_pnmj_j2k(n, j);
    printf("Fourier coefficient for degree %lu, order %lu and "
           "wave-number %lu = %0.16e\n", n, m, k, pnmj->pnmj[m][n - m][j]);


    /* Free the heap memory */
    charm_leg_pnmj_free(pnmj);
    charm_err_free(err);


    printf("\nGreat, all done!\n");


    exit(CHARM_SUCCESS);
}
