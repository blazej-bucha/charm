#include <stdio.h>
#include <stdlib.h>
#include <charm/charm.h>


int main(void)
{
    /* Type of the first spherical harmonic function, its harmonic degree and
     * order, respectively */
    _Bool i1 = 0;  /* The "cos" spherical harmonic function */
    unsigned long n1 = 287;
    unsigned long m1 = 122;


    /* Type of the second spherical harmonic function, its harmonic degree and
     * order, respectively */
    _Bool i2 = 1;  /* The "sin" spherical harmonic function */
    unsigned long n2 = 34;
    unsigned long m2 = 9;


    /* Minimum and maximum co-latitudes of the integration domain */
    double cltmin = 0.1;
    double cltmax = 0.9;


    /* Minimum and maximum longitudes of the integration domain */
    double lonmin = 0.5;
    double lonmax = 0.6;


    /* Initialize the Fourier coefficients of Legendre functions */
    unsigned long nmax = CHARM_MAX(n1, n2);
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


    /* Compute the integral of a product of two spherical harmonics */
    double iy = charm_integ_yi1n1m1yi2n2m2(cltmin, cltmax, lonmin, lonmax,
                                           i1, n1, m1,
                                           i2, n2, m2, pnmj, err);
    charm_err_handler(err, 1);


    /* Print the value of the integral */
    printf("The integral of the product of two spherical harmonics "
           "i1 = %d, n1 = %lu, m1 = %lu, i2 = %d, n2 = %lu, m2 = %lu "
           "is %0.16e\n", i1, n1, m1, i2, n2, m2, iy);


    /* Now compute the integral of a product of two Legendre functions over a
     * restricted domain */
    double ip = charm_integ_pn1m1pn2m2(cltmin, cltmax, n1, m1, n2, m2, pnmj,
                                       err);
    charm_err_handler(err, 1);


    printf("The integral of the product of two Legendre functions "
           "n1 = %lu, m1 = %lu, n2 = %lu, m2 = %lu "
           "is %0.16e\n", n1, m1, n2, m2, ip);


    /* Free the heap memory */
    charm_err_free(err);
    charm_leg_pnmj_free(pnmj);


    exit(CHARM_SUCCESS);
}
