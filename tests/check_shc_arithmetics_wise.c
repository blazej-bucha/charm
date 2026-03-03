/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../src/prec.h"
#include "cmp_vals.h"
#include "cmp_arrays.h"
#include "parameters.h"
#include "array2file.h"
#include "validate.h"
#include "error_messages.h"
#include "check_shc_arithmetics_wise.h"
/* ------------------------------------------------------------------------- */






/* This file checks arithmetics with "charm_shc".
 *
 * * "shc_add", "shc_sub", "shc_mul", "shc_div" -- Arithmetics functions to
 *                                                 check.
 *
 * */
long int check_shc_arithmetics_wise(void (*shc_arithmetic_wise)(CHARM(shc) *,
                                                                const REAL *,
                                                                unsigned long,
                                                                unsigned long,
                                                                CHARM(err) *))
{
    /* --------------------------------------------------------------------- */
    CHARM(err) *err = CHARM(err_init)();
    if (err == NULL)
    {
        fprintf(stderr, "%s", ERR_MSG_ERR);
        exit(CHARM_FAILURE);
    }


    int long e = 0;
    unsigned long N = 4;  /* Must be larger than zero in this routine */


    /* Get name of the routine that is being validated */
    char func[NSTR_SHORT];
    if (shc_arithmetic_wise == CHARM(shc_mul_degree_wise))
        snprintf(func, NSTR_SHORT, "shc_mul_degree_wise");
    else if (shc_arithmetic_wise == CHARM(shc_mul_order_wise))
        snprintf(func, NSTR_SHORT, "shc_mul_order_wise");
    else if (shc_arithmetic_wise == CHARM(shc_div_degree_wise))
        snprintf(func, NSTR_SHORT, "shc_div_degree_wise");
    else if (shc_arithmetic_wise == CHARM(shc_div_order_wise))
        snprintf(func, NSTR_SHORT, "shc_div_order_wise");
    /* --------------------------------------------------------------------- */


    /* Create structures to hold spherical harmonic coefficients */
    /* --------------------------------------------------------------------- */
    CHARM(shc) *shcs = CHARM(shc_calloc)(N, PREC(1.0), PREC(1.0));
    if (shcs == NULL)
    {
        fprintf(stderr, "%s", ERR_MSG_SHC);
        exit(CHARM_FAILURE);
    }
    /* --------------------------------------------------------------------- */


    /* Initialize coefficients of "shcs" */
    /* --------------------------------------------------------------------- */
    for (unsigned long m = 0; m <= shcs->nmax; m++)
    {
        for (unsigned long n = m; n <= shcs->nmax; n++)
        {
            shcs->c[m][n - m] = PREC(1.0);
            if (m > 0)
                shcs->s[m][n - m] = PREC(2.0);
        }
    }
    /* --------------------------------------------------------------------- */


    /* Initialize the degree- and order-dependent arrays */
    /* --------------------------------------------------------------------- */
    REAL *a = (REAL *)malloc((shcs->nmax + 1) * sizeof(REAL));
    if (a == NULL)
    {
        fprintf(stderr, "%s", CHARM_ERR_MALLOC_FAILURE"\n");
        exit(CHARM_FAILURE);
    }


    for (size_t i = 0; i <= shcs->nmax; i++)
        a[i] = (REAL)(i + 2);
    /* --------------------------------------------------------------------- */


    /* --------------------------------------------------------------------- */
    unsigned long nmin = 0;
    unsigned long nmax = shcs->nmax;
    shc_arithmetic_wise(shcs, a, nmin, nmax, err);
    CHARM(err_handler)(err, 1);


    char filec[NSTR_LONG];
    char files[NSTR_LONG];
    snprintf(filec, NSTR_LONG, "%s/%s_c_nmin%lu_nmax%lu%s",
             FOLDER, func, nmin, nmax, FTYPE);
    snprintf(files, NSTR_LONG, "%s/%s_s_nmin%lu_nmax%lu%s",
             FOLDER, func, nmin, nmax, FTYPE);
#ifdef GENREF
    e += array2file(filec, shcs->c[0], shcs->nc);
    e += array2file(files, shcs->s[0], shcs->ns);
#else
    e += validate(filec, shcs->c[0], shcs->nc,
                  PREC(10.0) * CHARM(glob_threshold));
    e += validate(files, shcs->s[0], shcs->ns,
                  PREC(10.0) * CHARM(glob_threshold));
#endif


    /* Before the next test, we have to reset " shcs" */
    for (unsigned long m = 0; m <= shcs->nmax; m++)
    {
        for (unsigned long n = m; n <= shcs->nmax; n++)
        {
            shcs->c[m][n - m] = PREC(1.0);
            if (m == 0)
                shcs->s[m][n - m] = PREC(0.0);
            else
                shcs->s[m][n - m] = PREC(2.0);
        }
    }
    nmin = 1;
    nmax = 2;
    if ((shc_arithmetic_wise == CHARM(shc_mul_degree_wise)) ||
        (shc_arithmetic_wise == CHARM(shc_div_degree_wise)))
        shc_arithmetic_wise(shcs, &a[nmin], nmin, nmax, err);
   else
        shc_arithmetic_wise(shcs, a, nmin, nmax, err);

    CHARM(err_handler)(err, 1);


    snprintf(filec, NSTR_LONG, "%s/%s_c_nmin%lu_nmax%lu%s",
             FOLDER, func, nmin, nmax, FTYPE);
    snprintf(files, NSTR_LONG, "%s/%s_s_nmin%lu_nmax%lu%s",
             FOLDER, func, nmin, nmax, FTYPE);
#ifdef GENREF
    e += array2file(filec, shcs->c[0], shcs->nc);
    e += array2file(files, shcs->s[0], shcs->ns);
#else
    e += validate(filec, shcs->c[0], shcs->nc,
                  PREC(10.0) * CHARM(glob_threshold));
    e += validate(files, shcs->s[0], shcs->ns,
                  PREC(10.0) * CHARM(glob_threshold));
#endif
    /* --------------------------------------------------------------------- */


    /* Check whether some common errors in the input parameters are treated
     * properly */
    /* --------------------------------------------------------------------- */
    shc_arithmetic_wise(shcs, a, 0, shcs->nmax + 1, err);
    if (CHARM(err_isempty)(err))
    {
        printf("\n\n            Didn't throw an error for "
               "\"nmax > shcs->nmax\".");
        e += 1;
    }
    CHARM(err_reset)(err);


    shc_arithmetic_wise(shcs, a, shcs->nmax + 1, shcs->nmax + 2, err);
    if (CHARM(err_isempty)(err))
    {
        printf("\n\n            Didn't throw an error for \"nmin > "
               " shcs->nmax\".");
        e += 1;
    }
    CHARM(err_reset)(err);


    shc_arithmetic_wise(shcs, a, shcs->nmax, shcs->nmax - 1, err);
    if (CHARM(err_isempty)(err))
    {
        printf("\n\n            Didn't throw an error for \"nmin > nmax\".");
        e += 1;
    }
    CHARM(err_reset)(err);


    if ((shc_arithmetic_wise == CHARM(shc_div_degree_wise)) ||
        (shc_arithmetic_wise == CHARM(shc_div_order_wise)))
    {
        /* Check division by zero */
        for (size_t i = 0; i <= shcs->nmax; i++)
            a[i] = PREC(0.0);


        shc_arithmetic_wise(shcs, a, 0, shcs->nmax, err);
        if (CHARM(err_isempty)(err))
        {
            printf("\n\n            Didn't throw an error for division by "
                   "zero.");
            e += 1;
        }
    }
    /* --------------------------------------------------------------------- */


    /* Free all memory allocated by this function */
    /* --------------------------------------------------------------------- */
    free(a);
    CHARM(shc_free)(shcs);
    CHARM(err_free)(err);
    /* --------------------------------------------------------------------- */


    return e;
}
