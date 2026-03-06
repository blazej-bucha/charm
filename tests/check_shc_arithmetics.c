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
#include "check_shc_arithmetics.h"
/* ------------------------------------------------------------------------- */






/* This file checks arithmetics with "charm_shc".
 *
 * * "shc_add", "shc_sub", "shc_mul", "shc_div" -- Arithmetics functions to
 *                                                 check.
 *
 * */
long int check_shc_arithmetics(void (*shc_arithmetic)(CHARM(shc) *,
                                                      CHARM(shc) *,
                                                      CHARM(shc) *,
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
    unsigned long Nrop = 4;  /* Must be larger than zero in this routine */
    unsigned long Nop1 = 4;
    unsigned long Nop2 = 3;
    REAL c = PREC(0.1);


    /* Get name of the routine that is being validated */
    char func[NSTR_SHORT];
    if (shc_arithmetic == CHARM(shc_add))
        snprintf(func, NSTR_SHORT, "shc_add");
    else if (shc_arithmetic == CHARM(shc_sub))
        snprintf(func, NSTR_SHORT, "shc_sub");
    else if (shc_arithmetic == CHARM(shc_mul))
        snprintf(func, NSTR_SHORT, "shc_mul");
    else if (shc_arithmetic == CHARM(shc_div))
    {
        snprintf(func, NSTR_SHORT, "shc_div");
        Nop2 += 2;  /* To avoid divisions by zero */
    }
    /* --------------------------------------------------------------------- */


    /* Create structures to hold spherical harmonic coefficients */
    /* --------------------------------------------------------------------- */
    CHARM(shc) *rop = CHARM(shc_calloc)(Nrop, PREC(1.0), PREC(1.0));
    if (rop == NULL)
    {
        fprintf(stderr, "%s", ERR_MSG_SHC);
        exit(CHARM_FAILURE);
    }


    CHARM(shc) *op1 = CHARM(shc_calloc)(Nop1, PREC(1.0), PREC(1.0));
    if (op1 == NULL)
    {
        fprintf(stderr, "%s", ERR_MSG_SHC);
        exit(CHARM_FAILURE);
    }


    CHARM(shc) *op2 = CHARM(shc_calloc)(Nop2, PREC(1.0), PREC(1.0));
    if (op2 == NULL)
    {
        fprintf(stderr, "%s", ERR_MSG_SHC);
        exit(CHARM_FAILURE);
    }
    /* --------------------------------------------------------------------- */


    /* Initialize coefficients of "op1" and "op2" */
    /* --------------------------------------------------------------------- */
    for (unsigned long m = 0; m <= op1->nmax; m++)
    {
        for (unsigned long n = m; n <= op1->nmax; n++)
        {
            op1->c[m][n - m] = PREC(1.0);
            if (m > 0)
                op1->s[m][n - m] = PREC(2.0);
        }
    }


    for (unsigned long m = 0; m <= op2->nmax; m++)
    {
        for (unsigned long n = m; n <= op2->nmax; n++)
        {
            op2->c[m][n - m] = PREC(10.0);
            if (m > 0)
                op2->s[m][n - m] = PREC(20.0);
        }
    }
    /* --------------------------------------------------------------------- */


    /* --------------------------------------------------------------------- */
    unsigned long nmin = 0;
    unsigned long nmax = rop->nmax;
    shc_arithmetic(rop, op1, op2, nmin, nmax, err);
    CHARM(err_handler)(err, 1);


    char filec[NSTR_LONG];
    char files[NSTR_LONG];
    snprintf(filec, NSTR_LONG, "%s/%s_c_nmin%lu_nmax%lu%s",
             FOLDER, func, nmin, nmax, FTYPE);
    snprintf(files, NSTR_LONG, "%s/%s_s_nmin%lu_nmax%lu%s",
             FOLDER, func, nmin, nmax, FTYPE);
#ifdef GENREF
    e += array2file(filec, rop->c[0], rop->nc);
    e += array2file(files, rop->s[0], rop->ns);
#else
    e += validate(filec, rop->c[0], rop->nc,
                  PREC(10.0) * CHARM(glob_threshold));
    e += validate(files, rop->s[0], rop->ns,
                  PREC(10.0) * CHARM(glob_threshold));
#endif


    /* Before the next test, we have to reset "rop" */
    for (unsigned long m = 0; m <= rop->nmax; m++)
    {
        for (unsigned long n = m; n <= rop->nmax; n++)
        {
            rop->c[m][n - m] = PREC(0.0);
            rop->s[m][n - m] = PREC(0.0);
        }
    }
    nmin = 1;
    nmax = 2;
    shc_arithmetic(rop, op1, op2, nmin, nmax, err);
    CHARM(err_handler)(err, 1);


    snprintf(filec, NSTR_LONG, "%s/%s_c_nmin%lu_nmax%lu%s",
             FOLDER, func, nmin, nmax, FTYPE);
    snprintf(files, NSTR_LONG, "%s/%s_s_nmin%lu_nmax%lu%s",
             FOLDER, func, nmin, nmax, FTYPE);
#ifdef GENREF
    e += array2file(filec, rop->c[0], rop->nc);
    e += array2file(files, rop->s[0], rop->ns);
#else
    e += validate(filec, rop->c[0], rop->nc,
                  PREC(10.0) * CHARM(glob_threshold));
    e += validate(files, rop->s[0], rop->ns,
                  PREC(10.0) * CHARM(glob_threshold));
#endif
    /* --------------------------------------------------------------------- */


    /* Check whether some common errors in the input parameters are treated
     * properly */
    /* --------------------------------------------------------------------- */
    shc_arithmetic(rop, op1, op2, 0, rop->nmax + 1, err);
    if (CHARM(err_isempty)(err))
    {
        printf("\n\n            Didn't throw an error for "
               "\"nmax > rop->nmax\".");
        e += 1;
    }
    CHARM(err_reset)(err);


    shc_arithmetic(rop, op1, op2, rop->nmax + 1, rop->nmax + 2, err);
    if (CHARM(err_isempty)(err))
    {
        printf("\n\n            Didn't throw an error for \"nmin > "
               "rop->nmax\".");
        e += 1;
    }
    CHARM(err_reset)(err);


    shc_arithmetic(rop, op1, op2, rop->nmax, rop->nmax - 1, err);
    if (CHARM(err_isempty)(err))
    {
        printf("\n\n            Didn't throw an error for \"nmin > nmax\".");
        e += 1;
    }
    CHARM(err_reset)(err);


    rop->mu += c;
    shc_arithmetic(rop, op1, op2, 0, rop->nmax, err);
    if (CHARM(err_isempty)(err))
    {
        printf("\n\n            Didn't throw an error for unmatching "
               "\"mu\".");
        e += 1;
    }
    CHARM(err_reset)(err);
    rop->mu -= c;


    op1->mu += c;
    shc_arithmetic(rop, op1, op2, 0, rop->nmax, err);
    if (CHARM(err_isempty)(err))
    {
        printf("\n\n            Didn't throw an error for unmatching "
               "\"mu\".");
        e += 1;
    }
    CHARM(err_reset)(err);
    op1->mu -= c;


    op2->mu += c;
    shc_arithmetic(rop, op1, op2, 0, rop->nmax, err);
    if (CHARM(err_isempty)(err))
    {
        printf("\n\n            Didn't throw an error for unmatching \"mu\".");
        e += 1;
    }
    CHARM(err_reset)(err);
    op2->mu -= c;


    rop->r += c;
    shc_arithmetic(rop, op1, op2, 0, rop->nmax, err);
    if (CHARM(err_isempty)(err))
    {
        printf("\n\n            Didn't throw an error for unmatching \"r\".");
        e += 1;
    }
    CHARM(err_reset)(err);
    rop->r -= c;


    op1->r += c;
    shc_arithmetic(rop, op1, op2, 0, rop->nmax, err);
    if (CHARM(err_isempty)(err))
    {
        printf("\n\n            Didn't throw an error for unmatching \"r\".");
        e += 1;
    }
    CHARM(err_reset)(err);
    op1->r -= c;


    op2->r += c;
    shc_arithmetic(rop, op1, op2, 0, rop->nmax, err);
    if (CHARM(err_isempty)(err))
    {
        printf("\n\n            Didn't throw an error for unmatching \"r\".");
        e += 1;
    }
    CHARM(err_reset)(err);
    op2->r -= c;


    if (shc_arithmetic == CHARM(shc_div))
    {
        /* Check division by zero */
        for (unsigned long m = 0; m <= rop->nmax; m++)
        {
            for (unsigned long n = m; n <= rop->nmax; n++)
            {
                op2->c[m][n - m] = PREC(0.0);
                op2->s[m][n - m] = PREC(0.0);
            }
        }


        shc_arithmetic(rop, op1, op2, 0, rop->nmax, err);
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
    CHARM(shc_free)(rop);
    CHARM(shc_free)(op1);
    CHARM(shc_free)(op2);
    CHARM(err_free)(err);
    /* --------------------------------------------------------------------- */


    return e;
}
