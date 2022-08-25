/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
#include "shc_reset_coeffs.h"
#include "shc_read_line.h"
#include "shc_read_mtdt.h"
#include "../err/err_set.h"
#include "../err/err_propagate.h"
/* ------------------------------------------------------------------------- */






void CHARM(shc_read_tbl)(FILE *stream, unsigned long nmax, CHARM(shc) *shcs,
                         CHARM(err) *err)
{
    /* Read the metadata of spherical harmonic coefficients */
    /* --------------------------------------------------------------------- */
    unsigned long nmax_file;
    CHARM(shc_read_mtdt)(stream, &nmax_file, &(shcs->mu), &(shcs->r), err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return;
    }
    /* --------------------------------------------------------------------- */







    /* Check maximum harmonic degrees */
    /* --------------------------------------------------------------------- */
    if (shcs->nmax < nmax)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Too low maximum degree \"shcs->nmax\" to read "
                       "coefficients up to degree \"nmax\".");
        return;
    }


    if (nmax_file < nmax)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Too low maximum degree inside the input file to read "
                       "coefficients up to degree \"nmax\".");
        return;
    }
    /* --------------------------------------------------------------------- */






    /* Read the table of spherical harmonic coefficients */
    /* --------------------------------------------------------------------- */
    /* At first, reset all coefficients in "shcs" to zero. */
    CHARM(shc_reset_coeffs)(shcs);


    unsigned long n, m;
    REAL cnm, snm;
    while (CHARM(shc_read_line)(stream, &n, &m, &cnm, &snm, SHC_READ_LINE_TBL,
                                err) != -1)
    {
        if (!CHARM(err_isempty)(err))
        {
            CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
            return;
        }


        if (n > nmax)
            continue;


        shcs->c[m][n - m] = cnm;
        shcs->s[m][n - m] = snm;
    }
    /* --------------------------------------------------------------------- */






    return;
}
