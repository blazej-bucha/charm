/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../src/prec.h"
#include "parameters.h"
#include "modify_low_degree_coefficients.h"
/* ------------------------------------------------------------------------- */






/* Modifies coefficients of degrees "0" and "1" to allow for an accurate
 * validation in all precisions and to make sure that the degree-1 coefficients
 * are non-zero, so that the tests are sensitive to all harmonic degrees and
 * orders. */
/* ------------------------------------------------------------------------- */
void modify_low_degree_coefficients(CHARM(shc) *shcs)
{
    shcs->c[0][0 - 0] = (REAL)C00;


    if (shcs->nmax > 0)
    {
        shcs->c[0][1 - 0] = (REAL)C10;
        shcs->c[1][1 - 1] = (REAL)C11;
        shcs->s[1][1 - 1] = (REAL)S11;
    }


    return;
}
/* ------------------------------------------------------------------------- */
