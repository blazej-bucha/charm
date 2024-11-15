/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../src/prec.h"
#include "../src/misc/misc_is_nearly_equal.h"
#include "parameters.h"
#include "cmp_vals.h"
/* ------------------------------------------------------------------------- */






/* Compare two unsigned long values */
/* ------------------------------------------------------------------------- */
long int cmp_vals_ulong(unsigned long x,
                        unsigned long xref)
{
    long int ret = 0;
    if (x != xref)
    {
        printf("\n");
        printf("            result:    %lu\n", x);
        printf("            reference: %lu\n", xref);
        ret = 1;
    }


    return ret;
}
/* ------------------------------------------------------------------------- */






/* Compares two real values up to some threshold */
/* ------------------------------------------------------------------------- */
long int cmp_vals_real(REAL x,
                       REAL xref,
                       REAL eps)
{
    long int ret = 0;
#if defined(CHARM_QUAD)
    int flt128_n;
    char buf1[BUF_QUAD];
    char buf2[BUF_QUAD];
#endif


    if (!CHARM(misc_is_nearly_equal)(x, xref, eps))
    {
        printf("\n");
#if defined(CHARM_FLOAT)
        printf("            result:    %14.7e\n", x);
        printf("            reference: %14.7e\n", xref);
#elif defined(CHARM_QUAD)
        flt128_n = quadmath_snprintf(buf1, sizeof(buf1), "%41.34Qe", 70,
                                     x);
        if ((size_t)flt128_n >= BUF_QUAD)
        {
            fprintf(stderr, "Failed to convert a \"__float128\" "
                            "number to a string.\n");
            exit(CHARM_FAILURE);
        }
        flt128_n = quadmath_snprintf(buf2, sizeof(buf2), "%41.34Qe", 70,
                                     xref);
        if ((size_t)flt128_n >= BUF_QUAD)
        {
            fprintf(stderr, "Failed to convert a \"__float128\""
                            "number to a string.\n");
            exit(CHARM_FAILURE);
        }


        printf("            result:    %s\n", buf1);
        printf("            reference: %s\n", buf2);
#else
        printf("            result:    %23.16e\n", x);
        printf("            reference: %23.16e\n", xref);
#endif
        printf("            WARNING: Possibly too large error!\n");


        ret = 1;
    }


    return ret;
}
/* ------------------------------------------------------------------------- */

