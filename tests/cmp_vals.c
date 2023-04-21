/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../src/prec.h"
#include "../src/misc/misc_is_nearly_equal.h"
#include "parameters.h"
/* ------------------------------------------------------------------------- */






/* Compares two values up to some threshold */
/* ------------------------------------------------------------------------- */
long int cmp_vals(REAL val1, REAL val2, REAL eps)
{
    long int ret = 0;
#if defined(CHARM_QUAD)
    int flt128_n;
    char buf1[BUF_QUAD];
    char buf2[BUF_QUAD];
#endif


    if (!CHARM(misc_is_nearly_equal)(val1, val2, eps))
    {
        printf("\n");
#if defined(CHARM_FLOAT)
        printf("            result:    %14.7e\n", val1);
        printf("            reference: %14.7e\n", val2);
#elif defined(CHARM_QUAD)
        flt128_n = quadmath_snprintf(buf1, sizeof(buf1), "%41.34Qe", 70,
                                     val1);
        if ((size_t)flt128_n >= BUF_QUAD)
        {
            fprintf(stderr, "Failed to convert a \"__float128\" "
                            "number to a string.\n");
            exit(CHARM_FAILURE);
        }
        flt128_n = quadmath_snprintf(buf2, sizeof(buf2), "%41.34Qe", 70,
                                     val2);
        if ((size_t)flt128_n >= BUF_QUAD)
        {
            fprintf(stderr, "Failed to convert a \"__float128\""
                            "number to a string.\n");
            exit(CHARM_FAILURE);
        }


        printf("            result:    %s\n", buf1);
        printf("            reference: %s\n", buf2);
#else
        printf("            result:    %23.16e\n", val1);
        printf("            reference: %23.16e\n", val2);
#endif
        printf("            WARNING: Possibly too large error!\n");


        ret = 1;
    }


    return ret;
}
/* ------------------------------------------------------------------------- */

