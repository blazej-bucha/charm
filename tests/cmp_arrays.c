/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef _MSC_VER
#   define _USE_MATH_DEFINES
#endif
#include <math.h>
#include "../src/prec.h"
/* ------------------------------------------------------------------------- */






/* Compares elements of two arrays up to some threshold */
/* ------------------------------------------------------------------------- */
int cmp_arrays(REAL *arr1, REAL *arr2, size_t n, REAL eps)
{
    int ret = 0;
#if defined(CHARM_QUAD)
    int flt128_n;
    char buf1[256];
    char buf2[256];
#endif


    for (size_t i = 0; i < n; i++)
    {
        if (!CHARM(misc_is_nearly_equal)(arr1[i], arr2[i], eps))
        {
            printf("\n");
#if defined(CHARM_FLOAT)
            printf("        result:    %14.7e\n", arr1[i]);
            printf("        reference: %14.7e\n", arr2[i]);
#elif defined(CHARM_QUAD)
            flt128_n = quadmath_snprintf(buf1, sizeof(buf1), "%41.34Qe", 70,
                                         arr1[i]);
            if ((size_t)flt128_n >= 256)
            {
                fprintf(stderr, "Failed to convert a \"__float128\" "
                                "number to a string.\n");
                exit(CHARM_FAILURE);
            }
            flt128_n = quadmath_snprintf(buf2, sizeof(buf2), "%41.34Qe", 70,
                                         arr2[i]);
            if ((size_t)flt128_n >= 256)
            {
                fprintf(stderr, "Failed to convert a \"__float128\""
                                "number to a string.\n");
                exit(CHARM_FAILURE);
            }


            printf("        result:    %s\n", buf1);
            printf("        reference: %s\n", buf2);
#else
            printf("        result:    %23.16e\n", arr1[i]);
            printf("        reference: %23.16e\n", arr2[i]);
#endif
            printf("        WARNING: Possibly too large error!\n");


            ret = 1;
        }
    }


    if (ret == 1)
        printf("\n");


    return ret;
}
/* ------------------------------------------------------------------------- */
