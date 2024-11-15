/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
/* ------------------------------------------------------------------------- */






/* ------------------------------------------------------------------------- */
void test_suite_end(long int e)
{
    if (e)
        printf("%ld unexpected %s in total\n",
               e, e == 1 ? "result" : "results");
    else
        printf("All tests passed\n");
    printf("\n");
    printf("End of the test suite\n");
    printf("\n");
    printf("==================================\n\n");


    return;
}
/* ------------------------------------------------------------------------- */
