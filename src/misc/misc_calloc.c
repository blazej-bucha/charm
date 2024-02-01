/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
#include "misc_calloc.h"
/* ------------------------------------------------------------------------- */






/* Defines a custom "calloc" function with the same interface as "malloc". */
void *CHARM(misc_calloc)(size_t size)
{
    return calloc(size, 1);
}
