/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <string.h>
#include "../prec.h"
#include "malloc_aligned.h"
#include "calloc_aligned.h"
/* ------------------------------------------------------------------------- */






/* The same as "malloc_aligned", but additionally sets all allocated bytes to
 * zero. */
void *CHARM(calloc_aligned)(size_t alignment, size_t nmemb, size_t size)
{
    void *p = CHARM(malloc_aligned)(alignment, nmemb * size);
    if (p == NULL)
        return NULL;


    memset(p, 0, nmemb * size);


    return p;
}
