/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#if HAVE_MPI
#   include <mpi.h>
#endif
#include "../prec.h"
/* ------------------------------------------------------------------------- */






CHARM(err) *CHARM(err_init)(void)
{
    /* Initialize the "CHARM(err)" structure */
    /* --------------------------------------------------------------------- */
    /* Allocate memory for the "CHARM(err)" data type */
    CHARM(err) *err = (CHARM(err) *)malloc(sizeof(CHARM(err)));
    if (err == NULL)
        return NULL;


    err->file = NULL;
    err->line = NULL;
    err->func = NULL;
    err->msg  = NULL;
    /* --------------------------------------------------------------------- */


    /* Allocate some members of "err" to be able to store "CHARM_ERR_MAX_LEVEL"
     * error levels. */
    /* --------------------------------------------------------------------- */
    err->file = malloc(CHARM_ERR_MAX_LEVEL * sizeof(char *));
    if (err->file == NULL)
        goto FAILURE_1;


    err->line = malloc(CHARM_ERR_MAX_LEVEL * sizeof(unsigned int));
    if (err->file == NULL)
        goto FAILURE_1;


    err->func = malloc(CHARM_ERR_MAX_LEVEL * sizeof(char *));
    if (err->func == NULL)
        goto FAILURE_1;
    /* --------------------------------------------------------------------- */


    /* Now allocate memory for each level, so that some of the member of "err"
     * will be able to store strings for file and function names within each
     * level. */
    /* --------------------------------------------------------------------- */
    for (size_t e = 0; e < CHARM_ERR_MAX_LEVEL; e++)
        err->file[e] = NULL;


    for (size_t e = 0; e < CHARM_ERR_MAX_LEVEL; e++)
        err->func[e] = NULL;


    for (size_t e = 0; e < CHARM_ERR_MAX_LEVEL; e++)
    {
        err->file[e] = malloc(CHARM_ERR_MAX_FILE * sizeof(char));
        if (err->file[e] == NULL)
            goto FAILURE_2;
    }


    for (size_t e = 0; e < CHARM_ERR_MAX_LEVEL; e++)
    {
        err->func[e] = malloc(CHARM_ERR_MAX_FUNC * sizeof(char));
        if (err->func[e] == NULL)
            goto FAILURE_2;
    }


    err->msg  = malloc(CHARM_ERR_MAX_MSG  * sizeof(char));
    if (err->msg == NULL)
        goto FAILURE_2;
    /* --------------------------------------------------------------------- */


    /* Now reset "err" to the default values */
    /* --------------------------------------------------------------------- */
    CHARM(err_reset)(err);
    /* --------------------------------------------------------------------- */


#if HAVE_MPI
    /* Set the MPI-related members */
    /* --------------------------------------------------------------------- */
    err->distributed = 0;
    err->comm        = MPI_COMM_NULL;
    /* --------------------------------------------------------------------- */
#endif


    /* Exit the function */
    /* --------------------------------------------------------------------- */
EXIT:
    return err;


FAILURE_2:
    for (size_t e = 0; e < CHARM_ERR_MAX_LEVEL; e++)
        free(err->file[e]);


    for (size_t e = 0; e < CHARM_ERR_MAX_LEVEL; e++)
        free(err->func[e]);


FAILURE_1:
    free(err->file);
    free(err->func);
    free(err->line);
    free(err->msg);
    free(err);


    err = NULL;


    goto EXIT;
    /* --------------------------------------------------------------------- */
}
