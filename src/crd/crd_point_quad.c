/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
/* ------------------------------------------------------------------------- */






CHARM(point) *CHARM(crd_point_quad)(unsigned long nmax,
                                    REAL r,
                                    void (*quad_shape)(unsigned long,
                                                       size_t *,
                                                       size_t *),
                                   CHARM(point) *(*quad_chunk)(unsigned long,
                                                               REAL,
                                                               size_t,
                                                               size_t,
                                                               CHARM(err) *))
{
    size_t nlat, nlon;
    quad_shape(nmax, &nlat, &nlon);


    CHARM(point) *pnt = NULL;
    CHARM(err)   *err = NULL;


    err = CHARM(err_init)();
    if (err == NULL)
        return NULL;


    /* "quad_chunk" takes an error structure as an input argument.  This is
     * because of the MPI variant of this function.  But in this non-MPI
     * variant, the error should be obvious, so we do not propagate it to the
     * caller. */
    pnt = quad_chunk(nmax, r, nlat, 0, err);
    if (!CHARM(err_isempty)(err))
    {
        /* Just in case there is an error in "err", but "pnt" is not "NULL". */
        CHARM(crd_point_free)(pnt);
        pnt = NULL;
    }


    CHARM(err_free)(err);


    return pnt;
}
