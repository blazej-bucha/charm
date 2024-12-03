/* This header file is not a part of API. */


#ifndef __ERR_IS_NULL_PTR_H__
#define __ERR_IS_NULL_PTR_H__


#include <config.h>
#include "../prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern _Bool CHARM(err_is_null_ptr)(CHARM(err) *
#if HAVE_MPI
                                    , _Bool distributed,
                                    MPI_Comm
#endif
                                   );


#ifdef __cplusplus
}
#endif


#endif
