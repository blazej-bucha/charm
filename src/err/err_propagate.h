/* This header file is not a part of API. */


#ifndef __ERR_PROPAGATE_H__
#define __ERR_PROPAGATE_H__


#include <config.h>
#include "../prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


/** Propagates ``err`` to the caller by appending ``file``, ``line`` and
 * ``func`` to ``err``, respectively.  ``err->level`` is updated automatically.
 * */
extern void CHARM(err_propagate)(CHARM(err) *err,
                                 const char *file,
                                 unsigned int line,
                                 const char *func);


#ifdef __cplusplus
}
#endif


#endif
