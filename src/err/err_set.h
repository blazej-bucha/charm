/* This header file is not a part of API. */


#ifndef __ERR_SET_H__
#define __ERR_SET_H__


#include <config.h>
#include "../prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


/** Sets members of ``err`` to ``file``, ``line``, ``func``, ``code`` and
 * ``msg``, respectively.  ``err->level`` is updated automatically.  */
extern void CHARM(err_set)(CHARM(err) *err,
                           const char *file,
                           unsigned int line,
                           const char *func,
                           int code,
                           const char *msg);


#ifdef __cplusplus
}
#endif


#endif
