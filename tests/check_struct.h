/* This header file is not a part of API. */


#ifndef __CHECK_STRUCT_H__
#define __CHECK_STRUCT_H__


#include <config.h>
#include "../src/prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern long int check_struct__Bool(_Bool,
                                   _Bool,
                                   _Bool,
                                   const char *,
                                   const char *,
                                   const char *);
extern long int check_struct_int(int,
                                 int,
                                 _Bool,
                                 const char *,
                                 const char *,
                                 const char *);
extern long int check_struct_size_t(size_t,
                                    size_t,
                                    _Bool,
                                    const char *,
                                    const char *,
                                    const char *);
extern long int check_struct_ulong(unsigned long,
                                   unsigned long,
                                   _Bool,
                                   const char *,
                                   const char *,
                                   const char *);
extern long int check_struct_REAL(REAL,
                                  REAL,
                                  _Bool,
                                  const char *,
                                  const char *,
                                  const char *);
extern long int check_struct_ptr(void *,
                                 void *,
                                 _Bool,
                                 const char *,
                                 const char *,
                                 const char *);


#ifdef __cplusplus
}
#endif


#endif
