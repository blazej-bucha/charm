/* This header file is not a part of API.
 *
 * Typical usage of "sscanf", "fscanf" with macros from this header file is as
 * follows,
 *
 *      char buf1[SCANF_BUFFER];
 *
 *      char buf2[SCANF_BUFFER];
 *
 *      sscanf(some_string,
 *             SCANF_SFS(SCANF_WIDTH) " " SCANF_SFS(SCANF_WIDTH), buf1, buf2);
 *
 * */


#ifndef __MISC_SCANF_H__
#define __MISC_SCANF_H__


#include <config.h>


#ifdef __cplusplus
extern "C"
{
#endif


/* Field width of strings parsed by "sscanf", "fscanf", etc. */
#undef SCANF_WIDTH
#define SCANF_WIDTH 255


/* All "char" arrays to store strings parsed from "sscanf", "fscanf",
 * etc. must be have this number of elements (including the terminating null
 * byte) */
#undef SCANF_BUFFER
#define SCANF_BUFFER (SCANF_WIDTH + 1)


/* When calling "sscanf", "fscanf", etc., use "SCANF_SFS" instead of "%s" */
#undef SCANF_SFS
#undef SCANF_SFS_CAT
#define SCANF_SFS_CAT(FIELD_WIDTH) "%" #FIELD_WIDTH "s"
#define SCANF_SFS(S) SCANF_SFS_CAT(S)


#ifdef __cplusplus
}
#endif


#endif
