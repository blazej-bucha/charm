/* This file defines functions designed to be especially useful when comparing
 * various members of custom structs.  In that case, we want to print errors
 * differently as done by "cmp_vals". */






/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../src/misc/misc_is_nearly_equal.h"
#include "parameters.h"
/* ------------------------------------------------------------------------- */






/* Prints a message when struct members are found to be wrong.
 *
 * "vi" should be either "VALID" or "INVALID",
 *
 * "func_call_str" should be a string representing the function call,
 * e.g. "shc_malloc(10, 1.0, 1.0)",
 *
 * "err_str" should be an error string explaining what happened but should not
 * happen, e.g. "returned a NULL pointer" or "didn't return a NULL pointer". */
void print_msg(char *vi,
               char *func_call_str,
               char *err_str)
{
    printf("\n            ERROR: %s function call %s %s",
           vi, func_call_str, err_str);
    fflush(stdout);


    return;
}






/* Macro to create functions comparing struct members of various data types.
 *
 * Checks whether "in"
 *
 *      is equal to (if "eq == 1")
 *
 * or
 *
 *      is not equal to (if "eq == 0")
 *
 * "ref" for a function call represented as a string "func_call_str".  If "vi"
 * is "VALID", the function call in "func_call_str" is assumed to be valid; if
 * "vi" is "INVALID", the call is assumed to be "invalid".  If the condition
 * "in" "eq" "ref" is not satisfied, an error message is printed based on "vi",
 * "func_call_str" and "err_str" and returned is "1"; otherwise, returned is
 * "0". */
#define CHECK_STRUCT(CHECK_STRUCT_TYPE, DATA_TYPE)                            \
                                                                              \
        long int CHECK_STRUCT_TYPE(DATA_TYPE in,                              \
                                   DATA_TYPE ref,                             \
                                   _Bool eq,                                  \
                                   char *vi,                                  \
                                   char *func_call_str,                       \
                                   char *err_str)                             \
        {                                                                     \
            long int e = 0;                                                   \
                                                                              \
                                                                              \
            if (eq)                                                           \
            {                                                                 \
                if (in == ref)                                                \
                    e += 1;                                                   \
            }                                                                 \
            else                                                              \
            {                                                                 \
                if (in != ref)                                                \
                    e += 1;                                                   \
            }                                                                 \
                                                                              \
                                                                              \
            if (e > 0)                                                        \
                print_msg(vi, func_call_str, err_str);                        \
                                                                              \
                                                                              \
            return e;                                                         \
        }






/* The same as "CHECK_STRUCT_TYPE" but for floating points, which need
 * a special treatment.  Technically, "CHECK_STRUCT_TYPE" could also be used
 * with "REAL", but this might perhaps not be safe enough. */
long int check_struct_REAL(REAL in,
                           REAL ref,
                           _Bool eq,
                           char *vi,
                           char *func_call_str,
                           char *err_str)
{
    long int e = 0;


    if (eq)
    {
        if (CHARM(misc_is_nearly_equal)(in, ref, CHARM(glob_threshold)))
            e += 1;
    }
    else
    {
        if (!CHARM(misc_is_nearly_equal)(in, ref, CHARM(glob_threshold)))
            e += 1;
    }


    if (e > 0)
        print_msg(vi, func_call_str, err_str);


    return e;
}






CHECK_STRUCT(check_struct__Bool, _Bool)
CHECK_STRUCT(check_struct_int, int)
CHECK_STRUCT(check_struct_size_t, size_t)
CHECK_STRUCT(check_struct_ulong, unsigned long)
CHECK_STRUCT(check_struct_ptr, void *)

