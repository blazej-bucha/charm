/* This header file is not a part of API.
 *
 * Defines some symbolic constants used throughout the tests. */


#ifndef __ERROR_MESSAGES_H__
#define __ERROR_MESSAGES_H__


/* Error messages to be printed if the CHarm's structures from the public API
 * couldn't be created.  This should provided better clues on the error than
 * a general message "CHARM_ERR_MALLOC_FAILURE" */


#include "../src/prec.h"


#undef ERR_MSG_SHC
#define ERR_MSG_SHC "Failed to create the \"charm" CHARM_SUFFIX \
                    "_shc\" structure.\n"


#undef ERR_MSG_ERR
#define ERR_MSG_ERR "Failed to create the \"charm" CHARM_SUFFIX \
                    "_err\" structure.\n"


#undef ERR_MSG_POINT
#define ERR_MSG_POINT "Failed to create the \"charm" CHARM_SUFFIX \
                      "_point\" structure.\n"


#undef ERR_MSG_CELL
#define ERR_MSG_CELL "Failed to create the \"charm" CHARM_SUFFIX \
                      "_cell\" structure.\n"


#undef ERR_MSG_PNMJ
#define ERR_MSG_PNMJ "Failed to create the \"charm" CHARM_SUFFIX \
                      "_pnmj\" structure.\n"


#endif

