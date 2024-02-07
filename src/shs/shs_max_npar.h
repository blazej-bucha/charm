/* This header file is not a part of API. */


#ifndef __SHS_MAX_NPAR_H__
#define __SHS_MAX_NPAR_H__


/* Maximum number of quantities that can be synthesised in a single run.
 * Presently, the value is "6", because CHarm can synthesised the second-order
 * gradients at maximum, that is, the gravitational tensor, which has 6 unique
 * elements */
#undef SHS_MAX_NPAR
#define SHS_MAX_NPAR (6)


#endif

