/* This header file is not a part of API.
 *
 * Defines some constants used throughout the tests. */


#ifndef __PARAMETERS_H__
#define __PARAMETERS_H__


#include <config.h>


/* Size of char arrays to store file names to be imported, etc. */
#undef NSTR
#define NSTR (2024)


/* Size of some small char arrays. */
#undef NSTR2
#define NSTR2 (20)


/* Suffix of the refenrece data files to be loaded. */
#undef FTYPE
#define FTYPE ".txt"


/* Maximum degree for harmonic analysis and synthesis (except for the area-mean
 * values on irregular surfaces) */
#undef NMAX
#define NMAX (3)


/* Auxiliary maximum degree to synthesize area-mean values on irregular
 * surfaces */
#undef NMAX2
#define NMAX2 (15)


/* Number of radial layers to test for the solid synthesis and analysis. */
#undef NDELTAR
#define NDELTAR (2)


/* Distances between the radial layers. */
#undef DELTAR
#define DELTAR (1000.0)


/* Number of various versions of custom grids to test SHA and SHS (except for
 * area-mean values on irregular surfaces). */
#undef NCUSTOM_GRD
#define NCUSTOM_GRD (4)


/* Number of custom grids to test SHS of area-mean values on irregular
 * surfaces. */
#undef NCUSTOM_GRD_ISURF
#define NCUSTOM_GRD_ISURF (3)


/* Number of scattered points/cells to test SHS. */
#undef NSCTR
#define NSCTR (9)


/* Path to the folder with the reference test data */
#undef FOLDER
#if CHARM_FLOAT
#   define FOLDER "../data/tests/single"
#elif CHARM_QUAD
#   define FOLDER "../data/tests/quad"
#else
#   define FOLDER "../data/tests/double"
#endif


/* In the tests, we want to use reference potential coefficients that are all
 * non-zero, but this is not the case with the degree-1 coefficients of our
 * reference model.  Therefore, in some of the tests, we artificially set these
 * coefficients to some more or less random non-zero value.  Also, the
 * zero-degree coefficient of the reference model is far larger than the other
 * coefficients.  For numerical reasons, it is more favourable to lower its
 * value by a few orders of magnitude. */
#undef  C00
#define C00 (0.0000001)
#undef  C10
#define C10 (0.0000002)
#undef  C11
#define C11 (0.0000003)
#undef  S11
#define S11 (0.0000004)


/* Radius to reference the topographic heights */
#undef RREF
#define RREF (6378136.3)


/* Constant used to intentionally break the latitudinal symmetry of cell
 * grids. */
#undef BREAK_SYMM
#define BREAK_SYMM (0.001)


/* Co-latitude to split the "[0, PI]" domain of co-latitudes into two
 * intervals, "[0, CLT0]" and "[CLT0, PI]".  This is useful to check, for
 * instance, the orthogonality of Legendre functions as "[0, CLT0] + [CLT0, PI]
 * = [0, PI]", where the right-hand side is computed using the orthogonality
 * property of Legendre functions. */
#undef CLT0
#define CLT0 (PREC(0.3234))


/* Longitude to split the "[0, 2 * PI]" domain of longitudes into two
 * intervals, similarly as "CLT0" splits the domain of co-latitudes. */
#undef LON0
#define LON0 (PREC(4.697))


/* Buffer size when converting "__float128" floating point numbers to strings.
 * */
#undef BUF_QUAD
#define BUF_QUAD (256)



#endif

