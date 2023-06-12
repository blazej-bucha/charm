/* This header file is not a part of API.
 *
 * Defines some symbolic constants used throughout the tests. */


#ifndef __PARAMETERS_H__
#define __PARAMETERS_H__


#include <config.h>


/* Maximum harmonic degree of the input topography in "SHCS_IN_PATH_TOPO_MTX" 
 * */
#undef SHCS_NMAX_TOPO
#define SHCS_NMAX_TOPO (4)


/* Path to input spherical harmonic coefficients of the topography in the 
 * matrix format */
#undef SHCS_IN_PATH_TOPO_MTX
#define SHCS_IN_PATH_TOPO_MTX "../data/input/Earth2014-RET2014-degree4-mtx.txt"


/* Path to output spherical harmonic coefficients of the topography in the bin 
 * format. */
#undef SHCS_OUT_PATH_TOPO_BIN
#define SHCS_OUT_PATH_TOPO_BIN "../data/output/Earth2014-RET2014-degree4" \
                               ".shcs"


/* Path to output spherical harmonic coefficients of the topography in the mtx 
 * format. */
#undef SHCS_OUT_PATH_TOPO_MTX
#define SHCS_OUT_PATH_TOPO_MTX "../data/output/Earth2014-RET2014-degree4" \
                               "-mtx.txt"


/* Maximum harmonic degree of the input potential in "SHCS_IN_PATH_POT_MTX", 
 * "SHCS_IN_PATH_POT_GFC" and "SHCS_IN_PATH_POT_TBL" */
#undef SHCS_NMAX_POT
#define SHCS_NMAX_POT (10)


/* Path to input spherical harmonic coefficients of the potential in the matrix 
 * format */
#undef SHCS_IN_PATH_POT_MTX
#define SHCS_IN_PATH_POT_MTX "../data/input/EGM96-degree10-mtx.txt"


/* Path to input spherical harmonic coefficients of the potential in the gfc 
 * format */
#undef SHCS_IN_PATH_POT_GFC
#define SHCS_IN_PATH_POT_GFC "../data/input/EGM96-degree10.gfc"


/* Path to input spherical harmonic coefficients of the potential in the tbl 
 * format */
#undef SHCS_IN_PATH_POT_TBL
#define SHCS_IN_PATH_POT_TBL "../data/input/EGM96-degree10-tbl.txt"


/* Path to input spherical harmonic coefficients of the potential in the binary 
 * format */
#undef SHCS_OUT_PATH_POT_BIN
#define SHCS_OUT_PATH_POT_BIN "../data/output/EGM96-degree10.shcs"


/* Path to input spherical harmonic coefficients of the potential in the mtx 
 * format */
#undef SHCS_OUT_PATH_POT_MTX
#define SHCS_OUT_PATH_POT_MTX "../data/output/EGM96-degree10-mtx.txt"


/* Path to output spherical harmonic coefficients of the potential in the tbl 
 * n format */
#undef SHCS_OUT_PATH_POT_TBL_N
#define SHCS_OUT_PATH_POT_TBL_N "../data/output/EGM96-degree10-tbl-n.txt"


/* Path to output spherical harmonic coefficients of the potential in the tbl 
 * m format */
#undef SHCS_OUT_PATH_POT_TBL_M
#define SHCS_OUT_PATH_POT_TBL_M "../data/output/EGM96-degree10-tbl-m.txt"


/* Size of char arrays for long string */
#undef NSTR_LONG
#define NSTR_LONG (2024)


/* Size of char arrays for short string */
#undef NSTR_SHORT
#define NSTR_SHORT (64)


/* Suffix of the reference data files to be loaded. */
#undef FTYPE
#define FTYPE ".txt"


/* Maximum degree for harmonic analysis and synthesis (except for the area-mean
 * values on irregular surfaces) */
#undef NMAX
#define NMAX (4)


/* Auxiliary maximum degree to synthesize area-mean values on irregular
 * surfaces */
#undef NMAX2
#define NMAX2 (15)


/* Maximum degrees to test dynamical switching and loop unrolling in point
 * synthesis and analysis */
#undef NMAX_DS_MIN
#define NMAX_DS_MIN (150)
#undef NMAX_DS_MAX
#define NMAX_DS_MAX (154)


/* Number of radial layers to test for the solid synthesis and analysis. */
#undef NDELTAR
#define NDELTAR (2)


/* Distances between the radial layers in metres. */
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
#if GENREF
#   if CHARM_FLOAT
#      define FOLDER "./genref-output/single"
#   elif CHARM_QUAD
#      define FOLDER "./genref-output/quad"
#   else
#      define FOLDER "./genref-output/double"
#   endif
#else
#   if CHARM_FLOAT
#      define FOLDER "../data/tests/single"
#   elif CHARM_QUAD
#      define FOLDER "../data/tests/quad"
#   else
#      define FOLDER "../data/tests/double"
#   endif
#endif


/* Format to print/write floating points */
#if CHARM_FLOAT
#   define FORMAT "%0.7e"
#elif CHARM_QUAD
#   define FORMAT "%0.34Qe"
#else
#   define FORMAT "%0.16e"
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


/* Logical value representing "equal" */
#undef EQ
#define EQ 1


/* Logical value representing "not equal" */
#undef NEQ
#define NEQ 0


/* A part of the string to be printed for valid function calls */
#undef VALID
#define VALID "valid"


/* A part of the string to be printed for invalid function calls */
#undef INVALID
#define INVALID "invalid"


#endif

