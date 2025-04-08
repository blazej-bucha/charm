/* This header file is not a part of API.
 *
 * Defines some symbolic constants used throughout the tests. */


#ifndef __PARAMETERS_H__
#define __PARAMETERS_H__


#include <config.h>


/* Maximum harmonic degree of the input topography in "SHCS_IN_PATH_TOPO_MTX" 
 * */
#undef SHCS_NMAX_TOPO
#define SHCS_NMAX_TOPO (4UL)


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
#define SHCS_NMAX_POT (10UL)


/* Path to input spherical harmonic coefficients of the potential in the matrix 
 * format */
#undef SHCS_IN_PATH_POT_MTX
#define SHCS_IN_PATH_POT_MTX "../data/input/EGM96-degree10-mtx.txt"


/* Path to input spherical harmonic coefficients of the potential in the gfc 
 * format */
#undef SHCS_IN_PATH_POT_GFC
#define SHCS_IN_PATH_POT_GFC "../data/input/EGM96-degree10.gfc"


/* Path to input time variable spherical harmonic coefficients of the potential
 * in the gfc format */
#undef SHCS_IN_PATH_POT_GFC_TVG
#define SHCS_IN_PATH_POT_GFC_TVG "../data/input/ITSG-Grace2018s.gfc"


/* Path to input spherical harmonic coefficients of the potential in the tbl 
 * format */
#undef SHCS_IN_PATH_POT_TBL
#define SHCS_IN_PATH_POT_TBL "../data/input/EGM96-degree10-tbl.txt"


/* Path to input spherical harmonic coefficients of the potential in the dov
 * format */
#undef SHCS_IN_PATH_POT_DOV
#define SHCS_IN_PATH_POT_DOV "../data/input/EGM96-degree10-dov.txt"


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


/* Path to output spherical harmonic coefficients of the potential in the dov
 * n format */
#undef SHCS_OUT_PATH_POT_DOV_N
#define SHCS_OUT_PATH_POT_DOV_N "../data/output/EGM96-degree10-dov-n.txt"


/* Path to output spherical harmonic coefficients of the potential in the dov
 * m format */
#undef SHCS_OUT_PATH_POT_DOV_M
#define SHCS_OUT_PATH_POT_DOV_M "../data/output/EGM96-degree10-dov-m.txt"


/* Multiplication factors to test rescaling of spherical harmonics coefficients
 * */
#undef SHCS_RESCALE_MU_FACTOR
#define SHCS_RESCALE_MU_FACTOR (PREC(1.1))
#undef SHCS_RESCALE_R_FACTOR
#define SHCS_RESCALE_R_FACTOR (PREC(1.5))


/* Size of char arrays for long string */
#undef NSTR_LONG
#define NSTR_LONG (1024)


/* Size of char arrays for short string */
#undef NSTR_SHORT
#define NSTR_SHORT (64)


/* Suffix of the reference data files to be loaded. */
#undef FTYPE
#define FTYPE ".txt"


/* Maximum degree for harmonic analysis and synthesis (except for the area-mean
 * values on irregular surfaces) */
#undef NMAX
#define NMAX (4UL)


/* Auxiliary maximum degree to synthesize area-mean values on irregular
 * surfaces */
#undef NMAX2
#define NMAX2 (15UL)


/* Maximum degree for quadrature tests */
#undef NMAX_QUAD
#define NMAX_QUAD (10UL)


/* Maximum degrees to test dynamical switching and loop unrolling in point
 * synthesis and analysis */
#undef NMAX_DS_MIN
#define NMAX_DS_MIN (150UL)
#undef NMAX_DS_MAX
#define NMAX_DS_MAX (154UL)


/* Maximum harmonic degree to check the distribution of "charm_shc" and
 * "charm_point".  Cannot be too low, because some of the algorithms to
 * distribute spherical harmonic coefficients across MPI processes are not
 * completely robust. */
#undef NMAX_MPI
#define NMAX_MPI (50UL)


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
#define NSCTR (12)


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


/* Epoch for time variable "gfc" models */
#undef TVG_EPOCH
#define TVG_EPOCH "20020305.1245"


#if HAVE_MPI
/* Maximum number of chunks to split spherical harmonic coefficients at
 * a single process */
#undef NCHUNK_MAX
#define NCHUNK_MAX 3
#endif


/* Path to input spherical harmonic coefficients of the lunar shape */
#undef SHCS_IN_PATH_MOON_SHAPE
#define SHCS_IN_PATH_MOON_SHAPE "../data/input/MoonTopo2600p_to10-tbl.txt"


/* Paths to input spherical harmonic coefficients of the lunar crust density */
#undef SHCS_IN_PATH_MOON_DENSITY0
#define SHCS_IN_PATH_MOON_DENSITY0 "../data/input/moon-rho0_to10.tbl"
#undef SHCS_IN_PATH_MOON_DENSITY1
#define SHCS_IN_PATH_MOON_DENSITY1 "../data/input/moon-rho1_to10.tbl"


/* Order of the polynomial density coefficients */
#undef GFM_MOON_IMAX
#define GFM_MOON_IMAX (1U)


/* Maximum harmonic degree of the lunar shape in "SHCS_IN_PATH_MOON_SHAPE" */
#undef SHCS_NMAX_MOON_SHAPE
#define SHCS_NMAX_MOON_SHAPE (10UL)


/* Maximum harmonic degree of the lunar crust density in
 * "SHCS_IN_PATH_MOON_DENSITY*" */
#undef SHCS_NMAX_MOON_DENSITY0
#define SHCS_NMAX_MOON_DENSITY0 (10UL)
#undef SHCS_NMAX_MOON_DENSITY1
#define SHCS_NMAX_MOON_DENSITY1 (10UL)


/* Moon's Bjerhammar sphere */
#undef MOON_RREF
#define MOON_RREF (PREC(1728000.0))


/* Newton's gravitational constant ("kg^-1 * m^3 * s^-2") */
#undef GRAV_CONST
#define GRAV_CONST (PREC(6.67430e-11))


/* Mass of the Moon ("kg") */
#undef MOON_MASS
#define MOON_MASS (PREC(7.346e22))


/* Minimum and maximum topography powers in GFM */
#undef GFM_MOON_PMIN
#define GFM_MOON_PMIN (1U)
#undef GFM_MOON_PMAX
#define GFM_MOON_PMAX (4U)


/* Maximum harmonic degree of the implied gravitational potential in SGFM */
#undef GFM_NMAX_POTENTIAL
#define GFM_NMAX_POTENTIAL (10UL)


/* Constant density of the lunar crust */
#undef GFM_MOON_DENSITY_CONST
#define GFM_MOON_DENSITY_CONST (PREC(2550.0))


/* Minimum and maximum orders of the radial derivatives of the output
 * gravitational quantity with cap-modified SGFM */
#undef GFM_CAP_KMIN
#define GFM_CAP_KMIN (0U)
#undef GFM_CAP_KMAX
#define GFM_CAP_KMAX (2U)


/* Integration radius with cap-modified SGFM */
#undef GFM_CAP_PSI0
#define GFM_CAP_PSI0 (PREC(1.0))


/* Height of the evaluation sphere above the Bjerhammar sphere with
 * cap-modified SGFM */
#undef GFM_CAP_HEIGHT
#define GFM_CAP_HEIGHT (PREC(25000.0))


/* Number of bits to represent significants of floating point numbers used to
 * compute truncation coefficients */
#undef GFM_CAP_NBITS
#define GFM_CAP_NBITS (256)


/* Maximum order of the potential derivative in cap-modified SGFM */
#undef GFM_CAP_UMAX
#define GFM_CAP_UMAX (2U)


/* Maximum order of polynomial density coefficients in tests of truncation
 * coefficients */
#undef GFM_Q_IMAX
#define GFM_Q_IMAX (3U)


/* Threshold to just whether two "mpfr_t" numbers are nearly equal (similar to
 * "CHARM(glob_threshold)").  The value must be proportional to
 * "GFM_CAP_NBITS", meaning that it must be neither too low, nor too large.  No
 * "PREC" macro should be used here. */
#undef GFM_Q_THRESHOLD
#define GFM_Q_THRESHOLD (1e-70)


#endif

