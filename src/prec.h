/*
 * Defines macros, symbolic constants, etc., that depend on the numerical
 * precision but are not supposed to be a part of the API.
 *
 * The numerical values of the symbolic constants related to the "X-numbers"
 * ("BIG", ...) were obtained with 50-digit arithmetics (in "Python" using
 * "mpmath").  The advantage of the literal definition is that the constants do
 * not need to be re-computed by the "pow" function over and over again when
 * the functions in which they occur are called multiple times.
 *
 * */


#ifndef __PREC_H__
#define __PREC_H__


/* Header files required by this module */
/* ------------------------------------------------------------------------- */
#include <config.h>
/* ------------------------------------------------------------------------- */


#ifndef DLL_EXPORT
#   define DLL_EXPORT
#endif


#undef CAT
#define CAT(x, y) x ## y
#undef CAT2
#define CAT2(x, y, z) x ## y ## z


#undef CHARM
#undef CHARM_SUFFIX
#undef REAL
#undef INT
#undef PREC
#undef FFTW
#undef FFTWC
#undef FFTW3_OMP


#undef SIN
#undef COS
#undef ASIN
#undef FABS
#undef POW
#undef SQRT
#undef PI
#undef PI_2
#undef ROOT3
#undef EPS
#undef BIG
#undef BIGI
#undef BIGS
#undef BIGSI


#undef STR2REAL


/* At first, let's check if one precision only is defined in "config.h".  The
 * configure script does not allow this, but if CHarm is compiled without the
 * autotools, this might perhaps happen. */
#if CHARM_FLOAT && CHARM_QUAD
#   error "One precision only must be defined in config.h."
#endif


#include <stdint.h>


/* ------------------------------------------------------------------------- */
/* Single precision */
/* ------------------------------------------------------------------------- */
#if CHARM_FLOAT


#   include <float.h>
#   include "../charm/charmf.h"


#   define CHARM(x)         CAT(charmf_, x)
#   define CHARM_SUFFIX     "f"
#   define FFTW(x)          CAT(fftwf_, x)
#   if MSVC_UNDERSCORE_PATCH
        /* Microsoft's MSVC compiler incorrectly expands "FFTW(complex)" to
         * "fftw__complex" if "FFTW" is defined as above and "math.h" is
         * included, so we need this stupid patch and to always use
         * "FFTWC(complex)" instead of "FFTW(complex)". */
#       define FFTWC(x)     CAT(fftwf, x)
#   else
#       define FFTWC        FFTW
#   endif
#   define REAL      float
#   define INT       int32_t
#   define PREC(x)   CAT(x, f)
#   if HAVE_LIBFFTW3F_OMP
#       define FFTW3_OMP 1
#   endif


#   define SIN        sinf
#   define COS        cosf
#   define ASIN       asinf
#   define FABS       fabsf
#   define POW        powf
#   define SQRT       sqrtf
#   define PI         (3.141592653589793238462643383279502884f) /* pi */
#   define PI_2       (1.570796326794896619231321691639751442f) /* pi / 2 */
#   define ROOT3      (1.732050807568877293527446341505872367f) /* sqrt(3) */
#   define EPS        FLT_EPSILON


    /* Constants related to "X-numbers" */
    /* --------------------------------------------------------------------- */

    /*
     * Radix of the ``X-numbers``.  The value can be obtained as
     *
     *      const float BIG = powf(2.0f, 120);
     *
     * */
#   define BIG (1329227995784915872903807060280344576.0f)


    /*
     * ``1 / BIG``.  The value can be obtained as
     *
     *      const float BIGI = powf(2.0f, -120);
     *
     * */
#   define BIGI (7.5231638452626400509999138382223723380394595633413601e-37f)


    /*
     * ``BIG^(1 / 2)``.  The value can be obtained as
     *
     *      const float BIGS = powf(2.0f, 120 / 2);
     *
     * */
#   define BIGS (1152921504606846976.0f)


    /*
     * ``1 / BIGS``.  The value can be obtained as
     *
     *      const float BIGSI = powf(2.0f, -120 / 2);
     *
     * */
#   define BIGSI (8.67361737988403547205962240695953369140625e-19f)
    /*---------------------------------------------------------------------- */


#   define STR2REAL strtof


/* ------------------------------------------------------------------------- */
/* Quadruple precision */
/* ------------------------------------------------------------------------- */
#elif CHARM_QUAD


#   include <quadmath.h>
#   include "../charm/charmq.h"


#   define CHARM(x)         CAT(charmq_, x)
#   define CHARM_SUFFIX     "q"
#   define FFTW(x)          CAT(fftwq_, x)
#   if MSVC_UNDERSCORE_PATCH
        /* See above for the explanation of this */
#       define FFTWC(x)     CAT(fftwq, x)
#   else
#       define FFTWC        FFTW
#   endif
#   define REAL      __float128
#   define INT       int64_t
#   define PREC(x)   CAT(x, q)
#   if HAVE_LIBFFTW3Q_OMP
#       define FFTW3_OMP 1
#   endif


#   define SIN        sinq
#   define COS        cosq
#   define ASIN       asinq
#   define FABS       fabsq
#   define POW        powq
#   define SQRT       sqrtq
#   define PI         (3.141592653589793238462643383279502884q) /* pi */
#   define PI_2       (1.570796326794896619231321691639751442q) /* pi / 2 */
#   define ROOT3      (1.732050807568877293527446341505872367q) /* sqrt(3) */
#   define EPS        FLT128_EPSILON


    /* Constants related to "X-numbers" */
    /* --------------------------------------------------------------------- */

    /*
     * Radix of the ``X-numbers``.  The value can be obtained as
     *
     *      const __float128 BIG = powq(2.0q, 16000);
     *
     * */
#   define BIG (3.0194693372392275795306584466152797092952625113753e+4816q)


    /*
     * ``1 / BIG``.  The value can be obtained as
     *
     *      const __float128 BIGI = powq(2.0q, -16000);
     *
     * */
#   define BIGI (3.3118402219455015713947284908357856923046345036674e-4817q)


    /*
     * ``BIG^(1 / 2)``.  The value can be obtained as
     *
     *      const __float128 BIGS = powq(2.0q, 16000 / 2);
     *
     * */
#   define BIGS (1.7376620319380945659998244594943562706193978610012e+2408q)


    /*
     * ``1 / BIGS``.  The value can be obtained as
     *
     *      const __float128 BIGSI = powq(2.0q, -16000 / 2);
     *
     * */
#   define BIGSI (5.7548590095201303475301830222419605375861964661775e-2409q)
    /*---------------------------------------------------------------------- */


#   define STR2REAL strtoflt128


/* ------------------------------------------------------------------------- */
/* Double precision */
/* ------------------------------------------------------------------------- */
#else


#   include <float.h>
#   include "../charm/charm.h"


#   define CHARM(x)         CAT(charm_, x)
#   define CHARM_SUFFIX     ""
#   define FFTW(x)          CAT(fftw_, x)
#   if MSVC_UNDERSCORE_PATCH
        /* See above for the explanation of this */
#       define FFTWC(x)     CAT(fftw, x)
#   else
#       define FFTWC        FFTW
#   endif
#   define REAL      double
#   define INT       int64_t
#   define PREC(x)   CAT(x,)
#   if HAVE_LIBFFTW3_OMP
#       define FFTW3_OMP 1
#   endif


#   define SIN        sin
#   define COS        cos
#   define ASIN       asin
#   define FABS       fabs
#   define POW        pow
#   define SQRT       sqrt
#   define PI         (3.14159265358979323846) /* pi */
#   define PI_2       (1.57079632679489661923) /* pi / 2 */
#   define ROOT3      (1.73205080756887729353) /* sqrt(3) */
#   define EPS        DBL_EPSILON


    /* Constants related to "X-numbers" */
    /* --------------------------------------------------------------------- */

    /*
     * Radix of the ``X-numbers``.  The value can be obtained as
     *
     *      const double BIG = pow(2.0, 960);
     *
     * */
#   define BIG (9.7453140113999990803533823878751883108762268575950075e+288)


    /*
     * ``1 / BIG``.  The value can be obtained as
     *
     *      const double BIGI = pow(2.0, -960);
     *
     * */
#   define BIGI (1.0261342003245940623340073222912829685219053068308029e-289)


    /*
     * ``BIG^(1 / 2)``.  The value can be obtained as
     *
     *      const double BIGS = pow(2.0, 960 / 2);
     *
     * */
#   define BIGS (3.1217485503159922313815972297931663057485981426649712e+144)


    /*
     * ``1 / BIGS``.  The value can be obtained as
     *
     *      const double BIGSI = pow(2.0, -960 / 2);
     *
     * */
#   define BIGSI (3.203332952292961479087336344218362393174372768098368e-145)
    /*---------------------------------------------------------------------- */


#   define STR2REAL strtod


#endif


#endif
