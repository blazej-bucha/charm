dnl                                            -*- Autoconf -*-
dnl Process this file with autoconf to produce a configure script.

dnl ===========================================================================
define(VERSION_SEMANTIC, m4_normalize(m4_include([VERSION])))
AC_INIT([CHarm], [VERSION_SEMANTIC], [blazej.bucha@stuba.sk], [charm], [https://github.com/blazej-bucha/charm])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([src/shc/shc_init.c])
AM_INIT_AUTOMAKE([foreign])
AM_MAINTAINER_MODE
AC_DISABLE_SHARED
AC_CANONICAL_HOST


dnl Get the current version number of CHarm
AX_SPLIT_VERSION
AC_SUBST([CHARM_MAJOR], ["$AX_MAJOR_VERSION"])
AC_SUBST([CHARM_MINOR], ["$AX_MINOR_VERSION"])
AC_SUBST([CHARM_PATCH], ["$AX_POINT_VERSION"])
dnl ===========================================================================






dnl Checks for programs
dnl ===========================================================================
AC_PROG_CC
AC_PROG_CPP
AC_LANG(C)
AC_PROG_CC_C_O
AC_PROG_CC_C99
if test "$ac_cv_prog_cc_c99" = "yes" ; then
    AC_MSG_ERROR([couldn't set $CC to accept C99])
fi
AC_PROG_INSTALL
AC_PROG_LIBTOOL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AX_COMPILER_VENDOR
dnl ===========================================================================






dnl Check for C keywords
dnl ===========================================================================
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T
dnl ===========================================================================





dnl Checks for libraries
dnl ===========================================================================
AC_CHECK_LIB([m], [sin], [], [AC_MSG_ERROR([libm not found])])


dnl On some older systems, librt needs to be linked separatelly
AC_SEARCH_LIBS(clock_gettime, rt, [AC_DEFINE([HAVE_CLOCK_GETTIME],1, [Define to 1 if you have the clock_gettime function.])], [AC_MSG_WARN([librt not found; the benchmark program will not be able to measure execution times])])
dnl ===========================================================================






dnl Checks for header files
dnl ===========================================================================
AC_CHECK_HEADER([math.h], [], [AC_MSG_ERROR([math.h not found])])
AC_CHECK_HEADER([stdio.h], [], [AC_MSG_ERROR([stdio.h not found])])
AC_CHECK_HEADER([stdlib.h], [], [AC_MSG_ERROR([stdlib.h not found])])
AC_CHECK_HEADER([string.h], [], [AC_MSG_ERROR([string.h not found])])
AC_CHECK_HEADER([time.h], [], [AC_MSG_ERROR([time.h not found])])
AC_CHECK_HEADER([errno.h], [], [AC_MSG_ERROR([errno.h not found])])
AC_CHECK_HEADER([float.h], [], [AC_MSG_ERROR([float.h not found])])
AC_CHECK_HEADER([fftw3.h], [], [AC_MSG_ERROR([fftw3.h not found])])
dnl ===========================================================================






dnl Check the presence of some C macros
dnl ===========================================================================
dnl Do we have the "isfinite" macro in "math.h"?
math_isfinite=no
AC_RUN_IFELSE([AC_LANG_PROGRAM([
#include <stdlib.h>
#include <math.h>
], [
#ifndef isfinite
exit(1);
#endif
])], [math_isfinite=yes])
dnl ===========================================================================






dnl Single vs. double vs. quadruple precision
dnl ===========================================================================
AC_ARG_ENABLE([single-precision], [AS_HELP_STRING([--enable-single-precision], [compile in single precision (float data type for floating point numbers) [default=no]])], sp=$enableval, sp=no)
AC_ARG_ENABLE([double-precision], [AS_HELP_STRING([--enable-double-precision], [compile in double precision (double data type for floating point numbers) [default=yes]])], dp=$enableval, dp=no)
AC_ARG_ENABLE([quad-precision], [AS_HELP_STRING([--enable-quad-precision], [compile in quadruple precision (__float128 data type for floating point numbers) [default=no]])], qp=$enableval, qp=no)


if test "$sp" = "no" -a "$dp" = "no" -a "$qp" = "no"; then
    AC_MSG_NOTICE([precision not specified, using default double precision])
    dp=yes
fi


dnl Single precision
if test "$sp" = "yes"; then
    if test "$dp" = "yes" -o "$qp" = "yes"; then
        AC_MSG_ERROR([more than one precision enabled])
    fi


    AC_MSG_NOTICE([compiling in single precision])
    AC_CHECK_LIB([fftw3f], [fftwf_malloc], [], [AC_MSG_ERROR([libfftw3f not found])])
    AC_DEFINE([CHARM_FLOAT], [1], [Define to 1 to compile in single precision.])
    AC_SUBST([PC], [F])
    AC_SUBST([P], [f])
    AC_SUBST([RDT], [float])
fi


dnl Double precision
if test "$dp" = "yes"; then
    if test "$sp" = "yes" -o "$qp" = "yes"; then
        AC_MSG_ERROR([more than one precision enabled])
    fi


    AC_MSG_NOTICE([compiling in double precision])
    AC_CHECK_LIB([fftw3], [fftw_malloc], [], [AC_MSG_ERROR([libfftw3 not found])])
    AC_SUBST([PC], [])
    AC_SUBST([P], [])
    AC_SUBST([RDT], [double])
fi


dnl Quadruple precision
if test "$qp" = "yes"; then
    if test "$sp" = "yes" -o "$dp" = "yes"; then
        AC_MSG_ERROR([more than one precision enabled])
    fi


    dnl gnu compiler is required to compile in quadruple
    if test "$ax_cv_c_compiler_vendor" != "gnu"; then
        AC_MSG_ERROR([gcc compiler required to compile in quadruple precision])
    fi


    dnl Get the version of gcc; 4.6 or newer is required.  The line that
    dnl follows, including the "AX_GCC_VERSION" macro, is taken from the FFTW
    dnl installation.
    AX_GCC_VERSION(4,6,0, [], [AC_MSG_ERROR([gcc 4.6 or later required for quad precision support])])



    AC_MSG_NOTICE([compiling in quadruple precision])
    AC_CHECK_LIB([quadmath], [sinq], [], [AC_MSG_ERROR([libquadmath not found])])
    AC_CHECK_LIB([fftw3q], [fftwq_malloc], [], [AC_MSG_ERROR([libfftw3q not found])])
    AC_CHECK_HEADER([quadmath.h], [], [AC_MSG_ERROR([quadmath.h not found])])
    AC_DEFINE([CHARM_QUAD], [1], [Define to 1 to compile in quadruple precision.])
    AC_SUBST([PC], [Q])
    AC_SUBST([P], [q])
    AC_SUBST([RDT], [__float128])


    dnl On FreeBSD, the "isfinite" macro does not seem to compile/work properly
    dnl with the "__float128" data type.  We therefore need to additionally
    dnl check whether an actual code that employes the macro can be compiled
    dnl and whether it returns the expected value.
    if test "$math_isfinite" = "yes" ; then
        dnl Reset the value of $math_isfinite
        math_isfinite=no
        AC_RUN_IFELSE([AC_LANG_PROGRAM([
        #include <stdlib.h>
        #include <math.h>
        ], [
        if (!isfinite(1.0q))
            exit(1);
        ])], [math_isfinite=yes])
    fi
fi


if test "$math_isfinite" != "no" ; then
   AC_DEFINE(HAVE_ISFINITE, [1], [Define to 1 if you have the isfinite macro defined in the <math.h> header file.])
fi
AC_MSG_NOTICE([isfinite macro in math.h found and works correctly... $math_isfinite])


dnl Documentation can only be built if compiling CHarm in double precision
AM_CONDITIONAL(DOC, test "$dp" = "yes")
dnl ===========================================================================






dnl SIMD
dnl ===========================================================================
AC_ARG_ENABLE([avx], [AS_HELP_STRING([--enable-avx], [enable AVX instructions to improve performance [default=no]])], avx=$enableval, avx=no)
AC_ARG_ENABLE([avx2], [AS_HELP_STRING([--enable-avx2], [enable AVX2 instructions to improve performance [default=no]])], avx2=$enableval, avx2=no)
AC_ARG_ENABLE([avx-512], [AS_HELP_STRING([--enable-avx-512], [enable AVX-512 instructions to improve performance [default=no]])], avx512=$enableval, avx512=no)


if test "$avx" = "yes" -a "$avx2" = "yes"; then
    AC_MSG_ERROR([more than one type of SIMD instructions enabled])
fi
if test "$avx2" = "yes" -a "$avx512" = "yes"; then
    AC_MSG_ERROR([more than one type of SIMD instructions enabled])
fi
if test "$avx" = "yes" -a "$avx512" = "yes"; then
    AC_MSG_ERROR([more than one type of SIMD instructions enabled])
fi
if test "$avx" = "yes" -o "$avx2" = "yes" -o "$avx512" = "yes"; then
    if test "$dp" != "yes" -a "$sp" != "yes"; then
        AC_MSG_ERROR([SIMD instructions can be enabled in single and double precision only])
    fi
fi


dnl --enable-avx
dnl ---------------------------------------------------------------------------
if test "$avx" = "yes"; then

    AC_MSG_NOTICE([compiling with AVX instructions enabled])
    AC_DEFINE([HAVE_AVX_INSTRUCTIONS], [1], [Define if avx instructions are supported])
    AX_CHECK_COMPILE_FLAG([-mavx], [CFLAGS_SIMD="-mavx"], [AC_MSG_ERROR([unrecognized -mavx compiler flag.  Try using the gcc compiler.])])

fi
dnl ---------------------------------------------------------------------------


dnl --enable-avx2
dnl ---------------------------------------------------------------------------
if test "$avx2" = "yes"; then

    AC_MSG_NOTICE([compiling with AVX2 instructions enabled])
    AC_DEFINE([HAVE_AVX2_INSTRUCTIONS], [1], [Define if avx2 instructions are supported])
    AX_CHECK_COMPILE_FLAG([-mavx2], [CFLAGS_SIMD="-mavx2"], [AC_MSG_ERROR([unrecognized -mavx2 compiler flag.  Try using the gcc compiler.])])

fi
dnl ---------------------------------------------------------------------------


dnl --enable-avx-512
dnl ---------------------------------------------------------------------------
if test "$avx512" = "yes"; then

    AC_MSG_NOTICE([compiling with AVX-512 instructions enabled])
    AC_DEFINE([HAVE_AVX512F_INSTRUCTIONS], [1], [Define if avx512f instructions are supported])
    AX_CHECK_COMPILE_FLAG([-mavx512f], [CFLAGS_SIMD="-mavx512f"], [AC_MSG_ERROR([unrecognized -mavx512f compiler flag.  Try using the gcc compiler.])])
    AX_CHECK_COMPILE_FLAG([-mavx512dq], [CFLAGS_SIMD="-mavx512dq"], [AC_MSG_ERROR([unrecognized -mavx512dq compiler flag.  Try using the gcc compiler.])])

fi
dnl ---------------------------------------------------------------------------


if test "$avx" = "yes" -o "$avx2" = "yes" -o "$avx512" = "yes"; then

    AC_CHECK_FUNCS([posix_memalign], [have_posix_memalign=yes], [have_posix_memalign=no])
    AC_CHECK_FUNCS([aligned_alloc], [have_aligned_alloc=yes], [have_aligned_alloc=no])
    AC_CHECK_HEADER([mm_malloc.h], [AC_DEFINE([HAVE_MM_MALLOC_H], [1], [Define to 1 if you have the <mm_malloc.h> header file.])], [have_mm_malloc_h=no])


    if test "$have_posix_memalign" = "no" -a "$have_aligned_alloc" = "no" -a "$have_mm_malloc_h" = "no"; then
        AC_MSG_ERROR([couldn't find any system function to allocate an aligned block of memory.  Try to compile with SIMD instructions disabled.])
    fi


    AC_CHECK_HEADER([immintrin.h], [], [AC_MSG_ERROR([immintrin.h not found])])


    CFLAGS="$CFLAGS $CFLAGS_SIMD"


    if test "$host_vendor" = "apple"; then
        if test "$ax_cv_c_compiler_vendor" = "gnu" -o "$ax_cv_c_compiler_vendor" = "intel"; then
            AX_CHECK_COMPILE_FLAG([-Wa,-q], [CFLAGS="$CFLAGS -Wa,-q"], [])
            AX_CHECK_COMPILE_FLAG([-Wl,-no_compact_unwind], [CFLAGS="$CFLAGS -Wl,-no_compact_unwind"], [])
        fi
    fi

fi


if test "$avx" = "no" -a "$avx2" = "no" -a "$avx512" = "no"; then
    AC_MSG_NOTICE([compiling with SIMD instructions disabled])
fi
dnl ===========================================================================






dnl Parallelization via OpenMP
dnl ===========================================================================
AC_ARG_ENABLE([openmp], [AS_HELP_STRING([--enable-openmp], [enable OpenMP parallelization [default=no]])])


if test "$enable_openmp" = "yes"; then
    AC_MSG_NOTICE([compiling with OpenMP parallelization enabled])
    AC_OPENMP
    AX_OPENMP([], [AC_MSG_ERROR([couldn't enable OpenMP parallelization])])
    CFLAGS="$CFLAGS $OPENMP_CFLAGS"


    AC_DEFINE([CHARM_PARALLEL], [1], [Define to 1 to parallelize CHarm.])
    AC_CHECK_HEADER([omp.h], [], [AC_MSG_ERROR([omp.h not found])])


    if test "$sp" = "yes"; then
        AC_CHECK_LIB([fftw3f_omp], [fftwf_init_threads], [], [AC_MSG_WARN([libfftw3f_omp not found; CHarm will be parallelized, but not the FFTW computations; compile FFTW with the OpenMP support to improve the speed])])
    fi
    if test "$dp" = "yes"; then
        AC_CHECK_LIB([fftw3_omp], [fftw_init_threads], [], [AC_MSG_WARN([libfftw3_omp not found; CHarm will be parallelized, but not the FFTW computations; compile FFTW with the OpenMP support to improve the speed])])
    fi
    if test "$qp" = "yes"; then
        AC_CHECK_LIB([fftw3q_omp], [fftwq_init_threads], [], [AC_MSG_WARN([libfftw3q_omp not found; CHarm will be parallelized, but not the FFTW computations; compile FFTW with the OpenMP support to improve the speed])])
    fi


    AC_SUBST([OMP], [_omp])
    AC_SUBST([CHARM_LIB], [-lcharm$P$OMP])
else
    AC_MSG_NOTICE([compiling with OpenMP parallelization disabled])


    AC_SUBST([OMP], [])
    AC_SUBST([CHARM_LIB], [-lcharm$P])
fi


AM_CONDITIONAL(OPENMP, test "$enable_openmp" = "yes")
dnl ===========================================================================






# Guards in "config.h"
dnl ===========================================================================
AH_TOP([
#ifndef __CONFIG_H__
#define __CONFIG_H__
])


AH_BOTTOM([
#endif
])
dnl ===========================================================================






AC_CONFIG_FILES([Makefile \
                 charm/charm${P}.h:src/charm.h.in \
                 charm/charm${P}_crd.h:src/crd/crd.h.in \
                 charm/charm${P}_err.h:src/err/err.h.in \
                 charm/charm${P}_integ.h:src/integ/integ.h.in \
                 charm/charm${P}_leg.h:src/leg/leg.h.in \
                 charm/charm${P}_misc.h:src/misc/misc.h.in \
                 charm/charm${P}_sha.h:src/sha/sha.h.in \
                 charm/charm${P}_shc.h:src/shc/shc.h.in \
                 charm/charm${P}_shs.h:src/shs/shs.h.in \
                 charm/charm${P}_xnum.h:src/xnum/xnum.h.in \
                 charm/charm${P}_glob.h:src/glob/glob.h.in \
                 charm/Makefile \
                 src/Makefile \
                 src/crd/Makefile \
                 src/err/Makefile \
                 src/integ/Makefile \
                 src/leg/Makefile \
                 src/misc/Makefile \
                 src/sha/Makefile \
                 src/shc/Makefile \
                 src/shs/Makefile \
                 src/xnum/Makefile \
                 src/glob/Makefile \
                 src/simd/Makefile \
                 tests/Makefile \
                 bench/Makefile \
                 docs/Makefile \
                 m4/Makefile \
                 cookbook/Makefile \
                 data/Makefile \
                 charm-tmpl.pc], [], [P=${P}])
AC_OUTPUT
