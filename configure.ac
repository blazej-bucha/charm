dnl                                            -*- Autoconf -*-
dnl Process this file with autoconf to produce a configure script.

dnl ===========================================================================
define(VERSION_SEMANTIC, m4_normalize(m4_include([VERSION])))
AC_INIT([CHarm], [VERSION_SEMANTIC], [blazej.bucha@stuba.sk], [charm], [https://github.com/blazej-bucha/charm])
AC_CONFIG_HEADER([config.h])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([src/shc/shc_init.c])
AM_INIT_AUTOMAKE([foreign])
AM_MAINTAINER_MODE
AC_DISABLE_SHARED


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
if test "$ac_cv_prog_cc_c99" == no ; then
    AC_MSG_ERROR([couldn't set $CC to accept C99])
fi
AC_PROG_INSTALL
AC_PROG_LIBTOOL
AC_PROG_LN_S
AC_PROG_MAKE_SET
dnl ===========================================================================






dnl Check for C types
dnl ===========================================================================
AC_C_CONST
AC_TYPE_SIZE_T
dnl ===========================================================================





dnl Checks for libraries
dnl ===========================================================================
AC_CHECK_LIB([m], [sin], [], [AC_MSG_ERROR([libm not found])])


dnl On some older systems, librt needs to be linked separatelly
AC_SEARCH_LIBS(clock_gettime, rt, [AC_DEFINE([HAVE_CLOCK_GETTIME],1,[Define to 1 if you have the clock_gettime function])],[AC_MSG_WARN([librt not found; the benchmark program will not be able to measure execution times])])
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
],[
#ifndef isfinite
exit(1);
#endif
])], [math_isfinite=yes])
dnl ===========================================================================






dnl Single vs. double vs. quadruple precision
dnl ===========================================================================
AC_ARG_ENABLE([single-precision], [AS_HELP_STRING([--enable-single-precision], [compile CHarm in single precision (float data type for floating point numbers) [default=no]])], sp=$enableval, sp=no)
AC_ARG_ENABLE([double-precision], [AS_HELP_STRING([--enable-double-precision], [compile CHarm in double precision (double data type for floating point numbers) [default=yes]])], dp=$enableval, dp=no)
AC_ARG_ENABLE([quad-precision], [AS_HELP_STRING([--enable-quad-precision], [compile CHarm in quadruple precision (__float128 data type for floating point numbers) [default=no]])], qp=$enableval, qp=no)


if test "$sp" = no -a "$dp" = no -a "$qp" = no; then
    AC_MSG_NOTICE([precision not specified, using default double precision])
    dp="yes"
fi


dnl Single precision
if test "$sp" = yes; then
    if test "$dp" = yes -o "$qp" = yes; then
        AC_MSG_ERROR([more than one precision enabled, please enable one precision only])
    fi


    AC_MSG_NOTICE([compiling in single precision])
    AC_CHECK_LIB([fftw3f], [fftwf_malloc], [], [AC_MSG_ERROR([libfftw3f not found])])
    AC_DEFINE([CHARM_FLOAT], [1], [Define to 1 to compile CHarm in single precision.])
    AC_SUBST([PC], [F])
    AC_SUBST([P], [f])
    AC_SUBST([RDT], [float])
fi


dnl Double precision
if test "$dp" = yes; then
    if test "$sp" = yes -o "$qp" = yes; then
        AC_MSG_ERROR([more than one precision enabled, please enable one precision only])
    fi


    AC_MSG_NOTICE([compiling in double precision])
    AC_CHECK_LIB([fftw3], [fftw_malloc], [], [AC_MSG_ERROR([libfftw3 not found])])
    AC_SUBST([PC], [])
    AC_SUBST([P], [])
    AC_SUBST([RDT], [double])
fi


dnl Quadruple precision
if test "$qp" = yes; then
    if test "$sp" = yes -o "$dp" = yes; then
        AC_MSG_ERROR([more than one precision enabled, please enable one precision only])
    fi


    dnl Get the compiler vendor; gnu is required to compile in quadruple
    dnl precision
    AX_COMPILER_VENDOR
    if test "$ax_cv_c_compiler_vendor" != "gnu"; then
        AC_MSG_ERROR([gcc compiler required to compile in quadruple precision])
    fi


    dnl Get the version of gcc; 4.6 or newer is required.  The line that
    dnl follows, including the "AX_GCC_VERSION" macro, is taken from the FFTW
    dnl installation.
    AX_GCC_VERSION(4,6,0,[],[AC_MSG_ERROR([gcc 4.6 or later required for quad precision support])])



    AC_MSG_NOTICE([compiling in quadruple precision])
    AC_CHECK_LIB([quadmath], [sinq], [], [AC_MSG_ERROR([libquadmath not found])])
    AC_CHECK_LIB([fftw3q], [fftwq_malloc], [], [AC_MSG_ERROR([libfftw3q not found])])
    AC_CHECK_HEADER([quadmath.h], [], [AC_MSG_ERROR([quadmath.h not found])])
    AC_DEFINE([CHARM_QUAD], [1], [Define to 1 to compile CHarm in quadruple precision.])
    AC_SUBST([PC], [Q])
    AC_SUBST([P], [q])
    AC_SUBST([RDT], [__float128])


    dnl On FreeBSD, the "isfinite" macro does not seem to compile/work properly
    dnl with the "__float128" data type.  We therefore need to additionally
    dnl check whether an actual code that employes the macro can be compiled
    dnl and whether it returns the expected value.
    if test "$math_isfinite" == yes ; then
        dnl Reset the value of $math_isfinite
        math_isfinite=no
        AC_RUN_IFELSE([AC_LANG_PROGRAM([
        #include <stdlib.h>
        #include <math.h>
        ],[
        if (!isfinite(1.0q))
            exit(1);
        ])], [math_isfinite=yes])
    fi
fi


if test "$math_isfinite" != no ; then
   AC_DEFINE(HAVE_ISFINITE,[1],[Define to 1 if you have the isfinite macro defined in the <math.h> header file.])
fi
AC_MSG_NOTICE([isfinite macro in math.h found and works correctly... $math_isfinite])


dnl Documentation can only be built if compiling CHarm in double precision
AM_CONDITIONAL(DOC, test "$dp" = "yes")
dnl ===========================================================================






dnl Parallelization via OpenMP
dnl ===========================================================================
AC_ARG_ENABLE([openmp], [AS_HELP_STRING([--enable-openmp], [parallelize CHarm with OpenMP [default=no]])])


if test "$enable_openmp" = yes; then
    AC_MSG_NOTICE([compiling with parallelization enabled])
    AC_OPENMP
    AX_OPENMP([], [AC_MSG_ERROR([couldn't enable OpenMP parallelization])])
    CFLAGS="$CFLAGS $OPENMP_CFLAGS"


    AC_DEFINE([CHARM_PARALLEL], [1], [Define to 1 to parallelize CHarm.])
    AC_CHECK_HEADER([omp.h], [], [AC_MSG_ERROR([omp.h not found])])


    if test "$sp" = yes; then
        AC_CHECK_LIB([fftw3f_omp], [fftwf_init_threads], [], [AC_MSG_WARN([libfftw3f_omp not found; CHarm will be parallelized, but not the FFTW computations; compile FFTW with the OpenMP support to improve the speed])])
    fi
    if test "$dp" = yes; then
        AC_CHECK_LIB([fftw3_omp], [fftw_init_threads], [], [AC_MSG_WARN([libfftw3_omp not found; CHarm will be parallelized, but not the FFTW computations; compile FFTW with the OpenMP support to improve the speed])])
    fi
    if test "$qp" = yes; then
        AC_CHECK_LIB([fftw3q_omp], [fftwq_init_threads], [], [AC_MSG_WARN([libfftw3q_omp not found; CHarm will be parallelized, but not the FFTW computations; compile FFTW with the OpenMP support to improve the speed])])
    fi


    AC_SUBST([OMP], [_omp])
    AC_SUBST([CHARM_LIB], [-lcharm$P$OMP])
else
    AC_MSG_NOTICE([compiling with parallelization disabled])


    AC_SUBST([OMP], [])
    AC_SUBST([CHARM_LIB], [-lcharm$P])
fi


AM_CONDITIONAL(OPENMP, test "$enable_openmp" = "yes")
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
                 tests/Makefile \
                 bench/Makefile \
                 docs/Makefile \
                 charm-tmpl.pc],[],[P=${P}])
AC_OUTPUT
