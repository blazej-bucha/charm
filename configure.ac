dnl                                            -*- Autoconf -*-
dnl Process this file with autoconf to produce a configure script.

dnl ===========================================================================
define(VERSION_SEMANTIC, m4_normalize(m4_include([VERSION])))
AC_INIT([CHarm],[VERSION_SEMANTIC],[blazej.bucha@stuba.sk],[charm],[https://github.com/blazej-bucha/charm])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([src/shc/shc_calloc.c])
dnl Automake version "1.16.5" or newer is needed to enable the "--with-python*"
dnl flags
AM_INIT_AUTOMAKE([1.16.5 foreign])
AM_MAINTAINER_MODE
AC_CANONICAL_HOST


dnl Get the current version number of CHarm
AC_PROG_SED
AX_SPLIT_VERSION
AC_SUBST([CHARM_MAJOR], ["$AX_MAJOR_VERSION"])
AC_SUBST([CHARM_MINOR], ["$AX_MINOR_VERSION"])
AC_SUBST([CHARM_PATCH], ["$AX_POINT_VERSION"])
dnl If "-" is found in "AX_POINT_VERSION", remove all characters that follow
dnl it, including the "-" character itself.  This is useful to get the
dnl patch version number without the "-dev" string.
AC_MSG_CHECKING([Point version number])
AX_POINT_VERSION_NUMBER=`echo "$AX_POINT_VERSION" | $SED 's/-.*//'`
AC_MSG_RESULT([$AX_POINT_VERSION_NUMBER])
AC_SUBST([CHARM_PATCH_NUMBER], ["$AX_POINT_VERSION_NUMBER"])
dnl ===========================================================================






dnl Checks for programs
dnl ===========================================================================
AC_PROG_CC
AC_PROG_CPP
AC_LANG(C)
AC_PROG_CC_C_O
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AX_COMPILER_VENDOR
LT_INIT([win32-dll])
dnl ===========================================================================






dnl Check for the C99 support
dnl ===========================================================================
dnl Microsoft's MSVC compiler doesn't seem to support C99.  So we do no check
dnl in this case and hope for the best.
if test "$ax_cv_c_compiler_vendor" != "microsoft"; then
    if test "$ac_cv_prog_cc_c99" = "no" ; then
        AC_MSG_ERROR([couldn't set $CC to accept C99])
    fi
fi
dnl ===========================================================================






dnl Check for C keywords
dnl ===========================================================================
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T
if test "$ac_cv_type_size_t" != "yes"; then
    AC_MSG_ERROR([couldn't define size_t data type.])
fi
AC_TYPE_INT32_T
AC_TYPE_INT64_T
dnl ===========================================================================





dnl Checks for libraries
dnl ===========================================================================
dnl Check the "math" library if not using MSVC, which apparently links this
dnl library automatically
if test "$ax_cv_c_compiler_vendor" != "microsoft"; then
    AC_CHECK_LIB([m], [sin], [], [AC_MSG_ERROR([libm not found])])
fi


dnl On some older systems, librt needs to be linked separatelly
AC_SEARCH_LIBS(clock_gettime, rt, [AC_DEFINE([HAVE_CLOCK_GETTIME],1, [Define to 1 if you have the clock_gettime function.])], [AC_MSG_WARN([librt not found; the benchmark program will not be able to measure execution times])])
dnl ===========================================================================






dnl Checks for header files
dnl ===========================================================================
AC_CHECK_HEADER([math.h], [], [AC_MSG_ERROR([math.h not found])])
AC_CHECK_HEADER([stdio.h], [], [AC_MSG_ERROR([stdio.h not found])])
AC_CHECK_HEADER([stdint.h], [], [AC_MSG_ERROR([stdint.h not found])])
AC_CHECK_HEADER([stdlib.h], [], [AC_MSG_ERROR([stdlib.h not found])])
AC_CHECK_HEADER([string.h], [], [AC_MSG_ERROR([string.h not found])])
AC_CHECK_HEADER([ctype.h], [], [AC_MSG_ERROR([ctype.h not found])])
AC_CHECK_HEADER([time.h], [], [AC_MSG_ERROR([time.h not found])])
AC_CHECK_HEADER([errno.h], [], [AC_MSG_ERROR([errno.h not found])])
AC_CHECK_HEADER([float.h], [], [AC_MSG_ERROR([float.h not found])])
AC_CHECK_HEADER([limits.h], [], [AC_MSG_ERROR([limits.h not found])])
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
    AC_CHECK_LIB([fftw3], [fftwf_version], [AC_DEFINE([HAVE_FFTW_VERSION], [1], [Define to 1 if fftwf_version is available in the fftw3f lib.])], [])
    AC_DEFINE([CHARM_FLOAT], [1], [Define to 1 to compile in single precision.])
    AC_SUBST([PC], [F])
    AC_SUBST([P], [f])
    AC_SUBST([RDT], [float])


    if test "$ac_cv_c_int32_t" != "yes"; then
        AC_MSG_ERROR([couldn't define a signed 32 bits wide integer.])
    fi
fi


dnl Double precision
if test "$dp" = "yes"; then
    if test "$sp" = "yes" -o "$qp" = "yes"; then
        AC_MSG_ERROR([more than one precision enabled])
    fi


    AC_MSG_NOTICE([compiling in double precision])
    AC_CHECK_LIB([fftw3], [fftw_malloc], [], [AC_MSG_ERROR([libfftw3 not found])])
    AC_CHECK_LIB([fftw3], [fftw_version], [AC_DEFINE([HAVE_FFTW_VERSION], [1], [Define to 1 if fftw_version is available in the fftw3 lib.])], [])
    AC_SUBST([PC], [])
    AC_SUBST([P], [])
    AC_SUBST([RDT], [double])


    if test "$ac_cv_c_int64_t" != "yes"; then
        AC_MSG_ERROR([couldn't define a signed 64 bits wide integer.])
    fi
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
    AC_CHECK_LIB([fftw3], [fftwq_version], [AC_DEFINE([HAVE_FFTW_VERSION], [1], [Define to 1 if fftwq_version is available in the fftw3q lib.])], [])
    AC_CHECK_HEADER([quadmath.h], [], [AC_MSG_ERROR([quadmath.h not found])])
    AC_DEFINE([CHARM_QUAD], [1], [Define to 1 to compile in quadruple precision.])
    AC_SUBST([PC], [Q])
    AC_SUBST([P], [q])
    AC_SUBST([RDT], [__float128])


    if test "$ac_cv_c_int64_t" != "yes"; then
        AC_MSG_ERROR([couldn't define a signed 64 bits wide integer.])
    fi


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


if test "$math_isfinite" != "no"; then
   AC_DEFINE(HAVE_ISFINITE, [1], [Define to 1 if you have the isfinite macro defined in the <math.h> header file.])
fi
AC_MSG_NOTICE([isfinite macro in math.h found and works correctly... $math_isfinite])


dnl Documentation can only be built if compiling CHarm in double precision
AM_CONDITIONAL(DOC, test "$dp" = "yes")
dnl ===========================================================================






dnl SIMD
dnl ===========================================================================
AC_ARG_ENABLE([avx], [AS_HELP_STRING([--enable-avx], [enable AVX instructions on x86_64 to improve performance [default=no]])], avx=$enableval, avx=no)
AC_ARG_ENABLE([avx2], [AS_HELP_STRING([--enable-avx2], [enable AVX2 instructions on x86_64 to improve performance [default=no]])], avx2=$enableval, avx2=no)
AC_ARG_ENABLE([avx-512], [AS_HELP_STRING([--enable-avx-512], [enable AVX-512 on x86_64 instructions to improve performance [default=no]])], avx512=$enableval, avx512=no)
AC_ARG_ENABLE([neon], [AS_HELP_STRING([--enable-neon], [enable NEON instructions on ARM64 to improve performance [default=no]])], neon=$enableval, neon=no)


if test "$avx" = "yes" -a "$avx2" = "yes"; then
    AC_MSG_ERROR([more than one type of SIMD instructions enabled])
fi
if test "$avx" = "yes" -a "$avx512" = "yes"; then
    AC_MSG_ERROR([more than one type of SIMD instructions enabled])
fi
if test "$avx" = "yes" -a "$neon" = "yes"; then
    AC_MSG_ERROR([more than one type of SIMD instructions enabled])
fi
if test "$avx2" = "yes" -a "$avx512" = "yes"; then
    AC_MSG_ERROR([more than one type of SIMD instructions enabled])
fi
if test "$avx2" = "yes" -a "$neon" = "yes"; then
    AC_MSG_ERROR([more than one type of SIMD instructions enabled])
fi
if test "$avx512" = "yes" -a "$neon" = "yes"; then
    AC_MSG_ERROR([more than one type of SIMD instructions enabled])
fi
if test "$avx" = "yes" -o "$avx2" = "yes" -o "$avx512" = "yes" -o "$neon" = "yes"; then
    if test "$qp" = "yes"; then
        AC_MSG_ERROR([SIMD instructions are not supported in quadruple precision])
    fi
fi


dnl ---------------------------------------------------------------------------
simd_flag1="no"
simd_flag2="no"


dnl --enable-avx
if test "$avx" = "yes"; then
    AC_MSG_NOTICE([compiling with AVX instructions enabled])
    AC_DEFINE([HAVE_AVX], [1], [Define to 1 to enable avx instructions])
    if test "$ax_cv_c_compiler_vendor" = "microsoft"; then
        simd_flag1=/arch:AVX
    else
        simd_flag1=-mavx
    fi
fi


dnl --enable-avx2
if test "$avx2" = "yes"; then
    AC_MSG_NOTICE([compiling with AVX2 instructions enabled])
    AC_DEFINE([HAVE_AVX2], [1], [Define to 1 to enable avx2 instructions])
    if test "$ax_cv_c_compiler_vendor" = "microsoft"; then
        simd_flag1=/arch:AVX2
    else
        simd_flag1=-mavx2
    fi
fi


dnl --enable-avx-512
if test "$avx512" = "yes"; then
    AC_MSG_NOTICE([compiling with AVX-512 instructions enabled])
    AC_DEFINE([HAVE_AVX512F], [1], [Define to 1 to enable avx512f instructions])
    if test "$ax_cv_c_compiler_vendor" = "microsoft"; then
        simd_flag1=/arch:AVX512
    else
        simd_flag1=-mavx512f
        simd_flag2=-mavx512dq
    fi
fi


dnl --enable-neon
if test "$neon" = "yes"; then
    AC_MSG_NOTICE([compiling with NEON instructions enabled])
    AC_DEFINE([HAVE_NEON], [1], [Define to 1 to enable neon instructions])
    dnl With NEON, we do not add compiler flags due to the large variety of the
    dnl flags.  The flags must be specified by the user using the "CFLAGS"
    dnl variable.
fi
dnl ---------------------------------------------------------------------------


if test "$avx" = "yes" -o "$avx2" = "yes" -o "$avx512" = "yes" -o "$neon" = "yes"; then

    if test "$neon" != "yes"; then
        AX_CHECK_COMPILE_FLAG([$simd_flag1], [CFLAGS_SIMD="$simd_flag1"], [AC_MSG_ERROR([unrecognized $simd_flag1 compiler flag.  Try using the gcc compiler.])])
        if test "$simd_flag2" != "no"; then
            AX_CHECK_COMPILE_FLAG([$simd_flag2], [CFLAGS_SIMD="$CFLAGS_SIMD $simd_flag2"], [AC_MSG_ERROR([unrecognized $simd_flag2 compiler flag.  Try using the gcc compiler.])])
        fi
    fi

    have_posix_memalign=no
    have_aligned_alloc=no
    have_mm_malloc_h=no
    have__aligned_malloc=no
    if test "$ax_cv_c_compiler_vendor" = "microsoft"; then
        AC_CHECK_FUNCS([_aligned_malloc], [have__aligned_malloc=yes], [have__aligned_malloc=no])
    else
        AC_CHECK_FUNCS([posix_memalign], [have_posix_memalign=yes], [have_posix_memalign=no])
        AC_CHECK_FUNCS([aligned_alloc], [have_aligned_alloc=yes], [have_aligned_alloc=no])
        AC_CHECK_HEADER([mm_malloc.h], [AC_DEFINE([HAVE_MM_MALLOC_H], [1], [Define to 1 if you have the <mm_malloc.h> header file.])], [have_mm_malloc_h=no])
    fi


    if test "$have_posix_memalign" = "no" -a "$have_aligned_alloc" = "no" -a "$have_mm_malloc_h" = "no" -a "$have__aligned_malloc" = "no"; then
        AC_MSG_ERROR([couldn't find any system function to allocate an aligned block of memory.  Try to compile with SIMD instructions disabled.])
    fi


    if test "$neon" = "yes"; then
        AC_CHECK_HEADER([arm_neon.h], [], [AC_MSG_ERROR([arm_neon.h not found])])
    else
        AC_CHECK_HEADER([immintrin.h], [], [AC_MSG_ERROR([immintrin.h not found])])
    fi
    if test "$ax_cv_c_compiler_vendor" = "microsoft"; then
        AC_CHECK_HEADER([malloc.h], [], [AC_MSG_ERROR([malloc.h not found])])
    fi


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
    if test "$ax_cv_c_compiler_vendor" = "microsoft"; then
        dnl Without "/openmp:llvm", the Microsoft's MSVC compiler does not
        dnl accept data type other than "int" as the main loop variable
        CFLAGS="$CFLAGS /openmp:llvm"
    fi


    AC_DEFINE([CHARM_OPENMP], [1], [Define to 1 to parallelize CHarm using OpenMP.])
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


    AC_SUBST([CHARM_LIB], [-lcharm])
else
    AC_MSG_NOTICE([compiling with OpenMP parallelization disabled])


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






dnl Python
dnl ===========================================================================
AC_ARG_ENABLE([python], [AS_HELP_STRING([--enable-python], [enable Python interface (PyHarm) [default=no]])], py=$enableval, py=no)


if test "$py" = "yes"; then
    AC_MSG_NOTICE([enabling Python interface])
    AM_PATH_PYTHON([3.6], [], [AC_MSG_ERROR([couldn't find suitable Python interpreter])])


    if test "$enable_shared" != yes; then
        AC_MSG_ERROR([shared library is required to build the Python interface; add the --enable-shared flag])
    fi


    if test "$sp" = "yes"; then
        AC_SUBST([NP_FLOAT],  [float32])
        AC_SUBST([CT_FLOAT],  [float])
        AC_SUBST([PY_FLOAT],  [])
    fi
    if test "$dp" = "yes"; then
        AC_SUBST([NP_FLOAT],  [float64])
        AC_SUBST([CT_FLOAT],  [double])
        AC_SUBST([PY_FLOAT],  [float, ]) dnl The comma is not a typo
    fi
    if test "$qp" = "yes"; then
        AC_MSG_ERROR([Python interface is not supported for quadruple precision])
    fi


    AC_SUBST([OBJDIR], ["$objdir"])
    eval shared_ext="$shrext_cmds"
    AC_SUBST([SHARED_EXT], ["$shared_ext"])

fi


AM_CONDITIONAL(PY, test "$py" = "yes")


dnl Python interface must be enabled to build the documentation
AM_CONDITIONAL(DOC_PY, test "$py" = "yes")


dnl For single precision version of PyHarm, a symlink "wrap/pyharmf" is create
dnl for "wrap/pyharm"
AM_CONDITIONAL(SINGLE_PY, test "$py" = "yes" -a "$sp" = "yes")
dnl ===========================================================================






dnl Fix bugs in MSVC
dnl ===========================================================================
AC_ARG_ENABLE([msvc-underscore-patch], [AS_HELP_STRING([--enable-msvc-underscore-patch], [enable on Windows if MSVC throws the "'fftw__complex': undeclared identifier" error on compilation; this is a bug in MSVC, which evidently cannot properly expand a valid C macro [default=no]])], msvc_underscore_patch=$enableval, msvc_underscore_patch=no)
if test "$msvc_underscore_patch" == "yes"; then
    AC_MSG_NOTICE([patching the MSVC underscore bug enabled])
    AC_DEFINE(MSVC_UNDERSCORE_PATCH, [1], [Define to 1 to patch the MSVC underscore bug on Windows])
fi
dnl ===========================================================================






dnl Export some useful user-defined options from the "./configure" call into
dnl "config.h".
dnl ===========================================================================
AC_DEFINE_UNQUOTED(_CHARM_CC,["$CC"],[Used compiler.])
AC_DEFINE_UNQUOTED(_CHARM_CFLAGS,["$CFLAGS"],[CFLAGS used to compile CHarm.])
dnl "CPPFLAGS" and "LDFLAGS" may contain "\".  In that case,
dnl "AC_DEFINE_UNQUOTED" perfoms backslash escaping which is not desired here.
dnl So before each "\" in these variables, we add yet another "\", so that
dnl after the backslash escaping, we get to the original point.  This is useful
dnl especially for Windows paths like "C:\Users\isaac\Documents\fftw".
CPPFLAGS_MOD=`echo "$CPPFLAGS" | $SED 's/\\\\/\\\\\\\\/g'`
LDFLAGS_MOD=`echo "$LDFLAGS" | $SED 's/\\\\/\\\\\\\\/g'`
AC_DEFINE_UNQUOTED(_CHARM_CPPFLAGS,["$CPPFLAGS_MOD"],[CPPFLAGS used to compile CHarm.])
AC_DEFINE_UNQUOTED(_CHARM_LDFLAGS,["$LDFLAGS_MOD"],[LDFLAGS used to compile CHarm.])
AC_DEFINE_UNQUOTED(_CHARM_LIBS,["$LIBS"],[LDFLAGS used to compile CHarm.])
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
                 wrap/pyharm/Makefile \
                 charm-tmpl.pc \
                 setup.py \
                 wrap/pyharm/__init__.py \
                 wrap/pyharm/_data_types.py], [], [P=${P}])


AC_OUTPUT
