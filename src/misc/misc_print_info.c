/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#if HAVE_MPI
#   include <mpi.h>
#endif
#include "../prec.h"
#include "misc_buildopt.h"
/* ------------------------------------------------------------------------- */






void CHARM(misc_print_info)(void)
{
    printf(


#ifdef PACKAGE_NAME
    "Name: " PACKAGE_NAME "\n"
#endif
    "Description: C library to work with spherical harmonics up to "
        "almost arbitrarily high degrees\n"
#ifdef PACKAGE_VERSION
    "Version: " PACKAGE_VERSION "\n"
#endif
    "Compiled: " __DATE__ " " __TIME__"\n"
#ifdef PACKAGE_URL
    /* Old versions of autotools do not seem to define "PACKAGE_URL". */
    "URL: " PACKAGE_URL "\n"
#endif
#ifdef PACKAGE_BUGREPORT
    "Bug-report: " PACKAGE_BUGREPORT "\n"
#endif
    "License of the source code: The 3-Clause BSD License\n"
    "License of the compiled library: GNU General Public License version 2 "
        "or any later\n\n"


    );


    printf("Precision: ");
    int ret = CHARM(misc_buildopt_precision)();
    if (ret == BUILDOPT_PRECISION_SINGLE)
        printf("single");
    else if (ret == BUILDOPT_PRECISION_DOUBLE)
        printf("double");
    else if (ret == BUILDOPT_PRECISION_QUAD)
        printf("quadruple");
    else
        printf("unsupported value, recompile CHarm");
    printf("\n");


    printf("SIMD: ");
    ret = CHARM(misc_buildopt_simd)();
    if (ret == BUILDOPT_SIMD_NONE)
        printf("disabled");
    else if (ret == BUILDOPT_SIMD_AVX)
        printf("avx");
    else if (ret == BUILDOPT_SIMD_AVX2)
        printf("avx2");
    else if (ret == BUILDOPT_SIMD_AVX512)
        printf("avx-512");
    else if (ret == BUILDOPT_SIMD_NEON)
        printf("neon");
    else if (ret == BUILDOPT_SIMD_SSE41)
        printf("sse4.1");
    else
        printf("unsupported value, recompile CHarm");
    printf(" (vector size %d)", CHARM(misc_buildopt_simd_vector_size)());
    printf("\n");


    printf("OpenMP in CHarm: ");
    if (CHARM(misc_buildopt_omp_charm)())
        printf("enabled");
    else
        printf("disabled");
    printf("\n");


    printf("OpenMP in FFTW: ");
    if (CHARM(misc_buildopt_omp_fftw)())
        printf("enabled");
    else
        printf("disabled");
    printf("\n");


    printf("MPI: ");
    if (CHARM(misc_buildopt_mpi)())
        printf("enabled");
    else
        printf("disabled");
    printf("\n");


    printf("MPFR: ");
    if (CHARM(misc_buildopt_mpfr)())
        printf("enabled");
    else
        printf("disabled");
    printf("\n");


    printf("FFTW version: %s\n", CHARM(misc_buildopt_version_fftw)());


    int mpi_major_h, mpi_minor_h, mpi_major_l, mpi_minor_l;
    CHARM(misc_buildopt_version_mpi)(&mpi_major_h, &mpi_minor_h,
                                     &mpi_major_l, &mpi_minor_l);
    printf("MPI standard version (header): ");
#if HAVE_MPI
    printf("%d.%d\n", mpi_major_h, mpi_minor_h);
#else
    printf("%s\n", LIB_NA_STR);
#endif
    printf("MPI standard version (library): ");
#if HAVE_MPI
    printf("%d.%d\n", mpi_major_l, mpi_minor_l);
#else
    printf("%s\n", LIB_NA_STR);
#endif


    printf("MPI implementation version (library): ");
#if HAVE_MPI
    char mpi_lib_version[MPI_MAX_LIBRARY_VERSION_STRING];
    int len;
    MPI_Get_library_version(mpi_lib_version, &len);
    printf("%s\n", mpi_lib_version);
#else
    printf("%s\n", LIB_NA_STR);
#endif


    int mpfr_major, mpfr_minor, mpfr_patch;
    const char *mpfr_lib = CHARM(misc_buildopt_version_mpfr)(&mpfr_major,
                                                             &mpfr_minor,
                                                             &mpfr_patch);
    printf("MPFR version (header): ");
#if HAVE_MPFR
    printf("%d.%d.%d\n", mpfr_major, mpfr_minor, mpfr_patch);
#else
    printf("%s\n", LIB_NA_STR);
#endif
    printf("MPFR version (library): %s\n", mpfr_lib);


    int gmp_major, gmp_minor, gmp_patch;
    const char *gmp_lib = CHARM(misc_buildopt_version_gmp)(&gmp_major,
                                                           &gmp_minor,
                                                           &gmp_patch);
    printf("GMP version (header): ");
#if HAVE_MPFR
    printf("%d.%d.%d\n", gmp_major, gmp_minor, gmp_patch);
#else
    printf("%s\n", LIB_NA_STR);
#endif
    printf("GMP version (library): %s\n", gmp_lib);


    printf("isfinite macro in math.h: ");
    if (CHARM(misc_buildopt_isfinite)())
        printf("yes");
    else
        printf("no");
    printf("\n\n");


#ifdef _CHARM_CC
    printf("Compiler (CC): %s\n", _CHARM_CC);
#endif
#ifdef _CHARM_CFLAGS
    printf("Debugging and optimization options (CFLAGS): %s\n",
           _CHARM_CFLAGS);
#endif
#ifdef _CHARM_CPPFLAGS
    printf("Preprocessor options (CPPFLAGS): %s\n", _CHARM_CPPFLAGS);
#endif
#ifdef _CHARM_LDFLAGS
    printf("Options for the linker (LDFLAGS): %s\n", _CHARM_LDFLAGS);
#endif
#ifdef _CHARM_LIBS
    printf("-l options passed to the linker (LIBS): %s\n", _CHARM_LIBS);
#endif


    return;
}
