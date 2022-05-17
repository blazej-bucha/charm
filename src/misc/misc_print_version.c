/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
/* ------------------------------------------------------------------------- */






void CHARM(misc_print_version)(void)
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
    "License: The 3-Clause BSD License\n"


    "\n"
    "Precision: "
#if CHARM_FLOAT
        "single"
#elif CHARM_QUAD
        "quadruple"
#else
        "double"
#endif


    "\n"
    "OpenMP in CHarm: "
#if CHARM_PARALLEL
        "enabled"
#else
        "disabled"
#endif


    "\n"
    "OpenMP in FFTW: "
#if FFTW3_OMP
        "enabled"
#else
        "disabled"
#endif


    "\n"
    "isfinite macro in math.h: "
#if HAVE_ISFINITE
        "yes"
#else
        "no"
#endif


    "\n\n"


    );


    return;
}
