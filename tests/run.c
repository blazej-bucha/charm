/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _MSC_VER
#   define _USE_MATH_DEFINES
#endif
#include <math.h>
#include "../src/prec.h"
/* ------------------------------------------------------------------------- */






/* Function prototypes */
/* ------------------------------------------------------------------------- */
int shc(unsigned long, unsigned long, char *, char *, char *, char *,
                                      char *, char *, char *);
int crd(unsigned long);
int shs(unsigned long, char *, unsigned long, char *);
int sha(unsigned long, char *);
int leg(unsigned long);
/* ------------------------------------------------------------------------- */






/*
 * This program tests the "CHarm" library.
 * */
int main(void)
{
    /* Inputs */
    /* --------------------------------------------------------------------- */
    /* Maximum harmonic degree of the topography */
    unsigned long nmax_topo = 4;


    /* Input file with spherical harmonic coefficients of the topography */
    char SHCs_in_topo_mtx_file[] =
        "../data/input/Earth2014-RET2014-degree4-mtx.txt";


    /* Maximum harmonic degree of the potential */
    unsigned long nmax_pot = 10;


    /* Input file with spherical harmonic coefficients of the potential in the
     * matrix text format */
    char SHCs_in_pot_mtx_file[] = "../data/input/EGM96-degree10-mtx.txt";


    /* Input file with spherical harmonic coefficients of the potential in the
     * gfc format */
    char SHCs_in_pot_gfc_file[] = "../data/input/EGM96-degree10.gfc";


    /* Input file with spherical harmonic coefficients of the potential in the
     * table txt format */
    char SHCs_in_pot_tbl_file[] = "../data/input/EGM96-degree10-tbl.txt";


    /* Output file with spherical harmonic coefficients.  This is to test
     * writing of the computed coefficients to a binary file.  The path in
     * "SHCs_out_topo_bin_file" (if any) must already exist. */
    char SHCs_out_topo_bin_file[] =
        "../data/output/Earth2014-RET2014-degree4-out.shcs";


    /* Output file with spherical harmonic coefficients.  This is to test
     * writing coefficients to a matrix text file.  The path in
     * "SHCs_out_topo_mtx_file" (if any) must already exist. */
    char SHCs_out_topo_mtx_file[] =
        "../data/output/Earth2014-RET2014-degree4-out-mtx.txt";


    /* Output file with spherical harmonic coefficients.  This is to test
     * writing coefficients to a table file following the "N" ordering scheme.
     * The path in "SHCs_out_tbl_n_file" (if any) must already exist. */
    char SHCs_out_pot_tbl_n_file[] =
        "../data/output/EGM96-degree10-out-tbl-n.txt";


    /* Output file with spherical harmonic coefficients.  This is to test
     * writing coefficients to a table file following the "M" ordering scheme.
     * The path in "SHCs_out_tbl_m_file" (if any) must already exist. */
    char SHCs_out_pot_tbl_m_file[] =
        "../data/output/EGM96-degree10-out-tbl-m.txt";
    /* --------------------------------------------------------------------- */






    /* Tests */
    /* --------------------------------------------------------------------- */
    printf("\n==================================\n");
    printf("\nStart of the validation.\n\n");
    int err = 0;
    int err_sum = 0;


    /* ..................................................................... */
    printf("Testing the \"shc\" module...\n");
    err = shc(nmax_topo,
              nmax_pot,
              SHCs_in_topo_mtx_file,
              SHCs_in_pot_gfc_file,
              SHCs_in_pot_tbl_file,
              SHCs_out_topo_bin_file,
              SHCs_out_topo_mtx_file,
              SHCs_out_pot_tbl_n_file,
              SHCs_out_pot_tbl_m_file);
    if (err != 0)
        printf("    The module didn't passed!\n");
    else
        printf("    The module passed.\n");
    err_sum += err;
    /* ..................................................................... */


    /* ..................................................................... */
    printf("Testing the \"crd\" module...\n");
    err = crd(nmax_topo);
    if (err != 0)
        printf("    The module didn't passed!\n");
    else
        printf("    The module passed.\n");
    err_sum += err;
    /* ..................................................................... */


    /* ..................................................................... */
    printf("Testing the \"shs\" module...\n");
    err = shs(nmax_topo, SHCs_in_topo_mtx_file, nmax_pot,
              SHCs_in_pot_mtx_file);
    if (err != 0)
        printf("    The module didn't passed!\n");
    else
        printf("    The module passed.\n");
    err_sum += err;
    /* ..................................................................... */


    /* ..................................................................... */
    printf("Testing the \"sha\" module...\n");
    err = sha(nmax_pot, SHCs_in_pot_mtx_file);
    if (err != 0)
        printf("    The module didn't passed!\n");
    else
        printf("    The module passed.\n");
    err_sum += err;
    /* ..................................................................... */


    /* ..................................................................... */
    printf("Testing the \"leg\" module...\n");
    err = leg(10);
    if (err != 0)
        printf("    The module didn't passed!\n");
    else
        printf("    The module passed.\n");
    err_sum += err;
    /* ..................................................................... */


    printf("\n");
    if (err_sum)
        printf("Oh, no...  %d test(s) didn't pass.\n", err_sum);
    else
        printf("All tests seem to pass.\n");
    printf("\n");
    printf("End of the validation.\n");
    printf("\n");
    printf("==================================\n\n");


    return CHARM_SUCCESS;
    /* --------------------------------------------------------------------- */


}
