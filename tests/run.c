/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../src/prec.h"
#include "shs.h"
#include "shc.h"
#include "crd.h"
#include "sha.h"
#include "leg.h"
/* ------------------------------------------------------------------------- */






/* Function prototypes */
/* ------------------------------------------------------------------------- */
static void print_test_outcome(int, char *);
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


    /* Print version number, etc. of the compiled CHarm library */
    printf("Printing info on CHarm...\n");
    printf("..................................\n");
    CHARM(misc_print_version());
    printf("..................................\n\n");


    /* ..................................................................... */
    printf("Testing the \"shc\" module...\n");


    int err = 0;
    int err_sum = 0;


    err = shc(nmax_topo,
              nmax_pot,
              SHCs_in_topo_mtx_file,
              SHCs_in_pot_gfc_file,
              SHCs_in_pot_tbl_file,
              SHCs_out_topo_bin_file,
              SHCs_out_topo_mtx_file,
              SHCs_out_pot_tbl_n_file,
              SHCs_out_pot_tbl_m_file);
    print_test_outcome(err, "shc");
    err_sum += err;
    /* ..................................................................... */


    /* ..................................................................... */
    printf("Testing the \"crd\" module...\n");
    err = crd();
    print_test_outcome(err, "crd");
    err_sum += err;
    /* ..................................................................... */


    /* ..................................................................... */
    printf("Testing the \"shs\" module...\n");
    err = shs(nmax_topo, SHCs_in_topo_mtx_file, nmax_pot,
              SHCs_in_pot_mtx_file);
    print_test_outcome(err, "shs");
    err_sum += err;
    /* ..................................................................... */


    /* ..................................................................... */
    printf("Testing the \"sha\" module...\n");
    err = sha(nmax_pot, SHCs_in_pot_mtx_file);
    print_test_outcome(err, "sha");
    err_sum += err;
    /* ..................................................................... */


    /* ..................................................................... */
    printf("Testing the \"leg\" module...\n");
    err = leg();
    print_test_outcome(err, "leg");
    err_sum += err;
    /* ..................................................................... */


    printf("\n");
    if (err_sum)
        printf("Oh, no...  %d %s didn't pass.\n", 
               err_sum, err_sum == 1 ? "test" : "tests");
    else
        printf("All tests seem to pass.\n");
    printf("\n");
    printf("End of the validation.\n");
    printf("\n");
    printf("==================================\n\n");


    return CHARM_SUCCESS;
    /* --------------------------------------------------------------------- */


}





static void print_test_outcome(int err, char module[])
{
    if (err)
        printf("    The module \"%s\" didn't pass some tests!\n", module);
    else
        printf("    The module \"%s\" passed all tests.\n", module);
}
