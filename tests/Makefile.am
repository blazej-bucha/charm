# Name of the test program (will not be installed)
check_PROGRAMS = test_run@P@
if MPI
check_PROGRAMS += mpi_test_run@P@
endif


# Name of the program to generate reference data for the test program (will not 
# be installed)
noinst_PROGRAMS = genref_run@P@


EXTRA_DIST = *.h \
			 check_shs_point_all.c


# Source files common to both the test and the genref program
SOURCES = test_run.c \
		  test_suite_start.c \
		  test_suite_end.c \
		  cmp_vals.c \
		  cmp_arrays.c \
		  write_array.c \
		  write_val.c \
		  validate.c \
		  array2file.c \
		  point_touch_array_elements.c \
		  cell_touch_array_elements.c \
		  shc_touch_array_elements.c \
		  modify_low_degree_coefficients.c \
		  generate_point.c \
		  generate_cell.c \
		  module_crd.c \
		  module_shc.c \
		  module_shs.c \
		  module_sha.c \
		  module_gfm.c \
		  module_leg.c \
		  module_integ.c \
		  misc.c \
		  check_struct.c \
		  check_func.c \
		  check_outcome.c \
		  check_integ_pn1m1pn2m2.c \
		  check_integ_yi1n1m1yi2n2m2.c \
		  check_leg_pnmj_coeffs.c \
		  check_leg_pnmj_alloc.c \
		  check_crd_point_alloc.c \
		  check_crd_point_init.c \
		  check_crd_point_quad.c \
		  check_crd_cell_alloc.c \
		  check_crd_cell_init.c \
		  check_shs_point.c \
		  check_shs_point_grad1.c \
		  check_shs_point_grad2.c \
		  check_shs_point_guru.c \
		  check_shs_cell.c \
		  check_shs_cell_isurf.c \
		  check_sha_point.c \
		  check_sha_cell.c \
		  check_gfm_global_density_3d.c \
		  check_gfm_global_density_lateral.c \
		  check_gfm_global_density_const.c \
		  check_shc_read_type.c \
		  check_shc_read_gfc.c \
		  check_shc_write_type.c \
		  check_shc_dav.c \
		  check_shc_ddav.c \
		  check_shc_rescale.c \
		  check_shc_alloc.c \
		  check_shc_init.c \
		  check_simd_abs_r.c \
		  check_simd_neg_r.c \
		  check_simd_sum_r.c \
		  check_simd_masks.c \
		  check_simd_blend_r.c
if MPFR
SOURCES += check_gfm_cap_density_3d.c \
           check_gfm_cap_density_lateral.c \
           check_gfm_cap_density_const.c \
		   check_gfm_cap_q.c \
		   check_gfm_cap_q_check_prec.c \
		   mpfr_write_val.c \
		   mpfr_write_array.c \
		   mpfr_array2file.c \
		   mpfr_cmp_vals.c \
		   mpfr_cmp_arrays.c \
		   mpfr_validate.c
endif


if MPI
SOURCES_MPI = mpi_test_run.c \
			  test_suite_start.c \
			  test_suite_end.c \
			  check_func.c \
			  check_outcome.c \
			  partition_interval.c \
			  cmp_arrays.c \
			  cmp_vals.c \
			  mpi_shc_distribute.c \
			  mpi_crd_point_distribute.c \
			  module_mpi.c \
			  check_mpi_shc_alloc.c \
			  check_mpi_shc_init.c \
			  check_mpi_shc_local_ncs.c \
			  check_mpi_shc.c \
			  check_mpi_crd_point_init.c \
			  check_mpi_crd_point_quad.c \
			  check_mpi_crd_point_alloc.c \
			  check_mpi_crd_point.c \
			  check_mpi_shs_point.c \
			  check_mpi_sha_point.c \
			  point_touch_array_elements.c \
			  shc_touch_array_elements.c \
			  generate_point.c \
			  check_struct.c
endif


# Source code of the test program
test_run@P@_SOURCES = $(SOURCES)
genref_run@P@_SOURCES = $(SOURCES)
if MPI
mpi_test_run@P@_SOURCES = $(SOURCES_MPI)
endif


# Path to header files
test_run@P@_CPPFLAGS = -I$(top_srcdir)/charm
genref_run@P@_CPPFLAGS = -I$(top_srcdir)/charm
if MPI
mpi_test_run@P@_CPPFLAGS = -I$(top_srcdir)/charm
endif


# CFLAGS
genref_run@P@_CFLAGS = -DGENREF
if OPENMP
test_run@P@_CFLAGS = $(OPENMP_CFLAGS)
genref_run@P@_CFLAGS += $(OPENMP_CFLAGS)
if MPI
mpi_test_run@P@_CFLAGS = $(OPENMP_CFLAGS)
endif
endif


# Link CHarm library
test_run@P@_LDADD = $(top_builddir)/libcharm@P@.la
genref_run@P@_LDADD = $(top_builddir)/libcharm@P@.la
if MPI
mpi_test_run@P@_LDADD = $(top_builddir)/libcharm@P@.la
endif


# Get the absolute path to the test programs
test_run@P@dir = $(abs_top_builddir)/tests
genref_run@P@dir = $(abs_top_builddir)/tests
if MPI
mpi_test_run@P@dir = $(abs_top_builddir)/tests
endif


check-local: test_run@P@$(EXEEXT)
	$(test_run@P@dir)/test_run@P@
if MPI
	$(MPIEXEC) -n 1 $(mpi_test_run@P@dir)/mpi_test_run@P@
	$(MPIEXEC) -n 2 $(mpi_test_run@P@dir)/mpi_test_run@P@
	$(MPIEXEC) -n 3 $(mpi_test_run@P@dir)/mpi_test_run@P@
	$(MPIEXEC) -n 4 $(mpi_test_run@P@dir)/mpi_test_run@P@
endif


generate-reference-data-local: genref_run@P@$(EXEEXT)
	mkdir -p ./genref-output
	mkdir -p ./genref-output/single
	mkdir -p ./genref-output/double
	mkdir -p ./genref-output/quad
	$(genref_run@P@dir)/genref_run@P@


remove-files:
	rm -rf test_run test_runf test_runq test_run-*o genref_run* genref-output


clean-local: remove-files


distclean-local: remove-files
