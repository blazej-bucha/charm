# Name of the test program (will not be installed)
check_PROGRAMS = test_run@P@@OMP@


EXTRA_DIST = *.h


# Source code of the test program
test_run@P@@OMP@_SOURCES = run.c \
						   cmp_vals.c \
						   cmp_arrays.c \
						   validate.c \
						   generate_point.c \
						   generate_cell.c \
						   module_crd.c \
						   module_shc.c \
						   module_shs.c \
						   module_sha.c \
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
						   check_shs_cell.c \
						   check_shs_cell_isurf.c \
						   check_sha_point.c \
						   check_sha_cell.c \
						   check_shc_read_type.c \
						   check_shc_write_type.c \
						   check_shc_dav.c \
						   check_shc_ddav.c \
						   check_shc_rescale.c \
						   check_shc_alloc.c \
						   check_shc_init.c \
						   check_simd_abs_r.c


# Path to header files
test_run@P@@OMP@_CPPFLAGS = -I$(top_srcdir)/charm


# CFLAGS
if OPENMP
test_run@P@@OMP@_CFLAGS = $(OPENMP_CFLAGS)
endif


# Link CHarm library
test_run@P@@OMP@_LDADD = $(top_builddir)/libcharm@P@@OMP@.la


# Get the absolute path to the test program
test_run@P@@OMP@dir = $(abs_top_builddir)/tests


check-local: test_run@P@@OMP@$(EXEEXT)
	$(test_run@P@@OMP@dir)/test_run@P@@OMP@


remove-files:
	rm -f test_run*


clean-local: remove-files


distclean-local: remove-files
