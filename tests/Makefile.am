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
						   crd.c \
						   shc.c \
						   shs.c \
						   sha.c \
						   leg.c \
						   integ.c \
						   misc.c \
						   check_integ_pn1m1pn2m2.c \
						   check_integ_yi1n1m1yi2n2m2.c \
						   abs_r.c


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
