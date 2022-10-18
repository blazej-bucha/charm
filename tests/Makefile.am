# Name of the test program (will not be installed)
check_PROGRAMS = test_run@P@@OMP@


EXTRA_DIST = *.h


# Source code of the test program
test_run@P@@OMP@_SOURCES = run.c \
						   cmp_arrays.c \
						   validate.c \
						   generate_point.c \
						   generate_cell.c \
						   crd.c \
						   shc.c \
						   shs.c \
						   sha.c \
						   leg.c


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
