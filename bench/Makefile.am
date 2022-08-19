# Name of the benchmark program (will not be installed)
check_PROGRAMS = bench_run@P@@OMP@


EXTRA_DIST = plot-bench.py \
			 plot-memory.py


# Source code of the benchmark program
bench_run@P@@OMP@_SOURCES = bench.c


# Path to header files
bench_run@P@@OMP@_CPPFLAGS = -I$(top_srcdir)/charm


# CFLAGS
if OPENMP
bench_run@P@@OMP@_CFLAGS = $(OPENMP_CFLAGS)
endif


# Link CHarm library
bench_run@P@@OMP@_LDADD = $(top_builddir)/libcharm@P@@OMP@.la


# Get the absolute path to the benchmark program
bench_run@P@@OMP@dir = $(abs_top_builddir)/bench


bench-local: bench_run@P@@OMP@$(EXEEXT)
	$(bench_run@P@@OMP@dir)/bench_run@P@@OMP@


remove-files:
	rm -f bench_run*


clean-local: remove-files


distclean-local: remove-files
