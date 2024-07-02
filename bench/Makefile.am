# Name of the benchmark program (will not be installed)
noinst_PROGRAMS = bench_run@P@


EXTRA_DIST = plot-bench.py \
			 plot-memory.py


# Source code of the benchmark program
bench_run@P@_SOURCES = bench.c


# Path to header files
bench_run@P@_CPPFLAGS = -I$(top_srcdir)/charm


# CFLAGS
if OPENMP
bench_run@P@_CFLAGS = $(OPENMP_CFLAGS)
endif


# Link CHarm library
bench_run@P@_LDADD = $(top_builddir)/libcharm@P@.la


# Get the absolute path to the benchmark program
bench_run@P@dir = $(abs_top_builddir)/bench


$(top_builddir)/libcharm@P@.la:
	cd $(top_builddir) && $(MAKE)


bench-local: bench_run@P@$(EXEEXT)
	$(bench_run@P@dir)/bench_run@P@


remove-files:
	rm -f bench_run*


clean-local: remove-files


distclean-local: remove-files
