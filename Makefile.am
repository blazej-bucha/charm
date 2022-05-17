SUBDIRS = charm \
		  src/crd \
		  src/err \
		  src/integ \
		  src/leg \
		  src/misc \
		  src/sha \
		  src/shc \
		  src/shs \
		  src/xnum \
		  src/glob \
		  tests \
		  bench \
		  docs


ACLOCAL_AMFLAGS = -I m4


CLEANFILES = charm.pc \
			 charm[fq].pc \
			 charm*_omp.pc \
			 config.h \
			 libcharm*.la \
			 data/output/* \
			 charm/charm*.h


EXTRA_DIST = data \
			 docs/build/html \
			 docs/img \
			 docs/logo \
			 docs/src/*.rst \
			 docs/src/*.py \
			 docs/Doxyfile \
			 docs/Makefile-sphinx \
			 docs/index.html \
			 docs/.nojekyll \
			 src/prec.h \
			 m4 \
			 cookbook \
			 bench/plot-bench.py \
			 bench/plot-memory.py \
			 AUTHORS \
			 INSTALL \
			 LICENSE \
			 NEWS \
			 README.rst \
			 TODO \
			 VERSION \
			 charm-tmpl.pc.in \
			 config.h.in \
			 configure.ac


dist-hook:
	rm -f $(top_distdir)/charm/charm*.h
	rm -f $(top_distdir)/data/output/*


SUBLIBS = src/crd/libcharm@P@_crd@OMP@.la \
		  src/err/libcharm@P@_err@OMP@.la \
		  src/integ/libcharm@P@_integ@OMP@.la \
		  src/leg/libcharm@P@_leg@OMP@.la \
		  src/misc/libcharm@P@_misc@OMP@.la \
		  src/sha/libcharm@P@_sha@OMP@.la \
		  src/shc/libcharm@P@_shc@OMP@.la \
		  src/shs/libcharm@P@_shs@OMP@.la \
		  src/xnum/libcharm@P@_xnum@OMP@.la \
		  src/glob/libcharm@P@_glob@OMP@.la


lib_LTLIBRARIES = libcharm@P@@OMP@.la
libcharm@P@@OMP@_la_SOURCES =
libcharm@P@@OMP@_la_LIBADD  = $(SUBLIBS)
libcharm@P@@OMP@_la_CFLAGS  = $(AM_CFLAGS)


charm@P@@OMP@.pc: charm-tmpl.pc
	cp -f charm-tmpl.pc charm@P@@OMP@.pc
pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = charm@P@@OMP@.pc


.PHONY: docs bench


docs:
	cd docs && $(MAKE) docs-local


bench:
	cd bench && $(MAKE) bench-local
