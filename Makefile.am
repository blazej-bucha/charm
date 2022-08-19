ACLOCAL_AMFLAGS = -I m4


SUBDIRS = charm \
		  src \
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
		  src/simd \
		  tests \
		  bench \
		  docs \
		  m4 \
		  cookbook \
		  data


SUBLIBS = src/crd/libcharm@P@_crd@OMP@.la \
		  src/err/libcharm@P@_err@OMP@.la \
		  src/integ/libcharm@P@_integ@OMP@.la \
		  src/leg/libcharm@P@_leg@OMP@.la \
		  src/misc/libcharm@P@_misc@OMP@.la \
		  src/sha/libcharm@P@_sha@OMP@.la \
		  src/shc/libcharm@P@_shc@OMP@.la \
		  src/shs/libcharm@P@_shs@OMP@.la \
		  src/xnum/libcharm@P@_xnum@OMP@.la \
		  src/glob/libcharm@P@_glob@OMP@.la \
		  src/simd/libcharm@P@_simd@OMP@.la


CLEANFILES = charm.pc \
			 charm[fq].pc \
			 charm*_omp.pc


EXTRA_DIST = AUTHORS \
			 INSTALL \
			 LICENSE \
			 NEWS \
			 README.rst \
			 TODO \
			 VERSION \
			 charm-tmpl.pc.in \
			 config.h.in \
			 configure.ac


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
