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
		  src/mpi \
		  src/glob \
		  src/gfm \
		  src/mpfr \
		  src/simd \
		  . \
		  tests \
		  bench \
		  docs \
		  m4 \
		  cookbook \
		  data \
		  wrap/pyharm


SUBLIBS = src/crd/libcharm@P@_crd.la \
		  src/err/libcharm@P@_err.la \
		  src/integ/libcharm@P@_integ.la \
		  src/leg/libcharm@P@_leg.la \
		  src/misc/libcharm@P@_misc.la \
		  src/sha/libcharm@P@_sha.la \
		  src/shc/libcharm@P@_shc.la \
		  src/shs/libcharm@P@_shs.la \
		  src/xnum/libcharm@P@_xnum.la \
		  src/glob/libcharm@P@_glob.la \
		  src/gfm/libcharm@P@_gfm.la \
		  src/simd/libcharm@P@_simd.la
if MPFR
SUBLIBS += src/mpfr/libcharm@P@_mpfr.la
endif
if MPI
SUBLIBS += src/mpi/libcharm@P@_mpi.la
endif


CLEANFILES = charm.pc \
			 charm[fq].pc


EXTRA_DIST = AUTHORS \
			 INSTALL \
			 LICENSE \
			 COPYING \
			 NEWS \
			 README.rst \
			 TODO \
			 VERSION \
			 configure.ac \
			 pyharm-parser.py


lib_LTLIBRARIES = libcharm@P@.la
libcharm@P@_la_SOURCES =
libcharm@P@_la_LIBADD  = $(SUBLIBS)
libcharm@P@_la_CFLAGS  = $(AM_CFLAGS)


charm@P@.pc: charm-tmpl.pc
	cp -f charm-tmpl.pc charm@P@.pc
pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = charm@P@.pc


.PHONY: bench


bench:
	cd bench && $(MAKE) bench-local


pyharm:
if PY
	cd wrap/pyharm && $(MAKE) pyharm-local
else
	@echo !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	@echo
	@echo ERROR: --enable-python option is required to build PyHarm
	@echo
	@echo !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
endif


pyharm-parser: all
	$(PYTHON) pyharm-parser.py


generate-reference-data:
	cd tests && $(MAKE) generate-reference-data-local
