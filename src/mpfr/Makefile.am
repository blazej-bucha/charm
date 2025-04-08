if MPFR
noinst_LTLIBRARIES = libcharm@P@_mpfr.la


AM_CPPFLAGS = -I$(top_srcdir)/charm


libcharm@P@_mpfr_la_SOURCES = mpfr_binomial.c \
							  mpfr_check_bits.c \
							  mpfr_double_fact.c \
							  mpfr_fact.c \
							  mpfr_legendre.c \
							  mpfr_ndarray_free.c \
							  mpfr_ndarray_malloc.c \
							  mpfr_ndarray_check.c \
							  mpfr_is_nearly_equal.c
endif


EXTRA_DIST = *.h
