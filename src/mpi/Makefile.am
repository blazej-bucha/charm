if MPI
noinst_LTLIBRARIES = libcharm@P@_mpi.la


AM_CPPFLAGS = -I$(top_srcdir)/charm


libcharm@P@_mpi_la_SOURCES = mpi_crd_point_init.c \
							 mpi_crd_point_init_base.c \
							 mpi_crd_point_alloc.c \
							 mpi_crd_point_malloc.c \
							 mpi_crd_point_calloc.c \
							 mpi_crd_point_gl.c \
							 mpi_crd_point_dh1.c \
							 mpi_crd_point_dh2.c \
							 mpi_crd_point_quad.c \
							 mpi_crd_point_local2distributed.c \
							 mpi_crd_point_check_struct.c \
							 mpi_crd_point_issymm.c \
							 mpi_shc_init.c \
							 mpi_shc_alloc.c \
							 mpi_shc_malloc.c \
							 mpi_shc_calloc.c \
							 mpi_shc_local_ncs.c \
							 mpi_shc_local2distributed.c \
							 mpi_shc_check_struct.c \
							 mpi_size_t2charm_mpi_count.c \
							 mpi_err_init.c \
							 mpi_err_isempty.c \
							 mpi_err_isdistributed.c \
							 mpi_err_gather.c \
							 mpi_allequal.c \
							 mpi_two_comm_match.c \
							 mpi_three_comm_match.c \
							 mpi_min_int_max.c \
							 mpi_check_point_shc_err.c
endif


EXTRA_DIST = *.h
