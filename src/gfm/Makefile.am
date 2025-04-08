noinst_LTLIBRARIES = libcharm@P@_gfm.la


AM_CPPFLAGS = -I$(top_srcdir)/charm


libcharm@P@_gfm_la_SOURCES = gfm_global_density_3d.c \
							 gfm_global_density_lateral.c \
							 gfm_global_density_const.c \
							 gfm_check_p.c \
							 gfm_check_kminkmax.c
if MPFR
libcharm@P@_gfm_la_SOURCES += gfm_cap_density_3d.c \
							  gfm_cap_density_3d_pywrap.c \
							  gfm_cap_density_lateral.c \
							  gfm_cap_density_lateral_pywrap.c \
							  gfm_cap_density_const.c \
							  gfm_cap_density_const_pywrap.c \
							  gfm_cap_q_aps.c \
							  gfm_cap_q_check_psi.c \
							  gfm_cap_q_check_radius.c \
							  gfm_cap_q_check_type.c \
							  gfm_cap_q_norm_dist.c \
							  gfm_cap_q_dnorm_dist.c \
							  gfm_cap_q_ddist.c \
							  gfm_cap_q_dkernel.c \
							  gfm_cap_q_rpows.c \
							  gfm_cap_q_dr.c \
							  gfm_cap_q_ref.c \
							  gfm_cap_qu0.c \
							  gfm_cap_quu.c \
							  gfm_cap_q_check_prec.c \
							  gfm_cap_q.c \
							  gfm_cap_nq.c
if !CHARM_QUAD
libcharm@P@_gfm_la_SOURCES += gfm_cap_q_check_prec_pywrap.c \
							  gfm_cap_q_pywrap.c
endif

endif


EXTRA_DIST = *.h \
			 gfm_sgfm.c \
			 gfm_sgfm_cap_lateral.c
