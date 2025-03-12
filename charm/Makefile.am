pkginclude_HEADERS = charm@P@.h \
					 charm@P@_crd.h \
					 charm@P@_err.h \
					 charm@P@_integ.h \
					 charm@P@_leg.h \
					 charm@P@_misc.h \
					 charm@P@_sha.h \
					 charm@P@_shc.h \
					 charm@P@_shs.h \
					 charm@P@_glob.h \
					 charm@P@_gfm.h
if MPI
pkginclude_HEADERS += charm@P@_mpi.h
endif


remove-header-files:
	rm -f charm*.h


clean-local: remove-header-files


distclean-local: remove-header-files


dist-hook:
	rm -f $(top_distdir)/charm/charm*.h

