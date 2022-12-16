# Note that a few files from EXTRA_DIST generated during the compilation and 
# installation will still *not* be intentionally included in distributions (see 
# the "DEL_FILES" variable)
EXTRA_DIST = *.py \
			 lib


# Files to be delete when cleaning or when packing a distribution
DEL_FILES = __init__.py \
			_data_types.py \
			lib/*


cp-lib:
if PY
	cp -f $(top_srcdir)/@OBJDIR@/libcharm@P@@OMP@@SHARED_EXT@ lib
endif


install-exec-local: cp-lib
if PY
if SINGLE_PY
	cd $(top_srcdir)/wrap && rm -rf pyharm@P@ && $(LN_S) pyharm pyharm@P@
endif
	cd $(top_srcdir) && $(PYTHON) -m pip install . --target $(pythondir)
endif


uninstall-local:
if PY
	rm -rf $(pythondir)/pyharm@P@
	MAJOR=`echo "@CHARM_MAJOR@" | tr - .`; \
	MINOR=`echo "@CHARM_MINOR@" | tr - .`; \
	PATCH=`echo "@CHARM_PATCH@" | tr - .`; \
	rm -rf $(pythondir)/pyharm@P@-$$MAJOR.$$MINOR.$$PATCH*
endif


# This target should be executed regardless of whether or not PyHarm is build.
remove-files:
	rm -f $(DEL_FILES)
	rm -rf __pycache__
	rm -rf $(top_srcdir)/build
	rm -rf $(top_srcdir)/wrap/pyharm.egg-info
	rm -rf $(top_srcdir)/wrap/pyharmf
	rm -rf $(top_srcdir)/wrap/pyharmf.egg-info


clean-local: remove-files


distclean-local: remove-files


# This target should be executed regardless of whether or not PyHarm is build.
dist-hook:
	cd $(top_distdir)/wrap/pyharm && rm -f $(DEL_FILES)

