# Note that a few files from EXTRA_DIST generated during the build process will 
# intentionally still *not* be included in the distribution (see 
# "del_in_derived_files")
EXTRA_DIST = *.py


pyharm-build-local:
if PY
if SINGLE_PY
	cd $(top_srcdir)/wrap && rm -rf pyharm@P@ && $(LN_S) pyharm pyharm@P@
endif
	cd $(top_srcdir) && CC=$(CC) $(PYTHON) -m pip install . --verbose \
		--target $(pythondir) --upgrade
endif


pyharm-local: remove-py-build
	$(MAKE) pyharm-build-local


install-exec-local: remove-py-build
	rm -rf $(top_srcdir)/build
	$(MAKE) pyharm-local


uninstall-local:
if PY
	rm -rf $(pythondir)/pyharm@P@
	MAJOR=`echo "@CHARM_MAJOR@" | tr - .`; \
	MINOR=`echo "@CHARM_MINOR@" | tr - .`; \
	PATCH=`echo "@CHARM_PATCH@" | tr - .`; \
	rm -rf $(pythondir)/pyharm@P@-$$MAJOR.$$MINOR.$$PATCH*
endif


# Removes all "*.py" files derived from "*.py.in"
del_in_derived_files=`ls *.py.in | sed 's/\.in//g' | xargs rm -f`


# This target should be executed regardless of whether or not PyHarm is build.
remove-py-build:
	rm -f *.so
	rm -rf __pycache__
	rm -rf $(top_srcdir)/build
	rm -rf $(top_srcdir)/wrap/pyharm.egg-info
	rm -rf $(top_srcdir)/wrap/pyharmf
	rm -rf $(top_srcdir)/wrap/pyharmf.egg-info


clean-local: remove-py-build
	${del_in_derived_files}


distclean-local: remove-py-build
	${del_in_derived_files}


# This target should be executed regardless of whether or not PyHarm is build.
dist-hook:
	cd $(top_distdir)/wrap/pyharm && ${del_in_derived_files}
