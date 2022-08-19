EXTRA_DIST = build/html \
			 img \
			 logo \
			 src \
			 Doxyfile \
			 Makefile-sphinx \
			 index.html \
			 .nojekyll


# Rule to generate the HTML documentation.  Requires the "doxygen" 
# documentation system and Python modules "sphinx", "sphinx_rtd_theme" and 
# "breathe".
docs-local:
if DOC
	rm -rf build/html
	cd .. && ./config.status charm/charm_*.h && cd docs
	doxygen Doxyfile
	$(MAKE) -f Makefile-sphinx html
else
	@echo ERROR: Documentation can only be built if CHarm is compiled in \
		  double precision
endif


remove-doc-files:
	rm -rf build/doctrees xml


clean-local: remove-doc-files


distclean-local: remove-doc-files


maintainer-clean-local: remove-doc-files
	rm -rf build
