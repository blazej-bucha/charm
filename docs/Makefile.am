EXTRA_DIST = build/html \
			 img \
			 logo \
			 src \
			 Doxyfile \
			 index.html \
			 .nojekyll


# Documentation formats that are deleted from the "./build" directory *only* 
# when calling "make maintainer-clean"
SPHINXTARGETS = html


# Documentation formats that are deleted also with "clean" and "distclean"
SPHINXTARGETS_DEL = dirhtml \
					singlehtml \
					pickle \
					json \
					htmlhelp \
					qthelp \
					devhelp \
					epub \
					latex \
					latexpdf \
					latexpdfja \
					text \
					man \
					texinfo \
					gettext \
					changes \
					xml \
					pseudoxml \
					linkcheck \
					doctest \
					coverage
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
SOURCEDIR     = src
BUILDDIR      = build


# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)


$(top_srcdir)/@OBJDIR@/libcharm@P@@OMP@@SHARED_EXT@:
	cd $(top_srcdir) && $(MAKE) -j


lib-file: $(top_srcdir)/@OBJDIR@/libcharm@P@@OMP@@SHARED_EXT@
	cd $(top_srcdir)/wrap/pyharm && $(MAKE) cp-lib


doxygen:
	cd $(top_srcdir) && ./config.status
	doxygen Doxyfile


$(SPHINXTARGETS) $(SPHINXTARGETS_DEL): doxygen lib-file
if DOC
if DOC_PY
	rm -rf build/$@
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
else
	@echo ERROR: Python interface is required to build the documentation
endif
else
	@echo ERROR: Documentation can only be built if CHarm is compiled in \
		  double precision
endif


remove-doc-files:
	rm -rf xml
	rm -rf build/doctrees
	cd build && rm -rf $(SPHINXTARGETS_DEL)


clean-local: remove-doc-files


distclean-local: remove-doc-files


maintainer-clean-local: remove-doc-files
	rm -rf build

