EXTRA_DIST = build/html \
			 img \
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


doxygen:
	cd $(top_srcdir) && ./config.status
if !DOC
	@echo !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	@echo
	@echo ERROR: Double precision, Python interface, MPI and MPFR are \
		  required to build the documentation
	@echo
	@echo !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	@exit 1
endif
	doxygen Doxyfile


$(SPHINXTARGETS) $(SPHINXTARGETS_DEL): doxygen
	cd $(top_srcdir)/wrap/pyharm@P@ && $(MAKE) pyharm-build-local
	cd $(top_srcdir) && rm -f ./wrap/pyharm@P@/libcharm@P@*.so
	cd $(top_srcdir) && \
		cp ./build/lib*/pyharm/libcharm@P@*.so ./wrap/pyharm@P@
	rm -rf build/$@
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)


remove-doc-files:
	rm -rf xml
	rm -rf build/doctrees
	cd build && rm -rf $(SPHINXTARGETS_DEL)


clean-local: remove-doc-files


distclean-local: remove-doc-files


maintainer-clean-local: remove-doc-files
	rm -rf build

