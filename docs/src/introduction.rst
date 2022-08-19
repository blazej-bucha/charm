============
Introduction
============

CHarm is a C library to work with spherical harmonics up to almost arbitrarily
high degrees.


Features
========

* Supports real-valued fully-normalized surface and solid spherical harmonics
  (the geodetic norm).

* Performs FFT-based surface spherical harmonic analysis and solid spherical
  harmonic synthesis with minimized memory requirements.

* Stable up to high degrees and orders (tens of thousands and beyond).

* Available in single, double and quadruple precision.

* Supports point and mean data values (both analysis and synthesis).

* Supports synthesis at grids and at scattered points/cells.  Synthesis at
  grids is done by efficient FFT-based algorithms if possible.  For grids,
  where FFT cannot be used, applied are the Chebyshev recurrences along the
  latitude parallels.

* Supports the Gauss--Legendre and Driscoll--Healy quadratures.

* Integrates solid spherical harmonic expansions (e.g., of the gravitational
  potential) on band-limited irregular surfaces (e.g., on the Earth's
  surface). [#f1]_

* Computes Fourier coefficients of fully-normalized associated Legendre
  functions of the first kind up to ultra-high harmonic degrees.

* Supports `OpenMP <https://www.openmp.org/>`_ parallelization for
  shared-memory architectures.

* Supports AVX, AVX2 and AVX-512 SIMD CPU intructions to improve the 
  performance.

* Performs discrete FFT by `FFTW <http://www.fftw.org/>`_.

.. [#f1] This routine is unique to CHarm.


.. _download:

Source code
===========

GitHub: `https://github.com/blazej-bucha/charm
<https://github.com/blazej-bucha/charm>`_

* Releases are pushed to ``master`` and the development happens in
  ``develop``.  If you prefer the most up-to-date version, switch to the
  ``develop`` branch.


Tarball and zip files of releases:
`https://github.com/blazej-bucha/charm/tags
<https://github.com/blazej-bucha/charm/tags>`_


Documentation
=============

The documentation of the latest version from the ``develop`` branch is 
available at `https://blazej-bucha.github.io/charm/index.html  
<https://blazej-bucha.github.io/charm/index.html>`_.

A pre-compiled HTML documentation is also available in ``docs/build/html`` or 
it can be built by executing ``make docs`` (requires ``doxygen`` and Python 
modules ``sphinx``, ``sphinx_rtd_theme`` and ``breathe``).  PDF and other 
versions of the documentation can also be built (outside the installation 
system) as ``cd docs && make -f Makefile-sphinx latexpdf``, etc.


.. _contact:

Contact
=======

Should you have any comments, questions, bug report or criticism, please feel
free to contact the author, Blažej Bucha, at blazej.bucha@stuba.sk.  Further
products developed by the author can be found at `https://blazejbucha.com
<https://blazejbucha.com>`_.


Other spherical-harmonic-based libraries
========================================

Many other libraries for working with spherical harmonics are available, each
having its pros and cons.  Explore!  A few examples are:

* `SHTOOLS <https://github.com/SHTOOLS>`_: Fortran95 library with Python API,

* `SHTns <https://bitbucket.org/nschaeff/shtns>`_: a C library for spherical
  harmonic transforms,

* `ISPACK <https://www.gfd-dennou.org/arch/ispack/>`_: a Fortran library for
  spherical harmonic transforms,

* `Libsharp <https://github.com/Libsharp/libsharp>`_: a C99 library for
  spherical harmonic transforms,

* `healpy <https://healpy.readthedocs.io/en/latest/index.html>`_: a Python
  package to handle pixelated data on the sphere building on the `HEALPix
  <https://healpix.jpl.nasa.gov/>`_ C++ library,

* `HARMONIC_SYNTH
  <https://earth-info.nga.mil/index.php?dir=wgs84&action=wgs84>`_: a Fortran
  code for spherical harmonic synthesis written by the EGM2008 development
  team.

* `SPHEREPACK
  <https://github.com/NCAR/NCAR-Classic-Libraries-for-Geophysics>`_: a Fortran
  library of spherical harmonic transforms,

* `SHAVEL <https://doi.org/10.1016/j.cpc.2018.06.015>`_: a program for the
  spherical harmonic analysis of a horizontal vector field sampled in an
  equiangular grid on a sphere

* `ICGEM <http://icgem.gfz-potsdam.de/home>`_: Online calculation service for
  working with Earth and celestial gravitational models,

* `FaVeST <https://github.com/mingli-ai/FaVeST>`_: Fast Vector Spherical
  Harmonic Transforms in MATLAB.

* `SHBundle
  <https://www.gis.uni-stuttgart.de/en/research/downloads/shbundle/>`_:
  Spherical harmonic analysis and synthesis in MATLAB up to high degrees and
  orders,

* `Spherical Harmonics Manipulator
  <https://sourceforge.net/projects/hmanipulator/>`_: Spherical harmonic
  synthesis in sparse points and grids (no longer maintained),

* `GrafLab <https://blazejbucha.com/#GrafLab>`_ and `isGrafLab
  <https://blazejbucha.com/#isGrafLab>`_: MATLAB-based software packages for
  spherical harmonic synthesis of gravity field functionals up to high degrees
  and orders (tens of thousands and well beyond).
