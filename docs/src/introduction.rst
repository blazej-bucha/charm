============
Introduction
============

CHarm is a C library to work with spherical harmonics up to almost arbitrarily 
high degrees.  The library is accompanied by a Python wrapper called PyHarm.


Features
========

* Supports real-valued fully-normalized surface and solid spherical harmonics
  (the geodetic norm).

* Performs FFT-based surface spherical harmonic analysis and solid spherical
  harmonic synthesis with minimized memory requirements.

* Stable up to high degrees and orders (tens of thousands and beyond).

* Available in single, double and quadruple precision.

* Supports point and mean data values (both analysis and synthesis).

* Supports synthesis at grids and at scattered points/cells.  Grid-wise 
  computations are done by FFT whenever possible.  If FFT cannot be applied, 
  the less efficient Chebyshev recurrences are used along the latitude 
  parallels instead.

* Computes the full first- and second-order gradients at evaluation points 
  (e.g., the gravitational vector and the gravitational tensor).

* Supports the Gauss--Legendre and Driscoll--Healy quadratures.

* Integrates solid spherical harmonic expansions (e.g., of the gravitational
  potential) on band-limited irregular surfaces (e.g., on the Earth's
  surface). [#f1]_

* Computes Fourier coefficients of fully-normalized associated Legendre
  functions of the first kind up to ultra-high harmonic degrees.

* Supports `OpenMP <https://www.openmp.org/>`_ parallelization for
  shared-memory architectures.

* Supports AVX, AVX2 and AVX-512 SIMD CPU instructions to improve the 
  performance.

* Performs discrete FFT by `FFTW <http://www.fftw.org/>`_.

* Ships with a Python wrapper to enable high-level programming while retaining 
  the efficiency of the C language.  The wrapper, called PyHarm, wraps CHarm 
  using `ctypes <https://docs.python.org/3/library/ctypes.html>`_ and is fully 
  integrated with `numpy <https://numpy.org/>`_.

.. [#f1] This routine is unique to CHarm.


Installation
============

On Linux (x86_64), macOS (x86_64, ARM64) and Windows (x86_64), install the 
Python wrapper PyHarm using ``pip``:

.. code-block:: bash

   pip install pyharm

The C library CHarm has to be build from source.

Further installation details at 
`https://www.charmlib.org/build/html/install.html 
<https://www.charmlib.org/build/html/install.html>`_.


.. _download:

Source code
===========

GitHub: `https://github.com/blazej-bucha/charm 
<https://github.com/blazej-bucha/charm>`_

* Releases are pushed to ``master`` and the development happens in
  ``develop``.

* Tarball and zip files of releases: 
  `https://github.com/blazej-bucha/charm/releases 
  <https://github.com/blazej-bucha/charm/releases>`_


Documentation
=============

The documentation of the latest version from the ``master`` branch is available 
at `https://www.charmlib.org <https://www.charmlib.org>`_.

A pre-compiled HTML documentation is also available in ``docs/build/html``.  
Alternatively, it can be built by executing ``make html`` after the 
``configure`` call (requires ``doxygen`` and Python modules ``sphinx``, 
``sphinx_book_theme`` and ``breathe``).  Other formats of the documentation, 
for instance, a PDF file, can be built with ``cd docs && make latexpdf``, etc.  
To list all available formats, execute ``cd docs && make help``.


.. _contact:


Contact
=======

Should you have any comments, questions, bug report or criticism, please feel
free to contact the author, Blažej Bucha, at blazej.bucha@stuba.sk.  Further
products developed by the author can be found at `https://www.blazejbucha.com
<https://www.blazejbucha.com>`_.


Pronunciation
=============

We prefer to pronounce CHarm and PyHarm like the words ``see harm`` and ``pie 
harm``.  But it is indeed quite charming to pronounce CHarm like the word 
``charm``, especially when the library works like a charm.


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

* `GrafLab <https://github.com/blazej-bucha/graflab>`_ and `isGrafLab 
  <https://github.com/blazej-bucha/isgraflab>`_: MATLAB-based software packages 
  for spherical harmonic synthesis of gravity field functionals up to high 
  degrees and orders (tens of thousands and well beyond).
