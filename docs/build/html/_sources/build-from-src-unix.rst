.. _build_from_src_unix:

Building from source on Unix-based operating systems
====================================================

This section provides a guide through the build of CHarm and/or PyHarm from 
source on Unix-based operating systems, particularly :ref:`Linux 
<installation_linux>` and :ref:`macOS <installation_mac>`.  First discussed are 
the steps to compile :ref:`CHarm (the C library) <charm_installation>` and then 
follow instructions on how to (optionally) build :ref:`PyHarm (the Python 
wrapper) <pyharm_installation>`.


.. _charm_installation:

CHarm (the C library)
---------------------


.. _charm_requirements:

Dependencies
~~~~~~~~~~~~

Mandatory
"""""""""

* C compiler supporting C99.  `GCC <https://gcc.gnu.org/>`_ is recommended 
  (available by default on most Linux distributions).  Other compilers known to 
  successfully compile CHarm are listed in the :ref:`tested_platforms` chapter.

  *To compile in quadruple precision, GCC is mandatory (v4.6 or later)* and no
  other compiler is allowed.

* `Make <https://www.gnu.org/software/make/>`_ to build the library (available 
  by default on most Linux distributions).

* `pkg-config <https://www.freedesktop.org/wiki/Software/pkg-config/>`_ for
  managing library compile and link flags (available by default on most Linux
  distributions).

* `FFTW <http://www.fftw.org/>`_ for discrete fast Fourier transforms.

Optional
""""""""

* `MPI <https://www.mpi-forum.org/>`_ for parallelization on distributed-memory 
  systems (e.g., high-performance computing clusters).

* `MPFR <https://mpfr.org>`_ for arbitrary precision arithmetic on floating 
  point numbers with correct rounding.  CHarm can be compiled without the MPFR 
  support, though at the cost of loosing the ability to spatially limit the 
  integration radius in spectral gravity forward modelling.


.. _installation_linux:

Installation on Linux
~~~~~~~~~~~~~~~~~~~~~

At first, we will install FFTW (mandatory), MPFR (optional) and MPI (optional), 
after which we will proceed with the installation of CHarm.

.. _installation_FFTW_linux:

FFTW installation
"""""""""""""""""

FFTW can either be installed via your package manager or it can be built from 
the source.

* *Using your package manager*

  * Debian-based distributions:

    .. code-block:: console

         sudo apt install libfftw3-dev

  * Arch Linux: You may want to install, for instance, `fftw-mpi
    <https://aur.archlinux.org/packages/fftw-mpi/>`_ using your package
    manager.

  After the installation, make sure the following libraries were successfully
  installed, depending on whether you want to use CHarm in
  single/double/quadruple precision and with OpenMP enabled/disabled:

  * ``libfftw3f.so`` for single precision version of CHarm with OpenMP
    disabled,
  * ``libfftw3f_omp.so`` for single precision version of CHarm with OpenMP
    enabled,
  * ``libfftw3.so`` for double precision version of CHarm with OpenMP disabled,
  * ``libfftw3_omp.so`` for double precision version of CHarm with OpenMP
    enabled,
  * ``libfftw3q.so`` for quadruple precision version of CHarm with OpenMP
    disabled,
  * ``libfftw3q_omp.so`` for quadruple precision version of CHarm with OpenMP
    enabled.

* *Compilation from source*

  Please read the manual on the homepage of FFTW on how to install the library.
  Example FFTW installation could be:

   .. code-block:: console

      ./configure --enable-openmp --enable-shared
      make
      make check
      sudo make install

  Do not forget to add the ``--enable-shared`` flag to compile FFTW also as
  a shared library (mandatory to install CHarm).  To install CHarm in single or
  quadruple precision, add also the ``--enable-single`` or
  ``--enable-quad-precision`` flag, respectively, when calling the
  ``configure`` script.  If you do not want to parallelize FFTW, you may omit
  the ``--enable-openmp`` flag.

.. _installation_MPI_linux:

MPI installation
""""""""""""""""

This step is optional.  Skip this chapter if you do not intend to use CHarm on 
distributed-memory system like HPC clusters.

Choose an MPI implementation (e.g., `OpenMPI <https://www.open-mpi.org/>`_, 
`MPICH <https://www.mpich.org/>`_) and either build the MPI library from source 
or install it using your package manager.

Importantly, the MPI implementation must support the MPI standard version 3.0 
(September 2012) or newer.  Since the MPI standard does not support the 
``__float128`` data type, CHarm in quadruple precision cannot be built with the 
MPI support.


.. _installation_MPFR_linux:

MPFR installation
"""""""""""""""""

This step is optional.  If you don't care about spectral gravity forward 
modelling with spatially limited integration radius, you may skip this section.

MPFR can either be installed using your package manager or it can be built from 
the source.

* *Using your package manager*

  This is the simplest and for most users also the recommended way of 
  installing MPFR.

  * Debian-based distributions:

    .. code-block:: console

         sudo apt install libmpfr-dev

  * Arch Linux: You may try to install, for instance, this package 
    `<https://archlinux.org/packages/core/x86_64/mpfr/>`_ using your package 
    manager.

  This should automatically install also `GMP <https://gmplib.org>`_, the GNU 
  multiple precision arithmetic library, which is required by MPFR.

* *Compilation from source*

  * Download, build and install the `GMP <https://gmplib.org>`_ library, which 
    is required to build MPFR.  The compilation and installation of GMP may 
    look like:

     .. code-block:: console

        ./configure
        make
        make check
        sudo make install

    You may want to specify a custom installation path by adding the 
    ``--prefix=/your/path/to/gmp`` flag to the ``./configure`` call (see the 
    installation manual of GMP).

  * Then, download, build and install `MPFR <https://mpfr.org>`_:

     .. code-block:: console

        ./configure
        make
        make check
        sudo make install

    If you used a custom installation path for GMP, use:

     .. code-block:: console

        ./configure LDFLAGS=-L/your/path/to/gmp/lib \
                    CPPFLAGS=-I/your/path/to/gmp/include
        make
        make check
        sudo make install

    Optionally, you may want to use a custom installation path also for MPFR 
    (again by using the ``--prefix`` flag).


.. _default_installation_charm_linux:

Default CHarm installation
""""""""""""""""""""""""""

If you

* want to install CHarm to your default installation path,
* have installed FFTW (version ``3.X.X``) to the default path available to the
  compiler,
* want a double precision version of CHarm,
* do not want the MPFR-related features,
* do not want SIMD parallelization,
* do not want OpenMP parallelization,
* do not want MPI parallelization,
* have root privileges,

you may simply execute the following commands:

.. code-block:: console

   ./configure
   make
   make check
   sudo make install

Briefly, ``./configure`` checks the availability of all components necessary to
build CHarm and prepares makefiles and a few other files.  ``make`` compiles
the library.  ``make check`` compiles and executes a test program.  ``make
install`` installs the library.


.. _customized_installation_charm_linux:

Customized CHarm installation
"""""""""""""""""""""""""""""

The installation process can be tailored by appending one or more of the
following flags to the ``./configure`` call.

* ``--enable-single-precision`` or ``--enable-double-precision`` or 
  ``--enable-quad-precision`` to compile CHarm in single, double or quadruple 
  precision, respectively (``float``, ``double`` and ``__float128`` data types 
  for floating point numbers, respectively).  If not specified, double 
  precision is used as default.

* ``--enable-sse4.1`` or ``--enable-avx`` or ``--enable-avx2`` or 
  ``--enable-avx-512`` or ``--enable-neon`` to enable SSE4.1, AVX, AVX2, 
  AVX-512 or NEON CPU instructions, respectively (all disabled by default).

  SSE4.1, AVX, AVX2 and AVX-512 are SIMD instruction sets introduced by Intel 
  in 2006, 2011, 2013 and 2017, respectively.  At least one of them is almost 
  surely available on any reasonably old x86_64 CPU.  The most critical number 
  crunching parts of CHarm are hand-written to take advantage of these 
  instructions in order to significantly improve the performance.  As a general 
  rule, it is strongly recommended to enable the newest instruction set that is 
  supported by your processor.  On many Linux distributions, the ``lscpu`` 
  utility prints all SIMD instruction sets supported by your processor.

  NEON is a SIMD instruction set available on ARM CPUs.  CHarm requires 64-bit 
  ARMv8 CPUs or newer to employ NEON.  In addition to attaching 
  ``--enable-neon``, you may also need to manually specify proper compiler 
  flags to enable NEON (there is a large number of the associated flags, so the 
  ``configure`` script does not guess the proper one(s)).  To specify the 
  compiler flags, use the ``CFLAGS`` environment variable (see below; with the 
  SSE4.1 and AVX family of CPUs, no additional compiler flags are needed).

  On the hardware level, SIMD instructions are not supported in quadruple 
  precision, therefore they can be enabled only when compiling in single or 
  double precision.

* ``--enable-openmp`` to enable OpenMP parallelization for shared-memory 
  architectures (no parallelization by default).

  The number of threads can be set either in your code by 
  ``omp_set_num_threads(N)`` or by using the ``OMP_NUM_THREADS`` environment 
  variable.

* ``--enable-mpi`` to enable MPI parallelization for shared- and 
  distributed-memory architectures (no parallelization by default).

  MPI parallelization combined with distributed-memory systems like HPC 
  clusters allows you to conduct spherical harmonic transforms up to a few 
  hundred thousands.  The basic idea is to distribute the signal and spherical 
  harmonic coefficients over several computing nodes, because these data may 
  consume hundreds of GBs of RAM.  Such an amount of memory is only rarely 
  available on a single shared-memory system.

  In addition to HPC clusters, you can also take advantage of MPI if you have 
  a few ordinary PCs connected through some network protocol (e.g., SSH).  This 
  will allow you to distribute your large data sets across your PCs.

  Finally, MPI works on shared-memory architectures, too.  However, in the case 
  of CHarm, there are generally no or little advantages over OpenMP.  For 
  shared-memory systems, most users should therefore prefer OpenMP over MPI.

  For best performance with high-degree spherical harmonic transforms, you can 
  (and in fact should) combine MPI with OpenMP and SIMD.

* ``--enable-mpfr`` to compile CHarm with the MPFR support enabling spectral 
  gravity forward modelling with spatially limited integration radius (disabled 
  by default).  GMP and MPFR libraries are mandatory if ``--enable-mpfr`` is 
  used.

* ``--prefix=/your/custom/path`` to specify a custom installation path for
  CHarm.

* There are two ways to specify your custom installation paths for FFTW and, 
  optionally, also for MPFR, GMP and MPI libraries.  If you installed FFTW 
  (possibly also MPFR, GMP and MPI) to default installation paths available to 
  the compiler, you may skip this bullet point.

  * The easiest way is to use ``--with-build-path=/path/to/fftw`` or, if you 
    want to compile, say, with the MPFR and MPI support, 
    ``--with-build-path=/path/to/fftw:/path/to/mpfr:/path/to/gmp:/path/to/mpi``.  
    Specify only the top directories of the libraries, that is, do not add the 
    ``include`` or ``lib`` subdirectories.

  * The other way is to use the ``LDFLAGS`` and ``CPPFLAGS`` variables.

    * To specify a custom path to your FFTW libs, use 
      ``LDFLAGS=-L/your/path/to/FFTW/lib``.

      Optionally, if you use, say, the ``--enable-mpfr`` installation flag, you 
      may need to provide also the path(s) for the GMP and MPFR libraries if 
      you installed these to custom paths, e.g., 
      ``LDFLAGS="-L/your/path/to/FFTW/lib -L/your/path/to/MPFR/lib 
      -L/your/path/to/GMP/lib"``.

      Specify only the path(s) to the library (libraries); the lib files 
      themselves are linked automatically.

    * To specify a custom path to your FFTW header file, use 
      ``CPPFLAGS=-I/your/path/to/FFTW/include``.

      Optionally, if you compile with, say, the ``--enable-mpfr`` installation 
      flag, you may need to provide also the path(s) for the GMP and MPFR 
      header files if you installed these to custom paths, e.g., 
      ``LDFLAGS="-L/your/path/to/FFTW/include -L/your/path/to/MPFR/include 
      -L/your/path/to/GMP/include"``.

    Do not forget to add the ``lib`` and ``include`` subdirectories to each 
    entry in ``LDFLAGS`` and ``CPPFLAGS``, respectively.

* ``--disable-shared`` to not compile CHarm as a shared library.

* Other useful variables:

  * ``CC`` selects other than your system's default C compiler,
    e.g. ``CC=clang`` for Clang.

  * ``CFLAGS`` defines user-defined compiler flags.  For instance, for solid 
    performance with GCC, you may want to use ``CFLAGS="-O3"`` or, if you know 
    what you are doing, ``CFLAGS="-O3 -ffast-math"``.

To get a summary of all supported flags, execute ``./configure --help``.


Customized CHarm installation (examples)
""""""""""""""""""""""""""""""""""""""""

* *Example 1*

  * Custom CHarm installation directory

  * Custom FFTW installation directory

  * Quadruple precision

  * OpenMP parallelization enabled

  * SIMD instructions disabled

  * MPI parallelization disabled

  * Spatially restricted spectral gravity forward modelling (MPFR) disabled

  .. code-block:: console

     ./configure --prefix=/opt/charm --enable-openmp --enable-quad-precision \
          LDFLAGS=-L/opt/fftwq-3.3.10/lib CPPFLAGS=-I/opt/fftwq-3.3.10/include
     make
     make check
     sudo make install


* *Example 2*

  * Custom CHarm installation directory

  * Custom FFTW installation directory

  * Custom GMP installation directory

  * Custom MPFR installation directory

  * Double precision

  * AVX SIMD instructions enabled

  * OpenMP parallelization enabled

  * MPI parallelization disabled

  * Aggressive optimizations during the compile time

  .. code-block:: console

     ./configure --prefix=/opt/charm \
     --with-build-path=/opt/fftw-3.3.10:/opt/mpfr-4.2.0:/opt/gmp-6.2.1 \
     --enable-openmp --enable-avx --enable-mpfr CFLAGS="-O3 -ffast-math"
     make
     make check
     sudo make install


.. _installation_mac:


Installation on macOS
~~~~~~~~~~~~~~~~~~~~~

At first, we will install FFTW (optionally also MPFR and MPI) and then we will 
proceed with the installation of CHarm.

FFTW installation
"""""""""""""""""

FFTW can either be installed via your package manager or built from the source,
preferably with GCC.  The latter is strongly recommended on macOS.

* *Using your package manager*

  You can use one of the following commands, depending on the package manager
  you use:

  .. code-block:: console

     sudo port install fftw-3


  .. code-block:: console

     brew install fftw

  This, however, most likely does not install FFTW in quadruple precision
  and/or with OpenMP support.  You may therefore be able to compile CHarm only
  in single or double precision with OpenMP disabled.

* *Compilation from source*

  It is recommended to compile FFTW using GCC.  If you do not have GCC
  installed yet, you may execute one of the following commands:

  .. code-block:: console

     sudo port install gcc10


  .. code-block:: console

     brew install gcc@10

  Now, you should be ready to build FFTW by following the instructions in the
  :ref:`installation_FFTW_linux` chapter (Linux, compilation from source).
  There is, however, one **important** additional remark.  When calling the
  FFTW's ``./configure`` script, specify also your GCC compiler, including its
  version number, e.g.:

  .. code-block:: console

      ./configure --enable-openmp CC=gcc-10

  Without the ``CC`` flag, Clang will most likely be used which may cause an 
  installation failure when using the ``--enable-openmp`` and/or 
  ``--enable-quad-precision`` flag(s).  It may **not** be sufficient to add 
  ``CC=gcc`` (GCC version number omitted), as this will still likely call 
  Clang.

MPFR installation
"""""""""""""""""

This step is optional.  If you don't care about spectral gravity forward 
modelling with spatially limited integration radius, you may skip this section.

You can install `MPFR <https://mpfr.org>`_ and `GMP <https://gmplib.org>`_, 
which is required by MPFR, either by your package manager or they can be built 
from the source.

* *Using your package manager*

  You may try to install MPFR using one of the following commands:

   .. code-block:: console

     sudo port install mpfr

   .. code-block:: console

     brew install mpfr

* *Compilation from source*

  See *Compilation from source* in :ref:`installation_MPFR_linux`.

MPI installation
""""""""""""""""

See :ref:`MPI installation on Linux <installation_MPI_linux>`.

CHarm installation
""""""""""""""""""

Having installed FFTW (and optionally GMP, MPFR, MPI), you may proceed with the 
same instructions as given in the :ref:`default_installation_charm_linux` and 
:ref:`customized_installation_charm_linux` chapters for Linux.  Similarly as 
when installing FFTW, it is recommended to use the GCC compiler via the ``CC`` 
variable when calling the ``./configure`` script from the CHarm installation.


A few installation notes
~~~~~~~~~~~~~~~~~~~~~~~~

* The output lib names depend on the precision used to compile CHarm:

   * ``libcharmf`` -- single precision,

   * ``libcharm`` -- double precision,

   * ``libcharmq`` -- quadruple precision.

* You may install CHarm in single, double and quadruple precision to the same 
  installation path.  You don't have to worry about overwriting the header and 
  lib files.


.. _charm_uninstallation:

Uninstallation
~~~~~~~~~~~~~~

Execute ``sudo make uninstall``.


.. _pyharm_installation:

PyHarm (the Python wrapper)
---------------------------

Before reading this chapter, make sure you know how to compile :ref:`CHarm 
<charm_installation>`.  Otherwise, you won't be able to build PyHarm.


Requirements
~~~~~~~~~~~~

*Additional* prerequisites when compared with :ref:`dependencies 
<charm_requirements>`:

* Python interpreter 3.6 or newer,

* Python module `pip <https://docs.python.org/3/installing/index.html>`_,

* Python module `numpy <https://numpy.org/>`_ (reasonably old version),

* Python module `ctypes <https://docs.python.org/3/library/ctypes.html>`_ 
  (reasonably old version).

PyHarm does not support quadruple precision and the MPI parallelization.


Building PyHarm
~~~~~~~~~~~~~~~

Installation of PyHarm is disabled by default.  To enable it, you have to add 
the ``--enable-python`` flag to the ``configure`` call *in addition* to the 
flags discussed in the :ref:`CHarm (the C library) <charm_installation>` 
chapter.

The following flags may be used in addition to ``--enable-python``.

* The ``PYTHON`` variable specifies the Python interpreter you want to use.  
  For instance, ``PYTHON=python3.9`` will ensure that the build is done 
  with/for Python version 3.9.  Use the appropriate version (depends on your 
  machine).

* By default, PyHarm is built to the ``${prefix}/lib`` directory.  The path in 
  ``${prefix}`` is taken from the ``--prefix`` flag (see :ref:`CHarm (the 
  C library) <charm_installation>`).  The default installation path can be 
  replaced by a custom one using the ``--with-python_prefix`` flag, for 
  instance, ``--with-python_prefix=/home/isaac/pyharm``.

  Using the correct path in ``--with-python_prefix`` is crucial for Python to 
  find PyHarm.  Otherwise, when calling

  .. code-block:: python

    >>> import pyharm

  from within the Python shell or a Python script, ``ModuleNotFoundError`` will 
  be thrown.

  There are several strategies to choose the installation path.

  * If you are not really confident with all this, create and activate a Python 
    virtual environment:

    .. code-block:: bash

        python3 -m venv /path/to/your/virtual/environment/
        source /path/to/your/virtual/environment/bin/activate

    Then use ``--with-python_prefix=/path/to/your/virtual/environment`` when 
    calling the ``configure`` script.  After executing ``make`` and ``make 
    install``, you are ready to import PyHarm in a Python shell or a Python 
    script:

    .. code-block:: python

        >>> import pyharm

  * If you want to install PyHarm as a user, find the lib path of your Python 
    user packages, for instance,

    .. code-block:: bash

        python3 -m site --user-site

    The output might like ``/home/isaac/.local/lib/python3.9/site-packages``, 
    depending on the version of your Python and on your OS.  Based on this 
    path, you can specify your installation path; in this case it is
    ``--with-python_prefix=/home/isaac/.local``.  Note that the 
    ``lib/python3.9/site-packages`` directories have to be omitted, as they are 
    added to the installation path automatically.

  * If you want to install PyHarm next to your system Python packages, you must 
    specify neither ``--prefix`` nor ``--with-python_prefix``.

  Note that if you use a custom path in ``--prefix`` but do not specify 
  ``--with-python_prefix``, you will most likely not be able to (easily) import 
  PyHarm.

Example installation
""""""""""""""""""""

A typical installation of PyHarm to a Python virtual environment looks like 
this:

.. code-block:: bash

  python3 -m venv /tmp/python-venv
  source /tmp/python-venv/bin/activate
  ./configure --prefix=/tmp/charm --enable-openmp \
     --with-build-path=/opt/fftw-3.3.10:/opt/mpfr-4.2.1:/opt/gmp-6.3.0 \
     --enable-python --enable-mpfr PYTHON=python3.9 \
     --with-python_prefix=/tmp/python-venv
  make
  make check
  make install

Then, open Python:

.. code-block:: bash

  python3

From within Python, you can now work with PyHarm:

.. code-block:: python

  >>> import pyharm as ph
  >>> ph.misc.print_info()
  >>> quit()

Deactivate the virtual environment from the shell:

.. code-block:: bash

  deactivate


Uninstallation
~~~~~~~~~~~~~~

See :ref:`charm_uninstallation`.
