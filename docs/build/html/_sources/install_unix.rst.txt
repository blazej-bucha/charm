.. _installation_unix:

Installation on Unix-based operating systems
============================================

This section provides an installation guide for Unix-based operating systems,
particularly :ref:`Linux <installation_linux>` and :ref:`macOS
<installation_mac>`.  First, discussed is the installation of :ref:`CHarm (the 
C library) <charm_installation>` and then follow instructions on how to build 
:ref:`PyHarm (the Python wrapper) <pyharm_installation>`.


.. _charm_installation:

CHarm (the C library)
---------------------


.. _charm_requirements:

Requirements
~~~~~~~~~~~~

* C compiler supporting C99.  `GCC <https://gcc.gnu.org/>`_ is recommended 
  (available by default on most Linux distributions).  Other compilers that are 
  known to successfully compile CHarm are listed in the :ref:`tested_platforms` 
  chapter.

  *To compile in quadruple precision, GCC is mandatory (v4.6 or later)* and no
  other compiler is allowed.

* `Make <https://www.gnu.org/software/make/>`_ to build the library (available 
  by default on most Linux distributions).

* `pkg-config <https://www.freedesktop.org/wiki/Software/pkg-config/>`_ for
  managing library compile and link flags (available by default on most Linux
  distributions).

* `FFTW <http://www.fftw.org/>`_ for discrete fast Fourier transforms.


.. _installation_linux:

Installation on Linux
~~~~~~~~~~~~~~~~~~~~~

At first, we will install FFTW and then we will proceed with the installation
of CHarm.

.. _installation_FFTW_linux:

FFTW installation
"""""""""""""""""

FFTW can either be installed via you package manager or built from the source.

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
  ``configure`` script.  If you do not want to parallelize CHarm, you may omit
  the ``--enable-openmp`` flag.


.. _default_installation_charm_linux:

Default CHarm installation
""""""""""""""""""""""""""

If you

* want to install CHarm to ``/usr/local``,
* have installed FFTW (version ``3.X.X``) to the default path available to the
  compiler,
* want a double precision version of CHarm,
* do not want OpenMP parallelization,
* do not want to enable SIMD CPU instructions, and
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

* ``--enable-openmp`` to enable OpenMP parallelization (no parallelization by 
  default).

  The number of threads can be set either in your code by 
  ``omp_set_num_threads(N)`` or by using the ``OMP_NUM_THREADS`` environment 
  variable.

* ``--enable-avx`` or ``--enable-avx2`` or ``--enable-avx-512`` to enable AVX, 
  AVX2 or AVX-512 CPU instructions, respectively (all disabled by default).

  AVX, AVX2 and AVX-512 are SIMD instructions introduced by Intel in 2011, 2013 
  and 2017, respectively.  The most critical number crunching parts of CHarm 
  are hand-written to take advantage of these instructions in order to 
  significantly improve the performance.  As a general rule, it is strongly 
  recommended to enable the latest set of AVX instructions that are supported 
  by your processor.  On many Linux distributions, you can find all the 
  supported CPU instructions by executing ``lscpu``.  On the hardware level, 
  SIMD instructions are not supported in quadruple precision, thus can be 
  enabled only when compiling in single or double precision.

* ``--prefix=/your/custom/path`` to specify a custom installation path for
  CHarm (default is ``--prefix=/usr/local``).

* ``LDFLAGS=-L/your/path/to/FFTW/lib`` to specify a custom path to your FFTW
  libs (empty by default, that is, default is to assume that FFTW is accessible
  to the compiler).

  You only need to specify the path to the FFTW library; the lib files
  themselves are linked automatically.

* ``CPPFLAGS=-I/your/path/to/FFTW/include`` to specify a custom path to your
  FFTW header file (empty by default, that is, default is to assume that FFTW
  is accessible to the compiler).

* ``--disable-shared`` to not compile CHarm as a shared library.

* Other useful variables:

  * ``CC`` selects other than your system's default C compiler,
    e.g. ``CC=clang`` for the ``Clang`` compiler, and

  * ``CFLAGS`` defines user-defined compiler flags, e.g.,  ``CFLAGS="-O3 
    -ffast-math"``
    (GCC).

* To get a summary of all the supported flags, execute ``./configure --help``.

An example installation

* with a custom CHarm installation directory,

* with a custom FFTW installation directory,

* in quadruple precision,

* with OpenMP parallelization enabled, and

* with SIMD instructions disabled

looks like:

.. code-block:: console

   ./configure --prefix=/opt/charm --enable-openmp --enable-quad-precision \
        LDFLAGS=-L/opt/fftwq-3.3.9/lib \
        CPPFLAGS=-I/opt/fftwq-3.3.9/include
   make
   make check
   sudo make install


.. _installation_mac:


Installation on macOS
~~~~~~~~~~~~~~~~~~~~~

At first, we will install FFTW and then we will proceed with the installation
of CHarm.

FFTW installation
"""""""""""""""""

FFTW can either be installed via you package manager or built from the source,
preferably with GCC.  The latter is strongly recommended on macOS.

* *Using your package manager*

  You can use one of the following commands, depending on the package manager
  you use:

  .. code-block:: console

     sudo port install fftw-3
     brew install fftw

  This, however, most likely does not install FFTW in quadruple precision
  and/or with OpenMP support.  You may therefore be able to compile CHarm only
  in single or double precision with OpenMP disabled.

* *Compilation from source*

  It is recommended to compile FFTW using GCC.  If you do not have GCC
  installed yet, you may execute one of the following commands:

  .. code-block:: console

     sudo port install gcc10
     brew install gcc@10

  Now, you should be ready to build FFTW by following the instructions in the
  :ref:`installation_FFTW_linux` chapter (Linux, compilation from source).
  There is, however, one **important** additional remark.  When calling the
  FFTW's ``./configure`` script, specify also your GCC compiler, including its
  version number, e.g.:

  .. code-block:: console

      ./configure --enable-openmp CC=gcc-10

  Without the ``CC`` flag, the ``Clang`` compiler will most likely be used
  which may cause an installation failure when using the ``--enable-openmp``
  and/or ``--enable-quad-precision`` flag(s).  It may **not** be sufficient to
  add ``CC=gcc`` (GCC version number omitted), as this will still likely call
  the ``Clang`` compiler.

CHarm installation
""""""""""""""""""

Having installed FFTW, you may proceed with the same instructions as given in
the :ref:`default_installation_charm_linux` and
:ref:`customized_installation_charm_linux` chapters for Linux.  Similarly as
when installing FFTW, it is recommended to use the GCC compiler via the ``CC``
variable when calling the ``./configure`` script from the CHarm installation.


A few installation notes
~~~~~~~~~~~~~~~~~~~~~~~~

* The output lib names depend on the user-defined compilation settings and
  follow the pattern:

   * ``libcharmf`` -- single precision with OpenMP disabled,

   * ``libcharmf_omp`` -- single precision with OpenMP enabled,

   * ``libcharm`` -- double precision with OpenMP disabled,

   * ``libcharm_omp`` --  double precision with OpenMP enabled,


   * ``libcharmq`` -- quadruple precision with OpenMP disabled,

   * ``libcharmq_omp`` -- quadruple precision with OpenMP enabled.

* You may install CHarm in single, double and quadruple precision, each with
  OpenMP enabled and disabled, to the same installation path.  You don't have
  to worry about overwriting the header and lib files.


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

*Additional* prerequisites when compared with :ref:`requirements 
<charm_requirements>`:

* Python interpreter 3.6 or newer,

* Python module `pip <https://docs.python.org/3/installing/index.html>`_,

* Python module `numpy <https://numpy.org/>`_ (reasonably old version),

* Python module `ctypes <https://docs.python.org/3/library/ctypes.html>`_ 
  (reasonably old version).


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
     LDFLAGS=-L/opt/fftw-3.3.9/lib CPPFLAGS=-I/opt/fftw-3.3.9/include \
     --enable-python PYTHON=python3.9 --with-python_prefix=/tmp/python-venv
  make
  make check
  make install


Open Python:

.. code-block:: bash

  python3

From within Python, you can now work with PyHarm:

.. code-block:: python

  >>> import pyharm as ph
  >>> ph.misc.print_version()
  >>> quit()

Deactivate the virtual environment from the shell:

.. code-block:: bash

  deactivate


Uninstallation
~~~~~~~~~~~~~~

See :ref:`charm_uninstallation`.
