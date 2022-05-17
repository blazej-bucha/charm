.. _installation_unix:

Installation on Unix-based operating systems
============================================

This section provides an installation guide for Unix-based operating systems,
particularly :ref:`Linux <installation_linux>` and :ref:`macOS
<installation_mac>`.  The instructions that follow will compile CHarm as
a static library and, if specifically asked for, also as a shared library.


.. _requirements:

Requirements
------------

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
---------------------

At first, we will install FFTW and then we will proceed with the installation
of CHarm.

.. _installation_FFTW_linux:

FFTW installation
^^^^^^^^^^^^^^^^^

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
^^^^^^^^^^^^^^^^^^^^^^^^^^

If you

* want to install CHarm to ``/usr/local``,
* have installed FFTW (version ``3.X.X``) to the default path available to the
  compiler,
* want a double precision version of CHarm,
* do not want parallelization via OpenMP, and
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
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The installation process can be tailored by appending one or more of the
following flags to the ``./configure`` call:

* ``--prefix=/your/custom/path`` to specify a custom installation path for
  CHarm (default is ``--prefix=/usr/local``).

* ``--enable-single-precision`` to compile CHarm in single precision (``float``
  data type for floating point numbers), **or**

  ``--enable-double-precision`` for double precision (``double`` data type),
  **or**

  ``--enable-quad-precision`` for quadruple precision (``__float128`` data
  type).

  If not specified, double precision is used as default.

* ``--enable-openmp`` to enable parallelization via OpenMP (no parallelization
  by default).

  The number of threads can be set either in your code by
  ``omp_set_num_threads(N)`` or via the ``OMP_NUM_THREADS`` environment
  variable.

* ``--enable-shared`` to compile CHarm as a shared library *in addition* to the
  static library.

* ``LDFLAGS=-L/your/path/to/FFTW/lib`` to specify a custom path to your FFTW
  libs (empty by default, that is, default is to assume that FFTW is accessible
  to the compiler).

  You only need to specify the path to the FFTW library; the lib files
  themselves are linked automatically.

* ``CPPFLAGS=-I/your/path/to/FFTW/include`` to specify a custom path to your
  FFTW header file (empty by default, that is, default is to assume that FFTW
  is accessible to the compiler).

* Other useful variables (see the Autotools documentation for further options):

  * ``CC`` selects other than your system's default C compiler,
    e.g. ``CC=clang`` for the ``Clang`` compiler, and

  * ``CFLAGS`` defines user-defined compiler flags, e.g.,  ``CFLAGS="-O3"``
    (GCC).

A possible installation

* with a custom CHarm installation directory,

* with a custom FFTW installation directory,

* in quadruple precision,

* with parallelization via OpenMP enabled, and

* with the shared library, too,

could look like

.. code-block:: console

   ./configure --prefix=/opt/charm --enable-openmp --enable-shared \
        --enable-quad-precision LDFLAGS=-L/opt/fftwq-3.3.9/lib \
        CPPFLAGS=-I/opt/fftwq-3.3.9/include
   make
   make check
   sudo make install


.. _installation_mac:


Installation on macOS
---------------------

At first, we will install FFTW and then we will proceed with the installation
of CHarm.

FFTW installation
^^^^^^^^^^^^^^^^^

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

      ./configure --enable-openmp --enable-shared CC=gcc-10

  Without the ``CC`` flag, the ``Clang`` compiler will most likely be used
  which may cause an installation failure when using the ``--enable-openmp``
  and/or ``--enable-quad-precision`` flag(s).  It may **not** be sufficient to
  add ``CC=gcc`` (GCC version number omitted), as this will still likely call
  the ``Clang`` compiler.

CHarm installation
^^^^^^^^^^^^^^^^^^

Having installed FFTW, you may proceed with the same instructions as given in
the :ref:`default_installation_charm_linux` and
:ref:`customized_installation_charm_linux` chapters for Linux.  Similarly as
when installing FFTW, it is recommended to use the GCC compiler via the ``CC``
variable when calling the ``./configure`` script from the CHarm installation.


A few installation notes
------------------------

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


Uninstallation
--------------

Execute ``sudo make uninstall``.
