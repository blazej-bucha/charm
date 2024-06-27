.. _install_pyharm_pip:

===============================
Installation of PyHarm with pip
===============================

On Linux, macOS and Windows running on the x86_64 CPU architecture, install 
PyHarm using the command:

.. code-block:: bash

    pip install pyharm

This will install the wrapper Python code of PyHarm and a pre-compiled CHarm 
library, the latter of which is called by PyHarm routines.

The CHarm binaries were compiled using the following settings.

* ``gcc`` compiler on Linux and macOS and ``MSVC`` on Windows.
* ``--enable-double-precision`` to compile in double precision.
* ``--enable-openmp`` to enable parallelization among CPU cores.
* ``--enable-avx`` to enable SIMD parallelism using the AVX instruction sets of 
  x86_64 CPUs.  AVX has been with us for more than 10 years now, so the 
  binaries should run on any reasonably old CPU.  Should you get an error like 
  ``Illegal instruction``, this means your CPU does not support AVX.  In that 
  case, you have to build PyHarm from source and avoid using any set of AVX 
  instructions (see :ref:`build_from_src`).
* ``-O3 -ffast-math -Wall -Wpedantic`` compiler flags on Linux and macOS and 
  ``/O2 /fp:fast /FS /GL`` on Windows to optimize the library for best 
  performance.  Importantly, you should be aware that ``-ffast-math`` and 
  ``/fp:fast`` relax some rules on floating point arithmetics in favour of 
  gaining some additional computing speed.  As a result, the accuracy may be 
  slightly *worse* than when not using these flags.  Normally, however, one 
  does not need to bother by this.  Should you need compliance with the IEEE 
  specifications, compile the library on your own without these flags.
* ``--enable-python`` to build the PyHarm wrapper.

.. tip::

   All this information is also available directly from PyHarm by calling:

    .. code-block:: python

        >>> import pyharm as ph
        >>> ph.misc.print_info()

   This will print some useful compilation details about your particular PyHarm 
   build.

An out-of-the-box support for ARM CPUs with NEON instructions is on the TODO 
list.
