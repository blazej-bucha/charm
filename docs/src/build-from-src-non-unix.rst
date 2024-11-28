.. _build_from_src_non_unix:

Building from source on Windows 10 and later
============================================

Perhaps the easiest way is to install the `Windows Subsystem for Linux
<https://en.wikipedia.org/wiki/Windows_Subsystem_for_Linux>`_ (WSL).  WSL
offers a fully functional Linux distribution that you can run from within
Windows 10 (no dual boot).  WSL is officially supported by Microsoft and its
installation is very easy.  Having installed WSL, you may proceed with
:ref:`build_from_src_unix`.

Other solution is, for instance, `Cygwin 
<https://en.wikipedia.org/wiki/Cygwin>`_.  However, to properly compile CHarm 
under Cygwin with SIMD instructions enabled (``--enable-avx``, 
``--enable-avx2`` or ``--enable-avx-512``), do *not* use GCC.  Under Cygwin, 
GCC doesn't seem to properly align the memory which is required to execute SIMD 
CPU instructions.  You won't get any compilation error, but when running your 
program or the check program (``make check``), you will most likely get 
a segfault error.  Clang, for instance, does its job properly, even under 
Cygwin, so is recommended in this case.

The last path you can (but do not want to) take is to build the library using 
the Microsoft's MSVC compiler.  To make it work, install some layer that will 
add the standard installation utilities (e.g. make) that are not found on 
Windows normally.  Next, you may want to use the `CCCL 
<https://github.com/swig/cccl>`_ wrapper for MSVC.  Finally, follow the 
instructions found in :ref:`build_from_src_unix`.  In fact, this is how we 
built the Python wheels for Windows that are available from PyPI.
