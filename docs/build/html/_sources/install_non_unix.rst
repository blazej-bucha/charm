.. _installation_win:

Installation on Windows 10
==========================

Perhaps the easiest way is to install the `Windows Subsystem for Linux
<https://en.wikipedia.org/wiki/Windows_Subsystem_for_Linux>`_ (WSL).  WSL
offers a fully functional Linux distribution that you can run from within
Windows 10 (no dual boot).  WSL is officially supported by Microsoft and its
installation is very easy.  Having installed WSL, you may proceed with
:ref:`installation_unix`.

Other solution is, for instance, `Cygwin 
<https://en.wikipedia.org/wiki/Cygwin>`_.  A note on enabling SIMD instructions 
(``--enable-avx``, ``--enable-avx2`` or ``--enable-avx-512``) under Cygwin is 
necessary here.  To properly compile CHarm under Cygwin with SIMD instructions 
enabled, do *not* use GCC.  Under Cygwin, GCC doesn't seem to properly align 
the memory (required to execute SIMD CPU instructions).  You won't get any 
compilation error, but when running your program or the check program (``make 
check``), you will most likely get a segfault error.  Clang, for instance, does 
its job properly, even under Cygwin, so is recommended in this case.
