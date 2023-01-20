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
<https://en.wikipedia.org/wiki/Cygwin>`_.  We have noticed, however, that when 
combined with SIMD instructions, Cygwin throws segmentation fault errors, even 
for some clearly valid parts of the code.  With SIMD instructions disabled, the 
library seems to work fine, but the performance is suboptimal.  OpenMP 
parallelization also seems to work properly.
