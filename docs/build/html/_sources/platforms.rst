.. _tested_platforms:

================
Tested platforms
================

.. _platforms:
.. list-table:: Known supported configurations
   :header-rows: 1
   :widths: 20 40 40

   * - Architecture
     - Operating System
     - C Compiler

   * - x86_64
     - Debian GNU/Linux 11 (bullseye)
     - GCC 10.2.1, Clang 11.0.1-2, ICC 2021.5.0

   * - x86_64
     - Scientific Linux release 6.4 (Carbon)
     - GCC 7.2, GCC 6.4, GCC 5.4, GCC 4.9.3, GCC 4.8.4, GCC 4.4.7

   * - x86_64
     - Manjaro Linux
     - GCC 11.2.0

   * - x86_64
     - Alpine Linux 3.17.0 (using Docker)
     - GCC 12.2.1 with musl 1.2.3-r4 as libc

   * - x86_64
     - FreeBSD 13.0-RELEASE
     - Clang 11.0.1, GCC 10.3.0

   * - x86_64
     - macOS Big Sur 11.4
     - Clang 12.0.5, GCC 10.3.0

   * - x86_64
     - Windows 10 + Windows Subsystem for Linux (Debian GNU/Linux, "bullseye" 
       release)
     - GCC 10.2.1

   * - x86_64
     - Windows 10 + Cygwin 3.4
     - GCC 11.4.0-1 [#f1]_

   * - x86_64
     - Windows 10 + Cygwin 3.4
     - Clang 8.0.1-1 [#f2]_

   * - x86_64
     - Windows 11
     - MSVC 14.40.33807

.. [#f1] SIMD instructions disabled (see :ref:`build_from_src_non_unix`).
.. [#f2] SIMD instructions enabled (see :ref:`build_from_src_non_unix`).


We would love to hear your own experience, especially with platforms not listed 
in the table!
