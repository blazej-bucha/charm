.. _cookbook:

========
Cookbook
========

This chapter offers a brief cookbook on how to perform some basic tasks with
CHarm.  For a detailed description of the functions to be used, see the
:ref:`Modules` chapter.  The source codes of the examples can be found in the
``cookbook`` directory.



CHarm in double precision
=========================

In this section, we assume you have successfully installed CHarm (see
:ref:`Installing`) in double precision with OpenMP enabled into the default
installation paths, and that you know how to compile a C code and link external
libraries.  A brief example of the compilation under Linux with GCC is provided
in the :ref:`Compilation_Linux` section below.


.. _Compilation_Linux:

Compilation on Linux
--------------------

Assuming your working directory is ``cookbook``, an example compilation
with the static library might look like this:

.. code-block:: shell

   gcc -fopenmp shcs.c -l:libcharm_omp.a -lfftw3 -lfftw3_omp -lm

Compilation with the shared library:

.. code-block:: shell

   gcc -fopenmp -Wl,-rpath -Wl,/usr/local/lib shcs.c \
        -lcharm_omp -lfftw3 -lfftw3_omp -lm

After a successful compilation, you are ready to execute the compiled binary as

.. code-block:: shell

   ./a.out

You should see the following output:

.. code-block:: none

    C(  9,  0) = 2.7671430085300001e-08
    S(  9,  0) = 0.0000000000000000e+00
    C(  9,  4) = -9.0017922533600003e-09
    S(  9,  4) = 1.9466677947499999e-08
    C(  9,  9) = -4.7747538613200003e-08
    S(  9,  9) = 9.6641284771399995e-08


    Degree variance for harmonic degree 0 = 1.0000000000000000e+00
    Degree variance for harmonic degree 4 = 2.5183167531865940e-12
    Degree variance for harmonic degree 10 = 1.2631494966213168e-13


    Difference degree variance for harmonic degree 0 = 0.0000000000000000e+00
    Difference degree variance for harmonic degree 4 = 0.0000000000000000e+00
    Difference degree variance for harmonic degree 10 = 0.0000000000000000e+00

    Great, all done!

.. note::

   Various approaches to compile the code with the static/shared library are
   possible.  Use whatever works best for you.


Working with spherical harmonic coefficients
--------------------------------------------

This example shows how to read/write spherical harmonic coefficients from/to
files and how to compute (difference) degree variances.

.. literalinclude:: ../../cookbook/shcs.c
   :language: c


Spherical harmonic synthesis and analysis
-----------------------------------------

This section shows several examples.

* A closed-loop test of analysis and synthesis with point data values.  First,
  spherical harmonic coefficients are loaded from a text file.  Then, they are
  used to synthesize a signal, from which a new set of coefficients is finally
  computed by surface spherical harmonic analysis.  The two coefficients sets
  are then finally compared (should be equal, that is, the difference degree
  amplitudes should be negligibly small).

* A solid spherical harmonic synthesis of point values at scattered points.

* A solid spherical harmonic synthesis of point values at a custom grid of
  points.

* A solid spherical harmonic synthesis of block-mean values at scattered cells.

* A solid spherical harmonic synthesis of block mean values at a grid of cells.

* Surface spherical harmonic analysis with block-mean values in cells.

.. literalinclude:: ../../cookbook/sha_shs.c
   :language: c


Fourier coefficients of Legendre functions
------------------------------------------

.. literalinclude:: ../../cookbook/fourier_legendre.c
   :language: c


Integrals
---------

Finally, let's compute an integral of a product of two spherical harmonics over
a restricted domain on the unit sphere.  Then, we do the same, but with the
product of Legendre functions only.

.. literalinclude:: ../../cookbook/integrals.c
   :language: c


.. _CHarm_in_single_and_quad_precision:

CHarm in single and quadruple precision
=======================================

A few simple rules need to be obeyed in single and quadruple precision.

* Make sure you have compiled CHarm in single or quadruple precision (see
  :ref:`Installing`), whichever you need.

* Replace every occurrence of the ``charm_*`` prefix with ``charmf_*`` for
  single precision and ``charmq_*`` for quadruple precision.  This applies to
  your own code, including the header files, and the compilation of your
  program.

  If the prefix is written in capital letters, do *not* alter it.  The
  ``CHARM_*`` symbols are common across all precisions of CHarm.

* In your code, replace every ``double`` variable entering CHarm functions by
  ``float`` for single precision or by ``__float128`` for quadruple precision.

* To compile your code in quadruple precision, link the ``libquadmath``
  library, that is, add ``-lquadmath`` before ``-lm`` (see
  :ref:`Compilation_Linux`).

* When compiling your program, use the ``fftwf*`` and ``fftwq*`` prefixes to
  link FFTW in single and quadruple precision, respectively.


Example compilation in single precision
---------------------------------------

Assuming your working directory is ``cookbook``, an example compilation
with the static library might look like this:

.. code-block:: shell

   gcc -fopenmp shcsf.c -l:libcharmf_omp.a -lfftw3f -lfftw3f_omp -lm

Compilation with the shared library:

.. code-block:: shell

   gcc -fopenmp -Wl,-rpath -Wl,/usr/local/lib shcsf.c \
        -lcharmf_omp -lfftw3f -lfftw3f_omp -lm


Example compilation in quadruple precision
------------------------------------------

Static library:

.. code-block:: shell

   gcc -fopenmp shcsq.c -l:libcharmq_omp.a -lfftw3q -lfftw3q_omp \
        -lquadmath -lm

Shared library:

.. code-block:: shell

   gcc -fopenmp -Wl,-rpath -Wl,/usr/local/lib shcsq.c \
        -lcharmq_omp -lfftw3q -lfftw3q_omp -lquadmath -lm
