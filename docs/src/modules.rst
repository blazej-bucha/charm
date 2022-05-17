.. _modules:

=======
Modules
=======

The source code of CHarm is divided into several modules, comprising functions
with conceptually similar tasks, such as spherical harmonic analysis or
spherical harmonic synthesis.  Each module has its own subdirectory inside the
``src`` directory, where the source files are stored.

A single header file of the entire library is available as:

* ``charm/charmf.h`` in single precision,

* ``charm/charm.h`` in double precision, and

* ``charm/charmq.h`` in quadruple precision.

In your code, you may either include the single header file to import the 
entire library or you may include only the header file(s) of the module(s), 
with which you need to interact.  You may combine various precisions in 
a single program, provided that you follow the rules from :ref:`cookbook`.  It 
is not recommended though to combine various CHarm releases in your program 
(e.g., single precision from v0.0.0 and double precision form v0.1.0).

Next follows the documentation of the CHarm modules in double precision.  There
are only a few easy-to-remember systematic differences in CHarm's API when
working in single or quadruple precision, so separate documentations are
omitted.  The list of the differences is given in
:ref:`CHarm_in_single_and_quad_precision`.

.. toctree::

   shc
   crd
   shs
   sha
   leg
   integ
   misc
   xnum
   err
   glob
   charm
