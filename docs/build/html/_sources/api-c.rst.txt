.. _c_api:

=============
CHarm (C API)
=============

This chapter documents the C interface.

CHarm is divided into several modules, comprising functions
with conceptually similar tasks, such as spherical harmonic analysis or
spherical harmonic synthesis.  Each module has its own subdirectory inside the
``src`` directory, where the source files are stored.

Next follows the documentation of the CHarm modules in double precision.  There
are only a few easy-to-remember systematic differences in CHarm's API when
working in single or quadruple precision, so separate documentations are
omitted.  The list of the differences is given in
:ref:`CHarm_in_single_and_quad_precision`.

.. toctree::

   api-c-shc
   api-c-crd
   api-c-shs
   api-c-sha
   api-c-leg
   api-c-integ
   api-c-misc
   api-c-xnum
   api-c-err
   api-c-glob
   api-c-charm
