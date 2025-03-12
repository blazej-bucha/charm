.. _py_api:

===================
PyHarm (Python API)
===================

This chapter documents the Python interface.

PyHarm is a Python wrapper for CHarm.  It relies on the Python's `ctypes 
<https://docs.python.org/3/library/ctypes.html>`_ module and is fully 
integrated with `NumPy <https://numpy.org/>`_ to work easily and efficiently 
with arrays.  Similarly as CHarm, PyHarm is divided into several modules, 
comprising classes and functions with conceptually similar tasks, such as 
spherical harmonic analysis or spherical harmonic synthesis.  Source files of 
the wrapper are located in ``wrap/pyharm``.  The wrapper is hand-written and 
follows the object-oriented paradigm of Python in order to provide a Pythonic 
feeling (sort of).  Therefore, there is no simple one-to-one relation between 
CHarm and PyHarm.

Next follows the documentation of the PyHarm modules in double precision. There 
are only a few easy-to-remember systematic differences in PyHarm’s API when 
working in single precision, so separate documentation is omitted. The list of 
the differences is given in :ref:`PyHarm_in_single_precision`.

.. toctree::
   :titlesonly:

   pyharm.shc <api-py-shc>
   pyharm.crd <api-py-crd>
   pyharm.shs <api-py-shs>
   pyharm.sha <api-py-sha>
   pyharm.gfm <api-py-gfm>
   pyharm.leg <api-py-leg>
   pyharm.integ <api-py-integ>
   pyharm.misc <api-py-misc>
   pyharm.glob <api-py-glob>
