.. _pyharm_cookbook:

======
PyHarm
======

This chapter demonstrates how to perform some basic tasks with
PyHarm.  For a detailed description of the classes, methods, functions, etc. to 
be used, see the :ref:`py_api` chapter.  The source codes of the examples can 
be found in the ``cookbook/python`` directory.

In your Python code, import PyHarm in double precision by:

.. code-block:: python

   >>> import pyharm as ph

To import PyHarm in single precision, use:

.. code-block:: python

   >>> import pyharmf as phf

Both precisions can be combined in a single Python file, too:

.. code-block:: python

   >>> import pyharm as ph
   >>> import pyharmf as phf

The chapters that follow show how to work with PyHarm in double and single 
precision, respectively.  Quadruple precision is not supported in PyHarm.


PyHarm in double precision
==========================

PyHarm attempts to follow the rules of the object-oriented programming 
paradigm.  Several classes are introduced to represent spherical harmonic 
coefficients, coordinates of evaluation points/cells and Fourier coefficients 
of Legendre functions.  Most of the classes are accompanied by factory methods 
for instantiation, so users should always instantiate the classes using these 
methods.  The factory methods are summarized within the class documentation.  
Alternatively, they can be found by searching for the *classmethod* flag that 
precedes the class method names in the documentation.

Integration with NumPy
----------------------

PyHarm supports NumPy arrays.  All arrays consisting of floating point data 
must

* be instances of the ``numpy.ndarray`` class,

* have the ``dtype`` flag set to ``numpy.float64`` (or to ``numpy.float32`` if 
  the library was compiled in single precision), and

* have the ``flags.c_contiguous`` attribute set to ``True``.

As an example, let's create a 1D numpy array ``x`` that meets all of these 
conditions.

.. code-block:: python

    >>> import numpy as np
    >>> x = np.zeros((5,))
    >>> x.dtype
    dtype('float64')
    >>> x.flags.c_contiguous
    True

In this code snippet, we did not specify the ``dtype`` and ``order`` flags when 
calling ``numpy.zeros``,  so numpy automatically used their default values, 
namely ``numpy.float64`` and ``'C'``.  In some cases, you might need to 
explicitly specify the two flags:

.. code-block:: python

    >>> import numpy as np
    >>> x = np.zeros((5,), dtype=np.float64, order='C')
    >>> x.dtype
    dtype('float64')
    >>> x.flags.c_contiguous
    True

For further details, see the `NumPy <https://numpy.org/>`_ documentation.

If an array passed to PyHarm does not meet all the required criteria, PyHarm 
will throw an error.

Working with spherical harmonic coefficients
--------------------------------------------

This example shows how to read/write spherical harmonic coefficients from/to
files and how to compute (difference) degree variances.

.. literalinclude:: ../../cookbook/python/shcs.py
   :language: python


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

.. literalinclude:: ../../cookbook/python/sha_shs.py
   :language: python


First- and second-order gradients (gravitational vector and tensor)
-------------------------------------------------------------------

This example shows how to compute the full gravitational vector and the full
gravitational tensor in the :ref:`lnof`.

.. literalinclude:: ../../cookbook/python/gradients.py
   :language: python


Gravity forward modelling
-------------------------

Given spherical harmonic coefficients of the Moon's topography and of its 
density (constant, lateral and 3D density), this example computes the 
gravitational field implied by the lunar topographic masses.

The lunar shape is here defined by the ``MoonTopo2600p.shape`` model due to

 * Wieczorek, M. A., Gravity and topography of the terrestrial planets, 
   Treatise on Geophysics, 10, 153-193, doi:10.1016/B978-0-444-53802-4.00169-X, 
   2015.

The density model is due to

 * Goossens, S., Sabaka, T.J., Wieczorek, M.A., Neumann, G.A., Mazarico, E., 
   Lemoine, F., Nicholas, J.B., Smith, D.E., Zuber, 
   M.T. (2020). High-resolution gravity field models from GRAIL data and 
   implications for models of the density structure of the Moon's crust, 
   Journal of Geophysical Research: Planets, 125, e2019JE006086, 
   doi:10.1029/2019JE006086.

The density model of Goossens et al. (2020) was transformed into the form 
required by our forward modelling method.  While the transformed models that we 
use in this example are realistic, it is not recommended to use them in any 
real task.  They were create for learning purposes only.

.. literalinclude:: ../../cookbook/python/gfm.py
   :language: python


Fourier coefficients of Legendre functions
------------------------------------------

.. literalinclude:: ../../cookbook/python/fourier_legendre.py
   :language: python


Integrals
---------

Finally, let's compute an integral of a product of two spherical harmonics over
a restricted domain on the unit sphere.  Then, we do the same, but with the
product of Legendre functions only.

.. literalinclude:: ../../cookbook/python/integrals.py
   :language: python


.. _PyHarm_in_single_precision:

PyHarm in single precision
==========================

A few simple rules need to be obeyed when working with PyHarm in single 
precision.

* Make sure you have compiled the library in single precision (see 
  :ref:`Installing`).

* In Python, import the PyHarm package as

  .. code:: python

    >>> import pyharmf as phf

* All floating point constants must be of the ``numpy.float32`` data type.  For 
  instance, while in double precision you can write this:

  .. code:: python

    >>> import pyharm as ph
    >>> shcs = ph.shc.Shc.from_zeros(10, 1.0, 1.0)

  in single precision, you have to write:

  .. code:: python

    >>> import numpy as np
    >>> import pyharmf as phf
    >>> shcs = phf.shc.Shc.from_zeros(10, np.float32(1.0), np.float32(1.0))

* All floating point numpy arrays must be of the ``numpy.float32`` data type.  
  So when creating a numpy array to be passed to PyHarm, you must set the 
  ``dtype`` flag to ``numpy.float32``:

  .. code:: python

    >>> import numpy as np
    >>> import pyharmf as phf
    >>> x = np.zeros((5,), dtype=np.float32)
    >>> phf.crd.PointSctr.from_arrays(x, x, x)

