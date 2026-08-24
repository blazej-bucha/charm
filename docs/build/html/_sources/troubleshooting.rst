===============
Troubleshooting
===============

This chapter discusses some relevant issues reported by users.

OMP: Error #15:
===============

**Error message:**

.. code-block:: none

   OMP: Error #15: Initializing libiomp5md.dll, but found libomp140.x86_64.dll 
   already initialized.

   OMP: Hint This means that multiple copies of the OpenMP runtime have been 
   linked into the program. That is dangerous, since it can degrade performance 
   or cause incorrect results. The best thing to do is to ensure that only 
   a single OpenMP runtime is linked into the process, e.g. by avoiding static 
   linking of the OpenMP runtime in any library. As an unsafe, unsupported, 
   undocumented workaround you can set the environment variable 
   KMP_DUPLICATE_LIB_OK=TRUE to allow the program to continue to execute, but 
   that may cause crashes or silently produce incorrect results. For more 
   information, please see http://www.intel.com/software/products/support/.

   Fatal Python error: Aborted

**Root cause:** Two or more of your Python modules are linked against different 
OpenMP implementations (e.g., PyHarm and Matplotlib).

**Solution:** You have to ensure that only one OpenMP library is loaded during 
the runtime.  There are several ways to achieve this, each with its pros and 
cons.

* *Solution 1*: Follow the advice from the error message and put this at the 
  top of your code, that is, before any ``import`` statements:

  .. code:: python

     >>> import os
     >>> os.environ['KMP_DUPLICATE_LIB_OK'] = 'True'

  This, however, is only a workaround and not a solution, as is clear from the 
  error message.  Having said this, users report that Python behaves as they 
  expect after adding this to their code.

* *Solution 2*: If you have installed ``numpy`` using the default ``conda``'s 
  channel, the `MKL 
  <https://anaconda.org/channels/anaconda/packages/mkl/overview>`_ module may 
  have been installed on your machine during the installation.  ``mkl`` is the 
  Intel's optimized math library that includes the Intel's OpenMP 
  implementation.  On the other hand, the PyHarm wheels distributed through 
  ``pip`` are linked against other than the Intel's OpenMP library.  Therefore, 
  if you import a Python module that is linked against the Intel's OpenMP 
  library from the ``mkl`` module and PyHarm, which relies on a different 
  OpenMP implementation, you may encounter the error.

  If you do not insist on using the Intel's math library and are willing to 
  replace it by open source implementations, install the ``nomkl`` module via 
  ``conda`` using the ``conda-forge`` channel:

  .. code:: bash

     conda install -c conda-forge nomkl

  This may resolve the issue.

* *Solution 3*: Compile PyHarm from source (see :ref:`installing`) using an 
  OpenMP implementation that is compatible with your other modules.
