"""
Module defining some miscellaneous functions, constants, etc.

Note
----
All functions that deal with numerics are written in double precision.
"""


from . import _libcharm, _CHARM
from ._data_types import _ct_int
import ctypes as _ct


def get_version():
    """
    Returns a string specifying the CHarm version number determined on
    compilation time.
    """

    func          = _libcharm[_CHARM + 'misc_get_version']
    func.restype  = _ct.c_char_p
    func.argtypes = None

    return func().decode()


def print_info():
    """
    Prints library info (library name, version, compilation date, precision,
    etc).
    """

    func          = _libcharm[_CHARM + 'misc_print_info']
    func.restype  = None
    func.argtypes = None

    func()

    return


def buildopt_precision():
    """
    Returns:

    * ``1`` if CHarm was compiled in single precision
      (``--enable-single-precision``),

    * ``2`` if CHarm was compiled in double precision
      (``--enable-double-precision`` or no precision flag),

    * ``3`` if CHarm was compiled quadruple precision
      (``--enable-quad-precision``).
    """

    func          = _libcharm[_CHARM + 'misc_buildopt_precision']
    func.argtypes = None
    func.restype  = _ct_int

    return func()


def buildopt_omp_charm():
    """
    Returns a non-zero value if CHarm was compiled with the OpenMP
    parallelization enabled (``--enable-openmp``).  Otherwise, zero is
    returned.
    """

    func          = _libcharm[_CHARM + 'misc_buildopt_omp_charm']
    func.argtypes = None
    func.restype  = _ct_int

    return func()


def buildopt_omp_fftw():
    """
    Returns a non-zero value if the host's FFTW library supports OpenMP
    parallelization.  Otherwise, zero is returned.

    If non-zero, all FFTW computations are parallelized.
    """

    func          = _libcharm[_CHARM + 'misc_buildopt_omp_fftw']
    func.argtypes = None
    func.restype  = _ct_int

    return func()


def buildopt_simd():
    """
    Returns:

    * ``0`` if CHarm was compiled with SIMD instructions disabled,

    * ``1`` if CHarm was compiled with AVX instructions enabled
      (``--enable-avx``),

    * ``2`` if CHarm was compiled with AVX2 instructions enabled
      (``--enable-avx2``),

    * ``3`` if CHarm was compiled with AVX-512 instructions enabled
      (``--enable-avx-512``),

    * ``4`` if CHarm was compiled with NEON instructions enabled
      (``--enable-neon``).
    """

    func          = _libcharm[_CHARM + 'misc_buildopt_simd']
    func.argtypes = None
    func.restype  = _ct_int

    return func()


def buildopt_version_fftw():
    """
    Returns a string specifying the FFTW version used to compile CHarm.
    """

    func          = _libcharm[_CHARM + 'misc_buildopt_version_fftw']
    func.restype  = _ct.c_char_p
    func.argtypes = None

    return func().decode()


def buildopt_isfinite():
    """
    Returns a non-zero value if correctly working ``isfinite`` macro was found
    in the system's ``math.h`` header file before the compilation.  Otherwise,
    zero is returned.

    On some systems, ``isfinite`` is available in ``math.h``, but it is not
    working correctly with ``__float128`` floating point data type (quadruple
    precision).  The macro may also not work correctly with the ``-ffast-math``
    compiler flag (e.g., ``gcc`` and ``clang``).  In these cases, zero is
    returned.
    """

    func          = _libcharm[_CHARM + 'misc_buildopt_isfinite']
    func.argtypes = None
    func.restype  = _ct_int

    return func()

