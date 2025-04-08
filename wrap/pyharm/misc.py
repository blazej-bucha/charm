"""
Module defining some miscellaneous functions, constants, etc.

Note
----
All functions that deal with numerics are written in double precision.
"""


from . import _libcharm, _CHARM
from ._data_types import _ct_int
import ctypes as _ct
from .gfm import _WITH_MPFR


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
      (``--enable-neon``),

    * ``5`` if CHarm was compiled with SSE4.1 instructions enabled
      (``--enable-sse4.1``).
    """

    func          = _libcharm[_CHARM + 'misc_buildopt_simd']
    func.argtypes = None
    func.restype  = _ct_int

    return func()


def buildopt_simd_vector_size():
    """
    Returns the size of SIMD vectors if CHarm was compiled with SIMD
    instructions enabled and ``1`` otherwise (SIMD instructions disabled).
    """

    func          = _libcharm[_CHARM + 'misc_buildopt_simd_vector_size']
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


def buildopt_mpi():
    """
    Returns a non-zero value if CHarm was compiled with the MPI
    parallelization enabled (``--enable-mpi``).  Otherwise, zero is
    returned.
    """

    func          = _libcharm[_CHARM + 'misc_buildopt_mpi']
    func.argtypes = None
    func.restype  = _ct_int

    return func()


def buildopt_mpfr():
    """
    Returns a non-zero value if CHarm was compiled with the MPFR enabled
    (``--enable-mpfr``).  Otherwise, zero is returned.
    """

    func          = _libcharm[_CHARM + 'misc_buildopt_mpfr']
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


def buildopt_version_mpi():
    """
    If CHarm was compiled with the MPI parallelization enabled
    (``--enable-mpi``), returns the version and subversion numbers of the MPI
    standard that is being supported by the linked MPI implementation.  If MPI
    parallelization is disabled, all four values are set to ``-1``.

    Note
    ----

    All four version numbers returned are related to the MPI standard.  To see
    the version of the linked MPI implementation, call :meth:`print_info()`.

    Returns
    -------
    major_header : integer
        Version of the MPI standard determined on compile time
    minor_header : integer
        Subversion of the MPI standard determined on compile time
    major_lib : integer
        Version of the MPI standard determined on runtime
    minor_lib : integer
        Subversion of the MPI standard determined on runtine
    """

    func          = _libcharm[_CHARM + 'misc_buildopt_version_mpi']
    func.restype  = None
    func.argtypes = [_ct.POINTER(_ct_int),
                     _ct.POINTER(_ct_int),
                     _ct.POINTER(_ct_int),
                     _ct.POINTER(_ct_int)]

    major_h, minor_h = _ct_int(0), _ct_int(0)
    major_l, minor_l = _ct_int(0), _ct_int(0)
    func(_ct.pointer(major_h), _ct.pointer(minor_h),
         _ct.pointer(major_l), _ct.pointer(minor_l))
    return major_h.value, minor_h.value, major_l.value, minor_l.value


def buildopt_version_mpfr():
    """
    If CHarm was compiled with MPFR enabled (``--enable-mpfr``), returns the
    version of the MPFR library.  If CHarm was compiled with MPFR disabled,
    returns ``n/a``, ``-1``, ``-1``, ``-1``, respectively.

    Returns
    -------
    ver_str : str
        Version of the MPFR library determined on runtime
    ver_major : integer
        Symbolic constant ``MPFR_VERSION_MAJOR`` determined on compile time
    ver_minor : integer
        Symbolic constant ``MPFR_VERSION_MINOR`` determined on compile time
    ver_patch : integer
        Symbolic constant ``MPFR_VERSION_PATCHLEVEL`` determined on compile
        time
    """

    func          = _libcharm[_CHARM + 'misc_buildopt_version_mpfr']
    func.restype  = _ct.c_char_p
    func.argtypes = [_ct.POINTER(_ct_int),
                     _ct.POINTER(_ct_int),
                     _ct.POINTER(_ct_int)]

    major, minor, patch = _ct_int(0), _ct_int(0), _ct_int(0)
    ret = func(_ct.pointer(major), _ct.pointer(minor), _ct.pointer(patch))
    return ret.decode(), major.value, minor.value, patch.value


def buildopt_version_gmp():
    """
    The same as :meth:`buildopt_version_mpfr()` but for the GMP library
    (note that ``--enable-gmp`` is not a valid installation flag, though).
    """

    func          = _libcharm[_CHARM + 'misc_buildopt_version_gmp']
    func.restype  = _ct.c_char_p
    func.argtypes = [_ct.POINTER(_ct_int),
                     _ct.POINTER(_ct_int),
                     _ct.POINTER(_ct_int)]

    major, minor, patch = _ct_int(0), _ct_int(0), _ct_int(0)
    ret = func(_ct.pointer(major), _ct.pointer(minor), _ct.pointer(patch))
    return ret.decode(), major.value, minor.value, patch.value


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

