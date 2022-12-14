"""
Module to work with X-numbers after Fukushima (2012).

Note
----
This documentation is written for double precision version of CHarm.


**References:**

* Fukushima, T. (2012) Numerical computation of spherical harmonics of 
  arbitrary degree and order by extending exponent of floating point 
  numbers. Journal of Geodesy 86:271-285.  """

import ctypes as _ct
import numpy as _np
from . import _libcharm, _CHARM
from ._check_types import _check_flt_scalar, _check_int_scalar
from ._data_types import _ct_flt, _ct_int, _pyharm_flt


def xnorm(x, ix):
    """
    Weakly normalizes an X-number. The function is due to Fukushima (2012),
    Table 7.

    Parameters
    ----------
    x : floating point
        Significand of the input X-number
    ix : integer
        Exponent of the input X-number

    Returns
    -------
    x : floating point
        Significand of the output X-number
    ix : integer
        Exponent of the output X-number
    """

    _check_flt_scalar(x,   'x')
    _check_int_scalar(ix, 'ix')

    func          = _libcharm[_CHARM + 'xnum_xnorm']
    func.restype  = None
    func.argtypes = [_ct.POINTER(_ct_flt),
                     _ct.POINTER(_ct_int)]

    x_array  = _np.array([x], dtype=_pyharm_flt)
    ix_array = _np.array([ix], dtype=_np.integer)

    func(x_array.ctypes.data_as(_ct.POINTER(_ct_flt)),
         ix_array.ctypes.data_as(_ct.POINTER(_ct_int)))

    return x_array[0], ix_array[0]


def xlsum2(f, x, g, y, ix, iy):
    """
    Computes a two-term linear sum of X-numbers with F-number coefficients.
    The function is due to Fukushima (2012), Table 8.

    Parameters
    ----------
    f : floating point
        F-number coefficient of the first linear term
    x : floating point
        Significand of the first X-number term
    g : floating point
        F-number coefficient of the second linear term
    y : floating point
        Significand of the second X-number term
    ix : integer
        Exponent of the first X-number term
    iy : integer
        Exponent of the second X-number term

    Returns
    -------
    z : floating point
        Pointer to the significand of the output X-number
    iz : floating point
        Pointer to the exponent of the output X-number
    """

    _check_flt_scalar(f,   'f')
    _check_flt_scalar(x,   'x')
    _check_flt_scalar(g,   'g')
    _check_flt_scalar(y,   'y')
    _check_int_scalar(ix, 'ix')
    _check_int_scalar(iy, 'iy')

    func          = _libcharm[_CHARM + 'xnum_xlsum2']
    func.restype  = None
    func.argtypes = [_ct_flt,
                     _ct_flt,
                     _ct_flt,
                     _ct_flt,
                     _ct.POINTER(_ct_flt),
                     _ct_int,
                     _ct_int,
                     _ct.POINTER(_ct_int)]

    z_array  = _np.array([0.0], dtype=_pyharm_flt)
    iz_array = _np.array([0], dtype=_np.integer)

    func(_ct_flt(f), _ct_flt(x), _ct_flt(g), _ct_flt(y),
         z_array.ctypes.data_as(_ct.POINTER(_ct_flt)), _ct_int(ix),
         _ct_int(iy), iz_array.ctypes.data_as(_ct.POINTER(_ct_int)))

    return z_array[0], iz_array[0]


def x2f(x, ix):
    """
    Translates an X-number into an F-number. The function is due to Fukushima
    (2012), Table 6.

    Parameters
    ----------
    x : floating point
        Significand of the X-number
    ix : integer
        Exponent of the X-number

    Returns
    -------
    out : floating point
        The F-number representation of an input X-number
    """

    _check_flt_scalar(x,   'x')
    _check_int_scalar(ix, 'ix')

    func          = _libcharm[_CHARM + 'xnum_x2f']
    func.restype  = _ct_flt
    func.argtypes = [_ct_flt,
                     _ct_int]

    return func(_ct_flt(x), _ct_int(ix))

