"""
Module to compute the following integrals:

* a product of two fully-normalized associated Legendre functions and a sine
  of a co-latitude,

  .. math::

      \int\limits_{\\theta_1}^{\\theta_2} \\bar{P}_{n_1,m_1}(\cos\\theta) \,
      \\bar{P}_{n_2,m_2}(\cos\\theta) \, \sin(\\theta) \,
      \mathrm{d}\\theta{,} \quad \\theta_1 \leq \\theta_2{,}


* a product of two ``4pi`` fully-normalized surface spherical harmonic
  functions over a rectangular cell on the unit sphere,

  .. math::

      \int\limits_{\\theta_1}^{\\theta_2}
      \int\limits_{\lambda_1}^{\lambda_2} \\bar{Y}_{n_1,m_1}(\\theta,
      \lambda) \, \\bar{Y}_{n_2,m_2}(\\theta, \lambda) \, \mathrm{d}\lambda
      \, \sin(\\theta) \, \mathrm{d}\\theta{,} \quad \\theta_1 \leq
      \\theta_2{,} \quad \lambda_1 \leq \lambda_2{.}


Note
----
This documentation is written for double precision version of PyHarm.
"""

import ctypes as _ct
import numpy as _np
from . import leg as _ph_leg
from . import _err as _ph_err
from . import _libcharm, _CHARM
from ._check_types import _check_deg_ord, _check_flt_scalar, _check_int_scalar
from ._data_types import _ct_flt, _ct_ulong


def pn1m1pn2m2(cltmin, cltmax, n1, m1, n2, m2, pnmj):
    """
    Analytically computes the integral

    .. math::

        \mathrm{ip} = \int\limits_{\\theta_{\mathrm{min}}}^{\\theta_{\mathrm{max}}}
        \\bar{P}_{n_1, m_1}(\cos\\theta) \, \\bar{P}_{n_2, m_2}(\cos\\theta) \,
        \sin\\theta \, \mathrm{d}\\theta{,} \quad \\theta_{\mathrm{min}} \leq
        \\theta_{\mathrm{max}}{.}

    The computation is based on the Fourier coefficients of the associated
    Legendre functions (see Eq. 33 of Pail and Plank, 2001).

    **References**:

    * Pail, R., Plank, G., Schuh, W.-D. (2001) Spatially restricted data
      distributions on the sphere: the method of orthonormalized functions and
      applications. Journal of Geodesy 75:44--56

    Parameters
    ----------
    cltmin : floating point
        Minimum co-latitude in radians
    cltmax : floating point
        Maximum co-latitude in radians
    n1 : integer
        Harmonic degree of the first Legendre function
    m1 : integer
        Harmonic order of the first Legendre function
    n2 : integer
        Harmonic degree of the second Legendre function
    m2 : integer
        Harmonic order of the second Legendre function
    pnmj : Pnmj
        Fourier coefficients of Legendre functions up to degree ``max(n1,
        n2)``.

    Returns
    -------
    out : floating point
        The output value ``ip`` of the integral.
    """

    _check_pn1m1pn2m2_inputs(cltmin, cltmax, n1, m1, n2, m2, pnmj)

    func          = _libcharm[_CHARM + 'integ_pn1m1pn2m2']
    func.restype  = _ct_flt
    func.argtypes = [_ct_flt,
                     _ct_flt,
                     _ct_ulong,
                     _ct_ulong,
                     _ct_ulong,
                     _ct_ulong,
                     _ct.POINTER(_ph_leg._Pnmj),
                     _ct.POINTER(_ph_err._Err)]

    err = _ph_err.init()
    ip = func(_ct_flt(cltmin), _ct_flt(cltmax), _ct_ulong(n1), _ct_ulong(m1),
              _ct_ulong(n2), _ct_ulong(m2), pnmj._Pnmj, err)
    _ph_err.handler(err, 1)
    _ph_err.free(err)

    return ip


def yi1n1m1yi2n2m2(cltmin, cltmax, lonmin, lonmax, i1, n1, m1, i2, n2, m2,
                   pnmj):
    """
    Analytically computes the integral

    .. math::

        \mathrm{iy} =
             \int\limits_{\\theta_{\mathrm{min}}}^{\\theta_{\mathrm{max}}}
             \int\limits_{\lambda_{\mathrm{min}}}^{\lambda_{\mathrm{max}}}
             \\bar{Y}_{i_1,n_1,m_1}(\\theta, \lambda) \,
             \\bar{Y}_{i_2,n_2,m_2}(\\theta, \lambda) \, \mathrm{d}
             \lambda \, \sin(\\theta) \, \mathrm{d}\\theta{,} \quad
             \\theta_{\mathrm{min}} \leq \\theta_{\mathrm{max}}{,} \quad
             \lambda_{\mathrm{min}} \leq \lambda_{\mathrm{max}}{,}

    where

    .. math::

        \\bar{Y}_{i,n,m}(\\theta, \lambda) = \\begin{cases}
        \\bar{P}_{nm}(\cos\\theta) \, \cos(m \, \lambda) \quad &\\textrm{if}
        \quad i = 0{,}\\\\ \\bar{P}_{nm}(\cos\\theta) \, \sin(m \, \lambda)
        \quad &\\textrm{if} \quad i = 1{.} \end{cases}

    Parameters
    ----------
    cltmin : floating point
        Minimum co-latitude in radians
    cltmax : floating point
        Maximum co-latitude in radians
    lonmin : floating point
        Minimum longitude in radians
    lonmax : floating point
        Maximum longitude in radians
    i1 : integer
        ``0`` if the first spherical harmonic function is of the ``cos``
        type, ``1`` for the ``sin`` type
    n1 : integer
        Harmonic degree of the first spherical harmonic function
    m1 : integer
        Harmonic order of the first spherical harmonic function
    i2 : integer
        The same as ``i1`` but for the second spherical harmonic function
    n2 : integer
        Harmonic degree of the second spherical harmonic function
    m2 : integer
        Harmonic order of the second spherical harmonic function
    pnmj : Pnmj
        Fourier coefficients of Legendre functions at least up to degree
        ``max(n1, n2)``.

    Returns
    -------
    out : floating point
        The value ``iy`` of the integral.
    """

    _check_pn1m1pn2m2_inputs(cltmin, cltmax, n1, m1, n2, m2, pnmj)
    _check_flt_scalar(lonmin, 'The minimum longitude')
    _check_flt_scalar(lonmax, 'The maximum longitude')
    _check_int_scalar(i1, '\'i1\'')
    _check_int_scalar(i2, '\'i2\'')
    if i1 != 0 and i1 != 1:
        raise ValueError('\'i1\' must be either \'0\' for spherical harmonic '
                         'with the \'cos\' term or \'1\' for the '
                         '\'sin\' type of spherical harmonic.')
    if i2 != 0 and i2 != 1:
        raise ValueError('\'i2\' must be either \'0\' for spherical harmonic '
                         'with the \'cos\' term or \'1\' for the '
                         '\'sin\' type of spherical harmonic.')

    func          = _libcharm[_CHARM + 'integ_yi1n1m1yi2n2m2']
    func.restype  = _ct_flt
    func.argtypes = [_ct_flt,
                     _ct_flt,
                     _ct_flt,
                     _ct_flt,
                     _ct.c_bool,
                     _ct_ulong,
                     _ct_ulong,
                     _ct.c_bool,
                     _ct_ulong,
                     _ct_ulong,
                     _ct.POINTER(_ph_leg._Pnmj),
                     _ct.POINTER(_ph_err._Err)]

    err = _ph_err.init()
    iy = func(_ct_flt(cltmin), _ct_flt(cltmax), _ct_flt(lonmin),
              _ct_flt(lonmax), _ct.c_bool(i1), _ct_ulong(n1), _ct_ulong(m1),
              _ct.c_bool(i2), _ct_ulong(n2), _ct_ulong(m2), pnmj._Pnmj, err)
    _ph_err.handler(err, 1)
    _ph_err.free(err)

    return iy


def _check_pn1m1pn2m2_inputs(cltmin, cltmax, n1, m1, n2, m2, pnmj):
    """
    A private function to check inputs to "pn1m1pn2m2" and some input to
    "yi1n1m1yi2n2m2".

    Parameters
    ----------
    cltmin : floating point
        Minimum co-latitude in radians
    cltmax : floating point
        Maximum co-latitude in radians
    n1 : integer
        Harmonic degree of the first spherical harmonic function
    m1 : integer
        Harmonic order of the first spherical harmonic function
    n2 : integer
        Harmonic degree of the second spherical harmonic function
    m2 : integer
        Harmonic order of the second spherical harmonic function
    pnmj : Pnmj
        A `pyharm.leg.Pnmj` structure with the Fourier coefficients of
        associated Legendre functions at least up to degree ``max(n1, n2)``.
        It is assumed that the instance is prepared beforehand.  Note that
        ``j`` is related to wave-numbers, but is *not* a wave-number (refer
        to `charm_leg <./api-c-leg.html>`_ for the full documentation).
    """

    _check_flt_scalar(cltmin, 'The minimum co-latitude')
    _check_flt_scalar(cltmax, 'The maximum co-latitude')
    _check_deg_ord(n1, 'degree')
    _check_deg_ord(m1, 'order')
    _check_deg_ord(n2, 'degree')
    _check_deg_ord(m2, 'order')
    if not isinstance(pnmj, _ph_leg.Pnmj):
        raise TypeError('\'Pnmj\' must be an instance of the '
                        '\'pyharm.leg.Pnmj\' class.')

    if m1 > n1:
        raise ValueError(f'\'m1 = {m1}\' cannot be larger than \'n1 = {n1}\'.')

    if m2 > n2:
        raise ValueError(f'\'m2 = {m2}\' cannot be larger than \'n2 = {n2}\'.')

    if max(n1, n2) > pnmj.nmax:
        raise ValueError('The \'Pnmj\' instance is initialize up to degree '
                         '%d, but must be initialized at least up to degree '
                         '\'max(n1, n2) = %d\'.' % (pnmj.nmax, max(n1, n2)))

