"""
Module to perform surface spherical harmonic analysis using:

    * point data values,

    * mean data values.

Note
----
This documentation is written for double precision version of PyHarm.
"""

import ctypes as _ct
import numpy as _np
from . import _libcharm, _CHARM
from . import crd as _ph_crd
from . import shc as _ph_shc
from . import _err as _ph_err
from ._get_module_constants import _get_module_constants
from ._data_types import _ct_ulong, _ct_flt, _ct_int
from ._check_types import _check_flt_ndarray, _check_deg_ord, \
                          _check_int_scalar, _check_flt_scalar, _check_radius
from .shc import _MU, _R


# Get the module constants from "_constants.py" and add them to the module's
# namespace
_get_module_constants('CHARM_SHA_')


#: int: Spherical harmonic analysis of cell data using the `approximate`
#: quadrature method.
CELL_AQ: int


def point(pnt, f, nmax, mu=_MU, r=_R):
    """
    Performs surface spherical harmonic analysis of point values ``f`` at
    ``pnt`` up to maximum degree ``nmax``.  Refer to `charm_sha
    <./api-c-sha.html>`_ for the full documentation.

    Parameters
    ----------
    pnt : PointGridDH1, PointGridDH2, PointGridGL
        Quadrature grid points at which ``f`` is sampled
    f : numpy array of floating points
        Signal to be harmonically analysed
    nmax : integer
        Maximum degree of the analysis
    mu : floating point
        Scaling parameter to be associated with the output spherical harmonic
        coefficients, optional.  Default is ``1.0``.
    r : floating point
        Radius of the reference sphere to be associated with the output
        spherical harmonic coefficients, optional.  Default is ``1.0``.

    Returns
    -------
    out : Shc
        Spherical harmonic coefficients of ``f``
    """

    if not isinstance(pnt, _ph_crd.PointGridGL) and \
        not isinstance(pnt, _ph_crd.PointGridDH1) and \
        not isinstance(pnt, _ph_crd.PointGridDH2):
        msg  = f'\'pnt\' must be an instance of one of the following classes: '
        msg += f'{_ph_crd.PointGridGL}, {_ph_crd.PointGridDH1}, '
        msg += f'{_ph_crd.PointGridDH2}.'
        raise TypeError(msg)

    _check_flt_ndarray(f, 2, 'The \'f\' variable')
    _check_deg_ord(nmax, 'degree')
    _check_flt_scalar(mu, 'Scaling parameter')
    _check_radius(r)

    if f.shape[0] != pnt.lat.size:
        msg  = f'The number of elements in \'f\' along the axis \'0\' does '
        msg += f'not match the number of latitudes in \'pnt\' '
        msg += f'({f.shape[0]} vs. {pnt.lat.size}).'
        raise ValueError(msg)

    if f.shape[1] != pnt.lon.size:
        msg  = f'The number of elements in \'f\' along the axis \'1\' does '
        msg += f'not match the number of latitudes in \'pnt\' '
        msg += f'({f.shape[1]} vs. {pnt.lon.size}).'
        raise ValueError(msg)

    func          = _libcharm[_CHARM + 'sha_point']
    func.restype  = None
    func.argtypes = [_ct.POINTER(_ph_crd._Point),
                     _ct.POINTER(_ct_flt),
                     _ct_ulong,
                     _ct.POINTER(_ph_shc._Shc),
                     _ct.POINTER(_ph_err._Err)]

    shcs = _ph_shc.Shc.from_garbage(nmax, mu, r)
    err = _ph_err.init()
    func(pnt._Point, f.ctypes.data_as(_ct.POINTER(_ct_flt)), _ct_ulong(nmax),
         shcs._Shc, err)
    _ph_err.handler(err, 1)
    _ph_err.free(err)

    return shcs


def cell(cell, f, nmax, method, mu=_MU, r=_R):
    """
    Performs surface spherical harmonic analysis of point values ``f`` at
    ``cell`` up to maximum degree ``nmax``.  Refer to `charm_sha
    <./api-c-sha.html>`_ for the full documentation.

    Parameters
    ----------
    cell : CellGrid
        Grid cells at which ``f`` is sampled
    f : numpy array of floating points
        Signal to be harmonically analysed
    nmax : integer
        Maximum degree of the analysis
    method : integer
        Method of the spherical harmonic analysis of area-mean values.
        Currently the only accepted value is :obj:`pyharm.sha.CELL_AQ`.
    mu : floating point
        Scaling parameter to be associated with the output spherical harmonic
        coefficients, optional.  Default is ``1.0``.
    r : floating point
        Radius of the reference sphere to be associated with the output
        spherical harmonic coefficients, optional.  Default is ``1.0``.

    Returns
    -------
    out : Shc
        Spherical harmonic coefficients of ``f``
    """

    if not isinstance(cell, _ph_crd.CellGrid):
        raise TypeError('\'cell\' must be a \'CellGrid\' class instance.')

    _check_flt_ndarray(f, 2, 'The \'f\' variable')
    _check_deg_ord(nmax, 'degree')
    _check_flt_scalar(mu, 'Scaling parameter')
    _check_radius(r)
    _check_int_scalar(method, 'The \'method\' variable')

    if method != CELL_AQ:
        raise ValueError('Unsupported value of \'method\'.')

    if f.shape[0] != cell.latmin.size:
        msg  = f'The number of elements in \'f\' along the axis \'0\' does '
        msg += f'not match the number of cells in the latitudinal direction '
        msg += f'in \'cell\' ({f.shape[0]} vs. {cell.latmin.size}).'
        raise ValueError(msg)

    if f.shape[1] != cell.lonmin.size:
        msg  = f'The number of elements in \'f\' along the axis \'1\' does '
        msg += f'not match the number of cells in the longitudinal direction '
        msg += f'in \'cell\' ({f.shape[1]} vs. {cell.lonmin.size}).'
        raise ValueError(msg)

    func          = _libcharm[_CHARM + 'sha_cell']
    func.restype  = None
    func.argtypes = [_ct.POINTER(_ph_crd._Cell),
                     _ct.POINTER(_ct_flt),
                     _ct_ulong,
                     _ct_int,
                     _ct.POINTER(_ph_shc._Shc),
                     _ct.POINTER(_ph_err._Err)]

    shcs = _ph_shc.Shc.from_garbage(nmax, mu, r)
    err = _ph_err.init()
    func(cell._Cell, f.ctypes.data_as(_ct.POINTER(_ct_flt)), _ct_ulong(nmax),
         _ct_int(method), shcs._Shc, err)
    _ph_err.handler(err, 1)
    _ph_err.free(err)

    return shcs

