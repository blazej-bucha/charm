"""
Module to perform solid spherical harmonic synthesis of point and mean values.

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
from ._check_types import _check_deg_ord
from ._data_types import _ct_ulong, _ct_flt, _pyharm_flt


def point(pnt, shcs, nmax):
    """
    Performs spherical harmonic synthesis of point values from :obj:`shcs` at
    :obj:`pnt` up to maximum degree :obj:`nmax`.  Refer to
    `https://blazej-bucha.github.io/charm/build/html/api-c-shs.html
    <https://blazej-bucha.github.io/charm/build/html/api-c-shs.html>`_ for the 
    full documentation.

    Parameters
    ----------
    pnt : PointGrid, PointGridDH1, PointGridDH2, PointGridGL or PointSctr
        Evaluation points
    shcs : Shc
        Spherical harmonic coefficients
    nmax : integer
        Maximum degree of the synthesis

    Returns
    -------
    out : numpy array of floating points
        Point values synthesized from :obj:`shcs` at :obj:`pnt`
    """

    if not isinstance(pnt, _ph_crd._PointBase):
        msg  = f'\'pnt\' must be an instance of one of following classes:'
        for i, cl in enumerate(_ph_crd._PointBase.__subclasses__()):
            msg += f' {_ph_crd._PointBase.__subclasses__()[i]}'
            if i < len(_ph_crd._PointBase.__subclasses__()) - 1:
                msg += f','
        msg += f'.'
        raise TypeError(msg)

    if not isinstance(shcs, _ph_shc.Shc):
        raise TypeError(f'\'shcs\' must be a \'{_ph_shc.Shc}\' '
                        f'class instance.')

    _check_deg_ord(nmax, 'degree')

    if nmax > shcs.nmax:
        msg  = f'The maximum degree of the synthesis \'nmax = {nmax}\' '
        msg += f'cannot be larger than the maximum degree of the '
        msg += f'coefficients \'shcs.nmax = {shcs.nmax}\'.'
        raise ValueError(msg)

    func          = _libcharm[_CHARM + 'shs_point']
    func.restype  = None
    func.argtypes = [_ct.POINTER(_ph_crd._Point),
                     _ct.POINTER(_ph_shc._Shc),
                     _ct_ulong,
                     _ct.POINTER(_ct_flt),
                     _ct.POINTER(_ph_err._Err)]

    if isinstance(pnt, _ph_crd.PointSctr):
        f_shape = (pnt.lat.size,)
    else:
        f_shape = (pnt.lat.size, pnt.lon.size)
    f = _np.zeros(f_shape, dtype=_pyharm_flt, order='C')

    err = _ph_err.init()
    func(pnt._Point, shcs._Shc, _ct_ulong(nmax),
         f.ctypes.data_as(_ct.POINTER(_ct_flt)), err)
    _ph_err.handler(err, 1)
    _ph_err.free(err)

    return f


def cell(cell, shcs, nmax):
    """
    Performs spherical harmonic synthesis of area-mean values from :obj:`shcs`
    at :obj:`cell` up to maximum degree :obj:`nmax`.  Refer to
    `https://blazej-bucha.github.io/charm/build/html/api-c-shs.html
    <https://blazej-bucha.github.io/charm/build/html/api-c-shs.html>`_ for the 
    full documentation.

    Parameters
    ----------
    cell : CellSctr, CellGrid
        Evaluation cells
    shcs : Shc
        Spherical harmonic coefficients
    nmax : integer
        Maximum degree of the synthesis

    Returns
    -------
    out : numpy array of floating points
        Area-mean values synthesized from :obj:`shcs` at :obj:`cell`
    """

    if not isinstance(cell, _ph_crd._CellBase):
        msg  = f'\'cell\' must be an instance of one of the following '
        msg += f'classes: {_ph_crd._CellBase.__subclasses__()}.'
        raise TypeError(msg)

    if not isinstance(shcs, _ph_shc.Shc):
        raise TypeError(f'\'shcs\' must be a \'{_ph_shc.Shc}\' '
                        f'class instance.')

    _check_deg_ord(nmax, 'degree')

    if nmax > shcs.nmax:
        msg  = f'The maximum degree of the synthesis \'nmax = {nmax}\' '
        msg += f'cannot be larger than the maximum degree of the '
        msg += f'coefficients \'shcs.nmax = {shcs.nmax}\'.'
        raise ValueError(msg)

    func          = _libcharm[_CHARM + 'shs_cell']
    func.restype  = None
    func.argtypes = [_ct.POINTER(_ph_crd._Cell),
                     _ct.POINTER(_ph_shc._Shc),
                     _ct_ulong,
                     _ct.POINTER(_ct_flt),
                     _ct.POINTER(_ph_err._Err)]

    if isinstance(cell, _ph_crd.CellSctr):
        f_shape = (cell.latmin.size,)
    elif isinstance(cell, _ph_crd.CellGrid):
        f_shape = (cell.latmin.size, cell.lonmin.size)
    f = _np.zeros(f_shape, dtype=_pyharm_flt, order='C')

    err = _ph_err.init()
    func(cell._Cell, shcs._Shc, _ct_ulong(nmax),
         f.ctypes.data_as(_ct.POINTER(_ct_flt)), err)
    _ph_err.handler(err, 1)
    _ph_err.free(err)

    return f


def cell_isurf(cell, shcs1, nmax1, shcs2, nmax2, nmax3, nmax4):
    """
    Performs spherical harmonic synthesis of area-mean values from :obj:`shcs1` 
    at :obj:`cell` residing on an irregular surface defined by :obj:`shcs2`.  
    The synthesis of area-mean values is done up to degree :obj:`nmax1` 
    and the irregular surface is expanded up to degree :obj:`nmax2`.  
    :obj:`nmax3` and :obj:`nmax4` represent the maximum harmonic degrees to 
    synthesize and analyze the :obj:`(shcs1.r / r)^(n + 1)` terms, where 
    :obj:`r` stands for the spherical radius of the irregular surface defined 
    by :obj:`shcs2`.  Refer to 
    `https://blazej-bucha.github.io/charm/build/html/api-c-shs.html 
    <https://blazej-bucha.github.io/charm/build/html/api-c-shs.html>`_ for the 
    full documentation.

    Parameters
    ----------
    cell : CellGrid
        Evaluation cells
    shcs1 : Shc
        Spherical harmonic coefficients of the function, the area-mean values
        of which are synthesized
    nmax1 : integer
        Maximum degree of the synthesis of the area-mean values
    shcs2 : Shc
        Spherical harmonic coefficients of the irregular surface, on which the
        area-mean values are synthesized
    nmax2 : integer
        Maximum degree of the synthesis of the irregular surface
    nmax3 : integer
        Maximum degree of the synthesis of :obj:`(shcs1.r / r)^(n + 1)`
    nmax4 : integer
        Maximum degree of the analysis of :obj:`(shcs1.r / r)^(n + 1)`

    Returns
    -------
    out : numpy array of floating points
        Area-mean values synthesized from :obj:`shcs1` at :obj:`cell` residing
        on the surface defined by :obj:`shcs2`
    """

    if not isinstance(cell, _ph_crd.CellGrid):
        msg  = f'\'cell\' must be a {_ph_crd.CellGrid} class instance.'
        raise TypeError(msg)

    if not isinstance(shcs1, _ph_shc.Shc):
        raise TypeError(f'\'shcs1\' must be a \'{_ph_shc.Shc}\' '
                        f'class instance.')

    if not isinstance(shcs2, _ph_shc.Shc):
        raise TypeError(f'\'shcs2\' must be a \'{_ph_shc.Shc}\' '
                        f'class instance.')

    _check_deg_ord(nmax1, 'degree')
    _check_deg_ord(nmax2, 'degree')
    _check_deg_ord(nmax3, 'degree')
    _check_deg_ord(nmax4, 'degree')

    if nmax1 > shcs1.nmax:
        msg  = f'The maximum degree of the synthesis \'nmax1 = {nmax1}\' '
        msg += f'cannot be larger than the maximum degree of the '
        msg += f'coefficients \'shcs1.nmax = {shcs1.nmax}\'.'
        raise ValueError(msg)

    if nmax2 > shcs2.nmax:
        msg  = f'The maximum degree of the synthesis \'nmax2 = {nmax2}\' '
        msg += f'cannot be larger than the maximum degree of the '
        msg += f'coefficients \'shcs2.nmax = {shcs2.nmax}\'.'
        raise ValueError(msg)

    if nmax3 > nmax4:
        raise ValueError('\'nmax3 = %d\' cannot be larger than '
                         '\'nmax4 = %d\'.' % (nmax3, nmax4))

    func          = _libcharm[_CHARM + 'shs_cell_isurf']
    func.restype  = None
    func.argtypes = [_ct.POINTER(_ph_crd._Cell),
                     _ct.POINTER(_ph_shc._Shc),
                     _ct_ulong,
                     _ct.POINTER(_ph_shc._Shc),
                     _ct_ulong,
                     _ct_ulong,
                     _ct_ulong,
                     _ct.POINTER(_ct_flt),
                     _ct.POINTER(_ph_err._Err)]

    f       = _np.zeros((cell.latmin.size, cell.lonmin.size),
                        dtype=_pyharm_flt, order='C')

    err = _ph_err.init()
    func(cell._Cell, shcs1._Shc, _ct_ulong(nmax1), shcs2._Shc,
         _ct_ulong(nmax2), _ct_ulong(nmax3), _ct_ulong(nmax4),
         f.ctypes.data_as(_ct.POINTER(_ct_flt)), err)
    _ph_err.handler(err, 1)
    _ph_err.free(err)

    return f

