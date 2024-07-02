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
from ._get_module_constants import _get_module_constants
from ._check_types import _check_deg_ord, _check_int_scalar
from ._data_types import _ct_uint, _ct_ulong, _ct_flt, _pyharm_flt


# Get the module constants from "_constants.py" and add them to the module's
# namespace.  Here, we intentional use "GRAD" instead of "GRAD_", so that the
# constants imported will have the "_" suffix, meaning they will be hidden for
# users.
_get_module_constants('GRAD')


def point(pnt, shcs, nmax):
    """
    Performs spherical harmonic synthesis of point values from ``shcs`` at
    ``pnt`` up to maximum degree ``nmax``.  Refer to `charm_shs
    <./api-c-shs.html>`_ for the full documentation.

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
        Point values synthesized from ``shcs`` at ``pnt``
    """

    _check_crd_instance(pnt, _ph_crd._PointBase, 'pnt')
    _check_shcs(shcs, 'shcs')
    _check_deg_ord(nmax, 'degree')
    _check_nmax(nmax, shcs.nmax, 'nmax', 'shcs.nmax')

    func          = _libcharm[_CHARM + 'shs_point']
    func.restype  = None
    func.argtypes = [_ct.POINTER(_ph_crd._Point),
                     _ct.POINTER(_ph_shc._Shc),
                     _ct_ulong,
                     _ct.POINTER(_ct_flt),
                     _ct.POINTER(_ph_err._Err)]

    if isinstance(pnt, _ph_crd.PointSctr):
        f_shape = (pnt.npoint,)
    else:
        f_shape = (pnt.lat.size, pnt.lon.size)
    f = _np.zeros(f_shape, dtype=_pyharm_flt, order='C')

    err = _ph_err.init()
    func(pnt._Point, shcs._Shc, _ct_ulong(nmax),
         f.ctypes.data_as(_ct.POINTER(_ct_flt)), err)
    _ph_err.handler(err, 1)
    _ph_err.free(err)

    return f


def point_grad1(pnt, shcs, nmax):
    """
    Performs the synthesis of point values of the first-order gradient in `LNOF
    <./definitions.html#lnof>`_ from ``shcs`` at ``pnt`` up to maximum
    degree ``nmax``.  Refer to `charm_shs <./api-c-shs.html>`_ for the full
    documentation.

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
    fx : numpy array of floating points
        The `x` element
    fy : numpy array of floating points
        The `y` element
    fz : numpy array of floating points
        The `z` element
    """

    return _point_gradn(pnt, shcs, nmax, 1)


def point_grad2(pnt, shcs, nmax):
    """
    Performs the synthesis of point values of the second-order gradient in
    `LNOF <./definitions.html#lnof>`_ from ``shcs`` at ``pnt`` up to
    maximum degree ``nmax``.  Refer to `charm_shs <./api-c-shs.html>`_ for
    the full documentation.

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
    fxx : numpy array of floating points
        The `xx` element
    fxy : numpy array of floating points
        The `xy` element
    fxz : numpy array of floating points
        The `xz` element
    fyy : numpy array of floating points
        The `yy` element
    fyz : numpy array of floating points
        The `yz` element
    fzz : numpy array of floating points
        The `zz` element
    """

    return _point_gradn(pnt, shcs, nmax, 2)


def _point_gradn(pnt, shcs, nmax, gradn):

    _check_crd_instance(pnt, _ph_crd._PointBase, 'pnt')
    _check_shcs(shcs, 'shcs')
    _check_deg_ord(nmax, 'degree')
    _check_nmax(nmax, shcs.nmax, 'nmax', 'shcs.nmax')

    if gradn == 1:
        nelem = 3
    elif gradn == 2:
        nelem = 6

    func          = _libcharm[_CHARM + 'shs_point_grad' + f'{gradn}']
    func.restype  = None
    func.argtypes = [_ct.POINTER(_ph_crd._Point),
                     _ct.POINTER(_ph_shc._Shc),
                     _ct_ulong,
                     _ct.POINTER(_ct.POINTER(_ct_flt) * nelem),
                     _ct.POINTER(_ph_err._Err)]

    fptr = (_ct.POINTER(_ct_flt) * nelem)(None)
    if isinstance(pnt, _ph_crd.PointSctr):
        f_shape = (pnt.npoint,)
    else:
        f_shape = (pnt.lat.size, pnt.lon.size)
    f = nelem * ['']
    for i in range(nelem):
        f[i] = _np.zeros(f_shape, dtype=_pyharm_flt, order='C')
        fptr[i] = f[i].ctypes.data_as(_ct.POINTER(_ct_flt))

    err = _ph_err.init()
    func(pnt._Point, shcs._Shc, _ct_ulong(nmax),
         _ct.pointer(fptr), err)

    _ph_err.handler(err, 1)
    _ph_err.free(err)

    # This is tricky.  In "fptr", the pointers are swapped after leaving the
    # CHarm function "func", but this of course does not affect the numpy
    # arrays in "f".  This means that, for instance, "fptr[2][0] == f[1][0,
    # 0]".  The only safe way of achieving nice ordering of vector and tensor
    # elements on output (as in "fptr") I could come up with is to change the
    # order of elements in the "f" list on return.  Importantly, such obtained
    # numpy arrays own their memory, which I was not able to achieve in any
    # other way.
    if gradn == 1:
        return f[_P], f[_L], f[_R]
    elif gradn == 2:
        return f[_PP], f[_LP], f[_RP], f[_LL], f[_LR], f[_RR]


def point_guru(pnt, shcs, nmax, dr, dlat, dlon):
    """
    Performs the synthesis of point values of

    .. math::

        \\frac{1}{r^{j + k} \, \cos^k\\varphi} \, \\frac{\partial^{i + j + k}
        f}{\partial r^i \, \partial \\varphi^j \, \lambda^k}

    for :math:`i = 0, 1, 2` (``dr``), :math:`j = 0, 1, 2` (``dlat``) and
    :math:`k = 0, 1, 2` (``dlon``) satisfying :math:`i + j + k \leq 2` from
    ``shcs`` at ``pnt`` up to maximum degree ``nmax``.  Refer to
    `charm_shs_point_guru` in `charm_shs <./api-c-shs.html>`_ for the full
    documentation.

    Parameters
    ----------
    pnt : PointGrid, PointGridDH1, PointGridDH2, PointGridGL or PointSctr
        Evaluation points
    shcs : Shc
        Spherical harmonic coefficients
    nmax : integer
        Maximum degree of the synthesis
    dr : integer
        Order of the radial derivative (variable :math:`i` in the equation
        above)
    dlat : integer
        Order of the latitudinal derivative (variable :math:`j` in the equation
        above)
    dlon : integer
        Order of the longitudinal derivative (variable :math:`k` in the
        equation above)

    Returns
    -------
    f : numpy array of floating points
        Output quantity depending on ``dr``, ``dlat`` and ``dlon``
    """

    _check_crd_instance(pnt, _ph_crd._PointBase, 'pnt')
    _check_shcs(shcs, 'shcs')
    _check_deg_ord(nmax, 'degree')
    _check_nmax(nmax, shcs.nmax, 'nmax', 'shcs.nmax')
    _check_int_scalar(dr, '\'dr\'')
    _check_int_scalar(dlat, '\'dlat\'')
    _check_int_scalar(dlon, '\'dlon\'')

    def _check_derivative(d, var):

        if d < 0:
            msg  = f'The derivative order \'{var}\' is \'{d}\', but '
            msg += f'it cannot be negative.'
            raise ValueError(msg)

        return

    _check_derivative(dr, 'dr')
    _check_derivative(dlat, 'dlat')
    _check_derivative(dlon, 'dlon')

    func          = _libcharm[_CHARM + 'shs_point_guru']
    func.restype  = None
    func.argtypes = [_ct.POINTER(_ph_crd._Point),
                     _ct.POINTER(_ph_shc._Shc),
                     _ct_ulong,
                     _ct_uint,
                     _ct_uint,
                     _ct_uint,
                     _ct.POINTER(_ct_flt),
                     _ct.POINTER(_ph_err._Err)]

    if isinstance(pnt, _ph_crd.PointSctr):
        f_shape = (pnt.npoint,)
    else:
        f_shape = (pnt.lat.size, pnt.lon.size)
    f = _np.zeros(f_shape, dtype=_pyharm_flt, order='C')

    err = _ph_err.init()
    func(pnt._Point,
         shcs._Shc,
         _ct_ulong(nmax),
         _ct_uint(dr),
         _ct_uint(dlat),
         _ct_uint(dlon),
         f.ctypes.data_as(_ct.POINTER(_ct_flt)),
         err)

    _ph_err.handler(err, 1)
    _ph_err.free(err)

    return f


def cell(cell, shcs, nmax):
    """
    Performs spherical harmonic synthesis of area-mean values from ``shcs``
    at ``cell`` up to maximum degree ``nmax``.  Refer to `charm_shs
    <./api-c-shs.html>`_ for the full documentation.

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
        Area-mean values synthesized from ``shcs`` at ``cell``
    """

    _check_crd_instance(cell, _ph_crd._CellBase, 'cell')
    _check_shcs(shcs, 'shcs')
    _check_deg_ord(nmax, 'degree')
    _check_nmax(nmax, shcs.nmax, 'nmax', 'shcs.nmax')

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
    Performs spherical harmonic synthesis of area-mean values from ``shcs1``
    at ``cell`` residing on an irregular surface defined by ``shcs2``.
    The synthesis of area-mean values is done up to degree ``nmax1``
    and the irregular surface is expanded up to degree ``nmax2``.
    ``nmax3`` and ``nmax4`` represent the maximum harmonic degrees to
    synthesize and analyze the ``(shcs1.r / r)^(n + 1)`` terms, where
    ``r`` stands for the spherical radius of the irregular surface defined
    by ``shcs2``.  Refer to `charm_shs <./api-c-shs.html>`_ for the full
    documentation.

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
        Maximum degree of the synthesis of ``(shcs1.r / r)^(n + 1)``
    nmax4 : integer
        Maximum degree of the analysis of ``(shcs1.r / r)^(n + 1)``

    Returns
    -------
    out : numpy array of floating points
        Area-mean values synthesized from ``shcs1`` at ``cell`` residing
        on the surface defined by ``shcs2``
    """

    _check_crd_instance(cell, _ph_crd.CellGrid, 'cell')
    _check_shcs(shcs1, 'shcs1')
    _check_shcs(shcs2, 'shcs2')
    _check_deg_ord(nmax1, 'degree')
    _check_deg_ord(nmax2, 'degree')
    _check_deg_ord(nmax3, 'degree')
    _check_deg_ord(nmax4, 'degree')
    _check_nmax(nmax1, shcs1.nmax, 'nmax1', 'shcs1.nmax')
    _check_nmax(nmax2, shcs2.nmax, 'nmax2', 'shcs2.nmax')

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


def _check_shcs(shcs, varname):
    """
    Private function to check whether 'shcs' is an 'Shc' class instance.

    Parameters
    ----------
    shcs : Shc
        A 'Shc' class instance to be checked
    varname : str
        Caller's name of the 'shcs' variable.
    """

    if not isinstance(shcs, _ph_shc.Shc):
        raise TypeError(f'\'{varname}\' must be a \'{_ph_shc.Shc}\' '
                        f'class instance.')

    return


def _check_crd_instance(instance, base, varname):
    """
    Private function to check whether a crd 'instance' of a 'base' class.

    Parameters
    ----------
    instance : some crd instance
        Some point or cell instance to be checked.
    base : base class
        Requried base class of 'instance'.
    varname : str
        Caller's name of the 'instance' variable.
    """

    if not isinstance(instance, base):
        msg = f'\'{varname}\' must be an instance of one of following classes:'
        for i, cl in enumerate(base.__subclasses__()):
            msg += f' {base.__subclasses__()[i]}'
            if i < len(base.__subclasses__()) - 1:
                msg += f','
        msg += f'.'
        raise TypeError(msg)

    return


def _check_nmax(nmax_synth, nmax_shcs, var_nmax_synth, var_nmax_shcs):
    """
    Private function to check whether 'nmax_synth <= nmax_shcs'.

    Parameters
    ----------
    nmax_synth : int
        Maximum harmonic degree of the synthesis.
    nmax_shcs : int
        Maximum harmonic degree of a "Shc" class instance.
    var_nmax_synth : str
        Caller's name of the 'nmax_synth' variable.
    var_nmax_shcs : str
        Caller's name of the 'nmax_shcs' variable.
    """

    if nmax_synth > nmax_shcs:
        msg  = f'The maximum degree of the synthesis '
        msg += f'\'{var_nmax_synth} = {nmax_synth}\' '
        msg += f'cannot be larger than the maximum degree of the '
        msg += f'coefficients \'{var_nmax_shcs} = {nmax_shcs}\'.'
        raise ValueError(msg)

    return
