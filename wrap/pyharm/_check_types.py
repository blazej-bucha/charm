"""
Internal checks of data types passed by users to PyHarm.
"""


import numpy as _np
from ._data_types import _pyharm_ints, _pyharm_flt_all, _pyharm_flt


def _check_int_scalar(x, s):
    """
    Checks whether "x" is an instance of one of the classes in `_pyharm_ints`
    and if so, whether it is a scalar in order to represent the quantity named
    `s`.

    Parameters
    ----------
    x : any data type
        Data to check
    s : str
        String with the name of the data that is being checked with the first
        letter being uppercase (e.g. 'The \'crd_type\' variable').
    """

    if not isinstance(x, _pyharm_ints):
        msg  = f'{s} must be an instance of one of the following classes: '
        msg += f'{_pyharm_ints}.'
        raise TypeError(msg)

    if not _np.isscalar(x):
        raise ValueError(f'{s} must be a scalar.')


def _check_flt_scalar(x, s):
    """
    Checks whether `x` is an instance of one of the classes in
    `_pyharm_flt_all` and if so, whether it is a scalar to represent
    the quantity named `s`.

    Parameters
    ----------
    x : any data type
        Data to check
    s : str
        String with the name of the data that is being checked with the first
        letter being uppercase (e.g. 'Scaling parameter').
    """

    if not isinstance(x, _pyharm_flt_all):
        msg  = f'{s} must be an instance of one of the following classes: '
        msg += f'{_pyharm_flt_all}.'
        raise TypeError(msg)

    if not _np.isscalar(x):
        raise ValueError(f'{s} must be a scalar.')


def _check_deg_ord(x, do):
    """
    Checks whether `x` satisfies all conditions to represent harmonic degree
    (`do == 'degree'`) or harmonic order (`do == 'order'`).

    Parameters
    ----------
    x : any data type
    do : str
        Either 'degree' or 'order'
    """

    # Check whether "x" is of a correct integer class and whether it is
    # a scalar
    _check_int_scalar(x, f'Spherical harmonic {do}')
    if do != 'degree' and do != 'order':
        msg  = f'The \'do\' variable must be set either to \'degree\' '
        msg += f'or to \'order\'.'
        raise ValueError(msg)

    if x < 0:
        msg  = f'Spherical harmonic {do} cannot be smaller than zero.'
        raise ValueError(msg)


def _check_radius(x):
    """
    Checks whether `x` satisfies all conditions to represent spherical radius.

    Parameters
    ----------
    x : any data type
    """

    _check_flt_scalar(x, 'Spherical radius')

    if x <= 0.0:
        raise ValueError('Spherical radius must be larger than zero.')


def _check_flt_ndarray(x, n, s):
    """
    Checks whether `x` satisfies all conditions to represent an nD floating
    point numpy array that can safely enter CHarm in order to represent
    quantity named `s`.

    Parameters
    ----------
    x : any data type
        Object to be checked
    n : integer
        Required number of dimensions of the nD floating point array
    s : str
        String with the name of the object that is being checked with the first
        letter being uppercase (e.g., 'Latitudes')
    """

    # Must be a "_np.ndarray"
    if not isinstance(x, _np.ndarray):
        msg  = f'{s} must be an instance of one of the following classes: '
        msg += f'{_np.ndarray}.'
        raise TypeError(msg)

    # Must be of the correct "dtype"
    if x.dtype != _pyharm_flt:
        msg = f'{s} must be of the following \'dtype\': {_pyharm_flt}.'
        raise ValueError(msg)

    # Must have the correct number of dimensions
    if x.ndim != n:
        msg = f'{s} must be of the following dimension: \'ndim={n}\'.'
        raise ValueError(msg)

    # Must be C-contiguous
    if not x.flags.c_contiguous:
        msg  = f'{s} must be C-contiguous in memory (see the \'order\' flag '
        msg += f'of numpy arrays).'
        raise ValueError(msg)


def _check_pointer(pointer, func_name, lib_name):

    if not pointer:
        msg  = f'Function call to \'{func_name}\' from the CHarm library '
        msg += f'\'{lib_name}\' returned a NULL pointer.'
        raise ValueError(msg)

