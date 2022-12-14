# Implements the "err" module of PyHarm.
#
# NOTE: Unlike in CHarm, users do *not* interact with the "err" module in
# PyHarm!


import ctypes as _ct
from . import _libcharm, _libcharmname, _CHARM
from ._get_module_constants import _get_module_constants
from ._constants import _globals
from ._data_types import _ct_int
from ._check_types import _check_pointer


# Get the module constants from "_constants.py" and add them to the module's
# namespace
_get_module_constants('CHARM_ERR_')


# Character to represent the C's termination character
NULL_CHAR = '\0'


class _Err(_ct.Structure):
    """
    Private class to represent the `charm_err` structure of CHarm.
    """

    _fields_ = [('level',       _ct.c_uint),
                ('file',        _ct.POINTER(_ct.POINTER(_ct.c_char))),
                ('line',        _ct.POINTER(_ct.c_uint)),
                ('func',        _ct.POINTER(_ct.POINTER(_ct.c_char))),
                ('code',        _ct_int),
                ('msg',         _ct.POINTER(_ct.c_char)),
                ('issaturated', _ct.c_bool)]


class CHarmError(Exception):
    pass


def init():
    """
    Creates an empty instance of the `_Err` class.

    Returns
    -------
    out : _Err
        An instance of the `_Err` class.
    """

    func          = _libcharm[_CHARM + 'err_init']
    func.restype  = _ct.POINTER(_Err)
    func.argtypes = None

    err = func()
    _check_pointer(err, _CHARM + 'err_init', _libcharmname)

    return err


def free(a):
    """
    Frees the memory associated with an instance `a` of the `_Err` class.

    Parameters
    ----------

    a : _Err
        An instance of the `_Err` class to be freed.
    """

    if a is not None:
        func          = _libcharm[_CHARM + 'err_free']
        func.restype  = None
        func.argtypes = [_ct.POINTER(_Err)]
        func(a)

    return


def isempty(a):
    """
    Returns `True` if an instance of the `_Err` class `a` is empty or `False` 
    otherwise.

    Parameters
    ----------
    a : _Err
        An instance of the `_Err` class to be checked

    Returns
    -------
    out : bool
        `True` if `a` is empty, `False` otherwise.
    """

    func          = _libcharm[_CHARM + 'err_isempty']
    func.restype  = _ct.c_bool
    func.argtypes = [_ct.POINTER(_Err)]

    return func(a)


def reset(a):
    """
    Resets an instance `a` of the `_Err` class.

    Parameters
    ----------
    a : _Err
        An instance of the `_Err` class to be reset.

    Returns
    -------
    out : _Err
        `a` reset to empty values
    """

    func          = _libcharm[_CHARM + 'err_reset']
    func.restype  = None
    func.argtypes = [_ct.POINTER(_Err)]

    return func(a)


def handler(a, terminate):
    """
    This is a Python implementation of the CHarm's `charm_err_handler`
    function.  This Python implementation allows better error handling in
    Python.

    Parameters
    ----------
    a : _Err
        An instance of the `_Err` class to be handled.
    terminate : bool
        If `True` and `a` is not empty, the execution of the program
        terminates.  If `False` and `a` is not empty, the error is printed and
        the program continues in code execution.  If `False` and `a` is empty,
        nothing happens.
    """

    if isempty(a):
        return

    e = a.contents

    err_msg  = f'\n'
    err_msg += f'Error code: {e.code}                   ' \
                'Traceback (most recent call last)\n'
    err_msg += f'\n'

    for l in range(e.level - 1, -1, -1):

        f    = e.file[l][:MAX_FILE].decode()
        idx1 = f.find(NULL_CHAR)

        func = e.func[l][:MAX_FUNC].decode()
        idx2 = func.find(NULL_CHAR)

        err_msg += f'   File \'{f[:idx1]}\', line: {e.line[l]}, '
        err_msg += f'function: \'{func[:idx2]}\'\n\n'

    if e.issaturated:
        err_msg += '   Warning: The levels of the error object are '
        err_msg += 'saturated.  Most recent function calls may therefore '
        err_msg += 'not be reported.\n'
        err_msg += f'\n'

    msg = e.msg[:MAX_MSG].decode()
    idx = msg.find(NULL_CHAR)
    err_msg += f'Error message: {msg[:idx]}\n'

    if terminate:
        raise CHarmError(err_msg)
    else:
        print(err_msg)

    return

