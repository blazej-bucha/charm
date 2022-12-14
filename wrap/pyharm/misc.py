"""
Module defining some miscellaneous functions, constants, etc.

Note
----
All functions that deal with numerics are written in double precision.
"""


from . import _libcharm, _CHARM


def print_version():
    """
    Prints library info (library name, version, compilation date, precision,
    etc).
    """

    func          = _libcharm[_CHARM + 'misc_print_version']
    func.restype  = None
    func.argtypes = None

    func()

    return

