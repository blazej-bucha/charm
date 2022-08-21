import os
import numpy as _np


def _load_lib(libdir, libname):
    """
    Private function to load a shared library 'libname' from directory 
    'libdir'.

    Parameters
    ----------
    libdir : str
        Path to the shared library
    libname : str
        Name of the shared library

    Returns
    -------
    out : ctypes.CDLL
        The CHarm shared library
    """

    try:
        return _np.ctypeslib.load_library(libname, libdir)
    except OSError:
        fullpath = os.path.join(libdir, libname)
        msg = f'Failed to load the shared library \'{fullpath}\'.'
        raise OSError(msg)

