import numpy as _np
from ._data_types import _pyharm_flt


def _get_empty_array():
    """
    Private function to return something that PyHarm considers to be an empty
    numpy array.

    Returns
    -------
    out : numpy floating point array
        Empty numpy ndarray
    """

    return _np.array([], dtype=_pyharm_flt, order='C')

