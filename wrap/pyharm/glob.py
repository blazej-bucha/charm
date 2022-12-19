"""
Module defining the :class:`pyharm.glob.Constants` class for global CHarm 
variables.
"""


import ctypes as _ct
import numpy as _np
from . import _libcharm, _CHARM
from ._check_types import _check_flt_scalar
from ._data_types import _charm_flt, _pyharm_flt


class Constants:
    """ Class for global CHarm variables.

    Warning
    -------

        This module is for experienced users only. Most users should not
        interact with it.
    """

    def __init__(self):
        self._threshold  = _pyharm_flt(_ct.c_double.in_dll(_libcharm,
                                       _CHARM + 'glob_threshold'))
        self._threshold2 = _pyharm_flt(_ct.c_double.in_dll(_libcharm,
                                       _CHARM + 'glob_threshold2'))


    @property
    def threshold(self):
        return self._threshold


    @threshold.getter
    def threshold(self):
        """
        The CHarm's :obj:`glob_threshold` variable.  To change its default
        value to, say, :obj:`1e-12`, use the following command

        >>> pyharm.glob.Constants().threshold = 1e-12

        """
        return self._threshold


    @threshold.setter
    def threshold(self, c):
        _check_flt_scalar(c, '\'threshold\'')
        if c < 0.0:
            raise ValueError('\'threshold\' must not be negative.')

        _ct.c_double.in_dll(_libcharm, _CHARM + 'glob_threshold').value = \
            _charm_flt(c)
        self._threshold = _charm_flt(c)


    @property
    def threshold2(self):
        return self._threshold2


    @threshold2.getter
    def threshold2(self):
        """
        The CHarm's :obj:`glob_threshold2` variable.  To change its default
        value to, say, :obj:`1e-11`, use the following command

        >>> pyharm.glob.Constants().threshold2 = 1e-11

        """
        return self._threshold2


    @threshold2.setter
    def threshold2(self, c):
        _check_flt_scalar(c, '\'threshold2\'')
        if c < 0.0:
            raise ValueError('\'threshold2\' must not be negative.')

        _ct.c_double.in_dll(_libcharm, _CHARM + 'glob_threshold2').value = \
            _charm_flt(c)
        self._threshold2 = _charm_flt(c)

