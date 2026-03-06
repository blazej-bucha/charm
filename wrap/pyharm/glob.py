"""
Module defining the :class:`pyharm.glob.Constants` class for global CHarm
variables.

Warning
-------

    This module is for experienced users only. Most users should not
    interact with it.
"""


import ctypes as _ct
import numpy as _np
from . import _libcharm, _CHARM
from ._check_types import _check_flt_scalar, _check_int_scalar
from ._data_types import _charm_flt, _pyharm_flt


class Constants:
    """
    Class for global CHarm variables.
    """


    @property
    def threshold(self):
        """
        The CHarm's ``charm_glob_threshold`` variable.  To change its
        default value to, say, ``1e-12``, use the following command:

        >>> pyharm.glob.Constants().threshold = 1e-12

        """
        return _pyharm_flt(_ct.c_double.in_dll(_libcharm,
                                               _CHARM + 'glob_threshold'))


    @threshold.setter
    def threshold(self, c):
        _check_flt_scalar(c, '\'threshold\'')
        if c < 0.0:
            raise ValueError('\'threshold\' must not be negative.')

        _ct.c_double.in_dll(_libcharm, _CHARM + 'glob_threshold').value = \
            _charm_flt(c)

        return


    @property
    def threshold2(self):
        """
        The CHarm's ``charm_glob_threshold2`` variable.  To change its
        default value to, say, ``1e-11``, use the following command:

        >>> pyharm.glob.Constants().threshold2 = 1e-11

        """
        return _pyharm_flt(_ct.c_double.in_dll(_libcharm,
                                               _CHARM + 'glob_threshold2'))


    @threshold2.setter
    def threshold2(self, c):
        _check_flt_scalar(c, '\'threshold2\'')
        if c < 0.0:
            raise ValueError('\'threshold2\' must not be negative.')

        _ct.c_double.in_dll(_libcharm, _CHARM + 'glob_threshold2').value = \
            _charm_flt(c)

        return


    @property
    def polar_optimization_a1(self):
        """
        The CHarm's ``charm_glob_polar_optimization_a1`` variable.  To
        change its default value to, say, ``50``, use the following command:

        >>> pyharm.glob.Constants().polar_optimization_a1 = 50

        """
        return _ct.c_ulong.in_dll(_libcharm,
                                  _CHARM + 'glob_polar_optimization_a1').value


    @polar_optimization_a1.setter
    def polar_optimization_a1(self, c):
        _check_int_scalar(c, '\'polar_optimization_a1\'')
        if c < 0:
            raise ValueError('\'polar_optimization_a1\' must not be negative.')

        _ct.c_ulong.in_dll(_libcharm,
                            _CHARM + 'glob_polar_optimization_a1').value = c

        return


    @property
    def polar_optimization_a2(self):
        """
        The CHarm's ``charm_glob_polar_optimization_a2`` variable.  To
        change its default value to, say, ``0.02``, use the following
        command:

        >>> pyharm.glob.Constants().polar_optimization_a2 = 0.02

        """
        return _pyharm_flt(_ct.c_double.in_dll(_libcharm,
                                               _CHARM + \
                                               'glob_polar_optimization_a2'))


    @polar_optimization_a2.setter
    def polar_optimization_a2(self, c):
        _check_flt_scalar(c, '\'polar_optimization_a2\'')

        _ct.c_double.in_dll(_libcharm,
                            _CHARM + 'glob_polar_optimization_a2').value = \
            _charm_flt(c)

        return

