"""
Module to work with spherical harmonic coefficients.  Defines the :class:`Shc` 
class, including its methods that are designed to:

    * read/write spherical harmonic coefficients,
    * rescale spherical harmonic coefficients, and
    * compute (difference) degree variances and amplitudes.

Note
----
This documentation is written for double precision version of PyHarm.
"""


import os as _os
import ctypes as _ct
import numpy as _np
from . import _libcharm, _libcharmname, _CHARM, _pyharm
from ._get_module_constants import _get_module_constants
from ._data_types import _ct_ulong, _ct_flt, _ct_int, _charm_flt, _pyharm_flt
from ._get_empty_array import _get_empty_array
from ._check_types import _check_deg_ord, _check_radius, _check_flt_scalar, \
                          _check_flt_ndarray, _check_pointer
from . import _err as _ph_err


# Get the module constants from "_constants.py" and add them to the module's
# namespace
_get_module_constants('CHARM_SHC_')


#: int: Ordering scheme to write spherical harmonic coefficients with
#: the :meth:`pyharm.shc.Shc.to_file_tbl` method: harmonic degree varies
#: fastest.
WRITE_TBL_N: int


#: int: Ordering scheme to write spherical harmonic coefficients with
#: the :meth:`pyharm.shc.Shc.to_file_tbl` method: harmonic order varies
#: fastest.
WRITE_TBL_M: int


# If radius of the reference sphere is not specified by the user when creating
# an `Shc` class instance, this value will be used by default (unit sphere).
_R = _pyharm_flt(1.0)

# If the scaling parameter is not specified by the user when creating an `Shc`
# class instance, this value will be used by default (implying no scaling
# parameter).
_MU = _pyharm_flt(1.0)


class _Shc(_ct.Structure):
    """ Private class to represent the `charm_shc` structure of CHarm. """

    _fields_ = [('nmax',  _ct_ulong),
                ('mu',    _ct_flt),
                ('r',     _ct_flt),
                ('nc',    _ct.c_size_t),
                ('ns',    _ct.c_size_t),
                ('c',     _ct.POINTER(_ct.POINTER(_ct_flt))),
                ('s',     _ct.POINTER(_ct.POINTER(_ct_flt))),
                ('owner', _ct.c_bool)]


class Shc:
    """
    Class for spherical harmonic coefficients.

    To create an :class:`Shc` class instance, always use one of the following
    factory methods:

        * :meth:`from_garbage`,
        * :meth:`from_zeros`,
        * :meth:`from_arrays`,
        * :meth:`from_file_gfc`,
        * :meth:`from_file_mtx`,
        * :meth:`from_file_tbl`,
        * :meth:`from_file_bin`.

    Examples
    --------
    >>> import pyharm as ph
    >>> shcs = ph.shc.Shc.from_garbage(10)

    Parameters
    ----------
    nmax : integer
        Maximum spherical harmonic degree
    mu : floating point
        Scaling parameter
    r : floating point
        Radius of the reference sphere
    coeffs : None, 0 or tuple
        Determines the way of initializing spherical harmonic coefficients:

            * :obj:`None` to not initialize spherical harmonic coefficients
              (`malloc` in C),

            * :obj:`0` to set all spherical harmonic coefficients to zero
              (`calloc` in C),

            * :obj:`(c, s)` to define spherical harmonic coefficients based on
              two numpy floating point arrays :obj:`c` and :obj:`s`, each of
              shape :obj:`(((nmax + 2) * (nmax + 1)) / 2,)`.  For ordering
              scheme of the coefficients in the `c` and `s` arrays,
              see the :attr:`c` and :attr:`s` properties.

    Note
    ----
    Once an :class:`Shc` class instance is created, its attributes are not 
    writable.  The :obj:`r` and :obj:`mu` attributes can properly be changed by 
    the :meth:`rescale` method.
    """


    @property
    def nmax(self):
        """ Maximum harmonic degree of the spherical harmonic coefficients. """
        return self._nmax


    @property
    def mu(self):
        """
        Scaling parameter :math:`\mu` associated with the spherical 
        harmonic coefficients, for instance, the geocentric gravitational 
        constant. In case the coefficients are not associated with any scaling 
        parameter (as it is, for instance, with planetary topographies), simply 
        set this variable to :obj:`1.0` (not to :obj:`0.0`!).
        """
        return self._mu


    @property
    def r(self):
        """
        Radius of the reference sphere :math:`R`, to which the spherical 
        harmonic coefficients refer (are scaled). The value must be greater 
        than zero. To get the unit sphere, as needed, for instance, when 
        working with planetary topographies, set this variable to :obj:`1.0`.
        """
        return self._r


    @property
    def c(self):
        """
        The :math:`\\bar{C}_{nm}` spherical harmonic coefficients.  A numpy
        array with the shape (((:attr:`nmax` + 2) * (:attr:`nmax` + 1)) / 2,)
        and the array structure `C0,0`, `C1,0`, ..., `Cnmax,0`, `C1,1`, `C2,1`,
        ..., `Cnmax,nmax`.
        """
        return self._c


    @property
    def s(self):
        """
        The same as :attr:`c` but for the :math:`\\bar{S}_{nm}` spherical
        harmonic coefficients.
        """
        return self._s


    @property
    def owner(self):
        """
        * If :obj:`True`, the memory associated with spherical harmonic
          coefficients is owned by CHarm and therefore it is automatically
          deallocated by CHarm when the user deletes :class:`Shc` class
          instances or when the instances get out of scope, etc.  The
          :attr:`owner` attribute is :obj:`True` for :class:`Shc` instances
          returned by all factory methods *except* for the :meth:`from_arrays`
          method.

          **Examples**

          >>> import pyharm as ph
          >>> shcs = ph.shc.Shc.from_garbage(10)
          >>> del shcs # Deletes "shcs" and properly deallocates all memory
          >>>          # that is associated with "shcs"

        * If :obj:`False`, neither CHarm nor PyHarm own the memory so neither
          CHarm nor PyHarm will ever deallocate it.  This is the case of
          :class:`Shc` class instances returned by the :meth:`from_arrays`
          factory method.  This method builds :class:`Shc` instances from
          user-owned external numpy arrays, so it is the responsibility of the
          user to delete the original arrays (if needed) once the :class:`Shc`
          class instance is deleted or got out of scope, etc.

          **Examples**

          >>> import numpy as np
          >>> import pyharm as ph
          >>> nmax = 10
          >>> ncs = ((nmax + 2) * (nmax + 1)) // 2
          >>> c = np.random.randn(ncs)
          >>> s = np.random.randn(ncs)
          >>> shcs = ph.shc.Shc.from_arrays(nmax, c, s)
          >>> del shcs # Deletes "shcs" but not the original "c" and "s" arrays
          >>> # Here, you can still work with "c" and "s"
          >>> del c, s # Finally, release the memory associated with "c" and
          >>>          # "s" if needed

        """
        return self._owner


    def __init__(self, nmax, mu, r, coeffs):

        self._nmax        = None
        self._mu          = None
        self._r           = None
        self._c           = _get_empty_array()
        self._nc          = None
        self._s           = _get_empty_array()
        self._ns          = None
        self._owner       = None
        self._from_method = None
        self._Shc         = None

        _check_deg_ord(nmax, 'degree')
        _check_flt_scalar(mu, 'Scaling parameter')
        _check_radius(r)

        f = ''
        if coeffs is None or coeffs == 0:

            if coeffs is None:
                f = _CHARM + 'shc_malloc'
            elif coeffs == 0:
                f = _CHARM + 'shc_calloc'

            func          = _libcharm[f]
            func.restype  = _ct.POINTER(_Shc)
            func.argtypes = [_ct_ulong,
                             _ct_flt,
                             _ct_flt]

            self._Shc = func(_ct_ulong(nmax), _ct_flt(mu), _ct_flt(r))

        elif isinstance(coeffs, tuple):

            if len(coeffs) != 2:
                raise ValueError('The length of the \'coeffs\' tuple must be '
                                 '2.')

            _check_flt_ndarray(coeffs[0], 1, 'The \'c\' item form the '
                                             '\'coeffs\' tuple')
            _check_flt_ndarray(coeffs[1], 1, 'The \'s\' item from the '
                                             '\'coeffs\' tuple')

            ncs = ((nmax + 2) * (nmax + 1)) // 2
            if coeffs[0].shape != (ncs,):
                msg  = f'Wrong shape of the input \'c\' array.  The required '
                msg += f'shape for \'nmax = {nmax}\' is \'({ncs},)\'.'
                raise ValueError(msg)
            if coeffs[1].shape != (ncs,):
                msg  = f'Wrong shape of the input \'s\' array.  The required '
                msg += f'shape for \'nmax = {nmax}\' is \'({ncs},)\'.'
                raise ValueError(msg)

            f = _CHARM + 'shc_init'
            func          = _libcharm[f]
            func.restype  = _ct.POINTER(_Shc)
            func.argtypes = [_ct_ulong,
                             _ct_flt,
                             _ct_flt,
                            _ct.POINTER(_ct_flt),
                            _ct.POINTER(_ct_flt)]

            self._Shc = func(_ct_ulong(nmax),
                             _ct_flt(mu),
                             _ct_flt(r),
                             coeffs[0].ctypes.data_as(_ct.POINTER(_ct_flt)),
                             coeffs[1].ctypes.data_as(_ct.POINTER(_ct_flt)))

        else:
            raise ValueError('Unsupported value of the \'coeffs\' input '
                             'parameter.')

        _check_pointer(self._Shc, f, _libcharmname)

        self._Shc2Shc()
        self._from_method = coeffs

        return


    def __str__(self):

        ret  = f'nmax = {self.nmax}\n\n'
        ret += f'mu = {self.mu}\n\n'
        ret += f'r = {self.r}\n\n'
        ret += f'c = {self.c}\n\n'
        ret += f's = {self.s}\n\n'
        ret += f'owner = {self.owner}\n'

        return ret


    def __repr__(self):

        return f'{_pyharm}.shc.Shc({self.nmax}, {self.mu}, {self.r}, ' \
               f'{self._from_method})'


    def __del__(self):

        self._free()

        return


    def __exit__(self):

        self._free()

        return


    @classmethod
    def from_garbage(cls, nmax, mu=_MU, r=_R):
        """
        Returns an :class:`Shc` class instance with uninitialized spherical
        harmonic coefficients (`malloc` in C).

        Parameters
        ----------
        nmax : integer
            Maximum spherical harmonic degree
        mu : floating point
            Scaling parameter, optional.  Default is :obj:`1.0`.
        r : floating point
            Radius of the reference sphere, optional.  Default is :obj:`1.0`.

        Returns
        -------
        out : Shc
            An :class:`Shc` class instance
        """

        return cls(nmax, mu, r, None)


    @classmethod
    def from_zeros(cls, nmax, mu=_MU, r=_R):
        """
        Returns an :class:`Shc` class instance with all spherical harmonic
        coefficients initialized to zero (`calloc` in C).

        Parameters
        ----------
        nmax : integer
            Maximum spherical harmonic degree
        mu : floating point
            Scaling parameter, optional.  Default is :obj:`1.0`.
        r : floating point
            Radius of the reference sphere, optional.  Default is :obj:`1.0`.

        Returns
        -------
        out : Shc
            An :class:`Shc` class instance
        """

        return cls(nmax, mu, r, 0)


    @classmethod
    def from_arrays(cls, nmax, c, s, mu=_MU, r=_R):
        """
        Returns an :class:`Shc` class instance with spherical harmonic
        coefficients copied from the `c` and `s` input arrays.  The copy is
        shallow, meaning that the :attr:`c` and :attr:`s` attributes of the
        returned :class:`Shc` class instance share the memory space with the
        input `c` and `s` arrays.

        Parameters
        ----------
        nmax : integer
            Maximum spherical harmonic degree
        c : numpy array of floating points
            Spherical harmonic coefficients `Cnm`.  For details on the
            structure of the array, see :attr:`c`.
        s : numpy array of floating points
            Spherical harmonic coefficients `Snm`.  For details on the
            structure of the array, see :attr:`s`.
        mu : floating point
            Scaling parameter, optional.  Default is :obj:`1.0`.
        r : floating point
            Radius of the reference sphere, optional.  Default is :obj:`1.0`.

        Returns
        -------
        out : Shc
            An :class:`Shc` class instance
        """

        return cls(nmax, mu, r, (c, s))


    @classmethod
    def from_file_gfc(cls, pathname, nmax):
        """
        Reads spherical harmonic coefficients up to degree `nmax` from a `gfc`
        file specified in `pathname`.  For the structure of the input file,
        refer to `charm_shc <./api-c-shc.html>`_.

        Parameters
        ----------
        pathname : str
            Input file path
        nmax : integer
            Maximum harmonic degree to read the spherical harmonic coefficients 
            from the input file

        Returns
        -------
        out : Shc
            An :class:`Shc` class instance
        """

        return Shc._read_shc('gfc', pathname, nmax)


    @classmethod
    def from_file_mtx(cls, pathname, nmax):
        """
        Reads spherical harmonic coefficients up to degree `nmax` from a `mtx`
        file specified in `pathname`.  For the structure of the input file,
        refer to `charm_shc <./api-c-shc.html>`_.

        Parameters
        ----------
        pathname : str
            Input file path
        nmax : integer
            Maximum harmonic degree to read the spherical harmonic coefficients 
            from the input file

        Returns
        -------
        out : Shc
            An :class:`Shc` class instance
        """

        return Shc._read_shc('mtx', pathname, nmax)


    @classmethod
    def from_file_tbl(cls, pathname, nmax):
        """
        Reads spherical harmonic coefficients up to degree `nmax` from a `tbl`
        file specified in `pathname`.  For the structure of the input file,
        refer to `charm_shc <./api-c-shc.html>`_.

        Parameters
        ----------
        pathname : str
            Input file path
        nmax : integer
            Maximum harmonic degree to read the spherical harmonic coefficients 
            from the input file

        Returns
        -------
        out : Shc
            An :class:`Shc` class instance
        """

        return Shc._read_shc('tbl', pathname, nmax)


    @classmethod
    def from_file_bin(cls, pathname, nmax):
        """
        Reads spherical harmonic coefficients up to degree `nmax` from a `bin`
        file specified in `pathname`.  For the structure of the input file,
        refer to `charm_shc <./api-c-shc.html>`_.

        Parameters
        ----------
        pathname : str
            Input file path
        nmax : integer
            Maximum harmonic degree to read the spherical harmonic coefficients 
            from the input file

        Returns
        -------
        out : Shc
            An :class:`Shc` class instance
        """

        return Shc._read_shc('bin', pathname, nmax)


    def to_file_mtx(self, nmax, format, pathname):
        """
        Writes an :class:`Shc` class instance up to degree `nmax` to a text
        file specified in `pathname` using formatting for floating point
        numbers `format`.  For the structure of the output file, refer to
        `charm_shc <./api-c-shc.html>`_.

        Parameters
        ----------
        nmax : integer
            Maximum harmonic degree to write the spherical harmonic coefficients
        format : str
            Format to write floating point numbers, e.g., :obj:`'%0.16e'`
        pathname : str
            Output file path
        """

        self._check_nmax(nmax)

        if not isinstance(format, str):
            raise TypeError('\'format\' must be a string.')

        if not isinstance(pathname, str):
            raise TypeError('\'pathname\' must be a string.')

        self._create_path(pathname)

        func_name    = _CHARM + 'shc_write_mtx'
        func         = _libcharm[func_name]
        func.restype = None
        func.argtype = [_ct.POINTER(_Shc),
                        _ct_ulong,
                        _ct.c_char_p,
                        _ct.c_char_p,
                        _ct.POINTER(_ph_err._Err)]

        err = _ph_err.init()
        func(self._Shc,
             _ct_ulong(nmax),
             _ct.create_string_buffer(format.encode()),
             _ct.create_string_buffer(pathname.encode()),
             err)
        _ph_err.handler(err, 1)
        _ph_err.free(err)

        return


    def to_file_tbl(self, nmax, format, ordering, pathname):
        """
        Writes an :class:`Shc` class instance up to degree `nmax` to a text
        file specified in `pathname` using formatting for floating point
        numbers `format` and ordering scheme `ordering`.  For the structure of
        the output file, refer to `charm_shc <./api-c-shc.html>`_.

        Parameters
        ----------
        nmax : integer
            Maximum harmonic degree to write the spherical harmonic coefficients
        format : str
            Format to write floating point numbers, e.g., :obj:`'%0.16e'`
        ordering : integer
            Scheme to sort spherical harmonic coefficients,
            either :obj:`pyharm.shc.WRITE_TBL_N` or
            :obj:`pyharm.shc.WRITE_TBL_M`
        pathname : str
            Output file path
        """

        self._check_nmax(nmax)

        if not isinstance(format, str):
            raise TypeError('\'format\' must be a string.')

        if not isinstance(ordering, int):
            raise TypeError('\'ordering\' must be an integer.')

        if ordering not in [WRITE_TBL_N, WRITE_TBL_M]:
            raise ValueError('Unsupported value of \'ordering\'.')

        if not isinstance(pathname, str):
            raise TypeError('\'pathname\' must be a string.')

        self._create_path(pathname)

        func_name    = _CHARM + 'shc_write_tbl'
        func         = _libcharm[func_name]
        func.restype = None
        func.argtype = [_ct.POINTER(_Shc),
                        _ct_ulong,
                        _ct.c_char_p,
                        _ct.c_int,
                        _ct.c_char_p,
                        _ct.POINTER(_ph_err._Err)]

        err = _ph_err.init()
        func(self._Shc,
             _ct_ulong(nmax),
             _ct.create_string_buffer(format.encode()),
             _ct_int(ordering),
             _ct.create_string_buffer(pathname.encode()),
             err)
        _ph_err.handler(err, 1)
        _ph_err.free(err)

        return


    def to_file_bin(self, nmax, pathname):
        """
        Writes an :class:`Shc` class instance up to degree `nmax` to a binary 
        file specified in `pathname`.  For the structure of the output file, 
        refer to `charm_shc <./api-c-shc.html>`_.

        Parameters
        ----------
        nmax : integer
            Maximum harmonic degree to write the spherical harmonic
            coefficients
        pathname : str
            Output file path
        """

        self._check_nmax(nmax)

        if not isinstance(pathname, str):
            raise TypeError('\'pathname\' must be a string.')

        self._create_path(pathname)

        func_name    = _CHARM + 'shc_write_bin'
        func         = _libcharm[func_name]
        func.restype = None
        func.argtype = [_ct.POINTER(_Shc),
                        _ct_ulong,
                        _ct.c_char_p,
                        _ct.POINTER(_ph_err._Err)]

        err = _ph_err.init()
        func(self._Shc,
             _ct_ulong(nmax),
             _ct.create_string_buffer(pathname.encode()),
             err)
        _ph_err.handler(err, 1)
        _ph_err.free(err)

        return


    def get_coeffs(self, n=None, m=None):
        """
        Returns spherical harmonic coefficients `Cnm` and `Snm` of degree 
        :obj:`n` and order :obj:`m`.  If the returned variables are arrays, 
        these are deep copies, meaning that the arrays have their own memory 
        space.

        The behaviour of the method depends on the type of the input variables 
        :obj:`n` and :obj:`m`.

        * If both :obj:`n` and :obj:`m` are integers, returned are two
          spherical harmonic coefficients `Cnm` and `Snm`, each of degree `n`
          and order `m`.  The two returned values are floating points.

        * If :obj:`n` is integer and :obj:`m` is :obj:`None`, returned are
          arrays spherical harmonic coefficients `Cnm` and `Snm` of degree 
          :obj:`n` and all corresponding orders :obj:`m = 0, 1, ..., n`.

        * If :obj:`n` is :obj:`None` and :obj:`m` is integer, returned are
          arrays spherical harmonic coefficients `Cnm` and `Snm` of order 
          :obj:`m` and all corresponding degrees :obj:`n = m, m + 1, ..., 
          self.nmax`.

        * If :obj:`n` and :obj:`m` are lists of equal size, returned are arrays 
          spherical harmonic coefficients `Cnm` and `Snm` of degrees taken from 
          the :obj:`n` list and orders taken from the :obj:`m` list.

        Note
        ----
        When the method returns numpy arrays, the output arrays are always deep 
        copies of the original spherical harmonic coefficients.

        Examples
        --------
        >>> import numpy as np
        >>> import pyharm as ph
        >>> nmax = 10
        >>> ncs = ((nmax + 2) * (nmax + 1)) // 2
        >>> c = np.random.randn(ncs)
        >>> s = np.random.randn(ncs)
        >>> shcs = ph.shc.Shc.from_arrays(nmax, c, s)
        >>> shcs.get_coeffs(3, 2) # Returns "C3,2" and "S3,2"
        >>> shcs.get_coeffs(n=3) # Returns "C3,0", "C3,1", "C3,2", "C3,3" and 
        >>>                      # "S3,0", "S3,1", "S3,2", "S3,3"
        >>> shcs.get_coeffs(m=8) # Returns "C8,8", "C9,8", "C10,8" and 
        >>>                      # "S8,8", "S9,8", "S10,8"
        >>> n = [3, 5, 6]
        >>> m = [2, 4, 6]
        >>> shcs.get_coeffs(n, m) # Returns "C3,2", "C5,4", "C6,6" and "S3,2", 
        >>>                       # "S5,4", "S6,6"


        Parameters
        ----------
        n : integer, list of integers
            Spherical harmonic degree, optional.
        m : integer, list of integers
            Spherical harmonic order, optional.

        Returns
        -------
        c : floating point, numpy array of floating points
            Spherical harmonic coefficient(s) `Cnm` of degree `n` and order `m`
        s : floating point, numpy array of floating points
            The same as `c` but with the `Snm` coefficient(s)
        """

        if n is None and m is None:
            return

        ls = (list, _np.ndarray)
        if isinstance(n, ls) and isinstance(m, ls):
            if len(n) != len(m):
                raise ValueError(f'The lengths of the \'n\' and \'m\' lists '
                                 f'do not match ({len(n)} vs. {len(m)}).')

            retc = _np.zeros((len(n),), dtype=_pyharm_flt, order='C')
            rets = _np.zeros((len(n),), dtype=_pyharm_flt, order='C')

            for i, (n1, m1) in enumerate(zip(n, m)):
                self._check_nm(n1, 'degree')
                self._check_nm(m1, 'order')
                self._check_mlen(n1, m1)
                retc[i] = self._Shc.contents.c[m1][n1 - m1]
                rets[i] = self._Shc.contents.s[m1][n1 - m1]
        elif not isinstance(n, ls) and not isinstance(m, ls):
            if n is not None:
                self._check_nm(n, 'degree')

            if m is not None:
                self._check_nm(m, 'order')

            if n is not None and m is not None:
                self._check_mlen(n, m)
                idx = self._get_m_idx(m)
                retc, rets = self.c[idx + (n - m)], self.s[idx + (n - m)]
            elif n is not None and m is None:
                idx    = (n + 1) * [None]
                idx[0] = n
                for m in range(1, n + 1):
                    idx[m] = idx[m - 1] + self.nmax + 1 - m

                retc, rets = self.c[idx].copy(), self.s[idx].copy()
            elif n is None and m is not None:
                idx_min = self._get_m_idx(m)
                idx_max = idx_min + self.nmax + 1 - m

                retc, rets = self.c[idx_min:idx_max].copy(), \
                             self.s[idx_min:idx_max].copy()
        else:
            msg  = f'Either both or none of \'n\' and \'m\' can be lists.'
            raise ValueError(msg)

        return retc, rets


    def set_coeffs(self, n=None, m=None, c=None, s=None):
        """
        Sets spherical harmonic coefficients `Cnm` and `Snm` of degree `n` and
        order `m` to the values of `c` and `s`, respectively.

        * If both :obj:`n` and :obj:`m` are integers, sets the spherical
          harmonic coefficients `Cnm` and/or `Snm` to the input parameters `c`
          and `s`, respectively.  `c` and `s` must be floating point scalars.

        * If :obj:`n` is integer and :obj:`m` is :obj:`None`, sets spherical
          harmonic coefficients of degree :obj:`n` and all corresponding
          harmonic orders :obj:`m = 0, 1, ..., n` to the input parameters `c`
          and/or `s`.  `c` and `s` must be numpy floating point arrays of
          shapes :obj:`(n + 1,)`.

        * If :obj:`n` is :obj:`None` and :obj:`m` is integer, sets spherical
          harmonic coefficients of order :obj:`m` and all corresponding
          harmonic degrees :obj:`n = m, m + 1, ..., self.nmax` to the input
          parameters `c` and/or `s`.  `c` and `s` must be numpy floating point
          arrays of shapes :obj:`(self.nmax + 1 - m,)`.

        * If :obj:`n` and :obj:`m` are lists of equal size, sets spherical 
          spherical harmonic coefficients of degrees and orders taken from the 
          :obj:`n` and :obj:`m` lists, respectively, to the corresponding 
          values taken from the input parameters :obj:`c` and/or :obj:`s`.  The 
          length of the input parameters must match.

        Note
        ----
        If the object's :obj:`owner` attribute is :obj:`True`, the copy of the 
        new coefficients is deep.  If :obj:`owner` is :obj:`False`, the copy is 
        shallow.

        Examples
        --------
        >>> import numpy as np
        >>> import pyharm as ph
        >>> nmax = 10
        >>> ncs = ((nmax + 2) * (nmax + 1)) // 2
        >>> c = np.random.randn(ncs)
        >>> s = np.random.randn(ncs)
        >>> shcs = ph.shc.Shc.from_arrays(nmax, c, s)
        >>> # Set "C3,2" and "S3,2"
        >>> shcs.set_coeffs(3, 2, 3.4, 1.3)
        >>> # Set "S3,2"
        >>> shcs.set_coeffs(3, 2, s=5.3)
        >>> # Set "C3,0", "C3,1", "C3,2", "C3,3" and "S3,0", "S3,1", "S3,2", 
        >>> #"S3,3"
        >>> shcs.set_coeffs(n=3, c=np.array([1.1, 1.2, 1.3, 1.4]),
        >>>                      s=np.array([0.0, 1.2, 1.3, 1.4]))
        >>>  # Set "C8,8", "C9,8", "C10,8" and "S8,8", "S9,8", "S10,8"
        >>> shcs.set_coeffs(m=8, c=np.array([1.1, 1.2, 1.3]),
        >>>                      s=np.array([1.1, 1.2, 1.3]))
        >>>  # Set "C3,2", "C5,4", "C6,6" and "S3,2", "S5,4", "S6,6"
        >>> shcs.set_coeffs(n=[3, 5, 6], m=[2, 4, 6],
        >>>                 c=np.array([1.1, 1.2, 1.3]),
        >>>                 s=np.array([1.1, 1.2, 1.3]))

        Parameters
        ----------
        n : integer, list of integers
            Spherical harmonic degree, optional.
        m : integer, list of integers
            Spherical harmonic order, optional.
        c : floating point or numpy array of floating points
            Spherical harmonic coefficient(s) `Cnm`, optional.
        s : floating point or numpy array of floating points
            Spherical harmonic coefficient(s) `Snm`, optional.
        """

        if n is None and m is None:
            return

        ls = (list, _np.ndarray)
        if isinstance(n, ls) and isinstance(m, ls):
            if len(n) != len(m):
                raise ValueError(f'The lengths of the \'n\' and \'m\' lists '
                                 f'do not match ({len(n)} vs. {len(m)}).')

            def set_cs(x, cs_str):
                _check_flt_ndarray(c, 1, 'The \'{cs_str}\' variable')

                if len(x) != len(n):
                    msg  = f'The length of \'{cs_str}\' is {len(x)}, but must '
                    msg += f'be {len(n)} to match the length of \'n\' and '
                    msg += f'\'m\'.'
                    raise ValueError(msg)

                for i, (n1, m1) in enumerate(zip(n, m)):
                    self._check_nm(n1, 'degree')
                    self._check_nm(m1, 'order')
                    self._check_mlen(n1, m1)
                    if cs_str == 'c':
                        self._Shc.contents.c[m1][n1 - m1] = x[i]
                    elif cs_str == 's':
                        self._Shc.contents.s[m1][n1 - m1] = x[i]

                return

            if c is not None:
                set_cs(c, 'c')

            if s is not None:
                set_cs(s, 's')

        elif not isinstance(n, ls) and not isinstance(m, ls):
            if n is not None:
                self._check_nm(n, 'degree')

            if m is not None:
                self._check_nm(m, 'order')

            if n is not None and m is not None:
                self._check_mlen(n, m)

                if c is not None:
                    _check_flt_scalar(c, 'The spherical harmonic coefficient '
                                         '\'c\'')
                    self._Shc.contents.c[m][n - m] = c

                if s is not None:
                    _check_flt_scalar(s, 'The spherical harmonic coefficient '
                                         '\'s\'')
                    self._Shc.contents.s[m][n - m] = s
            elif n is not None and m is None:
                idx    = (n + 1) * [None]
                idx[0] = n
                for m in range(1, n + 1):
                    idx[m] = idx[m - 1] + self.nmax + 1 - m

                nc = n + 1
                if c is not None:
                    _check_flt_ndarray(c, 1, 'The \'c\' variable')

                    if c.shape != (nc,):
                        msg  = f'The \'c\' array must have {nc} element(s) '
                        msg += f'when \'n = {n}\'.'
                        raise ValueError(msg)

                    self.c[idx] = c

                if s is not None:
                    _check_flt_ndarray(s, 1, 'The \'s\' variable')

                    if s.shape != (nc,):
                        msg  = f'The \'s\' array must have {nc} element(s) '
                        msg += f'when \'n = {n}\'.'
                        raise ValueError(msg)

                    self.s[idx] = s
            elif n is None and m is not None:
                idx_min = self._get_m_idx(m)
                idx_max = idx_min + self.nmax + 1 - m

                nc = self.nmax + 1 - m
                if c is not None:
                    _check_flt_ndarray(c, 1, 'The \'c\' variable')

                    if c.shape != (nc,):
                        msg  = f'The \'c\' array must have {nc} element(s) '
                        msg += f'when \'m = {m}\'.'
                        raise ValueError(msg)

                    self.c[idx_min:idx_max] = c

                if s is not None:
                    _check_flt_ndarray(s, 1, 'The \'s\' variable')

                    if s.shape != (nc,):
                        msg  = f'The \'s\' array must have {nc} element(s) '
                        msg += f'when \'m = {m}\'.'
                        raise ValueError(msg)

                    self.s[idx_min:idx_max] = s

        else:
            msg  = f'Either both or none of \'n\' and \'m\' can be lists.'
            raise ValueError(msg)

        return


    def rescale(self, mu=None, r=None):
        """
        Rescales spherical harmonic coefficients to a new scaling parameter
        `mu` and/or a new radius of the reference sphere `r`.

        Parameters
        ----------
        mu : floating point
            New scaling parameter, optional.  If not provided, the :attr:`mu`
            attribute is not changed.
        r : floating point
            New radius of the reference sphere, optional.  If not provided, the
            :attr:`r` attribute is not changed.
        """

        if mu is None:
            mu = self.mu
        else:
            _check_flt_scalar(mu, 'Scaling parameter')

        if r is None:
            r = self.r
        else:
            _check_radius(r)

        func_name    = _CHARM + 'shc_rescale'
        func         = _libcharm[func_name]
        func.restype = None
        func.argtype = [_ct.POINTER(_Shc),
                        _ct_flt,
                        _ct_flt,
                        _ct.POINTER(_ph_err._Err)]

        err = _ph_err.init()
        func(self._Shc, _ct_flt(mu), _ct_flt(r), err)
        _ph_err.handler(err, 1)
        _ph_err.free(err)

        self._mu = mu
        self._r  = r

        return


    def dv(self, nmax=None):
        """
        Computes degree variances up to degree `nmax`.

        Parameters
        ----------
        nmax : integer
            Maximum spherical harmonic degree to compute the degree variances,
            optional.  If not provided, it will be set to :attr:`nmax`.

        Returns
        -------
        out : numpy array of floating points
            A numpy floating point array with the degree variances
        """

        return self._da_dv(_CHARM + 'shc_dv', nmax)


    def da(self, nmax=None):
        """
        Computes degree amplitudes up to degree `nmax`.

        Parameters
        ----------
        nmax : integer
            Maximum spherical harmonic degree to compute the degree amplitudes,
            optional.  If not provided, it will be set to `self.nmax`.

        Returns
        -------
        out : numpy array of floating points
            A numpy floating point array with the degree amplitudes
        """

        return self._da_dv(_CHARM + 'shc_da', nmax)


    def ddv(self, shcs, nmax=None):
        """
        Computes difference degree variances with respect to `shcs`  up to
        degree `nmax`.

        Parameters
        ----------
        shcs : Shc
            An instance of the `Shc` class.
        nmax : integer
            Maximum spherical harmonic degree to compute the degree variances,
            optional.  If not provided, `nmax` is set to the smallest of
            :attr:`self.nmax` and :obj:`shcs.nmax`.

        Returns
        -------
        out : numpy array of floating points
            A numpy floating point array with the difference degree variances
        """

        return self._dda_ddv(_CHARM + 'shc_ddv', shcs, nmax)


    def dda(self, shcs, nmax=None):
        """
        Computes difference degree amplitudes with respect to `shcs` up to
        degree `nmax`.

        Parameters
        ----------
        shcs : Shc
            An instance of the `Shc` class.
        nmax : integer
            Maximum spherical harmonic degree to compute the degree amplitudes,
            optional.  If not provided, `nmax` is set to the smallest of
            `self.nmax` and `shcs.nmax`.

        Returns
        -------
        out : numpy array of floating points
            A numpy floating point array with the difference degree amplitudes
        """

        return self._dda_ddv(_CHARM + 'shc_dda', shcs, nmax)


    def _Shc2Shc(self):
        """
        Private function to convert an :class:`_Shc` class instance in
        `self._Shc.contents` to an :class:`Shc` class instance in `self`.  The
        :attr:`c` and :attr:`s` attributes of `self` share the same memory
        space as the corresponding attributes of `self._Shc.contents.c` and
        `self._Shc.contents.s`.
        """

        self._nmax = int(self._Shc.contents.nmax)
        self._mu   = _charm_flt(self._Shc.contents.mu)
        self._r    = _charm_flt(self._Shc.contents.r)
        self._nc   = int(self._Shc.contents.nc)
        self._ns   = int(self._Shc.contents.ns)

        if not self._Shc.contents.c:
            raise ValueError('\'self._Shc.contents.c\' is a \'NULL\' pointer.')
        self._c = _np.ctypeslib.as_array(self._Shc.contents.c[0],
                                         shape=(self._nc,))

        if not self._Shc.contents.s:
            raise ValueError('\'self._Shc.contents.s\' is a \'NULL\' pointer.')
        self._s = _np.ctypeslib.as_array(self._Shc.contents.s[0],
                                         shape=(self._ns,))

        self._owner = bool(self._Shc.contents.owner)

        return


    def _free(self):
        """
        Frees the memory associated with the object.
        """

        if self._Shc is not None:
            func         = _libcharm[_CHARM + 'shc_free']
            func.restype = None
            func.argtype = [_ct.POINTER(_Shc)]
            func(self._Shc)

        return


    @classmethod
    def _read_shc(cls, file_type, pathname, nmax):
        """
        Private function to call CHarm functions to load spherical harmonic
        coefficients up to degree `nmax` from a file of `file_type` specified
        by `pathname`.

        Parameters
        ----------
        file_type : str
            Either `gfc`, `tbl`, `mtx` or `bin
        pathname : str
            Input file path
        nmax : integer
            Maximum harmonic degree to read spherical harmonic coefficients

        Returns
        -------
        out : Shc
            An `Shc` class instance
        """

        if file_type == 'gfc':
            func_name = _CHARM + 'shc_read_gfc'
        elif file_type == 'tbl':
            func_name = _CHARM + 'shc_read_tbl'
        elif file_type == 'mtx':
            func_name = _CHARM + 'shc_read_mtx'
        elif file_type == 'bin':
            func_name = _CHARM + 'shc_read_bin'
        else:
            raise ValueError('Unknown CHarm function to read spherical '
                             'harmonic coefficients.')

        if not isinstance(pathname, str):
            raise TypeError('\'pathname\' must be a string.')

        _check_deg_ord(nmax, 'degree')

        func         = _libcharm[func_name]
        func.restype = None
        func.argtype = [_ct.c_char_p,
                        _ct_ulong,
                        _ct.POINTER(_Shc),
                        _ct.POINTER(_ph_err._Err)]

        shcs = cls.from_garbage(nmax)

        err = _ph_err.init()
        func(_ct.create_string_buffer(pathname.encode()),
             _ct_ulong(nmax),
             shcs._Shc,
             err)
        _ph_err.handler(err, 1)
        _ph_err.free(err)

        shcs._Shc2Shc()

        return shcs


    def _get_m_idx(self, m):
        """
        Private function to return index of spherical harmonic coefficients
        `Cmm` and `Smm` in `self.c` and `self.s` arrays.

        Parameters
        ----------
        m : integer
            Spherical harmonic order

        Returns
        -------
        out : Index of `Cmm` and `Smm` in `self.c` and `self.s`
        """

        return int(m * (self.nmax + 2) - (m**2 + m) / 2)


    def _check_nm(self, v, do_str):
        """
        Private function to check spherical harmonic degree `v` if `do_str ==
        'degree'` or spherical harmonic order `v` if `do_str == 'order'`.

        Parameters
        ----------
        v : integer
            Spherical harmonic degree or order
        do_str : str
            Either 'degree' or 'order'
        """

        if do_str == 'degree':
            do_var = 'n'
        elif do_str == 'order':
            do_var = 'm'
        else:
            raise ValueError('Unsupported value of \'do_str\'.')

        _check_deg_ord(v, do_str)

        if v > self.nmax:
            msg  = f'Spherical harmonic {do_str} \'{do_var} = {v}\' is '
            msg += f'larger than the maximum degree of the object '
            msg += f'\'self.nmax = {self.nmax}\'.'
            raise ValueError(msg)

        return


    def _check_nmax(self, nmax):
        """
        Private function to check maximum harmonic degree when an :class:`Shc` 
        class instance already exists.

        Parameters
        ----------
        nmax : integer
            Maximum harmonic degree
        """

        _check_deg_ord(nmax, 'degree')
        if nmax > self.nmax:
            msg =  f'\'nmax = {nmax}\' cannot be larger than '
            msg += f'\'self.nmax = {self.nmax}\'.'
            raise ValueError(msg)

        return


    @staticmethod
    def _check_mlen(n, m):
        """
        Private function to check whether spherical harmonic degree 'm' is 
        equal to or smaller than 'n'.

        Parameters
        ----------
        n : integer
            Spherical harmonic degree
        m : integer
            Spherical harmonic order
        """

        if n < m:
            msg  = f'Harmonic degree \'n = {n}\' cannot be smaller '
            msg += f'than harmonic order \'m = {m}\'.'
            raise ValueError(msg)

        return


    @staticmethod
    def _create_path(pathname):
        """
        Private function to create a path in 'pathname' if needed.

        Parameters
        ----------
        pathname : str
            Path to a file
        """

        dirpath = _os.path.dirname(pathname)

        if not _os.path.exists(dirpath):
            _os.makedirs(dirpath, exist_ok=True)

        return


    def _da_dv(self, f, nmax=None):
        """
        Private function to call CHarm functions 'charm_shc_da' or
        'charm_shc_dv' to compute degree amplitudes or degree variances,
        respectively.

        Parameters
        ----------
        f : str
            Either 'charm_shc_dv' or 'charm_shc_da'.
        nmax : integer
            Maximum spherical harmonic degree to compute the degree
            amplitudes/variances, optional.  If not provided, `nmax` will be
            set to `self.nmax`.

        Returns
        -------
        out : numpy array of floating points
            A numpy floating point array with the degree amplitudes or degree
            variances.
        """

        if nmax is None:
            nmax = self.nmax

        self._check_nmax(nmax)

        func         = _libcharm[f]
        func.restype = None
        func.argtype = [_ct.POINTER(_Shc),
                        _ct_ulong,
                        _ct.POINTER(_ct_flt),
                        _ct.POINTER(_ph_err._Err)]

        ret = _np.zeros((nmax + 1,), dtype=_pyharm_flt)

        err = _ph_err.init()
        func(self._Shc,
             _ct_ulong(nmax),
             ret.ctypes.data_as(_ct.POINTER(_ct_flt)),
             err)
        _ph_err.handler(err, 1)
        _ph_err.free(err)

        return ret


    def _dda_ddv(self, f, shcs, nmax=None):
        """
        Private function to call CHarm functions 'charm_shc_dda' or
        'charm_shc_ddv' to compute degree amplitudes or degree variances,
        respectively.

        Parameters
        ----------
        f : str
            Either 'charm_shc_ddv' or 'charm_shc_dda'.
        shcs : Shc
            An instance of the `Shc` class.
        nmax : integer
            Maximum spherical harmonic degree to compute the degree
            amplitudes/variances, optional.  If not provided, `nmax` will be
            set to `self.nmax`.

        Returns
        -------
        out : numpy array of floating points
            A numpy floating point array with the difference degree amplitudes
            or difference degree variances.
        """

        if not isinstance(shcs, Shc):
            raise TypeError('\'shcs\' must be instances of the \'Shc\' class.')

        if nmax is None:
            nmax = min(self.nmax, shcs.nmax)

        if nmax > self.nmax or nmax > shcs.nmax:
            msg  = f'\'nmax = {nmax}\' cannot be larger than '
            msg += f'\'self.nmax = {self.nmax}\' or '
            msg += f'\'shcs.nmax = {shcs.nmax}\'.'
            raise ValueError(msg)

        _check_deg_ord(nmax, 'degree')

        func         = _libcharm[f]
        func.restype = None
        func.argtype = [_ct.POINTER(_Shc),
                        _ct.POINTER(_Shc),
                        _ct_ulong,
                        _ct.POINTER(_ct_flt),
                        _ct.POINTER(_ph_err._Err)]

        ret = _np.zeros((nmax + 1,), dtype=_pyharm_flt)

        err = _ph_err.init()
        func(self._Shc,
             shcs._Shc,
             _ct_ulong(nmax),
             ret.ctypes.data_as(_ct.POINTER(_ct_flt)),
             err)
        _ph_err.handler(err, 1)
        _ph_err.free(err)

        return ret

