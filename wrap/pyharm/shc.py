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
from warnings import warn as _warn
from . import _libcharm, _libcharmname, _CHARM, _pyharm
from ._get_module_constants import _get_module_constants
from ._data_types import _ct_ulong, _ct_flt, _ct_int, _charm_flt, _pyharm_flt
from ._get_empty_array import _get_empty_array
from ._check_types import _check_deg_ord, _check_radius, _check_flt_scalar, \
                          _check_flt_ndarray, _check_pointer, _check_int_scalar
from . import _err as _ph_err


# Get the module constants from "_constants.py" and add them to the module's
# namespace
_get_module_constants('CHARM_SHC_')


#: int: Ordering scheme to write spherical harmonic coefficients with
#: the :meth:`pyharm.shc.Shc.to_file` method: harmonic degree varies
#: fastest.
WRITE_N: int


#: int: Ordering scheme to write spherical harmonic coefficients with
#: the :meth:`pyharm.shc.Shc.to_file` method: harmonic order varies
#: fastest.
WRITE_M: int


# If radius of the reference sphere is not specified by the user when creating
# an ``Shc`` class instance, this value will be used by default (unit sphere).
_R = _pyharm_flt(1.0)

# If the scaling parameter is not specified by the user when creating an
# ``Shc`` class instance, this value will be used by default (implying no
# scaling parameter).
_MU = _pyharm_flt(1.0)


# All possible file types to read and write spherical harmonic coefficients.
# An exception is that 'gfc' format is not supported to write spherical
# harmonic coefficients.
_FILE_TYPES = ['gfc', 'bin', 'mtx', 'tbl', 'dov']


class _Shc(_ct.Structure):
    """ Private class to represent the ``charm_shc`` structure of CHarm. """

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
        * :meth:`from_file`.

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

            * ``None`` to not initialize spherical harmonic coefficients
              (``malloc`` in C),

            * ``0`` to set all spherical harmonic coefficients to zero
              (``calloc`` in C),

            * ``(c, s)`` to define spherical harmonic coefficients based on
              two numpy floating point arrays ``c`` and ``s``, each of
              shape ``(((nmax + 2) * (nmax + 1)) // 2,)``.  For ordering
              scheme of the coefficients in the ``c`` and ``s`` arrays,
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
        set this variable to ``1.0`` (not to ``0.0``!).
        """
        return self._mu


    @property
    def r(self):
        """
        Radius of the reference sphere :math:`R`, to which the spherical
        harmonic coefficients refer (are scaled). The value must be greater
        than zero. To get the unit sphere, as needed, for instance, when
        working with planetary topographies, set this variable to ``1.0``.
        """
        return self._r


    @property
    def c(self):
        """
        The :math:`\\bar{C}_{nm}` spherical harmonic coefficients.  A numpy
        array with the shape (((:attr:`nmax` + 2) * (:attr:`nmax` + 1)) / 2,)
        and the array structure :math:`\\bar{C}_{00}`, :math:`\\bar{C}_{10}`,
        ..., :math:`\\bar{C}_{\mathrm{nmax},0}`, :math:`\\bar{C}_{1,1}`, ...,
        :math:`\\bar{C}_{\mathrm{nmax},\mathrm{nmax}}`.
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
        * If ``True``, the memory associated with spherical harmonic
          coefficients is owned by CHarm and therefore it is automatically
          deallocated by CHarm when the user deletes :class:`Shc` class
          instances or when the instances get out of scope, etc.  The
          :attr:`owner` attribute is ``True`` for :class:`Shc` instances
          returned by all factory methods *except* for the :meth:`from_arrays`
          method.

          **Examples**

          >>> import pyharm as ph
          >>> shcs = ph.shc.Shc.from_garbage(10)
          >>> del shcs # Deletes "shcs" and properly deallocates all memory
          >>>          # that is associated with "shcs"

        * If ``False``, neither CHarm nor PyHarm own the memory so neither
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
        harmonic coefficients (``malloc`` in C).

        Parameters
        ----------
        nmax : integer
            Maximum spherical harmonic degree
        mu : floating point
            Scaling parameter, optional.  Default is ``1.0``.
        r : floating point
            Radius of the reference sphere, optional.  Default is ``1.0``.

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
        coefficients initialized to zero (``calloc`` in C).

        Parameters
        ----------
        nmax : integer
            Maximum spherical harmonic degree
        mu : floating point
            Scaling parameter, optional.  Default is ``1.0``.
        r : floating point
            Radius of the reference sphere, optional.  Default is ``1.0``.

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
        coefficients copied from the ``c`` and ``s`` input arrays.  The copy is
        shallow, meaning that the :attr:`c` and :attr:`s` attributes of the
        returned :class:`Shc` class instance share the memory space with the
        input ``c`` and ``s`` arrays.

        Parameters
        ----------
        nmax : integer
            Maximum spherical harmonic degree
        c : numpy array of floating points
            Spherical harmonic coefficients :math:`\\bar{C}_{nm}`.  For details
            on the structure of the array, see :attr:`c`.
        s : numpy array of floating points
            Spherical harmonic coefficients :math:`\\bar{S}_{nm}`.  For details
            on the structure of the array, see :attr:`s`.
        mu : floating point
            Scaling parameter, optional.  Default is ``1.0``.
        r : floating point
            Radius of the reference sphere, optional.  Default is ``1.0``.

        Returns
        -------
        out : Shc
            An :class:`Shc` class instance
        """

        return cls(nmax, mu, r, (c, s))


    @classmethod
    def from_file(cls, file_type, pathname, nmax, epoch=None):
        """
        Reads spherical harmonic coefficients up to degree ``nmax`` from the
        ``pathname`` file of a given ``file_type``.  For time variable models
        stored in the ``gfc`` file type, ``epoch`` additionally transforms the
        coefficients to a given epoch.

        .. tip:: To get the maximum harmonic degree stored in ``pathname``, use
                 :meth:`nmax_from_file`.

        .. tip:: To print all supported file types, use :meth:`get_file_types`:

                 >>> import pyharm as ph
                 >>> ph.shc.Shc.get_file_types()

                 For the structure of the file types, refer to `charm_shc
                 <./api-c-shc.html>`_.

        Parameters
        ----------
        file_type : str
            Type of the input file
        pathname : str
            Input file path
        nmax : integer
            Maximum harmonic degree to read the spherical harmonic coefficients
        epoch : str
            Epoch, to which spherical harmonic coefficients from ``gfc`` file
            types are transformed in case of time variable gravity field
            models; optional.  For the structure of the string, refer to
            `charm_shc <./api-c-shc.html>`_.  This input parameter is used only
            if ``file_type`` is ``gfc``.  Default value is ``None``.

        Returns
        -------
        out : Shc
            An :class:`Shc` class instance
        """

        # "nmax" needs to be checked already here, because "_read_shc" accepts
        # also "_NMAX_MODEL"
        _check_deg_ord(nmax, 'degree')
        return Shc._read_shc(file_type, pathname, nmax, epoch)


    def to_file(self, file_type, nmax, pathname, formatting='%0.16e',
                ordering=WRITE_N):
        """
        Writes an :class:`Shc` class instance up to degree ``nmax`` to a file
        ``pathname`` of a given ``file_type``.  If ``file_type`` represents
        a text file (``dov``, ``tbl`` or ``mtx``), ``formatting`` specifies the
        formatting for all floating point numbers.  If ``file_type`` is ``dov``
        or ``tbl``, ``ordering`` defines the ordering of spherical harmonic
        coefficients in the output file.

        .. tip:: To print all supported file types, use :meth:`get_file_types`:

                 >>> import pyharm as ph
                 >>> ph.shc.Shc.get_file_types()

                 For the structure of the file types, refer to `charm_shc
                 <./api-c-shc.html>`_.

        Parameters
        ----------
        file_type : str
            Type of the input file
        nmax : integer
            Maximum harmonic degree to write the spherical harmonic coefficients
        pathname : str
            Output file path
        formatting : str
            Formatting to write floating point numbers to text formats,
            optional.  Default is ``'%0.16e'``.
        ordering : integer
            Scheme to sort spherical harmonic coefficients when ``file_type``
            is ``dov`` or ``tbl``, optional.  Accepted values are
            :obj:`pyharm.shc.WRITE_N` and :obj:`pyharm.shc.WRITE_M`.  Default
            is :obj:`pyharm.shc.WRITE_N`.
        """

        self._write_shc(file_type, nmax, pathname, formatting, ordering)
        return


    @staticmethod
    def get_file_types():
        """
        Prints all file types supported by the :meth:`to_file` and
        :meth:`from_file` methods of the :class:`Shc` class.
        """

        print(Shc._get_file_types())

        return


    @staticmethod
    def nmax_from_file(file_type, pathname):
        """
        Returns the maximum harmonic degree of coefficients stored in
        the ``pathname`` file that is of a given ``file_type``.

        .. tip:: Use this method to get the maximum degree of spherical
                 harmonic coefficients stored in a file and then load the file
                 up to its maximum harmonic degree:

                 >>> import pyharm as ph
                 >>> path = '/some/path/to/gfc'
                 >>> nmax = ph.shc.Shc.nmax_from_file('gfc', path)
                 >>> shcs = ph.shc.Shc.from_file('gfc', path, nmax)

        .. tip:: To print all supported file types, use :meth:`get_file_types`:

                 >>> import pyharm as ph
                 >>> ph.shc.Shc.get_file_types()

                 For the structure of the file types, refer to `charm_shc
                 <./api-c-shc.html>`_.

        Parameters
        ----------
        file_type : str
            Type of the input file
        pathname : str
            Input file path

        Returns
        -------
        out : integer
            Maximum harmonic degree of the model
        """

        return Shc._read_shc(file_type, pathname, _NMAX_MODEL)


    def get_coeffs(self, n=None, m=None):
        """
        Returns spherical harmonic coefficients :math:`\\bar{C}_{nm}` and
        :math:`\\bar{S}_{nm}` of degree ``n`` and order ``m``.  If the returned
        variables are arrays, these are deep copies, meaning that the arrays
        have their own memory space.

        The behaviour of the method depends on the type of the input variables
        ``n`` and ``m``.

        * If both ``n`` and ``m`` are integers, returned are two
          spherical harmonic coefficients :math:`\\bar{C}_{nm}` and
          :math:`\\bar{S}_{nm}`, each of degree ``n`` and order ``m``.  The two
          returned values are floating points.

        * If ``n`` is integer and ``m`` is ``None``, returned are
          two arrays of spherical harmonic coefficients
          :math:`\\bar{C}_{nm}` and :math:`\\bar{S}_{nm}` of degree
          ``n`` and all corresponding orders ``m = 0, 1, ..., n``.

        * If ``n`` is ``None`` and ``m`` is integer, returned are two arrays of
          spherical harmonic coefficients :math:`\\bar{C}_{nm}` and
          :math:`\\bar{S}_{nm}` of order ``m`` and all corresponding
          degrees ``n = m, m + 1, ..., self.nmax``.

        * If ``n`` and ``m`` are lists of equal size, returned are two
          arrays of spherical harmonic coefficients :math:`\\bar{C}_{nm}`
          and :math:`\\bar{S}_{nm}` of degrees taken from the ``n`` list
          and orders taken from the ``m`` list.

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
        n : integer, list of integers, None
            Spherical harmonic degree, optional.
        m : integer, list of integers, None
            Spherical harmonic order, optional.

        Returns
        -------
        c : floating point, numpy array of floating points
            Spherical harmonic coefficient(s) :math:`\\bar{C}_{nm}` of degree
            ``n`` and order ``m``
        s : floating point, numpy array of floating points
            The same as ``c`` but with the :math:`\\bar{S}_{nm}` coefficient(s)
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
        Sets spherical harmonic coefficients :math:`\\bar{C}_{nm}` and
        :math:`\\bar{S}_{nm}` of degree ``n`` and order ``m`` to the values of
        ``c`` and ``s``, respectively.

        * If both ``n`` and ``m`` are integers, sets the spherical
          harmonic coefficients :math:`\\bar{C}_{nm}` and/or
          :math:`\\bar{S}_{nm}` to the input parameters ``c`` and/or ``s``,
          respectively.  ``c`` and ``s`` must be floating point scalars.

        * If ``n`` is integer and ``m`` is ``None``, sets spherical
          harmonic coefficients of degree ``n`` and all corresponding
          harmonic orders ``m = 0, 1, ..., n`` to the input parameters ``c``
          and/or ``s``.  ``c`` and ``s`` must be numpy floating point arrays of
          shapes ``(n + 1,)``.

        * If ``n`` is ``None`` and ``m`` is integer, sets spherical
          harmonic coefficients of order ``m`` and all corresponding
          harmonic degrees ``n = m, m + 1, ..., self.nmax`` to the input
          parameters ``c`` and/or ``s``.  ``c`` and ``s`` must be numpy
          floating point arrays of shapes ``(self.nmax + 1 - m,)``.

        * If ``n`` and ``m`` are lists of equal size, sets spherical
          spherical harmonic coefficients of degrees and orders taken from the
          ``n`` and ``m`` lists, respectively, to the corresponding
          values taken from the input parameters ``c`` and/or ``s``.  The
          length of the input parameters must match.

        Note
        ----
        If the object's :obj:`owner` attribute is ``True``, the copy of the
        new coefficients is deep.  If :obj:`owner` is ``False``, the copy is
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
        n : integer, list of integers, None
            Spherical harmonic degree, optional.
        m : integer, list of integers, None
            Spherical harmonic order, optional.
        c : floating point, numpy array of floating points, None
            Spherical harmonic coefficient(s) :math:`\\bar{C}_{nm}`, optional.
        s : floating point, numpy array of floating points, None
            Spherical harmonic coefficient(s) :math:`\\bar{S}_{nm}`, optional.
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


    def get_degrees_orders(self):
        """
        Returns arrays of spherical harmonic degrees and orders matching the
        spherical harmonic coefficients stored in :obj:`c` and :obj:`s`.

        Returns
        -------
        degrees : numpy array of numpy.uint
            Array of spherical harmonic degrees matching the coefficients in
            :obj:`c` and :obj:`s`.
        orders : numpy array of numpy.uint
            Array of spherical harmonic orders matching the coefficients in
            :obj:`c` and :obj:`s`.
        """

        degrees = _np.zeros(self.c.shape, dtype=_np.uint)
        orders  = _np.zeros(self.c.shape, dtype=_np.uint)

        i = 0
        for m in range(self.nmax + 1):
            for n in range(m, self.nmax + 1):
                degrees[i] = n
                orders[i]  = m
                i += 1

        return degrees, orders


    def rescale(self, mu=None, r=None):
        """
        Rescales spherical harmonic coefficients to a new scaling parameter
        ``mu`` and/or a new radius of the reference sphere ``r``.

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

        func_name     = _CHARM + 'shc_rescale'
        func          = _libcharm[func_name]
        func.restype  = None
        func.argtypes = [_ct.POINTER(_Shc),
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
        Computes degree variances up to degree ``nmax``.

        Parameters
        ----------
        nmax : integer
            Maximum spherical harmonic degree to compute the degree variances,
            optional.  If not provided, it will be set to the object's
            :attr:`nmax`.

        Returns
        -------
        out : numpy array of floating points
            A numpy floating point array with the degree variances
        """

        return self._da_dv(_CHARM + 'shc_dv', nmax)


    def da(self, nmax=None):
        """
        Computes degree amplitudes up to degree ``nmax``.

        Parameters
        ----------
        nmax : integer
            Maximum spherical harmonic degree to compute the degree amplitudes,
            optional.  If not provided, it will be set to ``self.nmax``.

        Returns
        -------
        out : numpy array of floating points
            A numpy floating point array with the degree amplitudes
        """

        return self._da_dv(_CHARM + 'shc_da', nmax)


    def ddv(self, shcs, nmax=None):
        """
        Computes difference degree variances with respect to ``shcs``  up to
        degree ``nmax``.

        Parameters
        ----------
        shcs : Shc
            An instance of the :obj:`Shc` class.
        nmax : integer
            Maximum spherical harmonic degree to compute the degree variances,
            optional.  If not provided, ``nmax`` is set to the smallest of
            :attr:`self.nmax` and ``shcs.nmax``.

        Returns
        -------
        out : numpy array of floating points
            A numpy floating point array with the difference degree variances
        """

        return self._dda_ddv(_CHARM + 'shc_ddv', shcs, nmax)


    def dda(self, shcs, nmax=None):
        """
        Computes difference degree amplitudes with respect to ``shcs`` up to
        degree ``nmax``.

        Parameters
        ----------
        shcs : Shc
            An instance of the :obj:`Shc` class.
        nmax : integer
            Maximum spherical harmonic degree to compute the degree amplitudes,
            optional.  If not provided, ``nmax`` is set to the smallest of
            ``self.nmax`` and ``shcs.nmax``.

        Returns
        -------
        out : numpy array of floating points
            A numpy floating point array with the difference degree amplitudes
        """

        return self._dda_ddv(_CHARM + 'shc_dda', shcs, nmax)


    @staticmethod
    def _get_file_types():
        """
        Returns a help string on all supported values of the ``file_type``
        input parameter to :meth:`to_file` and :meth:`from_file` methods of the
        :class:`Shc` class.
        """

        ret  = f'The following strings are accepted as the \'file_type\' '
        ret += f'input parameter to \'to_file\' and \'from_file\' methods '
        ret += f'of the \'Shc\' class:'
        for ft in _FILE_TYPES:
            ret += '\n'
            ret += '\t\'%s\'' % ft
        ret += '\n'
        ret += 'Note that \'gfc\' is supported only by the \'from_file\' '
        ret += 'method.'

        return ret


    def _Shc2Shc(self):
        """
        Private function to convert an :class:`_Shc` class instance in
        ``self._Shc.contents`` to an :class:`Shc` class instance in ``self``.
        The :attr:`c` and :attr:`s` attributes of ``self`` share the same
        memory space as the corresponding attributes of
        ``self._Shc.contents.c`` and ``self._Shc.contents.s``.
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
            func          = _libcharm[_CHARM + 'shc_free']
            func.restype  = None
            func.argtypes = [_ct.POINTER(_Shc)]
            func(self._Shc)

        return


    @classmethod
    def _read_shc(cls, file_type, pathname, nmax, epoch=None):
        """
        Private function to call CHarm functions to load spherical harmonic
        coefficients up to degree ``nmax`` from a file of ``file_type``
        specified by ``pathname``.

        Parameters
        ----------
        file_type : str
            Either ``gfc``, ``tbl``, ``mtx`` or ``bin``
        pathname : str
            Input file path
        nmax : integer
            Maximum harmonic degree to read spherical harmonic coefficients
        epoch : str
            Epoch of the output spherical harmonic coefficients.  Relevant only
            for "gfc" files with time variable coefficients.  Default is
            'None'.

        Returns
        -------
        out : Shc, integer
            If ``nmax`` is a valid maximum harmonic degree, the function
        returns an :obj:`Shc` class instance.  If ``nmax`` is ``_NMAX_MODEL``,
        the function returns the maximum degree found in ``pathname`` without
        reading the coefficients inside the file.
        """

        Shc._check_file_type(file_type)

        if file_type != 'gfc' and epoch is not None:
            _warn(f'The input parameter \'epoch\' is ignored for '
                  f'\'file_type = {file_type}\'.')

        if not isinstance(pathname, str):
            raise TypeError('\'pathname\' must be a string.')

        if nmax == _NMAX_MODEL:
            # In this special case, "nmax" can be negative, so check for
            # integer only
            _check_int_scalar(nmax, '_NMAX_MODEL')
        else:
            _check_deg_ord(nmax, 'degree')

        if epoch is not None and not isinstance(epoch, str):
            raise TypeError('\'epoch\' must be a string or \'None\'.')

        if file_type == 'gfc':
            func_name = _CHARM + 'shc_read_gfc'
        elif file_type == 'tbl':
            func_name = _CHARM + 'shc_read_tbl'
        elif file_type == 'mtx':
            func_name = _CHARM + 'shc_read_mtx'
        elif file_type == 'bin':
            func_name = _CHARM + 'shc_read_bin'
        elif file_type == 'dov':
            func_name = _CHARM + 'shc_read_dov'

        if file_type == 'gfc':
            func_gfc          = _libcharm[func_name]
            func_gfc.restype  = _ct_ulong
            func_gfc.argtypes = [_ct.c_char_p,
                                 _ct_ulong,
                                 _ct.c_char_p,
                                 _ct.POINTER(_Shc),
                                 _ct.POINTER(_ph_err._Err)]
        else:
            func_rest          = _libcharm[func_name]
            func_rest.restype  = _ct_ulong
            func_rest.argtypes = [_ct.c_char_p,
                                  _ct_ulong,
                                  _ct.POINTER(_Shc),
                                  _ct.POINTER(_ph_err._Err)]

        if epoch is None:
            epoch_ptr = None
        else:
            epoch_ptr = _ct.create_string_buffer(epoch.encode())

        err = _ph_err.init()
        if nmax == _NMAX_MODEL:
            if file_type == 'gfc':
                ret = func_gfc(_ct.create_string_buffer(pathname.encode()),
                               _ct_ulong(_get_nmax_model()),
                               epoch_ptr,
                               None,
                               err)
            else:
                ret = func_rest(_ct.create_string_buffer(pathname.encode()),
                                _ct_ulong(_get_nmax_model()),
                                None,
                                err)
        else:
            shcs = cls.from_garbage(nmax)
            if file_type == 'gfc':
                ret = func_gfc(_ct.create_string_buffer(pathname.encode()),
                               _ct_ulong(nmax),
                               epoch_ptr,
                               shcs._Shc,
                               err)
            else:
                ret = func_rest(_ct.create_string_buffer(pathname.encode()),
                                _ct_ulong(nmax),
                                shcs._Shc,
                                err)
        _ph_err.handler(err, 1)
        _ph_err.free(err)

        if nmax == _NMAX_MODEL:
            # Returns maximum harmonic degree of the model
            return ret
        else:
            shcs._Shc2Shc()
            # Returns "Shc" class instance
            return shcs


    def _write_shc(self, file_type, nmax, pathname, formatting='%0.16e',
                   ordering=WRITE_N):
        """
        Private function to call CHarm functions to write spherical harmonic
        coefficients up to degree ``nmax`` to a file ``pathname`` of a given
        ``file_type``.  Default formatting for floating point numbers in text
        files is '%0.16e' and the default ordering scheme for ``dov`` and
        ``tbl`` file types is ``WRITE_N``.

        Parameters
        ----------
        file_type : str
            Either ``tbl``, ``mtx`` or ``bin``
        nmax : integer
            Maximum harmonic degree to write the spherical harmonic
            coefficients
        pathname : str
            Output file path
        formatting : str
            Formatting to write floating point numbers to text files.  Default
            is `'%0.16e'`,
        ordering : int
            Ordering scheme for ``dov`` and ``tbl`` file types.  Default is
            ``WRITE_N``
        """

        Shc._check_file_type(file_type)
        if file_type == 'gfc':
            raise ValueError(f'Writing to the \'gfc\' file type is not '
                             f'supported.')

        self._check_nmax(nmax)

        if not isinstance(pathname, str):
            raise TypeError('\'pathname\' must be a string.')

        if not isinstance(formatting, str):
            raise TypeError('\'formatting\' must be a string.')

        if file_type == 'dov' and ordering not in [WRITE_N, WRITE_M]:
            raise ValueError('Unsupported value of \'ordering\'.')
        elif file_type == 'tbl' and ordering not in [WRITE_N, WRITE_M]:
            raise ValueError('Unsupported value of \'ordering\'.')

        if file_type == 'tbl':
            func_name = _CHARM + 'shc_write_tbl'
        elif file_type == 'mtx':
            func_name = _CHARM + 'shc_write_mtx'
        elif file_type == 'bin':
            func_name = _CHARM + 'shc_write_bin'
        elif file_type == 'dov':
            func_name = _CHARM + 'shc_write_dov'

        self._create_path(pathname)

        func         = _libcharm[func_name]
        func.restype = None
        if file_type == 'tbl' or file_type == 'dov':
            func.argtypes = [_ct.POINTER(_Shc),
                             _ct_ulong,
                             _ct.c_char_p,
                             _ct_int,
                             _ct.c_char_p,
                             _ct.POINTER(_ph_err._Err)]
        elif file_type == 'bin':
            func.argtypes = [_ct.POINTER(_Shc),
                             _ct_ulong,
                             _ct.c_char_p,
                             _ct.POINTER(_ph_err._Err)]
        else:
            func.argtypes = [_ct.POINTER(_Shc),
                             _ct_ulong,
                             _ct.c_char_p,
                             _ct.c_char_p,
                             _ct.POINTER(_ph_err._Err)]

        err = _ph_err.init()
        if file_type == 'tbl' or file_type == 'dov':
            func(self._Shc,
                 _ct_ulong(nmax),
                 _ct.create_string_buffer(formatting.encode()),
                 _ct_int(ordering),
                 _ct.create_string_buffer(pathname.encode()),
                 err)
        elif file_type == 'bin':
            func(self._Shc,
                 _ct_ulong(nmax),
                 _ct.create_string_buffer(pathname.encode()),
                 err)
        else:
            func(self._Shc,
                 _ct_ulong(nmax),
                 _ct.create_string_buffer(formatting.encode()),
                 _ct.create_string_buffer(pathname.encode()),
                 err)
        _ph_err.handler(err, 1)
        _ph_err.free(err)

        return


    def _get_m_idx(self, m):
        """
        Private function to return index of spherical harmonic coefficients
        ``Cmm`` and ``Smm`` in ``self.c`` and ``self.s`` arrays.

        Parameters
        ----------
        m : integer
            Spherical harmonic order

        Returns
        -------
        out : Index of ``Cmm`` and ``Smm`` in ``self.c`` and ``self.s``
        """

        return int(m * (self.nmax + 2) - (m**2 + m) / 2)


    def _check_nm(self, v, do_str):
        """
        Private function to check spherical harmonic degree ``v`` if ``do_str
        == 'degree'` or spherical harmonic order ``v`` if ``do_str ==
        'order'``.

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
    def _check_file_type(file_type):
        """
        Private function to check file type to read/write spherical harmonic
        coefficients from/to a file.

        Parameters
        ----------
        file_type : str
            File type
        """

        if not isinstance(file_type, str):
            raise TypeError('\'file_type\' must be a string.')

        if file_type not in _FILE_TYPES:
            raise ValueError(f'Unsupported file type {file_type}.')

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
        if dirpath == '':
            # Saving the result to the current working directory, so no
            # directory needs to be created
            return

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
            amplitudes/variances, optional.  If not provided, ``nmax`` will be
            set to ``self.nmax``.

        Returns
        -------
        out : numpy array of floating points
            A numpy floating point array with the degree amplitudes or degree
            variances.
        """

        if nmax is None:
            nmax = self.nmax

        self._check_nmax(nmax)

        func          = _libcharm[f]
        func.restype  = None
        func.argtypes = [_ct.POINTER(_Shc),
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
            An instance of the :obj:`Shc` class.
        nmax : integer
            Maximum spherical harmonic degree to compute the degree
            amplitudes/variances, optional.  If not provided, ``nmax`` will be
            set to ``self.nmax``.

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

        func          = _libcharm[f]
        func.restype  = None
        func.argtypes = [_ct.POINTER(_Shc),
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


def _get_nmax_model():
    """
    Private function to return the CHarm's ``CHARM_SHC_NMAX_MODEL`` value.
    """

    func          = _libcharm[_CHARM + 'shc_get_nmax_model']
    func.argtypes = None
    func.restype  = _ct_ulong

    return func()

