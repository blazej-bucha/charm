"""
Module to work with the coordinates of evaluation points/cells:

    * defines the following classes :obj:`pyharm.crd.PointSctr`,
      :obj:`pyharm.crd.PointGrid`, :obj:`pyharm.crd.PointGridGL`,
      :obj:`pyharm.crd.PointGridDH1`, :obj:`pyharm.crd.PointGridDH2`,
      :obj:`pyharm.crd.CellSctr` and :obj:`pyharm.crd.CellGrid`,

    * computes Gauss--Legendre and Driscoll--Healy quadrature grids.

Note
----
This documentation is written for double precision version of PyHarm.
"""


import ctypes as _ct
import numpy as _np
from . import _libcharm, _libcharmname, _CHARM, _pyharm
from ._data_types import _ct_int, _ct_ulong, _ct_size_t, _ct_flt
from ._get_module_constants import _get_module_constants
from ._check_types import _check_deg_ord, _check_radius, _check_flt_ndarray, \
                          _check_int_scalar, _check_pointer, _check_pointer
from ._get_empty_array import _get_empty_array
from .shc import _R


# Get the module constants from "_constants.py" and add them to the module's
# namespace
_get_module_constants('CHARM_CRD_')


# List of all pre-defined quadrature grid types
_CRD_TYPES_QUADS = [POINT_GRID_GL, POINT_GRID_DH1, POINT_GRID_DH2]

# List of all point "Point" types
_CRD_TYPES_POINT = [POINT_SCATTERED, POINT_GRID] + _CRD_TYPES_QUADS

# List of all cell "crd" types
_CRD_TYPES_CELL = [CELL_SCATTERED, CELL_GRID]


class _Point(_ct.Structure):
    """
    Private class to represent the ``charm_point`` structure of CHarm.
    """

    _fields_ = [('type',   _ct_int),
                ('nlat',   _ct_size_t),
                ('nlon',   _ct_size_t),
                ('npoint', _ct_size_t),
                ('lat',    _ct.POINTER(_ct_flt)),
                ('lon',    _ct.POINTER(_ct_flt)),
                ('r',      _ct.POINTER(_ct_flt)),
                ('w',      _ct.POINTER(_ct_flt)),
                ('owner',  _ct.c_bool)]


class _PointBase:
    """
    Parent class of the derived :class:`PointSctr`, :class:`PointGrid`,
    :class:`PointGridGL`, :class:`PointGridDH1` and :class:`PointGridDH2`
    classes.

    Note
    ----
    It is not recommended for users to interact with this class.
    Always work with the derived classes.

    For details, refer to `charm_crd <./api-c-crd.html>`_.

    Parameters
    ----------
    crd_type : integer
        Point type, the ``type`` attribute of an instance of the ``_PointBase``
        class.  Accepted values are ``pyharm.crd.POINT_*``
    nlat : integer
        Number of latitudes
    nlon : integer
        Number of longitudes
    data : None, 0 or tuple
        Determines the way of initializing array elements:

            * ``None`` to not initialize array elements (``malloc`` in C),

            * ``0`` to set all array elements to zero (``calloc`` in C),

            * ``(lat, lon, r)`` to define the object's latitudes, longitudes
              and spherical radii, respectively.  The ``lat``, ``lon``
              and ``r`` variables must be numpy floating point arrays of
              dimension ``1``.  The number of array elements of the arrays
              depends on ``crd_type`` (for further detail, refer to
              `charm_crd <./api-c-crd.html>`_.)

    Note
    ----
    Once a :class:`_PointBase` class instance is created, its attributes are
    not writeable, but all array elements are writeable.
    """


    @property
    def lat(self):
        """ Spherical latitudes. """
        return self._lat


    @property
    def lon(self):
        """ Spherical longitudes. """
        return self._lon


    @property
    def r(self):
        """ Spherical radii. """
        return self._r


    @property
    def npoint(self):
        """ Total number of points. """
        return self._npoint


    @property
    def owner(self):
        """ Refer to :class:`pyharm.shc.Shc` for documentation. """
        return self._owner


    def __init__(self, crd_type=None, nlat=None, nlon=None, data=None):

        self._type         = None
        self._nlat         = None
        self._nlon         = None
        self._npoint       = None
        self._lat          = _get_empty_array()
        self._lon          = _get_empty_array()
        self._r            = _get_empty_array()
        self._w            = None
        self._owner        = None
        self._from_method  = None
        self._Point        = None

        if crd_type is None and nlat is None and nlon is None and data is None:
            return

        self._check_crd_type(crd_type)
        self._check_nlat_nlon(nlat, 'nlat')
        self._check_nlat_nlon(nlon, 'nlon')

        f = ''
        if data is None or data == 0:

            if data is None:
                f = _CHARM + 'crd_point_malloc'
            elif data == 0:
                f = _CHARM + 'crd_point_calloc'

            func          = _libcharm[f]
            func.restype  = _ct.POINTER(_Point)
            func.argtypes = [_ct_int,
                             _ct_size_t,
                             _ct_size_t]

            self._Point = func(crd_type, nlat, nlon)

        elif isinstance(data, tuple):

            self._check_data(crd_type, data)
            if nlat != data[0].shape[0]:
                msg  = f'The \'nlat = {nlat}\' input parameter does not match '
                msg += f'the size of the '
                msg += f'data[0].shape[0] = {data[0].shape[0]} value.'

            if nlon != data[1].shape[0]:
                msg  = f'The \'nlon = {nlon}\' input parameter does not match '
                msg += f'the size of the '
                msg += f'data[1].shape[0] = {data[1].shape[0]} value.'

            f             = _CHARM + 'crd_point_init'
            func          = _libcharm[f]
            func.restype  = _ct.POINTER(_Point)
            func.argtypes = [_ct_int,
                             _ct_size_t,
                             _ct_size_t,
                             _ct.POINTER(_ct_flt),
                             _ct.POINTER(_ct_flt),
                             _ct.POINTER(_ct_flt)]

            self._Point = func(crd_type,
                               nlat,
                               nlon,
                               data[0].ctypes.data_as(_ct.POINTER(_ct_flt)),
                               data[1].ctypes.data_as(_ct.POINTER(_ct_flt)),
                               data[2].ctypes.data_as(_ct.POINTER(_ct_flt)))

        else:
            raise ValueError('Unsupported value of the \'data\' input '
                             'parameter.')

        _check_pointer(self._Point, f, _libcharmname)
        self._Point2Point()
        self._from_method = data

        return


    def __str__(self):

        ret  = f'lat = {self.lat}\n\n'
        ret += f'lon = {self.lon}\n\n'
        ret += f'r = {self.r}\n\n'
        if self._w is not None:
            ret += f'w = {self.w}\n\n'
        ret += f'owner = {self.owner}\n'

        return ret


    def __repr__(self):

        return f'{_pyharm}.crd.{self.__class__.__name__}({self._type}, ' \
               f'{self._nlat}, {self._nlon}, {self._from_method})'


    def __del__(self):

        self._free()

        return


    def __exit__(self):

        self._free()

        return


    def _gl(self, nmax, r=_R):
        """
        Computes the Gauss--Legendre grid for harmonic degree ``nmax`` and
        spherical radius ``r``.
        """

        self._quadrature(_CHARM + 'crd_point_gl', nmax, r)

        return


    def _dh1(self, nmax, r=_R):
        """
        Computes the non-equiangular variant of the Driscoll--Healy grid for
        harmonic degree ``nmax`` and spherical radius ``r``.
        """

        self._quadrature(_CHARM + 'crd_point_dh1', nmax, r)

        return


    def _dh2(self, nmax, r=_R):
        """
        Computes the non-equiangular variant of the Driscoll--Healy grid for
        harmonic degree ``nmax`` and spherical radius ``r``.
        """

        self._quadrature(_CHARM + 'crd_point_dh2', nmax, r)

        return


    def _Point2Point(self):
        """
        Private function to transform the `_Point` class instance in
        ``self._Point`` to a :class:`_PointBase` class instance in ``self``.
        The array attributes of ``self._Point.contents.`` share the same memory
        space as the corresponding attributes of ``self``.
        """

        self._type   = int(self._Point.contents.type)
        self._nlat   = int(self._Point.contents.nlat)
        self._nlon   = int(self._Point.contents.nlon)
        self._npoint = int(self._Point.contents.npoint)


        def get_error_msg(arr):
            """
            Private function to return an error message.

            Parameters
            ----------
            arr : str
                Name of an array

            Returns
            -------
            out : str
                Error message
            """

            if not isinstance(arr, str):
                raise TypeError('\'arr\' must be a string.')

            ret  = f'Failed to convert a \'_Point\' class instance to a '
            ret += f'\'_PointBase\' class instance.  The \'{arr}\' attribute '
            ret += f'of the \'_Point\' class instance is a \'NULL\' pointer.'

            return ret

        if not self._Point.contents.lat:
            raise ValueError(get_error_msg('lat'))
        self._lat = _np.ctypeslib.as_array(self._Point.contents.lat,
                                           shape=(self._nlat,))

        if not self._Point.contents.lon:
            raise ValueError(get_error_msg('lon'))
        self._lon = _np.ctypeslib.as_array(self._Point.contents.lon,
                                           shape=(self._nlon,))

        if not self._Point.contents.r:
            raise ValueError(get_error_msg('r'))
        self._r = _np.ctypeslib.as_array(self._Point.contents.r,
                                         shape=(self._nlat,))

        if self._Point.contents.type in _CRD_TYPES_QUADS:
            if not self._Point.contents.w:
                raise ValueError(get_error_msg('w'))
            self._w = _np.ctypeslib.as_array(self._Point.contents.w,
                                             shape=(self._nlat,))
        else:
            # "NULL" pointer in "self._Point.contents.w" is perfectly valid for
            # non-quadrature grids
            self._w = None

        self._owner = bool(self._Point.contents.owner)

        return


    def _free(self):
        """
        Frees the memory associated with class instance.
        """

        if self._Point is not None:
            func          = _libcharm[_CHARM + 'crd_point_free']
            func.restype  = None
            func.argtypes = [_ct.POINTER(_Point)]
            func(self._Point)

        return


    def _quadrature(self, f, nmax, r=_R):
        """
        Private function to prepare quadrature grids.

        Parameters
        ----------
        f : str
            CHarm function to call, one of 'charm_crd_point_gl',
            'charm_crd_point_dh1' or 'charm_crd_point_dh2'.
        nmax : integer
            Maximum harmonic degree, for which CHarm will compute the
            quadrature grid specified by ``f``.
        r : floating point, optional
            Spherical radius of the grid points, default is the unit sphere.
        """

        if f != _CHARM + 'crd_point_gl' and f != _CHARM + 'crd_point_dh1' and \
            f != _CHARM + 'crd_point_dh2':
            f1 = _CHARM + 'crd_point_gl'
            f2 = _CHARM + 'crd_point_dh1'
            f3 = _CHARM + 'crd_point_dh2'
            raise ValueError(f'Unknown CHarm function \'{f}\'.  Supported '
                             f'functions are \'{f1}\', \'{f2}\' or \'{f3}\'')

        self._quad_grd_check_inputs(nmax, r)

        func          = _libcharm[f]
        func.restype  = _ct.POINTER(_Point)
        func.argtypes = [_ct_ulong,
                         _ct_flt]

        self._Point = func(_ct_ulong(nmax), _ct_flt(r))
        _check_pointer(self._Point, f, _libcharmname)

        self._Point2Point()

        return


    def _lock_arrays(self):
        """
        Private method to lock all array attributes.
        """

        self._lat.flags.writeable = False
        self._lon.flags.writeable = False
        self._r.flags.writeable   = False
        self._w.flags.writeable   = False

        return


    def _quad_grd_check_inputs(self, nmax, r):
        """
        Private function to check input parameters of the 'gl', 'dh1' and
        'dh2' functions that implement supported quadrature grids.

        Parameters
        ----------
        nmax : integer
            Maximum harmonic degree associated with the quadrature grid
        r : floating point
            Spherical radius of grid meridians
        """

        _check_deg_ord(nmax, 'degree')
        _check_radius(r)

        return


    def _check_nlat_nlon(self, x, nx):
        """
        Check whether ``x`` is a positive integer.

        Parameters
        ----------
        x : any data type
            Object to check
        nx : str
            Name of the value that ``x`` represents
        """

        _check_int_scalar(x, nx)

        if x <= 0:
            raise ValueError(f'\'{nx}\' must be positive.')

        return


    def _check_crd_type(self, crd_type):
        """
        Checks whether ``crd_type`` represents a supported point type.

        Parameters
        ----------
        crd_type : integer
            Point type, the ``type`` attribute of an instance of the
            ``_PointBase`` class
        """

        _check_int_scalar(crd_type, 'crd_type')
        if crd_type not in _CRD_TYPES_POINT:
            raise ValueError('Unsupported point type.')

        return


    @staticmethod
    def _check_data(crd_type, data):

        POINT_TYPES = [POINT_SCATTERED, POINT_GRID]
        if crd_type not in POINT_TYPES:
            msg  = f'The only supported \'crd_type\' values to '
            msg += f'create a class instance from arrays are '
            msg += f'\'{POINT_TYPES}\'.'
            raise ValueError(msg)

        if len(data) != 3:
            raise ValueError('The length of the \'data\' tuple must be '
                             '3.')

        _check_flt_ndarray(data[0], 1, 'The input parameter \'lat\' '
                                       '(\'data[0]\')')
        _check_flt_ndarray(data[1], 1, 'The input parameter \'lon\' '
                                       '(\'data[1]\')')
        _check_flt_ndarray(data[2], 1, 'The input parameter \'r\' '
                                       '(\'data[2]\')')

        if data[0].shape[0] != data[2].shape[0]:
            msg  = f'The number of spherical latitudes and spherical radii '
            msg += f'must match.'
            raise ValueError(msg)

        if crd_type == POINT_SCATTERED:
            if data[0].shape[0] != data[1].shape[0]:
                msg  = f'The number of spherical latitudes and spherical '
                msg += f'longitudes must match.'
                raise ValueError(msg)

        return


class PointSctr(_PointBase):
    """
    Class for scattered points.

    For details, refer to `charm_crd <./api-c-crd.html>`_.

    To create a :class:`PointSctr` class instance, always use one of the
    following factory methods:

        * :meth:`from_garbage`,
        * :meth:`from_zeros`,
        * :meth:`from_arrays`.

    Example
    --------
    * Create an instance with 5 scattered points and all array elements
      initialized to zero:

      >>> import pyharm as ph
      >>> grd = ph.crd.PointSctr.from_zeros(5)

    * Create an instance with 3 scattered points with latitudes, longitudes and
      spherical radii taken from numpy arrays:

      >>> import numpy as np
      >>> import pyharm as ph
      >>> lat  = np.array([0.1, 0.2, 0.3])
      >>> lon  = np.array([0.5, 0.6, 0.7])
      >>> r    = np.array([1.1, 1.2, 1.3])
      >>> sctr = ph.crd.PointGrid.from_arrays(lat, lon, r)

    Parameters
    ----------
    npoint : integer
        Number of scattered points
    data : None, 0 or tuple
        Determines the way of initializing array elements:

            * ``None`` to not initialize array elements (``malloc`` in C),

            * ``0`` to set all array elements to zero (``calloc`` in C),

            * ``(lat, lon, r)`` to define the object's latitudes, longitudes
              and spherical radii, respectively.  The ``lat``, ``lon``
              and ``r`` variables must be numpy floating point arrays of
              dimension ``1``.  The number of elements in ``lat``,
              ``lon`` and ``r`` must match.

    Note
    ----
    Once a :class:`PointSctr` class instance is created, its attributes are
    not writeable, but all array elements are writeable.
    """


    def __init__(self, npoint, data):

        super().__init__(POINT_SCATTERED, npoint, npoint, data)

        return


    def __repr__(self):

        return f'{_pyharm}.crd.{self.__class__.__name__}({self._nlat}, ' \
               f'{self._from_method})'


    @classmethod
    def from_garbage(cls, npoint):
        """
        Returns a :class:`PointSctr` class instance with uninitialized array
        elements (``malloc`` in C).

        Parameters
        ----------
        npoint : integer
            Number of scattered points

        Returns
        -------
        out : PointSctr
            :class:`PointSctr` class instance
        """

        return cls(npoint, None)

    @classmethod
    def from_zeros(cls, npoint):
        """
        Returns a :class:`PointSctr` class instance with all array elements
        initialized to zero (``calloc`` in C).

        Parameters
        ----------
        npoint : integer
            Number of scattered points

        Returns
        -------
        out : PointSctr
            :class:`PointSctr` class instance
        """

        return cls(npoint, 0)


    @classmethod
    def from_arrays(cls, lat, lon, r):
        """
        Returns a :class:`PointSctr` class instance with spherical latitudes,
        longitudes and spherical radii taken from ``lat``, ``lon`` and
        ``r``, respectively.

        Parameters
        ----------
        lat : numpy array of floating points
            Spherical latitudes
        lon : numpy array of floating points
            Spherical longitudes
        r : numpy array of floating points
            Spherical radii

        Returns
        -------
        out : PointSctr
            :class:`PointSctr` class instance
        """

        super()._check_data(POINT_SCATTERED, (lat, lon, r))
        return cls(lat.shape[0], (lat, lon, r))


class PointGrid(_PointBase):
    """
    Class for custom point grids.

    For details, refer to `charm_crd <./api-c-crd.html>`_.

    To create a :class:`PointGrid` class instance, always use one of the
    following factory methods:

        * :meth:`from_garbage`,
        * :meth:`from_zeros`,
        * :meth:`from_arrays`.

    Examples
    --------
    * Create an instance with 5 grid latitudes and 10 grid longitudes and
      initialize all array elements to zero:

      >>> import ph
      >>> grd = ph.crd.PointGrid.from_zeros(5, 10)

    * Create an instance with grid latitudes, longitudes and spherical radii
      taken from numpy arrays:

      >>> import numpy as np
      >>> import pyharm as ph
      >>> lat = np.array([0.1, 0.2, 0.3])
      >>> lon = np.array([0.5, 0.6, 0.7, 0.8, 0.9])
      >>> r   = np.array([1.1, 1.2, 1.3])
      >>> grd = ph.crd.PointGrid.from_arrays(lat, lon, r)


    Parameters
    ----------
    nlat : integer
        Number of grid latitudes
    nlon : integer
        Number of grid longitudes
    data : None, 0 or tuple
        Determines the way of initializing array elements:

            * ``None`` to not initialize array elements (``malloc`` in C),

            * ``0`` to set all array elements to zero (``calloc`` in C),

            * ``(lat, lon, r)`` to define the object's latitudes, longitudes
              and spherical radii, respectively.  The ``lat``, ``lon``
              and ``r`` variables must be numpy floating point arrays of the
              dimension ``1``.  The shape of ``lat`` and ``lon`` as
              well as the shape of ``lat`` and ``r`` must match.

    Note
    ----
    Once a :class:`PointGrid` class instance is created, its attributes are
    not writeable, but all array elements are writeable.
    """


    def __init__(self, nlat, nlon, data):

        super().__init__(POINT_GRID, nlat, nlon, data)

        return


    def __repr__(self):

        return f'{_pyharm}.crd.{self.__class__.__name__}({self._nlat}, ' \
               f'{self._nlon}, {self._from_method})'


    @classmethod
    def from_garbage(cls, nlat, nlon):
        """
        Returns a :class:`PointGrid` class instance with uninitialized array
        elements (``malloc`` in C).

        Parameters
        ----------
        nlat : integer
            Number of grid latitudes
        nlon : integer
            Number of grid longitudes

        Returns
        -------
        out : PointGrid
            :class:`PointGrid` class instance
        """

        return cls(nlat, nlon, None)


    @classmethod
    def from_zeros(cls, nlat, nlon):
        """
        Returns a :class:`PointGrid` class instance with all array
        elements initialized to zero (``calloc`` in C).

        Parameters
        ----------
        nlat : integer
            Number of grid latitudes
        nlon : integer
            Number of grid longitudes

        Returns
        -------
        out : PointGrid
            :class:`PointGrid` class instance
        """

        return cls(nlat, nlon, 0)


    @classmethod
    def from_arrays(cls, lat, lon, r):
        """
        Returns a :class:`PointGrid` class instance with spherical latitudes,
        longitudes and spherical radii taken from ``lat``, ``lon`` and
        ``r``, respectively.

        Parameters
        ----------
        lat : numpy array of floating points
            Spherical latitudes
        lon : numpy array of floating points
            Spherical longitudes
        r : numpy array of floating points
            Spherical radii

        Returns
        -------
        out : PointGrid
            :class:`PointGrid` class instance
        """

        super()._check_data(POINT_GRID, (lat, lon, r))
        return cls(lat.shape[0], lon.shape[0], (lat, lon, r))


class PointGridGL(_PointBase):
    """
    Class for Gauss--Legendre point grids.

    For details, refer to `charm_crd <./api-c-crd.html>`_.

    Create a :class:`PointGridGL` class instance like this:

    >>> import pyharm as ph
    >>> gl = ph.crd.PointGridGL(10)

    Parameters
    ----------
    nmax : integer
        Maximum harmonic degree associated with the quadrature grid
    r : floating point
        Spherical radius of grid meridians, optional.  Default value is
        ``1.0``.

    Note
    ----
    Once a :class:`PointGridGL` class instance is created, neither its
    attributes nor its array elements are writeable.
    """


    @property
    def w(self):
        """ Integration weights on the unit sphere. """
        return self._w


    @property
    def nmax(self):
        """
        Maximum harmonic degree associated with the quadrature grid points.
        """
        return self._nmax


    def __init__(self, nmax, r=_R):

        super().__init__()
        super()._gl(nmax, r)
        self._lock_arrays()
        self._nmax = nmax

        return


    def __repr__(self):

        return f'{_pyharm}.crd.{self.__class__.__name__}({self.nmax}, ' \
               f'{self.r[0]})'


class PointGridDH1(_PointBase):
    """
    Class for non-equiangular Driscoll--Healy point grids.

    For details, refer to `charm_crd <./api-c-crd.html>`_.

    Create a :class:`PointGridDH1` class instance like this:

    >>> import pyharm as ph
    >>> gl = ph.crd.PointGridDH1(10)

    Parameters
    ----------
    nmax : integer
        Maximum harmonic degree associated with the quadrature grid
    r : floating point
        Spherical radius of grid meridians, optional.  Default value is
        ``1.0``.

    Note
    ----
    Once a :class:`PointGridDH1` class instance is created, neither its
    attributes nor its array elements are writeable.
    """


    @property
    def w(self):
        """ Integration weights on the unit sphere. """
        return self._w


    @property
    def nmax(self):
        """
        Maximum harmonic degree associated with the quadrature grid points.
        """
        return self._nmax


    def __init__(self, nmax, r=_R):

        super().__init__()
        super()._dh1(nmax, r)
        self._lock_arrays()
        self._nmax = nmax

        return


    def __repr__(self):

        return f'{_pyharm}.crd.{self.__class__.__name__}({self.nmax}, ' \
               f'{self.r[0]})'


class PointGridDH2(_PointBase):
    """
    Class for equiangular Driscoll--Healy point grids.

    For details, refer to `charm_crd <./api-c-crd.html>`_.

    Create a :class:`PointGridDH2` class instance like this:

    >>> import pyharm as ph
    >>> gl = ph.crd.PointGridDH2(10)

    Parameters
    ----------
    nmax : integer
        Maximum harmonic degree associated with the quadrature grid
    r : floating point
        Spherical radius of grid meridians, optional.  Default value is
        ``1.0``.

    Note
    ----
    Once a :class:`PointGridDH2` class instance is created, neither its
    attributes nor its array elements are writeable.
    """


    @property
    def w(self):
        """
        Integration weights on the unit sphere.
        """
        return self._w


    @property
    def nmax(self):
        """
        Maximum harmonic degree associated with the quadrature grid points.
        """
        return self._nmax


    def __init__(self, nmax, r=_R):

        super().__init__()
        super()._dh2(nmax, r)
        self._lock_arrays()
        self._nmax = nmax

        return


    def __repr__(self):

        return f'{_pyharm}.crd.{self.__class__.__name__}({self.nmax}, ' \
               f'{self.r[0]})'


class _Cell(_ct.Structure):
    """
    Private class to represent the ``charm_cell`` structure of CHarm.
    """

    _fields_ = [('type',   _ct_int),
                ('nlat',   _ct_size_t),
                ('nlon',   _ct_size_t),
                ('ncell',  _ct_size_t),
                ('latmin', _ct.POINTER(_ct_flt)),
                ('latmax', _ct.POINTER(_ct_flt)),
                ('lonmin', _ct.POINTER(_ct_flt)),
                ('lonmax', _ct.POINTER(_ct_flt)),
                ('r',      _ct.POINTER(_ct_flt)),
                ('owner',  _ct.c_bool)]


class _CellBase:
    """
    Parent class of the derived :class:`CellSctr` and :class:`CellGrid`
    classes.

    Note
    ----
    It is not recommended for users to interact with this class.
    Always work with the derived classes.

    For details, refer to `charm_crd <./api-c-crd.html>`_.

    Parameters
    ----------
    crd_type : integer
        Cell type, the ``type`` attribute of an instance of the `_CellBase`
        class.  Accepted values are ``pyharm.crd.CELL_*``
    nlat : integer
        Number of latitudes
    nlon : integer
        Number of longitudes
    data : None, 0 or tuple
        Determines the way of initializing array elements:

            * ``None`` to not initialize array elements (``malloc`` in C),

            * ``0`` to set all array elements to zero (``calloc`` in C),

            * ``(latmin, latmax, lonmin, lonmax, r)`` to define the object's
              minimum and maximum latitudes, minimum and maximum longitudes and
              spherical radii, respectively.  The ``latmin``, ``latmax``,
              ``lonmin``, ``lonmax`` and ``r`` variables must be numpy
              floating point arrays of dimension ``1``.  The number of array
              elements of the arrays depends on ``crd_type`` (for further
              detail, refer to `charm_crd <./api-c-crd.html>`_.)

    Note
    ----
    Once a :class:`_CellBase` class instance is created, its attributes are
    not writeable, but all array elements are writeable.
    """


    @property
    def latmin(self):
        """ Minimum spherical cell latitudes. """
        return self._latmin


    @property
    def latmax(self):
        """ Maximum spherical cell latitudes. """
        return self._latmax


    @property
    def lonmin(self):
        """ Minimum spherical cell longitudes. """
        return self._lonmin


    @property
    def lonmax(self):
        """ Maximum spherical cell longitudes. """
        return self._lonmax


    @property
    def r(self):
        """ Spherical radii (constant over each cell). """
        return self._r


    @property
    def ncell(self):
        """ Total number of cells. """
        return self._ncell


    @property
    def owner(self):
        """ Refer to :class:`pyharm.shc.Shc` for documentation. """
        return self._owner


    def __init__(self, crd_type=None, nlat=None, nlon=None, data=None):

        self._type        = None
        self._nlat        = None
        self._nlon        = None
        self._ncell       = None
        self._latmin      = _get_empty_array()
        self._latmax      = _get_empty_array()
        self._lonmin      = _get_empty_array()
        self._lonmax      = _get_empty_array()
        self._r           = _get_empty_array()
        self._owner       = None
        self._from_method = None
        self._Cell        = None

        if crd_type is None and nlat is None and nlon is None and data is None:
            return

        self._check_crd_type(crd_type)
        self._check_nlat_nlon(nlat, 'nlat')
        self._check_nlat_nlon(nlon, 'nlon')

        f = ''
        if data is None or data == 0:

            if data is None:
                f = _CHARM + 'crd_cell_malloc'
            elif data == 0:
                f = _CHARM + 'crd_cell_calloc'

            func          = _libcharm[f]
            func.restype  = _ct.POINTER(_Cell)
            func.argtypes = [_ct_int,
                             _ct_size_t,
                             _ct_size_t]

            self._Cell = func(crd_type, nlat, nlon)

        elif isinstance(data, tuple):

            self._check_data(crd_type, data)
            for i in range(2):
                if nlat != data[i].shape[0]:
                    msg  = f'The \'nlat = {nlat}\' input parameter does not '
                    msg += f'match the size of the '
                    if i == 0:
                        msg += f'latmin.shape[0] = {data[i].shape[0]} '
                    elif i == 1:
                        msg += f'latmax.shape[0] = {data[i].shape[0]} '

                    msg += f'value.'

            for i in range(2, 4):
                if nlon != data[i].shape[0]:
                    msg  = f'The \'nlon = {nlon}\' input parameter does not '
                    msg += f'match the size of the '
                    if i == 2:
                        msg += f'lonmin.shape[0] = {data[i].shape[0]} '
                    elif i == 3:
                        msg += f'lonmax.shape[0] = {data[i].shape[0]} '
                    msg += f'value.'

            if nlat != data[4].shape[0]:
                msg  = f'The \'nlat = {nlat}\' input parameter does not '
                msg += f'match the size of the '
                msg += f'r.shape[0] = {data[4].shape[0]} value.'

            f             = _CHARM + 'crd_cell_init'
            func          = _libcharm[f]
            func.restype  = _ct.POINTER(_Cell)
            func.argtypes = [_ct_int,
                             _ct_size_t,
                             _ct_size_t,
                             _ct.POINTER(_ct_flt),
                             _ct.POINTER(_ct_flt),
                             _ct.POINTER(_ct_flt),
                             _ct.POINTER(_ct_flt),
                             _ct.POINTER(_ct_flt)]

            self._Cell = func(crd_type,
                              nlat,
                              nlon,
                              data[0].ctypes.data_as(_ct.POINTER(_ct_flt)),
                              data[1].ctypes.data_as(_ct.POINTER(_ct_flt)),
                              data[2].ctypes.data_as(_ct.POINTER(_ct_flt)),
                              data[3].ctypes.data_as(_ct.POINTER(_ct_flt)),
                              data[4].ctypes.data_as(_ct.POINTER(_ct_flt)))

        else:
            raise ValueError('Unsupported value of the \'data\' input '
                             'parameter.')

        _check_pointer(self._Cell, f, _libcharmname)
        self._Cell2Cell()
        self._from_method = data

        return


    def __str__(self):

        ret  = f'latmin = {self.latmin}\n\n'
        ret += f'latmax = {self.latmax}\n\n'
        ret += f'lonmin = {self.lonmin}\n\n'
        ret += f'lonmax = {self.lonmax}\n\n'
        ret += f'r = {self.r}\n\n'
        ret += f'owner = {self.owner}\n'

        return ret


    def __repr__(self):

        return f'{_pyharm}.crd.{self.__class__.__name__}({self._type}, ' \
               f'{self._nlat}, {self._nlon}, {self._from_method})'


    def __del__(self):

        self._free()

        return


    def __exit__(self):

        self._free()

        return


    def _Cell2Cell(self):
        """
        Private function to transform the `_Cell` class instance in
        ``self._Cell`` to a :class:`_CellBase` class instance in ``self``.  The
        array attributes of ``self._Cell.contents.`` share the same memory
        space as the corresponding attributes of ``self``.
        """

        self._type  = int(self._Cell.contents.type)
        self._nlat  = int(self._Cell.contents.nlat)
        self._nlon  = int(self._Cell.contents.nlon)
        self._ncell = int(self._Cell.contents.ncell)


        def get_error_msg(arr):
            """
            Private function to return an error message.

            Parameters
            ----------
            arr : str
                Name of an array

            Returns
            -------
            out : str
                Error message
            """

            if not isinstance(arr, str):
                raise TypeError('\'arr\' must be a string.')

            ret  = f'Failed to convert a \'_Cell\' class instance to a '
            ret += f'\'_CellBase\' class instance.  The \'{arr}\' attribute '
            ret += f'of the \'_Cell\' class instance is a \'NULL\' pointer.'

            return ret

        if not self._Cell.contents.latmin:
            raise ValueError(get_error_msg('latmin'))
        self._latmin = _np.ctypeslib.as_array(self._Cell.contents.latmin,
                                              shape=(self._nlat,))

        if not self._Cell.contents.latmax:
            raise ValueError(get_error_msg('latmax'))
        self._latmax = _np.ctypeslib.as_array(self._Cell.contents.latmax,
                                              shape=(self._nlat,))

        if not self._Cell.contents.lonmin:
            raise ValueError(get_error_msg('lonmin'))
        self._lonmin = _np.ctypeslib.as_array(self._Cell.contents.lonmin,
                                              shape=(self._nlon,))

        if not self._Cell.contents.lonmax:
            raise ValueError(get_error_msg('lonmax'))
        self._lonmax = _np.ctypeslib.as_array(self._Cell.contents.lonmax,
                                              shape=(self._nlon,))

        if not self._Cell.contents.r:
            raise ValueError(get_error_msg('r'))
        self._r = _np.ctypeslib.as_array(self._Cell.contents.r,
                                         shape=(self._nlat,))

        self._owner = bool(self._Cell.contents.owner)

        return


    def _free(self):
        """
        Frees the memory associated with class instance.
        """

        if self._Cell is not None:
            func          = _libcharm[_CHARM + 'crd_cell_free']
            func.restype  = None
            func.argtypes = [_ct.POINTER(_Cell)]
            func(self._Cell)

        return


    def _check_nlat_nlon(self, x, nx):
        """
        Check whether ``x`` is a positive integer.

        Parameters
        ----------
        x : any data type
            Object to check
        nx : str
            Name of the value that ``x`` represents
        """

        _check_int_scalar(x, nx)

        if x <= 0:
            raise ValueError(f'\'{nx}\' must be positive.')

        return


    def _check_crd_type(self, crd_type):
        """
        Checks whether ``crd_type`` represents a supported cell type.

        Parameters
        ----------
        crd_type : integer
            Cell type, the ``type`` attribute of an instance of the `_CellBase`
            class
        """

        _check_int_scalar(crd_type, 'crd_type')
        if crd_type not in _CRD_TYPES_CELL:
            raise ValueError('Unsupported cell type.')

        return


    @staticmethod
    def _check_data(crd_type, data):

        if crd_type not in _CRD_TYPES_CELL:
            msg  = f'The only supported \'crd_type\' values to '
            msg += f'create a class instance from arrays are '
            msg += f'\'{_CRD_TYPES_CELL}\'.'
            raise ValueError(msg)

        if len(data) != 5:
            raise ValueError('The length of the \'data\' tuple must be '
                             '5.')

        _check_flt_ndarray(data[0], 1, 'The \'latmin\' variable')
        _check_flt_ndarray(data[1], 1, 'The \'latmax\' variable')
        _check_flt_ndarray(data[2], 1, 'The \'lonmin\' variable')
        _check_flt_ndarray(data[3], 1, 'The \'lonmax\' variable')
        _check_flt_ndarray(data[4], 1, 'The \'r\' variable')

        if data[0].shape[0] != data[1].shape[0]:
            msg  = f'The number of minimum and maximum cell latitudes '
            msg += f'must match.'
            raise ValueError(msg)

        if data[2].shape[0] != data[3].shape[0]:
            msg  = f'The number of minimum and maximum cell longitudes '
            msg += f'must match.'
            raise ValueError(msg)

        if data[0].shape[0] != data[4].shape[0]:
            msg  = f'The number latitudes and radii must match.'
            raise ValueError(msg)

        if crd_type == CELL_SCATTERED:
            if data[0].shape[0] != data[2].shape[0]:
                msg  = f'The number latitudes and longitudes must match.'
                raise ValueError(msg)

        return


class CellSctr(_CellBase):
    """
    Class for scattered cells.

    For details, refer to `charm_crd <./api-c-crd.html>`_.

    To create a :class:`CellSctr` class instance, always use one of the
    following factory methods:

        * :meth:`from_garbage`,
        * :meth:`from_zeros`,
        * :meth:`from_arrays`.

    Examples
    --------
    * Create an instance with 5 scattered cells and all array elements
      initialized to zero:

      >>> import pyharm as ph
      >>> grd = ph.crd.CellSctr.from_zeros(5)

    * Create an instance with 3 scattered cells with latitudes, longitudes
      and spherical radii taken from numpy arrays:

      >>> import numpy as np
      >>> import pyharm as ph
      >>> latmin = np.array([0.1, 0.2, 0.3])
      >>> latmax = latmin + 0.1
      >>> lonmin = np.array([0.5, 0.6, 0.7])
      >>> lonmax = lonmin + 0.1
      >>> r      = np.array([1.1, 1.2, 1.3])
      >>> sctr   = ph.crd.CellSctr.from_arrays(latmin, latmax,
      >>>                                      lonmin, lonmax, r)

    Parameters
    ----------
    nlat : integer
        Number of grid latitudes
    nlon : integer
        Number of grid longitudes
    data : None, 0 or tuple
        Determines the way of initializing array elements:

            * ``None`` to not initialize array elements (``malloc`` in C),

            * ``0`` to set all array elements to zero (``calloc`` in C),

            * ``(latmin, latmax, lonmin, lonmax, r)`` to define the object's
              minimum and maximum cell latitudes, minimum and maximum cell
              longitudes and spherical radii, respectively.  The ``latmin``,
              ``latmax``, ``lonmin``, ``lonmax`` and ``r`` variables must be
              numpy floating point arrays of dimension ``1``.  All arrays must
              share the same number of elements.

    Note
    ----
    Once a :class:`CellSctr` class instance is created, its attributes are not
    writeable, but all array elements are writeable.
    """


    def __init__(self, ncell, data):

        super().__init__(CELL_SCATTERED, ncell, ncell, data)

        return


    def __repr__(self):

        return f'{_pyharm}.crd.{self.__class__.__name__}({self._nlat}, ' \
               f'{self._from_method})'


    @classmethod
    def from_garbage(cls, ncell):
        """
        Returns a :class:`CellSctr` class instance with uninitialized array
        elements (``malloc`` in C).

        Parameters
        ----------
        ncell : integer
            Number of scattered cells

        Returns
        -------
        out : CellSctr
            :class:`CellSctr` class instance
        """

        return cls(ncell, None)

    @classmethod
    def from_zeros(cls, ncell):
        """
        Returns a :class:`CellSctr` class instance with all array elements
        initialized to zero (``calloc`` in C).

        Parameters
        ----------
        ncell : integer
            Number of scattered cells

        Returns
        -------
        out : CellSctr
            :class:`CellSctr` class instance
        """

        return cls(ncell, 0)


    @classmethod
    def from_arrays(cls, latmin, latmax, lonmin, lonmax, r):
        """
        Returns a :class:`CellSctr` class instance with minimum and maximum
        spherical cell latitudes, minimum and maximum cell longitudes and
        spherical radii taken from ``latmin``, ``latmax``, ``lonmin``,
        ``lonmax`` and ``r``, respectively.

        Parameters
        ----------
        latmin : numpy array of floating points
            Minimum spherical cell latitudes
        latmax : numpy array of floating points
            Maximum spherical cell latitudes
        lonmin : numpy array of floating points
            Minimum spherical cell longitudes
        lonmax : numpy array of floating points
            Maximum spherical cell longitudes
        r : numpy array of floating points
            Spherical radii

        Returns
        -------
        out : CellSctr
            :class:`CellSctr` class instance
        """

        super()._check_data(CELL_SCATTERED,
                            (latmin, latmax, lonmin, lonmax, r))
        return cls(latmin.shape[0], (latmin, latmax, lonmin, lonmax, r))


class CellGrid(_CellBase):
    """
    Class for custom cell grids.

    For details, refer to `charm_crd <./api-c-crd.html>`_.

    To create a :class:`CellGrid` class instance, always use one of the
    following factory methods:

        * :meth:`from_garbage`,
        * :meth:`from_zeros`,
        * :meth:`from_arrays`.

    Examples
    --------
    * Create an instance with 5 grid latitudes and 10 grid longitudes and
      initialize all array elements to zero:

      >>> import pyharm as ph
      >>> grd = ph.crd.CellGrid.from_zeros(5, 10)

    * Create an instance with grid latitudes, longitudes and spherical radii
      taken from numpy arrays:

      >>> import numpy as np
      >>> import pyharm as ph
      >>> latmin = np.array([0.1, 0.2, 0.3])
      >>> latmax = latmin + 0.1
      >>> lonmin = np.array([0.5, 0.6, 0.7, 0.8, 0.9])
      >>> lonmax = lonmin + 0.1
      >>> r      = np.array([1.1, 1.2, 1.3])
      >>> grd    = ph.crd.CellGrid.from_arrays(latmin, latmax,
      >>>                                      lonmin, lonmax, r)

    Parameters
    ----------
    nlat : integer
        Number of grid latitudes
    nlon : integer
        Number of grid longitudes
    data : None, 0 or tuple
        Determines the way of initializing array elements:

            * ``None`` to not initialize array elements (``malloc`` in C),

            * ``0`` to set all array elements to zero (``calloc`` in C),

            * ``(latmin, latmax, lonmin, lonmax, r)`` to define the object's
              minimum and maximum cell latitudes, minimum and maximum cell
              longitudes and spherical radii, respectively.  The ``latmin``,
              ``latmax``, ``lonmin``, ``lonmax`` and ``r``
              variables must be numpy floating point arrays of dimension
              ``1``.  The ``latmin``, ``latmax`` and ``r`` arrays
              must have the same number of elements, the same is true for
              ``lonmin`` and ``lonmax``.

    Note
    ----
    Once a :class:`CellGrid` class instance is created, its attributes are not
    writeable, but all array elements are writeable.
    """


    def __init__(self, nlat, nlon, data):

        super().__init__(CELL_GRID, nlat, nlon, data)

        return


    def __repr__(self):

        return f'{_pyharm}.crd.{self.__class__.__name__}({self._nlat}, ' \
               f'{self._nlon}, {self._from_method})'


    @classmethod
    def from_garbage(cls, nlat, nlon):
        """
        Returns a :class:`CellGrid` class instance with uninitialized array
        elements (``malloc`` in C).

        Parameters
        ----------
        nlat : integer
            Number of grid latitudes
        nlon : integer
            Number of grid longitudes

        Returns
        -------
        out : CellGrid
            :class:`CellGrid` class instance
        """

        return cls(nlat, nlon, None)


    @classmethod
    def from_zeros(cls, nlat, nlon):
        """
        Returns a :class:`CellGrid` class instance with all array elements
        initialized to zero (``calloc`` in C).

        Parameters
        ----------
        nlat : integer
            Number of grid latitudes
        nlon : integer
            Number of grid longitudes

        Returns
        -------
        out : CellGrid
            :class:`CellGrid` class instance
        """

        return cls(nlat, nlon, 0)


    @classmethod
    def from_arrays(cls, latmin, latmax, lonmin, lonmax, r):
        """
        Returns a :class:`CellGrid` class instance with minimum and maximum
        spherical cell latitudes, minimum and maximum cell longitudes and
        spherical radii taken from ``latmin``, ``latmax``, ``lonmin``,
        ``lonmax`` and ``r``, respectively.

        Parameters
        ----------
        latmin : numpy array of floating points
            Minimum spherical cell latitudes
        latmax : numpy array of floating points
            Maximum spherical cell latitudes
        lonmin : numpy array of floating points
            Minimum spherical cell longitudes
        lonmax : numpy array of floating points
            Maximum spherical cell longitudes
        r : numpy array of floating points
            Spherical radii

        Returns
        -------
        out : CellGrid
            :class:`CellGrid` class instance
        """

        super()._check_data(CELL_GRID, (latmin, latmax, lonmin, lonmax, r))
        return cls(latmin.shape[0], lonmin.shape[0],
                   (latmin, latmax, lonmin, lonmax, r))

