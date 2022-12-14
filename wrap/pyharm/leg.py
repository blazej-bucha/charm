"""
Module to work with the fully-normalized associated Legendre functions of the 
first kind:

    * defines the :class:`pyharm.leg.Pnmj` class to store Fourier coefficients 
      of the Legendre functions,

    * computes Fourier coefficients of Legendre functions.

Note
----
This documentation is written for double precision version of PyHarm.
"""

import ctypes as _ct
import numpy as _np
from . import _libcharm, _libcharmname, _CHARM, _pyharm
from ._get_module_constants import _get_module_constants
from ._check_types import _check_deg_ord, _check_int_scalar, _check_pointer
from ._get_empty_array import _get_empty_array
from ._constants import _globals
from ._data_types import _ct_ulong, _ct_int, _ct_size_t, _ct_flt, _pyharm_flt
from . import _err as _ph_err


# Get the module constants from "_constants.py" and add them to the module's
# namespace
_get_module_constants('CHARM_LEG_')


# All possible ordering schemes of Fourier coefficients in "Pnmj" class
# instances
_ORDERING = [PNMJ_ORDER_MNJ, PNMJ_ORDER_MJN]


#: int: Ordering scheme of Fourier coefficients of Legendre functions: order
#: :obj:`m`, degree :obj:`n` and wave-number-related variable :obj:`j`.
#: For further details, refer to
#: `https://blazej-bucha.github.io/charm/build/html/api-c-leg.html
#: <https://blazej-bucha.github.io/charm/build/html/api-c-leg.html>`_.
PNMJ_ORDER_MNJ: int


#: int: Ordering scheme of Fourier coefficients of Legendre functions: order
#: :obj:`m`, wave-number-related variable :obj:`j` and degree :obj:`n`.
#: For further details, refer to
#: `https://blazej-bucha.github.io/charm/build/html/api-c-leg.html
#: <https://blazej-bucha.github.io/charm/build/html/api-c-leg.html>`_.
PNMJ_ORDER_MJN: int


class _Pnmj(_ct.Structure):
    """
    Private class to represent the `charm_leg` structure of CHarm.
    """

    _fields_ = [('nmax', _ct_ulong),
                ('ordering', _ct_int),
                ('npnmj', _ct_size_t),
                ('pnmj', _ct.POINTER(_ct.POINTER(_ct.POINTER(_ct_flt))))]


class Pnmj:
    """
    Class for Fourier coefficients of fully-normalized associated Legendre
    functions of the first kind.

    To create a :class:`Pnmj` class instance, always use one of the following
    factory methods:

        * :meth:`from_garbage`,
        * :meth:`from_zeros`.

    Examples
    --------
    >>> import pyharm
    >>> pnmj = pyharm.leg.Pnmj.from_garbage(10)

    Parameters
    ----------
    nmax : integer
        Maximum harmonic degree
    ordering : integer
        Ordering scheme of Fourier coefficients.  Use 
        :meth:`get_ordering_types` to get all supported ordering schemes.
    coeffs : None or :obj:`0`
        Determines the way of initializing Fourier coefficients:

            * :obj:`None` to not initialize Fourier coefficients (`malloc` in
              C),

            * :obj:`0` to set all Fourier coefficients to zero (`calloc` in C).

    Note
    ----
    Once a :class:`Pnmj` class instance is created, its :obj:`nmax`,
    :attr:`ordering` and :attr:`pnmj` attributes are not
    writable.

    For details, refer to
    `https://blazej-bucha.github.io/charm/build/html/api-c-leg.html
    <https://blazej-bucha.github.io/charm/build/html/api-c-leg.html>`_.

    """


    @property
    def nmax(self):
        """ Maximum harmonic degree of the Fourier coefficients. """
        return self._nmax


    @property
    def ordering(self):
        """
        Ordering scheme of the Fourier coefficients.  Use the
        :meth:`get_ordering_str` method to get all supported values.
        """
        return self._ordering


    @property
    def pnmj(self):
        """ Fourier coefficients. """
        return self._pnmj


    def __init__(self, nmax, ordering, coeffs):

        _check_deg_ord(nmax, 'degree')

        _check_int_scalar(ordering, 'ordering')
        if ordering not in _ORDERING:
            raise ValueError('Unsupported value of \'ordering\'.\n\n%s' % \
                             Pnmj._get_ordering_types())

        f = ''
        if coeffs == None:
            f = _CHARM + 'leg_pnmj_malloc'
        elif coeffs == 0:
            f = _CHARM + 'leg_pnmj_calloc'
        else:
            raise ValueError('Unsupported method \'%s\' to create a '
                             '\'Pnmj\' class instance.' % coeffs)

        self._nmax        = None
        self._ordering    = None
        self._npnmj       = None
        self._pnmj        = _get_empty_array()
        self._from_method = None
        self._Pnmj        = None

        func          = _libcharm[f]
        func.restype  = _ct.POINTER(_Pnmj)
        func.argtypes = [_ct_ulong, _ct_int]

        self._Pnmj = func(_ct_ulong(nmax), _ct_int(ordering))
        _check_pointer(self._Pnmj, f, _libcharmname)

        self._Pnmj2Pnmj()
        self._from_method = coeffs

        return


    def __str__(self):

        ordering = self.get_ordering_str()
        out  = f'nmax = {self.nmax}\n\n'
        out += f'ordering = {self.ordering}\n\n'
        out += f'pnmj = {self.pnmj}\n'

        return out


    def __repr__(self):

        return f'{_pyharm}.leg.Pnmj({self.nmax}, {self.ordering}, ' \
               f'\'{self._from_method}\')'


    def __del__(self):

        self._free()

        return


    def __exit__(self):

        self._free()

        return


    @classmethod
    def from_garbage(cls, nmax, ordering=PNMJ_ORDER_MNJ):
        """
        Returns a :class:`Pnmj` class instance with uninitialized Fourier
        coefficients (`malloc` in C).

        Parameters
        ----------
        nmax : integer
            Maximum harmonic degree
        ordering : integer
            Ordering scheme of Fourier coefficients, optional.  Default value
            is :obj:`PNMJ_ORDER_MNJ`.  Use :meth:`get_ordering_types` to
            get all supported ordering schemes.

        Returns
        -------
        out : Pnmj
            A :class:`Pnmj` class instance
        """

        return cls(nmax, ordering, None)


    @classmethod
    def from_zeros(cls, nmax, ordering=PNMJ_ORDER_MNJ):
        """
        Returns a :class:`Pnmj` class instance with all Fourier coefficients
        initialized to zero (`calloc` in C).

        Parameters
        ----------
        nmax : integer
            Maximum harmonic degree
        ordering : integer
            Ordering scheme of Fourier coefficients, optional.  Default value
            is :obj:`PNMJ_ORDER_MNJ`.  Use :meth:`get_ordering_types` to
            get all supported ordering schemes.

        Returns
        -------
        out : Pnmj
            A :class:`Pnmj` class instance
        """

        return cls(nmax, ordering, 0)


    def coeffs(self, nmax=None):
        """
        Computes Fourier coefficients of Legendre functions.

        Parameters
        ----------
        nmax : integer
            Maximum harmonic degree, optional.  If not provided, the object's
            :obj:`nmax` attribute is used.
        """

        if nmax is not None:
            _check_deg_ord(nmax, 'degree')
            if nmax > self.nmax:
                msg  = f'Couldn\'t compute Fourier coefficients up to degree '
                msg += f'{nmax}, because the \'Pnmj\' class instance is '
                msg += f'initialized only up to degree {self.nmax}.'
                raise ValueError(msg)
        else:
            nmax = self.nmax

        func         = _libcharm[_CHARM + 'leg_pnmj_coeffs']
        func.restype = None
        func.argtype = [_ct.POINTER(_Pnmj),
                        _ct_ulong,
                        _ct.POINTER(_ph_err._Err)]

        err = _ph_err.init()
        func(self._Pnmj, _ct_ulong(nmax), err)
        _ph_err.handler(err, 1)
        _ph_err.free(err)
        self._Pnmj2Pnmj()

        return


    def get_coeff(self, n, m, j):
        """
        Returns Fourier coefficient of degree :obj:`n`, order :obj:`m` and
        wave-number-related value :obj:`j`.

        Parameters
        ----------
        n : integer
            Spherical harmonic degree
        m : integer
            Spherical harmonic order
        j : integer
            Value related to the wave-number of Fourier coefficient (refer to
            :meth:`j2k` and :meth:`k2j` for further details)

        Returns
        -------
        out : floating point
            Fourier coefficient of degree :obj:`n`, order :obj:`m` and
            wave-number-related value :obj:`j`.
        """

        _check_deg_ord(n, 'degree')
        _check_deg_ord(m, 'order')
        Pnmj._check_j(j)
        if n > self.nmax:
            msg  = f'Couldn\'t get Fourier coefficient of degree {n}, because '
            msg += f'the object is initialized only up to degree {self.nmax}.'
            raise ValueError(msg)
        if m > self.nmax:
            msg  = f'Couldn\'t get Fourier coefficient of order {m}, because '
            msg += f'the object is initialized only up to degree {self.nmax}.'
            raise ValueError(msg)
        if n < m:
            msg  = f'Harmonic degree \'n = {n}\' cannot be smaller than '
            msg += f'harmonic order \'m = {m}\'.'
            raise ValueError(msg)


        if self.ordering == PNMJ_ORDER_MNJ:
            return _pyharm_flt(self._Pnmj.contents.pnmj[m][n - m][j])
        elif self.ordering == PNMJ_ORDER_MJN:
            return _pyharm_flt(\
                            self._Pnmj.contents.pnmj[m][j][n - max(m, 2 * j)])


    def get_ordering_str(self):
        """
        Transforms the object's :attr:`ordering` attribute to a pretty string.

        Returns
        -------
        out : str
            Pretty name of the :attr:`ordering` attribute
        """

        return _get_pnmj_ordering_str(self.ordering)


    @staticmethod
    def j2k(n, j):
        """
        Transforms a wave-number-related variable :obj:`j` to the wave-number
        :obj:`k` of a Fourier coefficient of fully-normalized associated
        Legendre function of degree :obj:`n`.

        For further details, see the references at
        `https://blazej-bucha.github.io/charm/build/html/api-c-leg.html
        <https://blazej-bucha.github.io/charm/build/html/api-c-leg.html>`_.

        Parameters
        ----------
        n : integer
            Spherical harmonic degree of the Legendre function
        j : integer
            Variable related to the wave-number :obj:`k`

        Returns
        -------
        out : integer
            Wave-number :obj:`k` of the Fourier coefficient of a Legendre function
        """

        _check_deg_ord(n, 'degree')
        Pnmj._check_j(j)

        func          = _libcharm[_CHARM + 'leg_pnmj_j2k']
        func.restype  = _ct_ulong
        func.argtypes = [_ct_ulong, _ct_ulong]

        return func(_ct_ulong(n), _ct_ulong(j))


    @staticmethod
    def k2j(k):
        """
        Transforms a wave-number :obj:`k` of a Fourier coefficient of
        fully-normalized associated Legendre functions to the
        wave-number-related variable :obj:`j`.

        For further details, see the references at
        `https://blazej-bucha.github.io/charm/build/html/api-c-leg.html
        <https://blazej-bucha.github.io/charm/build/html/api-c-leg.html>`_.

        Parameters
        ----------
        k : integer
            Wave-number

        Returns
        -------
        out : integer
            Wave-number-related variable :obj:`j` of the Fourier coefficient
            of a Legendre function
        """

        _check_int_scalar(k, 'k')
        if k < 0:
            raise ValueError('The wave-number \'k\' cannot be negative.')

        func          = _libcharm[_CHARM + 'leg_pnmj_k2j']
        func.restype  = _ct_ulong
        func.argtypes = [_ct_ulong]

        return func(_ct_ulong(k))


    @staticmethod
    def get_ordering_types():
        """
        Prints all names of symbolic constants that are accepted as the
        :attr:`ordering` attribute of the :class:`Pnmj` class.
        """

        print(Pnmj._get_ordering_types())

        return


    def _Pnmj2Pnmj(self):
        """
        Private function to transform the :class:`_Pnmj` class instance in
        `self._Pnmj.contents` to a :class:`Pnmj` class instance in `self`.  The
        :attr:`pnmj` attribute of `self` shares the same memory space as the
        corresponding attribute `self._Pnmj.contents.pnmj`.
        """

        if not isinstance(self._Pnmj.contents, _Pnmj):
            raise TypeError('\'self._Pnmj.contents\' must be a \'_Pnmj\' '
                            'class instance.')

        self._nmax     = int(self._Pnmj.contents.nmax)
        self._npnmj    = int(self._Pnmj.contents.npnmj)
        self._ordering = int(self._Pnmj.contents.ordering)

        if not self._Pnmj.contents.pnmj:
            raise ValueError('\'self._Pnmj.contents.pnmj\' is a \'NULL\' '
                             'pointer.')
        self._pnmj = _np.ctypeslib.as_array(self._Pnmj.contents.pnmj[0][0],
                                       shape=(int(self._Pnmj.contents.npnmj),))

        return


    def _free(self):
        """
        Frees the memory associated with the object.
        """

        if self._Pnmj is not None:
            func          = _libcharm[_CHARM + 'leg_pnmj_free']
            func.restype  = None
            func.argtypes = [_ct.POINTER(_Pnmj)]
            func(self._Pnmj)

        return


    @staticmethod
    def _get_ordering_types():
        """
        Returns a help string on all supported values of the :attr:`ordering`
        attribute of the :class:`Pnmj` class.
        """

        ret  = f'The following symbolic constants are accepted as the '
        ret += f'\'ordering\'\nattribute of the \'Pnmj\' class:\n'
        ret += f'\n'
        for ordering in _ORDERING:
            ret += '\t%s' % _get_pnmj_ordering_str(ordering)
            ret += '\n'
        ret += f'\n'
        ret += f'For further details, visit '
        ret += f'\'`https://blazej-bucha.github.io/charm/build/html/'
        ret += f'api-c-leg.html\'.'

        return ret


    @staticmethod
    def _check_j(j):
        """
        Private method to check whether the input :obj:`j` value can represent
        the wave-number-related variable.

        Parameters
        ----------
        j : integer
            Anything
        """

        _check_int_scalar(j, 'j')
        if j < 0:
            raise ValueError('The wave-number-related variable \'j\' must '
                             'not be negative.')


def fourier_coeffs(nmax, ordering=PNMJ_ORDER_MNJ):
    """
    Computes Fourier coefficients of Legendre functions up to degree
    :obj:`nmax` using the :obj:`ordering` ordering scheme of coefficients.

    Parameters
    ----------
    nmax : integer
        Maximum harmonic degree
    ordering : integer
        Ordering scheme of Fourier coefficients, optional.  Default value is
        :obj:`PNMJ_ORDER_MNJ`.  Use :meth:`Pnmj.get_ordering_types()` to get 
        all supported ordering schemes.

    Returns
    -------
    out : Pnmj
        Fourier coefficients
    """

    pnmj = Pnmj.from_garbage(nmax, ordering)
    pnmj.coeffs()

    return pnmj


def _get_pnmj_ordering_str(ordering):
    """
    Private function to transform :attr:`ordering` attribute of :class:`Pnmj`
    to a pretty string.

    Returns
    -------
    out : str
        Pretty name of the :attr:`ordering` attribute
    """

    if ordering == PNMJ_ORDER_MNJ:
        return f'{_pyharm}.leg.PNMJ_ORDER_MNJ'
    elif ordering == PNMJ_ORDER_MJN:
        return f'{_pyharm}.leg.PNMJ_ORDER_MJN'
    elif ordering is None:
        return None
