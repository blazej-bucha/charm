"""
Defines default encoding for strings passed to CHarm and a routine returing
ctype pointer to a Python string.
"""

# Default encoding
_default_encoding = 'utf-8'

# Default treatment of characters outside "_default_encoding" when calling the
# "decode" method on strings of the "bytes" class
_default_encoding_errors = 'replace'


import ctypes as _ct


def _str_ptr(string, encoding=_default_encoding):
    """
    Returns a ctypes pointer to `string` encoded with `encoding`.  If `string`
    is `None`, returned is `None`.

    Before entering this function, `string` should be checked whether it is
    indeed a string or `None`.
    """

    if string is None:
        return None
    else:
        return _ct.create_string_buffer(string.encode(encoding))


def _bytes_decode(string,
                  encoding=_default_encoding,
                  errors=_default_encoding_errors):
    """
    Returns `str` by encoding `string` of class `bytes` using
    `_default_encoding` and treating encoding errors by
    `_default_encoding_errors`.
    """

    return string.decode(encoding=_default_encoding,
                         errors=_default_encoding_errors)
