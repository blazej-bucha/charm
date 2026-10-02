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
    Returns a ctypes pointer to `string` encoded with `encoding`.
    """

    if isinstance(string, str):
        return _ct.create_string_buffer(string.encode(encoding))
    elif string is None:
        return None
    else:
        raise ValueError('\'string\' must be of \'str\' type or \'None\'.')


def _bytes_decode(string,
                  encoding=_default_encoding,
                  errors=_default_encoding_errors):
    """
    Returns `str` by encoding `string` of class `bytes` using
    `_default_encoding` and treating encoding errors by
    `_default_encoding_errors`.
    """

    if isinstance(string, bytes):
        return string.decode(encoding=_default_encoding,
                             errors=_default_encoding_errors)
    elif string is None:
        return None
    else:
        raise ValueError('\'string\' must be of \'bytes\' type.')
