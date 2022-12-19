import sys
from ._constants import _globals


def _get_module_constants(prefix):
    """
    Private function to dynamically add global constants from the `_globals`
    dictionary inside `_constants.py` to the module's namespace.  The variables
    are added under a modified name following the pattern:

        CHARM_CRD_CELLS_GRID --> CELLS_GRID

    where `prefix = 'CHARM_CRD_'` is the input variable.

    Parameters
    ----------
    prefix : str
        Prefix to be removed from the keys in `_globals` before adding to the
        module's namespace.  Should be, e.g., 'CHARM_CRD_', 'CHARM_SHA_', etc.
    """

    namespace = sys._getframe(1).f_globals

    for name, value in _globals.items():
        if name.startswith(prefix):
            namespace[name[len(prefix):]] = value

    return

