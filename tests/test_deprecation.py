import pytest
import numpy as np
from healpy import remove_monopole
from healpy.utils.deprecation import HealpyDeprecationWarning


def test_healpy_deprecation():
    with pytest.warns(HealpyDeprecationWarning):
        remove_monopole(np.zeros(12), verbose=True)
