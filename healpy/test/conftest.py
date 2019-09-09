import numpy as np
import pytest


def pytest_configure():
    """Set the Numpy print style to a fixed version to make doctest outputs
    reproducible."""
    try:
        np.set_printoptions(legacy="1.13")
    except TypeError:
        # On older versions of Numpy, the unrecognized 'legacy' option will
        # raise a TypeError.
        pass
