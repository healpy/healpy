import numpy as np


def pytest_configure():
    """Set the Numpy print style to a fixed version to make doctest outputs
    reproducible."""
    np.set_printoptions(legacy="1.13")
