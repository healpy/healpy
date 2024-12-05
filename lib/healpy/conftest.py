import numpy as np


def pytest_configure(config):
    np.set_printoptions(legacy="1.13")

    config.option.astropy_header = True
    try:
        from pytest_astropy_header.display import PYTEST_HEADER_MODULES

        PYTEST_HEADER_MODULES.pop("Pandas", None)
        PYTEST_HEADER_MODULES.pop("h5py", None)
        PYTEST_HEADER_MODULES["astropy"] = "astropy"
        PYTEST_HEADER_MODULES["cython"] = "cython"
    except ImportError:
        pass
