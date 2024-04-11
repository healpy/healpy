import numpy as np
import pytest
import requests

import healpy as hp
from astropy.utils import data


@pytest.fixture
def create_reference_alm():
    nside = 64
    lmax = 96
    alm_size = 4753
    np.random.seed(123)
    input_alm = np.ones(alm_size, dtype=complex)
    m = hp.alm2map(input_alm, nside=nside, lmax=lmax)
    return lmax, input_alm, m


def test_astropy_download_file():
    data.conf.dataurl = "https://healpy.github.io/healpy-data/"
    print(
        data.get_pkg_data_filename(
            "full_weights/healpix_full_weights_nside_0032.fits", package="healpy"
        )
    )


def test_pixelweights_local_datapath_missing():

    with pytest.raises(RuntimeError):
        hp.map2alm(np.zeros(12), use_pixel_weights=True, datapath="datapath/")


def test_pixelweights_local_datapath(tmp_path, create_reference_alm):
    lmax, input_alm, m = create_reference_alm
    datapath = tmp_path / "datapath" / "full_weights"
    datapath.mkdir(parents=True)
    pixel_weights_file = requests.get(
        "https://github.com/healpy/healpy-data/"
        "blob/master/full_weights/"
        "healpix_full_weights_nside_0064.fits?raw=true"
    )
    with open(datapath / "healpix_full_weights_nside_0064.fits", "wb") as f:
        f.write(pixel_weights_file.content)

    alm = hp.map2alm(
        m, use_pixel_weights=True, datapath=tmp_path / "datapath", lmax=lmax
    )
    np.testing.assert_allclose(input_alm, alm)
