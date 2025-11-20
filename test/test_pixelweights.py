from pathlib import Path

import numpy as np
import pytest

import healpy as hp
from astropy.utils import data

TEST_DATA = Path(__file__).resolve().parent / "data"
PIXEL_WEIGHTS_NSIDE = 64
PIXEL_WEIGHTS_FILE = TEST_DATA / f"healpix_full_weights_nside_{PIXEL_WEIGHTS_NSIDE:04d}.fits"


@pytest.fixture
def create_reference_alm():
    nside = PIXEL_WEIGHTS_NSIDE
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


import unittest.mock
from urllib.error import URLError


def test_map2alm_pixelweights_download_fail(caplog):
    nside = 16
    m = np.arange(hp.nside2npix(nside))

    with unittest.mock.patch('astropy.utils.data.get_pkg_data_filename') as mock_get_pkg_data_filename:
        mock_get_pkg_data_filename.side_effect = URLError('Simulated download error')
        
        alm = hp.map2alm(m, use_pixel_weights=True, lmax=3*nside-1)

        mock_get_pkg_data_filename.assert_called_once()

        assert "Could not download pixel weights" in caplog.text
        assert "Proceeding without pixel weights" in caplog.text

        assert isinstance(alm, np.ndarray)


def test_pixelweights_local_datapath_missing():

    with pytest.raises(RuntimeError):
        hp.map2alm(np.zeros(12), use_pixel_weights=True, datapath="datapath/")


def test_pixelweights_local_datapath(tmp_path, create_reference_alm):
    lmax, input_alm, m = create_reference_alm
    datapath = tmp_path / "datapath" / "full_weights"
    datapath.mkdir(parents=True)
    assert PIXEL_WEIGHTS_FILE.exists(), "bundled pixel weights file missing"
    destination = datapath / f"healpix_full_weights_nside_{PIXEL_WEIGHTS_NSIDE:04d}.fits"
    destination.write_bytes(PIXEL_WEIGHTS_FILE.read_bytes())

    alm = hp.map2alm(
        m, use_pixel_weights=True, datapath=tmp_path / "datapath", lmax=lmax
    )
    np.testing.assert_allclose(input_alm, alm)
