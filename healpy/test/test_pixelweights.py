import numpy as np

import unittest

import healpy as hp
from astropy.utils import data

import warnings

# disable new order warnings in tests
warnings.filterwarnings("ignore")


class TestMap2Alm(unittest.TestCase):

    def setUp(self):
        self.nside = 64
        self.lmax = 96
        alm_size = 4753
        np.random.seed(123)
        self.input_alm = np.ones(alm_size, dtype=np.complex)
        self.m = hp.alm2map(self.input_alm, nside=self.nside, lmax=self.lmax)

    def test_astropy_download_file(self):
        data.conf.dataurl = "https://healpy.github.io/healpy-data/"
        print(
            data.get_pkg_data_filename(
                "full_weights/healpix_full_weights_nside_0032.fits", package="healpy"
            )
        )


if __name__ == "__main__":
    unittest.main()
