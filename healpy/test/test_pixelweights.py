import numpy as np

import unittest

import healpy as hp

import warnings
# disable new order warnings in tests
warnings.filterwarnings('ignore')

class TestMap2Alm(unittest.TestCase):

    def setUp(self):
        self.nside = 64
        self.lmax = 96
        alm_size = 4753
        np.random.seed(123)
        self.input_alm = np.random.normal(size=alm_size) + 1j * np.random.normal(size=alm_size)
        self.m = hp.alm2map(self.input_alm, nside=self.nside, lmax=self.lmax)


    def test_pixelweights(self):
        alm = hp.map2alm(self.m, lmax=self.lmax, use_pixel_weights=True)
        np.testing.assert_allclose(self.input_alm, alm, rtol=1e-4)
