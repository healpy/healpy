import numpy as np

import unittest

import healpy as hp

import warnings
# disable new order warnings in tests
warnings.filterwarnings('ignore')

class TestMap2Alm(unittest.TestCase):

    def setUp(self):
        self.nside = 32
        self.lmax = 32
        self.input_alm = np.ones(561, dtype=np.complex)
        self.m = hp.alm2map(self.input_alm, nside=self.nside, lmax=self.lmax)


    def test_pixelweights(self):
        alm = hp.map2alm(self.m, lmax=self.lmax, use_pixel_weights=True)
        np.testing.assert_allclose(self.input_alm, alm, rtol=1e-4)
