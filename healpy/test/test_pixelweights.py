import numpy as np

import unittest

import healpy as hp

import warnings
# disable new order warnings in tests
warnings.filterwarnings('ignore')

class TestMap2Alm(unittest.TestCase):

    def setUp(self):
        self.lmax = 96
        self.nside = 32
        self.input_alm = np.ones(4753, dtype=np.complex)
        self.m = hp.alm2map(self.input_alm, nside=self.nside)


    def test_pixelweights(self):
        alm = hp.map2alm(self.m, lmax=self.lmax)
        np.testing.assert_allclose(self.input_alm, alm)
