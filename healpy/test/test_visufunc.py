import os
import pyfits
import unittest
import numpy as np

import healpy
from healpy.visufunc import *

class TestNoCrash(unittest.TestCase):
    
    def setUp(self):
        self.nside = 1
        self.m = np.arange(healpy.nside2npix(self.nside))

    def test_cartview_nocrash(self):
        cartview(self.m)

    def test_mollview_nocrash(self):
        mollview(self.m)

    def test_gnomview_nocrash(self):
        gnomview(self.m)
