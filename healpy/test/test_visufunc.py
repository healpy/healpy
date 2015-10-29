import matplotlib
matplotlib.use("agg")
import unittest
import numpy as np

import healpy as hp
from ..visufunc import *
from ..zoomtool import mollzoom

class TestNoCrash(unittest.TestCase):
    
    def setUp(self):
        self.nside = 1
        self.m = np.arange(hp.nside2npix(self.nside), dtype=np.double)
        self.ma = self.m.copy()
        self.ma[3] = hp.UNSEEN
        self.ma = hp.ma(self.ma)

    def test_cartview_nocrash(self):
        cartview(self.m)

    def test_mollview_nocrash(self):
        mollview(self.m)

    def test_gnomview_nocrash(self):
        gnomview(self.m)

    def test_orthview_nocrash(self):
        orthview(self.m)

    def test_azeqview_nocrash(self):
        azeqview(self.m)

    def test_mollzoom_nocrash(self):
        mollzoom(self.m)

    def test_cartview_ma_nocrash(self):
        cartview(self.ma)

    def test_mollview_ma_nocrash(self):
        mollview(self.ma)

    def test_gnomview_ma_nocrash(self):
        gnomview(self.ma)

    def test_orthview_ma_nocrash(self):
        orthview(self.ma)

    def test_mollzoom_ma_nocrash(self):
        mollzoom(self.ma)

    def test_azeqview_ma_nocrash(self):
        azeqview(self.ma)
