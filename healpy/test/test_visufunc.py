import matplotlib

matplotlib.use("agg")
import unittest
import numpy as np

import healpy as hp
from ..visufunc import *
from ..zoomtool import mollzoom


class TestNoCrash(unittest.TestCase):
    def setUp(self):
        self.nside = 16
        self.npix = hp.nside2npix(self.nside)
        self.m = np.arange(self.npix, dtype=np.double)
        self.m_wrong_npix = np.arange(self.npix + 1, dtype=np.double)
        self.ma = self.m.copy()
        self.ma[3] = hp.UNSEEN
        self.ma = hp.ma(self.ma)
        self.m2 = self.m.copy()
        self.m2[100:] = hp.UNSEEN
        self.ma2 = hp.ma(self.m2)

    def test_cartview_nocrash(self):
        cartview(self.m)

    def test_cartview_crash(self):
        with self.assertRaises(TypeError):
            cartview(self.m_wrong_npix)

    def test_mollview_nocrash(self):
        mollview(self.m)

    def test_mollview_crash(self):
        with self.assertRaises(TypeError):
            mollview(self.m_wrong_npix)

    def test_gnomview_nocrash(self):
        gnomview(self.m)

    def test_gnomview_crash(self):
        with self.assertRaises(TypeError):
            gnomview(self.m_wrong_npix)

    def test_orthview_nocrash(self):
        orthview(self.m)

    def test_orthview_crash(self):
        with self.assertRaises(TypeError):
            orthview(self.m_wrong_npix)

    def test_azeqview_nocrash(self):
        azeqview(self.m)

    def test_azeqview_crash(self):
        with self.assertRaises(TypeError):
            azeqview(self.m_wrong_npix)

    def test_mollzoom_nocrash(self):
        mollzoom(self.m)

    def test_mollzoom_histnocrash(self):
        mollzoom(self.m, norm="hist")

    def test_mollzoom_crash(self):
        with self.assertRaises(TypeError):
            mollzoom(self.m_wrong_npix)

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

    def test_mollzoom_ma_hist_nocrash(self):
        mollzoom(self.ma2, norm="hist")

    def test_azeqview_ma_nocrash(self):
        azeqview(self.ma)
