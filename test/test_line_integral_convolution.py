import os.path
import unittest
import numpy as np

from healpy import line_integral_convolution as lic
from healpy import read_map


class TestLIC(unittest.TestCase):
    def test_qu_shape_equal(self):
        self.assertRaises(ValueError, lic, np.empty(12), np.empty(48))

    def test_texture_shape_equal(self):
        self.assertRaises(
            ValueError, lic, np.empty(12), np.empty(12), texture=np.empty(48)
        )

    def test_qu_shape_dims(self):
        self.assertRaises(ValueError, lic, np.empty((12, 2)), np.empty((12, 2)))

    def test_kernel_steps_fewer_than_steps(self):
        self.assertRaises(
            ValueError, lic, np.empty(12), np.empty(12), kernel_steps=51, steps=50
        )

    def test_lic_no_crash(self):
        lic(np.empty(12), np.empty(12))
        lic(np.empty(12), np.empty(12), ell=-1)
        lic(np.empty(12), np.empty(12), ell=10)
        lic(np.empty(12), np.empty(12), texture=np.empty(12))
        lic(np.empty(12), np.empty(12), modulate=True)

    def test_lic_regression(self):
        path = os.path.dirname(os.path.realpath(__file__))
        Q, U = read_map(
            os.path.join(path, "data", "wmap_band_iqumap_r9_7yr_W_v4_udgraded32.fits"),
            (1, 2), dtype=np.float64
        )
        lic_result = lic(Q, U, step_radian=0.01)
        np.testing.assert_almost_equal(
            np.mean(np.abs(lic_result)), 0.54281382, decimal=8
        )
        np.testing.assert_almost_equal(
            np.std(np.abs(lic_result)), 0.13154804, decimal=8
        )
        lic_result = lic(Q, U, step_radian=0.01, ell=-1)
        np.testing.assert_almost_equal(
            np.mean(np.abs(lic_result)), 0.51307069, decimal=8
        )
        np.testing.assert_almost_equal(
            np.std(np.abs(lic_result)), 0.11628964, decimal=8
        )
