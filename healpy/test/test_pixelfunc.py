from ..pixelfunc import *
from .._query_disc import boundaries
from .._pixelfunc import ringinfo, pix2ring
import numpy as np
import unittest


class TestPixelFunc(unittest.TestCase):

    def setUp(self):
        # data fixture
        self.theta0 = [1.52911759, 0.78550497, 1.57079633, 0.05103658, 3.09055608]
        self.phi0 = [0., 0.78539816, 1.61988371, 0.78539816, 0.78539816]
        self.lon0 = np.degrees(self.phi0)
        self.lat0 = 90. - np.degrees(self.theta0)

    def test_nside2npix(self):
        self.assertEqual(nside2npix(512), 3145728)
        self.assertEqual(nside2npix(1024), 12582912)

    def test_nside2resol(self):
        self.assertAlmostEqual(nside2resol(512, arcmin=True), 6.87097282363)
        self.assertAlmostEqual(nside2resol(1024, arcmin=True), 3.43548641181)

    def test_max_pixrad(self):
        self.assertAlmostEqual(max_pixrad(512), 2.0870552355e-03)
        self.assertAlmostEqual(
            max_pixrad(512, degrees=True), np.rad2deg(2.0870552355e-03)
        )

    def test_nside2pixarea(self):
        self.assertAlmostEqual(nside2pixarea(512), 3.9947416351188569e-06)

    def test_ang2pix_ring(self):
        # ensure nside = 1 << 23 is correctly calculated
        # by comparing the original theta phi are restored.
        # NOTE: nside needs to be sufficiently large!
        id = ang2pix(1048576 * 8, self.theta0, self.phi0, nest=False)
        theta1, phi1 = pix2ang(1048576 * 8, id, nest=False)
        np.testing.assert_array_almost_equal(theta1, self.theta0)
        np.testing.assert_array_almost_equal(phi1, self.phi0)

    def test_ang2pix_ring_outofrange(self):
        # Healpy_Base2 works up to nside = 2**29.
        # Check that a ValueError is raised for nside = 2**30.
        self.assertRaises(
            ValueError, ang2pix, 1 << 30, self.theta0, self.phi0, nest=False
        )

    def test_ang2pix_nest(self):
        # ensure nside = 1 << 23 is correctly calculated
        # by comparing the original theta phi are restored.
        # NOTE: nside needs to be sufficiently large!
        # NOTE: with Healpy_Base this will fail because nside
        #       is limited to 1 << 13 with Healpy_Base.
        id = ang2pix(1048576 * 8, self.theta0, self.phi0, nest=True)
        theta1, phi1 = pix2ang(1048576 * 8, id, nest=True)
        np.testing.assert_array_almost_equal(theta1, self.theta0)
        np.testing.assert_array_almost_equal(phi1, self.phi0)

        self.assertTrue(np.allclose(theta1, self.theta0))
        self.assertTrue(np.allclose(phi1, self.phi0))

    def test_ang2pix_nest_outofrange_doesntcrash(self):
        # Healpy_Base2 works up to nside = 2**29.
        # Check that a ValueError is raised for nside = 2**30.
        self.assertRaises(
            ValueError, ang2pix, 1 << 30, self.theta0, self.phi0, nest=False
        )

    def test_ang2pix_negative_theta(self):
        self.assertRaises(ValueError, ang2pix, 32, -1, 0)

    def test_ang2pix_lonlat(self):
        # Need to decrease the precision of the check because deg not radians
        id = ang2pix(1048576 * 8, self.lon0, self.lat0, nest=False, lonlat=True)
        lon1, lat1 = pix2ang(1048576 * 8, id, nest=False, lonlat=True)
        np.testing.assert_array_almost_equal(lon1, self.lon0, decimal=5)
        np.testing.assert_array_almost_equal(lat1, self.lat0, decimal=5)

        # Now test nested
        id = ang2pix(1048576 * 8, self.theta0, self.phi0, nest=True)
        theta1, phi1 = pix2ang(1048576 * 8, id, nest=True)
        np.testing.assert_array_almost_equal(theta1, self.theta0)
        np.testing.assert_array_almost_equal(phi1, self.phi0)

    def test_vec2pix_lonlat(self):
        # Need to decrease the precision of the check because deg not radians
        vec = ang2vec(self.lon0, self.lat0, lonlat=True)
        lon1, lat1 = vec2ang(vec, lonlat=True)
        np.testing.assert_array_almost_equal(lon1, self.lon0, decimal=5)
        np.testing.assert_array_almost_equal(lat1, self.lat0, decimal=5)

    def test_get_interp_val_lonlat(self):
        m = np.arange(12.)
        val0 = get_interp_val(m, self.theta0, self.phi0)
        val1 = get_interp_val(m, self.lon0, self.lat0, lonlat=True)
        np.testing.assert_array_almost_equal(val0, val1)

    def test_get_interp_weights(self):
        p0, w0 = (np.array([0, 1, 4, 5]), np.array([1., 0., 0., 0.]))

        # phi not specified, theta assumed to be pixel
        p1, w1 = get_interp_weights(1, 0)
        np.testing.assert_array_almost_equal(p0, p1)
        np.testing.assert_array_almost_equal(w0, w1)

        # If phi is not specified, lonlat should do nothing
        p1, w1 = get_interp_weights(1, 0, lonlat=True)
        np.testing.assert_array_almost_equal(p0, p1)
        np.testing.assert_array_almost_equal(w0, w1)

        p0, w0 = (np.array([1, 2, 3, 0]), np.array([0.25, 0.25, 0.25, 0.25]))

        p1, w1 = get_interp_weights(1, 0, 0)
        np.testing.assert_array_almost_equal(p0, p1)
        np.testing.assert_array_almost_equal(w0, w1)

        p1, w1 = get_interp_weights(1, 0, 90, lonlat=True)
        np.testing.assert_array_almost_equal(p0, p1)
        np.testing.assert_array_almost_equal(w0, w1)

    def test_get_all_neighbours(self):
        ipix0 = np.array([8, 4, 0, -1, 1, 6, 9, -1])
        ipix1 = get_all_neighbours(1, np.pi / 2, np.pi / 2)
        ipix2 = get_all_neighbours(1, 90, 0, lonlat=True)
        np.testing.assert_array_almost_equal(ipix0, ipix1)
        np.testing.assert_array_almost_equal(ipix0, ipix2)

    def test_fit_dipole(self):
        nside = 32
        npix = nside2npix(nside)
        d = [0.3, 0.5, 0.2]
        vec = np.transpose(pix2vec(nside, np.arange(npix)))
        signal = np.dot(vec, d)
        mono, dipole = fit_dipole(signal)
        self.assertAlmostEqual(mono, 0.)
        self.assertAlmostEqual(d[0], dipole[0])
        self.assertAlmostEqual(d[1], dipole[1])
        self.assertAlmostEqual(d[2], dipole[2])

    def test_boundaries(self):
        """Test whether the boundary shapes look sane"""
        for lgNside in range(1, 5):
            nside = 1 << lgNside
            for pix in range(nside2npix(nside)):
                for res in range(1, 50, 7):
                    num = 4 * res  # Expected number of points
                    for nest in (True, False):
                        points = boundaries(nside, pix, res, nest=nest)
                        self.assertTrue(points.shape == (3, num))
                        dist = np.linalg.norm(
                            points[:, : num - 1] - points[:, 1:]
                        )  # distance between points
                        self.assertTrue((dist != 0).all())
                        dmin = np.min(dist)
                        dmax = np.max(dist)
                        self.assertTrue(dmax / dmin <= 2.0)

    def test_ring(self):
        for lgNside in range(1, 5):
            nside = 1 << lgNside
            numPix = nside2npix(nside)
            numRings = 4 * nside - 1  # Expected number of rings
            for nest in (True, False):
                pix = np.arange(numPix)
                ring = pix2ring(nside, pix, nest=nest)
                self.assertTrue(pix.shape == ring.shape)
                self.assertTrue(len(set(ring)) == numRings)
                if not nest:
                    first = ring[: numPix - 1]
                    second = ring[1:]
                    self.assertTrue(
                        np.logical_or(first == second, first == second - 1).all()
                    )

    def test_accept_ma_allows_only_keywords(self):
        """ Test whether the accept_ma wrapper accepts calls using only keywords."""
        ma = np.zeros(12 * 16 ** 2)
        try:
            ud_grade(map_in=ma, nside_out=32)
        except IndexError:
            self.fail("IndexError raised")


if __name__ == "__main__":
    unittest.main()
