from healpy.pixelfunc import *
from healpy._query_disc import boundaries
from healpy._pixelfunc import ringinfo, pix2ring
import exceptions
import numpy as np
import unittest

class TestPixelFunc(unittest.TestCase):
    
    def test_nside2npix(self):
        self.assertEqual(nside2npix(512), 3145728) 
        self.assertEqual(nside2npix(1024), 12582912) 

    def test_nside2resol(self):
        self.assertAlmostEqual(nside2resol(512,arcmin=True), 6.87097282363) 
        self.assertAlmostEqual(nside2resol(1024,arcmin=True), 3.43548641181) 

    def test_nside2pixarea(self):
        self.assertAlmostEqual(nside2pixarea(512), 3.9947416351188569e-06)

    def test_ang2pix_ring(self):
        theta0, phi0 = ([ 1.52911759,  0.78550497,  1.57079633,  0.05103658,  3.09055608], 
                      [ 0.        ,  0.78539816,  1.61988371,  0.78539816,  0.78539816])
        id = ang2pix(1048576 * 8, theta0, phi0, nest=False)
        theta1, phi1 = pix2ang(1048576 * 8, id, nest=False)
        np.testing.assert_array_almost_equal(theta1, theta0)
        np.testing.assert_array_almost_equal(phi1, phi0)

    def test_ang2pix_nest(self):
        theta0, phi0 = ([ 1.52911759,  0.78550497,  1.57079633,  0.05103658,  3.09055608], 
                      [ 0.        ,  0.78539816,  1.61988371,  0.78539816,  0.78539816])
        id = ang2pix(1048576 * 8, theta0, phi0, nest=True)
        theta1, phi1 = pix2ang(1048576 * 8, id, nest=True)
        np.testing.assert_array_almost_equal(theta1, theta0)
        np.testing.assert_array_almost_equal(phi1, phi0)

    def test_ang2pix_negative_theta(self):
        self.assertRaises(exceptions.AssertionError, ang2pix, 32, -1, 0)
      
    def test_fit_dipole(self):
        nside = 32
        npix = nside2npix(nside)
        d = [0.3, 0.5, 0.2]
        vec = np.transpose(pix2vec(nside, np.arange(npix)))
        signal = vec.dot(d)
        mono, dipole = fit_dipole(signal)
        self.assertAlmostEqual(mono, 0.)
        self.assertAlmostEqual(d[0], dipole[0])
        self.assertAlmostEqual(d[1], dipole[1])
        self.assertAlmostEqual(d[2], dipole[2])

    def test_boundaries(self):
        """Test whether the boundary shapes look sane"""
        for lgNside in range(1, 5):
            nside = 1<<lgNside
            for pix in range(nside2npix(nside)):
                for res in range(1, 50, 7):
                    num = 4*res # Expected number of points
                    for nest in (True, False):
                        points = boundaries(nside, pix, res, nest=nest)
                        self.assertTrue(points.shape == (3,num))
                        dist = np.linalg.norm(points[:,:num-1] - points[:,1:]) # distance between points
                        self.assertTrue((dist != 0).all())
                        dmin = np.min(dist)
                        dmax = np.max(dist)
                        self.assertTrue(dmax/dmin <= 2.0)

    def test_ring(self):
        for lgNside in range(1, 5):
            nside = 1<<lgNside
            numPix = nside2npix(nside)
            numRings = 4*nside - 1 # Expected number of rings
            for nest in (True, False):
                pix = np.arange(numPix)
                ring = pix2ring(nside, pix, nest=nest)
                self.assertTrue(pix.shape == ring.shape)
                self.assertTrue(len(set(ring)) == numRings)
                if not nest:
                    first = ring[:numPix-1]
                    second = ring[1:]
                    self.assertTrue(np.logical_or(first == second, first == second-1).all())


if __name__ == '__main__':
    unittest.main()
