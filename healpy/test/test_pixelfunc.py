from healpy.pixelfunc import *
import numpy
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
        self.assertTrue(numpy.allclose(theta1, theta0))
        self.assertTrue(numpy.allclose(phi1, phi0))

        id = ang2pix(1<<30, theta0, phi0, nest=False)
        theta1, phi1 = pix2ang(1<<30, id, nest=False)
        self.assertFalse(numpy.allclose(theta1, theta0))
        self.assertFalse(numpy.allclose(phi1, phi0))

    def test_ang2pix_nest(self):
        theta0, phi0 = ([ 1.52911759,  0.78550497,  1.57079633,  0.05103658,  3.09055608], 
                      [ 0.        ,  0.78539816,  1.61988371,  0.78539816,  0.78539816])
        id = ang2pix(1048576 * 8, theta0, phi0, nest=True)
        theta1, phi1 = pix2ang(1048576 * 8, id, nest=True)
        self.assertTrue(numpy.allclose(theta1, theta0))
        self.assertTrue(numpy.allclose(phi1, phi0))

        id = ang2pix(1<<30, theta0, phi0, nest=True)
        theta1, phi1 = pix2ang(1<<30, id, nest=True)
        self.assertFalse(numpy.allclose(theta1, theta0))
        self.assertFalse(numpy.allclose(phi1, phi0))
      
if __name__ == '__main__':
    unittest.main()
