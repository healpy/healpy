from healpy.pixelfunc import *
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

if __name__ == '__main__':
    unittest.main()
