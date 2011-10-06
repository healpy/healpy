from healpy.pixelfunc import *

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

if __name__ == '__main__':
    unittest.main()
