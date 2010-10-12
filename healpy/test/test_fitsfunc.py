import os
import pyfits
import unittest
import numpy as np

import healpy
from healpy.fitsfunc import *

class TestFitsFunc(unittest.TestCase):
    
    def setUp(self):
        self.nside = 512
        self.m = np.arange(healpy.nside2npix(self.nside))
        self.filename = 'testmap.fits'

    def test_write_map_IDL(self):
        write_map(self.filename, self.m, fits_IDL=True)
        read_m = pyfits.open(self.filename)[1].data.field(0)
        self.assertEqual(read_m.ndim, 2)
        self.assertEqual(read_m.shape[1], 1024)
        self.assertTrue(np.all(self.m == read_m.flatten()))

    def test_write_map_C(self):
        write_map(self.filename, self.m, fits_IDL=False)
        read_m = pyfits.open(self.filename)[1].data.field(0)
        self.assertEqual(read_m.ndim, 1)
        self.assertTrue(np.all(self.m == read_m))

    def tearDown(self):
        os.remove(self.filename)


if __name__ == '__main__':
    unittest.main()
