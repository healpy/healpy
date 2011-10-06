import warnings
import os
import numpy as np
import exceptions

import unittest

import healpy as hp

class TestSphtFunc(unittest.TestCase):

    def setUp(self):
        try:
            self.map = hp.ma(hp.read_map(os.path.join('data', 'wmap_band_imap_r9_7yr_W_v4.fits')))
            self.mask = hp.read_map(os.path.join('data', 'wmap_temperature_analysis_mask_r9_7yr_v4.fits')).astype(np.bool)
        except exceptions.IOError:
            warnings.warn("""Missing Wmap test maps from the data folder, please download them from Lambda and copy them in the test/data folder:
            http://lambda.gsfc.nasa.gov/data/map/dr4/skymaps/7yr/raw/wmap_band_imap_r9_7yr_W_v4.fits
            http://lambda.gsfc.nasa.gov/data/map/dr4/ancillary/masks/wmap_temperature_analysis_mask_r9_7yr_v4.fits
            on Mac or Linux you can run the bash script get_wmap_maps.sh from the same folder
            """)
            raise
        self.map.mask = np.logical_not(self.mask)
        self.cla = hp.read_cl(os.path.join('data', 'cl_wmap_fortran.fits'))
    
    def test_anafast(self):
        cl = hp.anafast(self.map.filled(), lmax = 1024)
        np.testing.assert_array_almost_equal(cl, self.cla, decimal=8)

if __name__ == '__main__':
    unittest.main()
