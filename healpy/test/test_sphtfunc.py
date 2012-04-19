import warnings
import pyfits
import os
import numpy as np
import exceptions
from itertools import chain

import unittest

import healpy as hp

class TestSphtFunc(unittest.TestCase):

    def setUp(self):
        self.path = os.path.dirname( os.path.realpath( __file__ ) )
        try:
            self.map1 = [hp.ma(m) for m in hp.read_map(os.path.join(self.path, 'data', 'wmap_band_iqumap_r9_7yr_W_v4.fits'), (0,1,2))]
            self.map2 = [hp.ma(m) for m in hp.read_map(os.path.join(self.path, 'data', 'wmap_band_iqumap_r9_7yr_V_v4.fits'), (0,1,2))]
            self.mask = hp.read_map(os.path.join(self.path, 'data', 'wmap_temperature_analysis_mask_r9_7yr_v4.fits')).astype(np.bool)
        except exceptions.IOError:
            warnings.warn("""Missing Wmap test maps from the data folder, please download them from Lambda and copy them in the test/data folder:
            http://lambda.gsfc.nasa.gov/data/map/dr4/skymaps/7yr/raw/wmap_band_iqumap_r9_7yr_W_v4.fits
            http://lambda.gsfc.nasa.gov/data/map/dr4/skymaps/7yr/raw/wmap_band_iqumap_r9_7yr_V_v4.fits
            http://lambda.gsfc.nasa.gov/data/map/dr4/ancillary/masks/wmap_temperature_analysis_mask_r9_7yr_v4.fits
            on Mac or Linux you can run the bash script get_wmap_maps.sh from the same folder
            """)
            raise
        for m in chain(self.map1, self.map2):
            m.mask = np.logical_not(self.mask)
        self.cla = hp.read_cl(os.path.join(self.path, 'data', 'cl_wmap_fortran.fits'))

    def test_anafast(self):
        cl = hp.anafast(self.map1[0].filled(), lmax = 1024)
        self.assertEqual(len(cl), 1025)
        np.testing.assert_array_almost_equal(cl, self.cla, decimal=8)

    def test_anafast_iqu(self):
        cl = hp.anafast([m.filled() for m in self.map1], lmax = 1024)
        cliqu = np.array(pyfits.open(os.path.join(self.path, 'data', 'cl_iqu_wmap_fortran.fits'))[1].data)
        self.assertEqual(len(cl[0]), 1025)
        self.assertEqual(len(cl), 6)
        for comp in range(6):
            np.testing.assert_array_almost_equal(cl[comp], np.double(cliqu[cliqu.dtype.names[comp]]), decimal=4)

    def test_anafast_xspectra(self):
        cl = hp.anafast(self.map1[0].filled(), self.map2[0].filled(), lmax = 1024)
        self.assertEqual(len(cl), 1025)
        clx = hp.read_cl(os.path.join(self.path, 'data', 'clx.fits'))
        np.testing.assert_array_almost_equal(cl, clx, decimal=8)

    def test_synfast(self):
        nside = 32
        lmax = 64
        fwhm_deg = 7.
        seed = 12345
        np.random.seed(seed)
        map_pregen = hp.read_map(os.path.join(self.path, 'data',
                                              'map_pregen_seed%d.fits' % seed))
        sim_map = hp.synfast(self.cla, nside, lmax = lmax, pixwin=False,
                             fwhm=np.radians(fwhm_deg), new=False, pol=False)
        np.testing.assert_array_almost_equal(sim_map, map_pregen,
                                             decimal=8)
        
if __name__ == '__main__':
    unittest.main()
