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
        for m in chain(self.map1, self.map2):
            m.mask = np.logical_not(self.mask)
        self.cla = hp.read_cl(os.path.join(self.path, 'data', 'cl_wmap_fortran.fits'))
        cls = pyfits.open(os.path.join(self.path, 'data',
                                       'cl_iqu_wmap_fortran.fits'))[1].data
        # order of HEALPIX is TB, EB while in healpy is EB, TB
        self.cliqu = [cls.field(i) for i in (0,1,2,3,5,4)]

    def test_anafast(self):
        cl = hp.anafast(self.map1[0].filled(), lmax = 1024)
        self.assertEqual(len(cl), 1025)
        np.testing.assert_array_almost_equal(cl, self.cla, decimal=8)

    def test_anafast_iqu(self):
        cl = hp.anafast([m.filled() for m in self.map1], lmax = 1024)
        self.assertEqual(len(cl[0]), 1025)
        self.assertEqual(len(cl), 6)
        for i in range(6):
            np.testing.assert_array_almost_equal(cl[i], self.cliqu[i], decimal=8)

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

    def test_smoothing(self):
        smoothed = hp.ud_grade(hp.smoothing(self.map1[0].data, fwhm=np.radians(1), lmax=1024), 32)
        smoothed_f90 = hp.read_map(os.path.join(self.path, 'data',
                  'wmap_band_imap_r9_7yr_W_v4_1deg_nside32.fits'))
        np.testing.assert_array_almost_equal(smoothed, smoothed_f90)

    def test_smoothing_masked(self):
        smoothed = hp.ud_grade(hp.smoothing(self.map1[0], fwhm=np.radians(1), lmax=1024), 32)
        smoothed_f90 = hp.ma(hp.read_map(os.path.join(self.path, 'data',
                  'wmap_band_imap_r9_7yr_W_v4_1deg_nside32_masked.fits')))
        np.testing.assert_array_almost_equal(smoothed.filled(), smoothed_f90.filled())
        np.testing.assert_array_almost_equal(smoothed, smoothed_f90)
        
    def test_gauss_beam(self):
        idl_gauss_beam = np.array(pyfits.open(os.path.join(self.path, 'data', 'gaussbeam_10arcmin_lmax512_pol.fits'))[0].data).T
        gauss_beam = hp.gauss_beam(np.radians(10./60.), lmax=512, pol=True)
        np.testing.assert_allclose(idl_gauss_beam, gauss_beam)

    def test_map2alm(self):
        nside = 32
        lmax = 64
        fwhm_deg = 7.
        seed = 12345
        np.random.seed(seed)
        orig = hp.synfast(self.cla, nside, lmax=lmax, pixwin=False,
                          fwhm=np.radians(fwhm_deg), new=False)
        tmp = np.empty(orig.size * 2)
        tmp[::2] = orig
        maps = [orig, orig.astype(np.float32), tmp[::2]]
        for input in maps:
            for regression in (False, True):
                alm = hp.map2alm(input, iter=10, regression=regression)
                output = hp.alm2map(alm, nside)
                np.testing.assert_allclose(input, output, atol=1e-4)

    def test_map2alm_pol(self):
        nside = 32
        lmax = 64
        fwhm_deg = 7.
        seed = 12345
        np.random.seed(seed)
        orig = hp.synfast(self.cliqu, nside, lmax=lmax, pixwin=False,
                          fwhm=np.radians(fwhm_deg), new=False)
        tmp = [np.empty(o.size*2) for o in orig]
        for t, o in zip(tmp, orig):
            t[::2] = o
        maps = [orig, [o.astype(np.float32) for o in orig],
                [t[::2] for t in tmp]]
        for input in maps:
            for regression in (False, True):
                alm = hp.map2alm(input, iter=10, regression=regression)
                output = hp.alm2map(alm, nside)
                for i, o in zip(input, output):
                    np.testing.assert_allclose(i, o, atol=1e-4)

if __name__ == '__main__':
    unittest.main()
