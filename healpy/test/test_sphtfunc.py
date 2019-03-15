import astropy.io.fits as pf
import os
import numpy as np
from copy import deepcopy
from itertools import chain

import unittest

import healpy as hp

import warnings

# disable new order warnings in tests
warnings.filterwarnings("ignore")


class TestSphtFunc(unittest.TestCase):
    def setUp(self):
        self.lmax = 64
        self.path = os.path.dirname(os.path.realpath(__file__))
        self.map1 = [
            hp.ma(m)
            for m in hp.read_map(
                os.path.join(
                    self.path, "data", "wmap_band_iqumap_r9_7yr_W_v4_udgraded32.fits"
                ),
                (0, 1, 2),
            )
        ]
        self.map2 = [
            hp.ma(m)
            for m in hp.read_map(
                os.path.join(
                    self.path, "data", "wmap_band_iqumap_r9_7yr_V_v4_udgraded32.fits"
                ),
                (0, 1, 2),
            )
        ]
        self.mask = hp.read_map(
            os.path.join(
                self.path,
                "data",
                "wmap_temperature_analysis_mask_r9_7yr_v4_udgraded32.fits",
            )
        ).astype(np.bool)
        for m in chain(self.map1, self.map2):
            m.mask = np.logical_not(self.mask)
        self.cla = hp.read_cl(
            os.path.join(
                self.path,
                "data",
                "cl_wmap_band_iqumap_r9_7yr_W_v4_udgraded32_II_lmax64_rmmono_3iter.fits",
            )
        )
        self.cl_fortran_nomask = hp.read_cl(
            os.path.join(
                self.path,
                "data",
                "cl_wmap_band_iqumap_r9_7yr_W_v4_udgraded32_II_lmax64_rmmono_3iter_nomask.fits",
            )
        )
        cls_file = pf.open(
            os.path.join(
                self.path,
                "data",
                "cl_wmap_band_iqumap_r9_7yr_W_v4_udgraded32_IQU_lmax64_rmmono_3iter.fits",
            )
        )
        # fix for pyfits to read the file with duplicate column names
        for i in range(2, 6):
            cls_file[1].header["TTYPE%d" % i] += "-%d" % i
        cls = cls_file[1].data
        # order of HEALPIX is TB, EB while in healpy is EB, TB
        self.cliqu = [np.array(cls.field(i)) for i in (0, 1, 2, 3, 5, 4)]
        nside = 32
        lmax = 64
        fwhm_deg = 7.0
        seed = 12345
        np.random.seed(seed)
        self.mapiqu = hp.synfast(
            self.cliqu,
            nside,
            lmax=lmax,
            pixwin=False,
            fwhm=np.radians(fwhm_deg),
            new=False,
        )

    def test_anafast(self):
        cl = hp.anafast(hp.remove_monopole(self.map1[0].filled()), lmax=self.lmax)
        self.assertEqual(len(cl), 65)
        np.testing.assert_array_almost_equal(cl, self.cla, decimal=8)

    def test_anafast_nomask(self):
        cl = hp.anafast(hp.remove_monopole(self.map1[0].data), lmax=self.lmax)
        self.assertEqual(len(cl), 65)
        np.testing.assert_array_almost_equal(cl, self.cl_fortran_nomask, decimal=8)

    def test_anafast_iqu(self):
        self.map1[0] = hp.remove_monopole(self.map1[0])
        cl = hp.anafast(self.map1, lmax=self.lmax)
        self.assertEqual(len(cl[0]), 65)
        self.assertEqual(len(cl), 6)
        for i in range(6):
            np.testing.assert_array_almost_equal(cl[i], self.cliqu[i], decimal=8)

    def test_anafast_xspectra(self):
        cl = hp.anafast(
            hp.remove_monopole(self.map1[0]),
            hp.remove_monopole(self.map2[0]),
            lmax=self.lmax,
        )
        self.assertEqual(len(cl), self.lmax + 1)
        clx = hp.read_cl(
            os.path.join(
                self.path,
                "data",
                "cl_wmap_band_iqumap_r9_7yr_WVxspec_v4_udgraded32_II_lmax64_rmmono_3iter.fits",
            )
        )
        np.testing.assert_array_almost_equal(cl, clx, decimal=8)

    def test_synfast(self):
        nside = 32
        lmax = 64
        fwhm_deg = 7.0
        seed = 12345
        np.random.seed(seed)
        map_pregen = hp.read_map(
            os.path.join(self.path, "data", "map_synfast_seed%d.fits" % seed), (0, 1, 2)
        )
        sim_map = hp.synfast(
            self.cliqu,
            nside,
            lmax=lmax,
            pixwin=False,
            fwhm=np.radians(fwhm_deg),
            new=False,
            pol=True,
        )
        np.testing.assert_array_almost_equal(sim_map, map_pregen, decimal=8)

    def test_smoothing_notmasked(self):
        smoothed = hp.smoothing(
            [m.data for m in self.map1], fwhm=np.radians(10), lmax=self.lmax
        )
        smoothed_f90 = hp.read_map(
            os.path.join(
                self.path,
                "data",
                "wmap_band_iqumap_r9_7yr_W_v4_udgraded32_smoothed10deg_fortran.fits",
            ),
            (0, 1, 2),
        )
        np.testing.assert_array_almost_equal(smoothed, smoothed_f90, decimal=6)

    def test_smoothing_masked(self):
        smoothed = hp.smoothing(self.map1, fwhm=np.radians(10), lmax=self.lmax)
        smoothed_f90 = hp.ma(
            hp.read_map(
                os.path.join(
                    self.path,
                    "data",
                    "wmap_band_iqumap_r9_7yr_W_v4_udgraded32_masked_smoothed10deg_fortran.fits",
                ),
                (0, 1, 2),
            )
        )
        # fortran does not restore the mask
        smoothed_f90.mask = smoothed.mask
        np.testing.assert_array_almost_equal(
            smoothed.filled(), smoothed_f90.filled(), decimal=6
        )

    def test_gauss_beam(self):
        idl_gauss_beam = np.array(
            pf.open(
                os.path.join(self.path, "data", "gaussbeam_10arcmin_lmax512_pol.fits")
            )[0].data
        ).T
        gauss_beam = hp.gauss_beam(np.radians(10.0 / 60.0), lmax=512, pol=True)
        np.testing.assert_allclose(idl_gauss_beam, gauss_beam)

    def test_alm2cl(self):
        nside = 32
        lmax = 64
        lmax_out = 100
        seed = 12345
        np.random.seed(seed)

        # Input power spectrum and alm
        alm_syn = hp.synalm(self.cla, lmax=lmax)

        cl_out = hp.alm2cl(alm_syn, lmax_out=lmax_out - 1)

        np.testing.assert_array_almost_equal(cl_out, self.cla[:lmax_out], decimal=4)

    def test_map2alm(self):
        nside = 32
        lmax = 64
        fwhm_deg = 7.0
        seed = 12345
        np.random.seed(seed)
        orig = hp.synfast(
            self.cla,
            nside,
            lmax=lmax,
            pixwin=False,
            fwhm=np.radians(fwhm_deg),
            new=False,
        )
        tmp = np.empty(orig.size * 2)
        tmp[::2] = orig
        maps = [orig, orig.astype(np.float32), tmp[::2]]
        for use_weights in [False, True]:
            for input in maps:
                alm = hp.map2alm(input, iter=10, use_weights=use_weights)
                output = hp.alm2map(alm, nside)
                np.testing.assert_allclose(input, output, atol=1e-4)

    def test_map2alm_pol(self):
        tmp = [np.empty(o.size * 2) for o in self.mapiqu]
        for t, o in zip(tmp, self.mapiqu):
            t[::2] = o
        maps = [
            self.mapiqu,
            [o.astype(np.float32) for o in self.mapiqu],
            [t[::2] for t in tmp],
        ]
        for use_weights in [False, True]:
            for input in maps:
                alm = hp.map2alm(input, iter=10, use_weights=use_weights)
                output = hp.alm2map(alm, 32)
                for i, o in zip(input, output):
                    np.testing.assert_allclose(i, o, atol=1e-4)

    def test_map2alm_pol_gal_cut(self):
        tmp = [np.empty(o.size * 2) for o in self.mapiqu]
        for t, o in zip(tmp, self.mapiqu):
            t[::2] = o
        maps = [
            self.mapiqu,
            [o.astype(np.float32) for o in self.mapiqu],
            [t[::2] for t in tmp],
        ]
        for use_weights in [False, True]:
            for input in maps:
                gal_cut = 30
                nside = hp.get_nside(input)
                npix = hp.nside2npix(nside)
                gal_mask = (
                    np.abs(hp.pix2ang(nside, np.arange(npix), lonlat=True)[1]) < gal_cut
                )
                alm = hp.map2alm(
                    input, iter=10, use_weights=use_weights, gal_cut=gal_cut
                )
                output = hp.alm2map(alm, 32)
                for i, o in zip(input, output):
                    # Testing requires low tolerances because of the
                    # mask boundary
                    i[gal_mask] = 0
                    np.testing.assert_allclose(i, o, atol=1e-2)

    def test_rotate_alm(self):
        almigc = hp.map2alm(self.mapiqu)
        alms = [almigc[0], almigc[0:2], almigc, np.vstack(almigc)]
        for i in alms:
            o = deepcopy(i)
            hp.rotate_alm(o, 0.1, 0.2, 0.3)
            hp.rotate_alm(o, -0.3, -0.2, -0.1)
            # FIXME: rtol=1e-6 works here, except on Debian with Python 3.4.
            np.testing.assert_allclose(i, o, rtol=1e-5)

    def test_rotate_alm_rotmatrix(self):
        """rotate_alm also support rotation matrix instead of angles"""
        lmax = 32
        nalm = hp.Alm.getsize(lmax)
        alm = np.zeros([3, nalm], dtype=np.complex)
        alm[0, 1] = 1
        alm[1, 2] = 1
        alm_rotated_angles = alm.copy()
        angles = hp.rotator.coordsys2euler_zyz(coord=["G", "E"])
        hp.rotate_alm(alm_rotated_angles, *angles)
        gal2ecl = hp.Rotator(coord=["G", "E"])
        hp.rotate_alm(alm, matrix=gal2ecl.mat)
        np.testing.assert_allclose(alm_rotated_angles, alm)

    def test_rotate_alm2(self):
        # Test rotate_alm against the Fortran library
        lmax = 64
        nalm = hp.Alm.getsize(lmax)
        alm = np.zeros([3, nalm], dtype=np.complex)
        for i in range(3):
            for ell in range(lmax + 1):
                for m in range(ell):
                    ind = hp.Alm.getidx(lmax, ell, m)
                    alm[i, ind] = (i + 1) * 10 + ell + 1j * m
        psi = np.pi / 3.0
        theta = 0.5
        phi = 0.01
        hp.rotate_alm(alm, psi, theta, phi)

        ref_0_0_0 = 0.00000000000 + 0.00000000000j
        ref_0_21_0 = -64.0056622444 + 0.00000000000j
        ref_0_21_21 = -3.19617408364 + 2.00219590117j
        ref_0_42_0 = 87.8201360825 + 0.00000000000j
        ref_0_42_21 = -6.57242309702 + 50.1128079361j
        ref_0_42_42 = 0.792592362074 - 0.928452597766j
        ref_0_63_0 = -49.6732554742 + 0.00000000000j
        ref_0_63_21 = -51.2812623888 - 61.6289129316j
        ref_0_63_42 = -9.32823219430 + 79.0787993482j
        ref_0_63_63 = -0.157204566965 + 0.324692958700j
        ref_1_0_0 = 0.00000000000 + 0.00000000000j
        ref_1_21_0 = -85.5520809077 + 0.00000000000j
        ref_1_21_21 = -3.57384285749 + 2.93255811219j
        ref_1_42_0 = 107.541172254 + 0.00000000000j
        ref_1_42_21 = -2.77944941833 + 57.1015322415j
        ref_1_42_42 = 0.794212854046 - 1.10982745343j
        ref_1_63_0 = -60.7153303746 + 0.00000000000j
        ref_1_63_21 = -61.0915123767 - 65.9943878923j
        ref_1_63_42 = -4.86354653261 + 86.5277253196j
        ref_1_63_63 = -0.147165377786 + 0.360474777237j

        ref = np.array(
            [
                ref_0_0_0,
                ref_0_21_0,
                ref_0_21_21,
                ref_0_42_0,
                ref_0_42_21,
                ref_0_42_42,
                ref_0_63_0,
                ref_0_63_21,
                ref_0_63_42,
                ref_0_63_63,
                ref_1_0_0,
                ref_1_21_0,
                ref_1_21_21,
                ref_1_42_0,
                ref_1_42_21,
                ref_1_42_42,
                ref_1_63_0,
                ref_1_63_21,
                ref_1_63_42,
                ref_1_63_63,
            ]
        )

        mine = []
        for i in [0, 1]:
            for ell in range(0, lmax + 1, 21):
                for m in range(0, ell + 1, 21):
                    ind = hp.Alm.getidx(lmax, ell, m)
                    mine.append(alm[i, ind])

        mine = np.array(ref)

        np.testing.assert_allclose(ref, mine, rtol=1e-10)

    def test_accept_ma_allows_only_keywords(self):
        """ Test whether 'smoothing' wrapped with accept_ma works with only
            keyword arguments. """

        ma = np.ones(12 * 16 ** 2)
        try:
            hp.smoothing(map_in=ma)
        except IndexError:
            self.fail()

    def test_beam2bl(self):
        """ Test beam2bl against analytical transform of Gaussian beam. """

        theta = np.linspace(0, np.radians(1.0), 1000)
        sigma = np.radians(10.0 / 60.0) / np.sqrt(8.0 * np.log(2.0))
        gaussian_beam = np.exp(-0.5 * (theta / sigma) ** 2) / (2 * np.pi * sigma ** 2)

        ell = np.arange(512 + 1.0)
        gaussian_window = np.exp(-0.5 * ell * (ell + 1) * sigma ** 2)

        bl = hp.beam2bl(gaussian_beam, theta, 512)
        np.testing.assert_allclose(gaussian_window, bl, rtol=1e-5)

    def test_bl2beam(self):
        """ Test bl2beam against analytical transform of Gaussian beam. """

        theta = np.linspace(0, np.radians(3.0), 1000)
        sigma = np.radians(1.0) / np.sqrt(8.0 * np.log(2.0))
        gaussian_beam = np.exp(-0.5 * (theta / sigma) ** 2) / (2 * np.pi * sigma ** 2)

        ell = np.arange(2048 + 1.0)
        gaussian_window = np.exp(-0.5 * ell * (ell + 1) * sigma ** 2)

        beam = hp.bl2beam(gaussian_window, theta)
        np.testing.assert_allclose(gaussian_beam, beam, rtol=1e-3)

    def test_max_nside_check(self):
        """ Test whether the max_nside_check correctly raises ValueErrors for nsides
        that are too large."""

        # Test an nside that is too large
        with self.assertRaises(ValueError):
            hp.check_max_nside(16384)

        # Test an nside that is valid
        # hp.check_max_nside will return 0 if no exceptions are raised
        self.assertEqual(hp.check_max_nside(1024), 0)

    def test_pixwin_base(self):
        # Base case
        nsides = [2 ** p for p in np.arange(1, 14)]
        [hp.pixwin(nside) for nside in nsides]

        # Test invalid nside
        with self.assertRaises(ValueError):
            hp.pixwin(15)

    def test_pixwin_pol(self):
        pixwin = hp.pixwin(128, pol=True)
        self.assertEqual(len(pixwin), 2)

    def test_pixwin_lmax(self):
        nside = 128
        pixwin = hp.pixwin(nside, lmax=None)
        self.assertEqual(len(pixwin), 3 * nside)

        lmax = 200
        pixwin = hp.pixwin(nside, lmax=lmax)
        self.assertEqual(len(pixwin)-1, lmax)


if __name__ == "__main__":
    unittest.main()
