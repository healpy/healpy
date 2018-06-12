import os
import astropy.io.fits as pf
import unittest
import numpy as np
import gzip

import healpy
from ..fitsfunc import *
from ..sphtfunc import *


class TestFitsFunc(unittest.TestCase):

    def setUp(self):
        self.nside = 512
        self.m = np.arange(healpy.nside2npix(self.nside))
        self.filename = "testmap.fits"

    def test_write_map_IDL(self):
        write_map(self.filename, self.m, fits_IDL=True)
        read_m = pf.open(self.filename)[1].data.field(0)
        self.assertEqual(read_m.ndim, 2)
        self.assertEqual(read_m.shape[1], 1024)
        self.assertTrue(np.all(self.m == read_m.flatten()))

    def test_write_map_units_string(self):
        write_map(self.filename, self.m, column_units="K")
        read_m = pf.open(self.filename)[1].data.field(0)

    def test_write_map_units_list(self):
        write_map(self.filename, [self.m, self.m], column_units=["K", "K"])
        read_m = pf.open(self.filename)[1].data.field(0)

    def test_write_map_C(self):
        write_map(self.filename, self.m, fits_IDL=False)
        read_m = pf.open(self.filename)[1].data.field(0)
        self.assertEqual(read_m.ndim, 1)
        self.assertTrue(np.all(self.m == read_m))

    def test_write_map_C_3comp(self):
        write_map(self.filename, [self.m, self.m, self.m], fits_IDL=False)
        read_m = pf.open(self.filename)[1].data
        for comp in range(3):
            self.assertTrue(np.all(self.m == read_m.field(comp)))

    def test_read_map_filename(self):
        write_map(self.filename, self.m)
        read_map(self.filename)

    def test_read_map_filename_with_header(self):
        write_map(self.filename, self.m)
        m, h = read_map(self.filename, h=True)
        header = dict(h)
        self.assertEqual(header["NSIDE"], 512)

    def test_read_map_hdulist(self):
        write_map(self.filename, self.m)
        hdulist = pf.open(self.filename)
        read_map(hdulist)

    def test_read_map_hdu(self):
        write_map(self.filename, self.m)
        hdu = pf.open(self.filename)[1]
        read_map(hdu)

    def test_read_map_all(self):
        write_map(self.filename, [self.m, self.m, self.m], overwrite=True)
        read_m = read_map(self.filename, None)
        for rm in read_m:
            np.testing.assert_array_almost_equal(self.m, rm)

    def test_read_map_multiple_dtype(self):
        dtypes = [np.int32, np.float32, np.float64]
        m = []
        for dtype in dtypes:
            m.append(self.m.astype(dtype))
        write_map(self.filename, m, overwrite=True)
        read_m = read_map(self.filename, None, dtype=dtypes)
        for rm, dtype in zip(read_m, dtypes):
            self.assertEqual(rm.dtype, dtype)

    def test_read_map_single_dtype(self):
        dtypes = [np.int32, np.float32, np.float64]
        m = []
        for dtype in dtypes:
            m.append(self.m.astype(dtype))
        write_map(self.filename, m, overwrite=True)
        dtype = np.float32
        read_m = read_map(self.filename, None, dtype=dtype)
        self.assertEqual(read_m.dtype, dtype)

    def test_read_map_all_with_header(self):
        write_map(self.filename, [self.m, self.m, self.m], overwrite=True)
        read_m, h = read_map(self.filename, None, h=True)
        header = dict(h)
        self.assertEqual(header["NSIDE"], 512)
        for rm in read_m:
            np.testing.assert_array_almost_equal(self.m, rm)

    def test_read_write_partial(self):
        m = self.m.astype(float)
        m[: 11 * self.nside * self.nside] = UNSEEN
        write_map(self.filename, m, partial=True)
        read_m = read_map(self.filename)
        np.testing.assert_array_almost_equal(m, read_m)

    def test_read_write_partial_3comp(self):
        m = self.m.astype(float)
        m[: 11 * self.nside * self.nside] = UNSEEN
        write_map(self.filename, [m, m, m], partial=True)
        read_m = read_map(self.filename, (0, 1, 2))
        for rm in read_m:
            np.testing.assert_array_almost_equal(m, rm)

    def test_read_write_dtype(self):
        write_map(self.filename, self.m, dtype=np.float64)
        read_m = read_map(self.filename, dtype=np.float32)
        np.testing.assert_almost_equal(np.float32(self.m), read_m)

    def test_read_write_dtype_multiarr(self):
        write_map(self.filename, [self.m, self.m, self.m], dtype=np.float64)
        read_m = read_map(self.filename, field=(0, 1, 2), dtype=np.float32)
        for rm in read_m:
            np.testing.assert_array_almost_equal(np.float32(self.m), rm)

    def test_read_write_iterable_dtype(self):
        write_map(
            self.filename,
            [self.m, self.m, self.m],
            dtype=[np.float64, np.float32, np.int32],
        )
        read_m = read_map(
            self.filename, field=(0, 1, 2), dtype=[np.float32, np.int32, np.float64]
        )
        for rm, dtype in zip(read_m, [np.float32, np.int32, np.float64]):
            np.testing.assert_almost_equal(dtype(self.m), rm)

    def tearDown(self):
        os.remove(self.filename)


class TestFitsFuncGzip(unittest.TestCase):

    def setUp(self):
        self.nside = 4
        self.m = np.arange(healpy.nside2npix(self.nside))
        self.filename = "testmap.fits.gz"

    def test_write_map(self):
        write_map(self.filename, self.m)
        # Make sure file is gzip-compressed
        gzfile = gzip.open(self.filename, "rb")
        gzfile.read()
        gzfile.close()
        read_m = pf.open(self.filename)[1].data.field(0)

    def tearDown(self):
        os.remove(self.filename)


class TestReadWriteAlm(unittest.TestCase):

    def setUp(self):

        s = Alm.getsize(256)
        self.alms = [
            np.arange(s, dtype=np.complex128),
            np.arange(s, dtype=np.complex128),
            np.arange(s, dtype=np.complex128),
        ]

    def tearDown(self):
        if os.path.exists("testalm_128.fits"):
            os.remove("testalm_128.fits")
        if os.path.exists("testalm_256_128.fits"):
            os.remove("testalm_256_128.fits")

    def test_write_alm(self):

        write_alm("testalm_128.fits", self.alms, lmax=128, mmax=128)
        a0 = read_alm("testalm_128.fits")
        # Sanity check of the file
        self.assertEqual(Alm.getlmax(len(a0)), 128)

        # Check the written data
        a0 = read_alm("testalm_128.fits")
        l0, m0 = Alm.getlm(128)
        # We extract 0 <= l <= 128 and 0 <= m <= 128 from self.alms
        idx = Alm.getidx(256, l0, m0)
        np.testing.assert_array_almost_equal(self.alms[0][idx], a0)

    def test_write_alm_256_128(self):
        write_alm("testalm_256_128.fits", self.alms, lmax=256, mmax=128)
        a0, mmax = read_alm("testalm_256_128.fits", return_mmax=True)
        self.assertEqual(mmax, 128)
        self.assertEqual(Alm.getlmax(len(a0), mmax=mmax), 256)

        # Check the written data
        a0 = read_alm("testalm_256_128.fits")

        l0, m0 = Alm.getlm(256)
        idx = Alm.getidx(256, l0, m0)
        # Extract 0 <= l <= 256 and 0 <= m <= 128
        idx_mmax = np.where(m0 <= mmax)
        idx = idx[idx_mmax]
        np.testing.assert_array_almost_equal(self.alms[0][idx], a0)

    def test_read_alm_filename(self):
        write_alm("testalm_128.fits", self.alms, lmax=128, mmax=128)
        read_alm("testalm_128.fits")

    def test_read_alm_filename_array(self):
        write_alm("testalm_256.fits", self.alms, overwrite=True)
        testalm1 = np.array(self.alms)
        testalm2 = read_alm("testalm_256.fits", [1, 2, 3])
        np.testing.assert_array_almost_equal(testalm1, testalm2)

    def test_read_alm_hdulist(self):
        write_alm("testalm_128.fits", self.alms, lmax=128, mmax=128)
        hdulist = pf.open("testalm_128.fits")
        read_alm(hdulist)

    def test_read_alm_hdu(self):
        write_alm("testalm_128.fits", self.alms, lmax=128, mmax=128)
        hdu = pf.open("testalm_128.fits")[1]
        read_alm(hdu)


class TestReadWriteCl(unittest.TestCase):

    def tearDown(self):
        os.remove("test_cl.fits")

    def test_write_read_cl_II(self):
        cl = np.arange(1025, dtype=np.double)
        write_cl("test_cl.fits", cl)
        cl_read = read_cl("test_cl.fits")
        np.testing.assert_array_almost_equal(cl, cl_read)

    def test_write_read_cl_4comp(self):
        cl = [np.arange(1025, dtype=np.double) for n in range(4)]
        write_cl("test_cl.fits", cl)
        cl_read = read_cl("test_cl.fits")
        for cl_column, cl_read_column in zip(cl, cl_read):
            np.testing.assert_array_almost_equal(cl_column, cl_read_column)

    def test_write_read_cl_6comp(self):
        cl = [np.arange(1025, dtype=np.double) for n in range(6)]
        write_cl("test_cl.fits", cl)
        cl_read = read_cl("test_cl.fits")
        for cl_column, cl_read_column in zip(cl, cl_read):
            np.testing.assert_array_almost_equal(cl_column, cl_read_column)

    def test_read_cl_filename(self):
        cl = np.arange(1025, dtype=np.double)
        write_cl("test_cl.fits", cl)
        read_cl("test_cl.fits")

    def test_read_cl_hdulist(self):
        cl = np.arange(1025, dtype=np.double)
        write_cl("test_cl.fits", cl)
        hdulist = pf.open("test_cl.fits")
        read_cl(hdulist)

    def test_read_cl_hdu(self):
        cl = np.arange(1025, dtype=np.double)
        write_cl("test_cl.fits", cl)
        hdu = pf.open("test_cl.fits")[1]
        read_cl(hdu)


if __name__ == "__main__":
    unittest.main()
