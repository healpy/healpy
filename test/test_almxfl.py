"""Tests for almxfl function."""
import numpy as np
import pytest
import healpy as hp


class TestAlmxfl:
    """Test suite for almxfl function."""

    def test_basic_multiplication(self):
        """Test basic multiplication of alm by fl."""
        lmax = 10
        mmax = 10
        n_alm = hp.Alm.getsize(lmax, mmax)
        alm = np.random.randn(n_alm) + 1j * np.random.randn(n_alm)
        fl = np.random.randn(lmax + 1)

        alm_out = hp.almxfl(alm, fl, mmax=mmax, inplace=False)

        # Verify a few elements
        for test_l in [0, 5, 10]:
            for test_m in range(min(test_l, mmax) + 1):
                idx = hp.Alm.getidx(lmax, test_l, test_m)
                expected = alm[idx] * fl[test_l]
                np.testing.assert_allclose(alm_out[idx], expected)

    def test_fl_shorter_than_lmax(self):
        """Test that fl shorter than lmax+1 pads with zeros."""
        lmax = 10
        mmax = 10
        n_alm = hp.Alm.getsize(lmax, mmax)
        alm = np.random.randn(n_alm) + 1j * np.random.randn(n_alm)
        fl_short = np.random.randn(5)  # Only goes up to l=4

        alm_out = hp.almxfl(alm, fl_short, mmax=mmax, inplace=False)

        # Check that l >= 5 are zeroed
        for test_l in range(5, lmax + 1):
            for test_m in range(min(test_l, mmax) + 1):
                idx = hp.Alm.getidx(lmax, test_l, test_m)
                np.testing.assert_allclose(alm_out[idx], 0.0)

        # Check that l < 5 are multiplied correctly
        for test_l in range(5):
            for test_m in range(min(test_l, mmax) + 1):
                idx = hp.Alm.getidx(lmax, test_l, test_m)
                expected = alm[idx] * fl_short[test_l]
                np.testing.assert_allclose(alm_out[idx], expected)

    def test_inplace_modification(self):
        """Test inplace modification."""
        lmax = 10
        mmax = 10
        n_alm = hp.Alm.getsize(lmax, mmax)
        alm = np.ascontiguousarray(
            np.random.randn(n_alm) + 1j * np.random.randn(n_alm),
            dtype=np.complex128
        )
        alm_copy = alm.copy()
        fl = np.random.randn(lmax + 1)

        alm_out = hp.almxfl(alm, fl, mmax=mmax, inplace=True)

        # Check that the modification was done in place
        assert alm_out is alm

        # Verify the result is correct
        for test_l in [0, 5, 10]:
            for test_m in range(min(test_l, mmax) + 1):
                idx = hp.Alm.getidx(lmax, test_l, test_m)
                expected = alm_copy[idx] * fl[test_l]
                np.testing.assert_allclose(alm_out[idx], expected)

    def test_not_inplace_modification(self):
        """Test that inplace=False creates a copy."""
        lmax = 10
        mmax = 10
        n_alm = hp.Alm.getsize(lmax, mmax)
        alm = np.random.randn(n_alm) + 1j * np.random.randn(n_alm)
        alm_copy = alm.copy()
        fl = np.random.randn(lmax + 1)

        alm_out = hp.almxfl(alm, fl, mmax=mmax, inplace=False)

        # Check that original is unchanged
        np.testing.assert_array_equal(alm, alm_copy)

        # Check that output is different object
        assert alm_out is not alm

    def test_different_lmax_mmax(self):
        """Test with lmax != mmax."""
        lmax = 100
        mmax = 50
        n_alm = hp.Alm.getsize(lmax, mmax)
        alm = np.random.randn(n_alm) + 1j * np.random.randn(n_alm)
        fl = np.random.randn(lmax + 1)

        alm_out = hp.almxfl(alm, fl, mmax=mmax, inplace=False)

        # Verify a few elements
        for test_l in [0, 25, 50, 75, 100]:
            for test_m in range(min(test_l, mmax) + 1):
                idx = hp.Alm.getidx(lmax, test_l, test_m)
                expected = alm[idx] * fl[test_l]
                np.testing.assert_allclose(alm_out[idx], expected)

    def test_large_dataset(self):
        """Test with larger dataset to ensure performance."""
        lmax = 500
        mmax = 500
        n_alm = hp.Alm.getsize(lmax, mmax)
        alm = np.random.randn(n_alm) + 1j * np.random.randn(n_alm)
        fl = np.random.randn(lmax + 1)

        alm_out = hp.almxfl(alm, fl, mmax=mmax, inplace=False)

        # Spot check a few elements
        for test_l in [0, 100, 250, 500]:
            for test_m in [0, min(test_l, min(50, mmax))]:
                idx = hp.Alm.getidx(lmax, test_l, test_m)
                expected = alm[idx] * fl[test_l]
                np.testing.assert_allclose(alm_out[idx], expected)

    def test_zeros_fl(self):
        """Test with fl containing zeros."""
        lmax = 10
        mmax = 10
        n_alm = hp.Alm.getsize(lmax, mmax)
        alm = np.random.randn(n_alm) + 1j * np.random.randn(n_alm)
        fl = np.zeros(lmax + 1)
        fl[0] = 1.0  # Only l=0 is non-zero

        alm_out = hp.almxfl(alm, fl, mmax=mmax, inplace=False)

        # Check that only m=0, l=0 is non-zero
        idx0 = hp.Alm.getidx(lmax, 0, 0)
        np.testing.assert_allclose(alm_out[idx0], alm[idx0])

        # Check all others are zero
        for test_l in range(1, lmax + 1):
            for test_m in range(min(test_l, mmax) + 1):
                idx = hp.Alm.getidx(lmax, test_l, test_m)
                np.testing.assert_allclose(alm_out[idx], 0.0)

    def test_complex_fl(self):
        """Test with complex fl values."""
        lmax = 10
        mmax = 10
        n_alm = hp.Alm.getsize(lmax, mmax)
        alm = np.random.randn(n_alm) + 1j * np.random.randn(n_alm)
        fl = np.random.randn(lmax + 1) + 1j * np.random.randn(lmax + 1)

        alm_out = hp.almxfl(alm, fl, mmax=mmax, inplace=False)

        # Verify a few elements
        for test_l in [0, 5, 10]:
            for test_m in range(min(test_l, mmax) + 1):
                idx = hp.Alm.getidx(lmax, test_l, test_m)
                expected = alm[idx] * fl[test_l]
                np.testing.assert_allclose(alm_out[idx], expected)
