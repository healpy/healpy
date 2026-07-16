"""Tests for the :func:`healpy.almxfl` function.

``almxfl`` multiplies every spherical-harmonic coefficient ``a_lm`` by a
factor ``f_l`` that depends only on the multipole ``l``.  The factor is
assumed to be zero where it is not defined (i.e. for ``l >= len(fl)``).

These tests pin the numerical contract of the function using a fixed
random seed, and exercise the supported code paths: real and complex
``fl``, ``fl`` shorter than / longer than ``lmax + 1``, the ``mmax=None``
default, ``inplace`` vs copy semantics and the ``mmax < lmax`` layout.
"""

import numpy as np
import pytest
import healpy as hp


# Fixed seed so the tests are reproducible.
SEED = 12345


def _random_alm(lmax, mmax, rng):
    n_alm = hp.Alm.getsize(lmax, mmax)
    return rng.standard_normal(n_alm) + 1j * rng.standard_normal(n_alm)


def _reference_almxfl(alm, fl, mmax=None):
    """Reference implementation using the old l-outer / m-inner order.

    This is the value the function is expected to return; the optimised
    implementation must reproduce it exactly.
    """
    lmax = hp.Alm.getlmax(alm.size, mmax)
    out = alm.copy()
    for l in range(lmax + 1):
        f = fl[l] if l < fl.size else 0.0
        for m in range(min(l, mmax) + 1):
            i = hp.Alm.getidx(lmax, l, m)
            out[i] *= f
    return out


@pytest.fixture
def rng():
    return np.random.default_rng(SEED)


@pytest.mark.parametrize("lmax,mmax", [(0, 0), (1, 1), (5, 5), (10, 3),
                                       (100, 50), (500, 500)])
def test_almxfl_matches_reference(rng, lmax, mmax):
    """almxfl reproduces the reference l-outer / m-inner computation."""
    alm = _random_alm(lmax, mmax, rng)
    fl = rng.standard_normal(lmax + 1)
    got = hp.almxfl(alm, fl, mmax=mmax)
    exp = _reference_almxfl(alm, fl, mmax)
    np.testing.assert_allclose(got, exp, rtol=1e-12, atol=0)


def test_almxfl_complex_fl(rng):
    """almxfl supports a complex ``fl``."""
    lmax = 64
    mmax = 64
    alm = _random_alm(lmax, mmax, rng)
    fl = rng.standard_normal(lmax + 1) + 1j * rng.standard_normal(lmax + 1)
    got = hp.almxfl(alm, fl, mmax=mmax)
    exp = _reference_almxfl(alm, fl, mmax)
    np.testing.assert_allclose(got, exp, rtol=1e-12, atol=0)


def test_almxfl_mmax_none_default(rng):
    """``mmax=None`` defaults to ``lmax`` (the full layout)."""
    lmax = 32
    mmax = None
    alm = _random_alm(lmax, lmax, rng)
    fl = rng.standard_normal(lmax + 1)
    got = hp.almxfl(alm, fl, mmax=mmax)
    exp = _reference_almxfl(alm, fl, mmax=lmax)
    np.testing.assert_allclose(got, exp, rtol=1e-12, atol=0)


def test_almxfl_short_fl_zeros_tail(rng):
    """``fl`` shorter than ``lmax + 1`` zero-pads the high-l modes."""
    lmax = 10
    mmax = 10
    alm = _random_alm(lmax, mmax, rng)
    fl_short = rng.standard_normal(5)  # defined for l = 0..4
    out = hp.almxfl(alm, fl_short, mmax=mmax)

    for l in range(5, lmax + 1):
        for m in range(min(l, mmax) + 1):
            idx = hp.Alm.getidx(lmax, l, m)
            np.testing.assert_allclose(out[idx], 0.0)

    for l in range(5):
        for m in range(min(l, mmax) + 1):
            idx = hp.Alm.getidx(lmax, l, m)
            np.testing.assert_allclose(out[idx], alm[idx] * fl_short[l])


def test_almxfl_long_fl_ignored_tail(rng):
    """``fl`` longer than ``lmax + 1`` ignores the trailing entries."""
    lmax = 10
    mmax = 10
    alm = _random_alm(lmax, mmax, rng)
    fl = rng.standard_normal(lmax + 1)
    fl_long = np.concatenate([fl, rng.standard_normal(20)])
    a = hp.almxfl(alm, fl, mmax=mmax)
    b = hp.almxfl(alm, fl_long, mmax=mmax)
    np.testing.assert_allclose(a, b, rtol=1e-12, atol=0)


def test_almxfl_inplace_modifies_input(rng):
    """``inplace=True`` modifies and returns the input array."""
    lmax = 32
    mmax = 32
    alm = _random_alm(lmax, mmax, rng)
    alm = np.ascontiguousarray(alm, dtype=np.complex128)
    alm_copy = alm.copy()
    fl = rng.standard_normal(lmax + 1)

    out = hp.almxfl(alm, fl, mmax=mmax, inplace=True)
    assert out is alm
    exp = _reference_almxfl(alm_copy, fl, mmax)
    np.testing.assert_allclose(out, exp, rtol=1e-12, atol=0)


def test_almxfl_default_returns_copy(rng):
    """By default ``almxfl`` returns a copy and leaves the input untouched."""
    lmax = 32
    mmax = 32
    alm = _random_alm(lmax, mmax, rng)
    alm_copy = alm.copy()
    fl = rng.standard_normal(lmax + 1)

    out = hp.almxfl(alm, fl, mmax=mmax)
    assert out is not alm
    np.testing.assert_array_equal(alm, alm_copy)
    exp = _reference_almxfl(alm_copy, fl, mmax)
    np.testing.assert_allclose(out, exp, rtol=1e-12, atol=0)


def test_almxfl_deterministic():
    """Two calls with the same seed return identical results."""
    rng1 = np.random.default_rng(SEED)
    rng2 = np.random.default_rng(SEED)
    lmax = 32
    mmax = 32
    alm1 = _random_alm(lmax, mmax, rng1)
    alm2 = _random_alm(lmax, mmax, rng2)
    fl = rng1.standard_normal(lmax + 1)
    a = hp.almxfl(alm1, fl, mmax=mmax)
    b = hp.almxfl(alm2, fl, mmax=mmax)
    np.testing.assert_array_equal(a, b)
