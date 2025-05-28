import numpy as np
import healpy as hp
import pytest


def test_fill_small_holes_pixel_threshold():
    nside = 8
    npix = hp.nside2npix(nside)
    mask = np.ones(npix, dtype=np.float64)
    # Create two holes: one small (single pixel), one large (10 pixels)
    mask[10] = 0
    mask[20:30] = 0
    filled = hp._masktools.fill_small_holes(mask, nside, min_size=2)
    # The single-pixel hole should be filled
    assert filled[10] == 1
    # The large hole should remain
    assert np.all(filled[20:30] == 0)
    # All other pixels should be unchanged
    assert np.all(filled[(mask == 1)] == 1)


def test_fill_small_holes_area_threshold():
    nside = 8
    npix = hp.nside2npix(nside)
    mask = np.ones(npix, dtype=np.float64)
    mask[42] = 0
    area_per_pix = 4 * np.pi / npix * (180 * 60 / np.pi) ** 2
    filled = hp._masktools.fill_small_holes(
        mask, nside, min_area_arcmin2=area_per_pix * 1.1
    )
    # The hole should be filled
    assert filled[42] == 1


def test_fill_small_holes_no_fill():
    nside = 8
    npix = hp.nside2npix(nside)
    mask = np.ones(npix, dtype=np.float64)
    mask[100:105] = 0
    # Threshold too small, nothing should be filled
    filled = hp._masktools.fill_small_holes(mask, nside, min_size=1)
    assert np.all(filled[100:105] == 0)
    # No change elsewhere
    assert np.all(filled[(mask == 1)] == 1)


def test_dist2holes_hole_min_size():
    nside = 8
    npix = hp.nside2npix(nside)
    mask = np.ones(npix, dtype=np.float64)
    # Create two holes: one small (single pixel), one large (10 pixels)
    mask[10] = 0
    mask[20:30] = 0
    # With hole_min_size=2, the single-pixel hole should be filled
    d1 = hp.dist2holes(mask, hole_min_size=2)
    # The pixel at 10 should now be treated as valid (filled)
    assert d1[10] > 0  # Should not be zero (not a hole anymore)
    # The large hole should remain
    assert np.any(d1[20:30] == 0)


def test_dist2holes_hole_min_surf_arcmin2():
    nside = 8
    npix = hp.nside2npix(nside)
    mask = np.ones(npix, dtype=np.float64)
    # Create a small hole (single pixel)
    mask[42] = 0
    # Compute area per pixel
    area_per_pix = 4 * np.pi / npix * (180 * 60 / np.pi) ** 2
    # Set threshold just above one pixel
    d2 = hp.dist2holes(mask, hole_min_surf_arcmin2=area_per_pix * 1.1)
    # The hole should be filled
    assert d2[42] > 0


def test_dist2holes_no_hole_filter():
    nside = 8
    npix = hp.nside2npix(nside)
    mask = np.ones(npix, dtype=np.float64)
    mask[100] = 0
    d = hp.dist2holes(mask)
    assert d[100] == 0
    assert np.all(d >= 0)


def test_dist2holes_all_valid():
    nside = 8
    npix = hp.nside2npix(nside)
    mask = np.ones(npix, dtype=np.float64)
    d = hp.dist2holes(mask)
    # All distances should be zero (no holes)
    assert np.all(d == 0)


def test_dist2holes_all_invalid():
    nside = 8
    npix = hp.nside2npix(nside)
    mask = np.zeros(npix, dtype=np.float64)
    d = hp.dist2holes(mask)
    # All distances should be zero (all holes)
    assert np.all(d == 0)
