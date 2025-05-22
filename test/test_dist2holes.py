import numpy as np
import healpy as hp
import pytest


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
