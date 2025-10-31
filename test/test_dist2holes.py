import numpy as np
import healpy as hp
import pytest


def make_small_hole_mask(nside):
    npix = hp.nside2npix(nside)
    mask = np.zeros(npix, dtype=bool)
    theta_small, phi_small = np.deg2rad(45.0), 0.0
    pix_small = hp.ang2pix(nside, theta_small, phi_small)
    mask[pix_small] = True
    return mask, pix_small


def make_mask_with_small_and_large_holes(nside):
    mask, pix_small = make_small_hole_mask(nside)
    vec_large = hp.ang2vec(np.deg2rad(120.0), np.deg2rad(60.0))
    large_pix = hp.query_disc(nside, vec_large, np.deg2rad(10.0))
    mask[large_pix] = True
    return mask, pix_small, large_pix


def test_fill_small_holes_pixel_threshold():
    nside = 8
    mask, pix_small, large_pix = make_mask_with_small_and_large_holes(nside)
    filled = hp._masktools.fill_small_holes(mask, nside, min_size=2)
    # The single-pixel hole should be filled
    assert not filled[pix_small]
    # The large hole should remain
    assert np.all(filled[large_pix])
    # All other pixels should be unchanged (remain False)
    assert not filled[~mask].any()


def test_fill_small_holes_area_threshold():
    nside = 8
    mask, pix_small = make_small_hole_mask(nside)
    npix = mask.size
    area_per_pix = 4 * np.pi / npix * (180 * 60 / np.pi) ** 2
    filled = hp._masktools.fill_small_holes(
        mask, nside, min_area_arcmin2=area_per_pix * 1.1
    )
    # The hole should be filled
    assert not filled[pix_small]


def test_fill_small_holes_no_fill():
    nside = 8
    mask, pix_small, large_pix = make_mask_with_small_and_large_holes(nside)
    # Threshold too small, nothing should be filled
    filled = hp._masktools.fill_small_holes(mask, nside, min_size=1)
    np.testing.assert_array_equal(filled, mask)


def test_dist2holes_hole_min_size():
    nside = 8
    mask, pix_small, large_pix = make_mask_with_small_and_large_holes(nside)
    # With hole_min_size=2, the single-pixel hole should be filled
    d1 = hp.dist2holes(mask, hole_min_size=2)
    # The pixel at 10 should now be treated as valid (filled)
    assert d1[pix_small] > 0  # Should not be zero (not a hole anymore)
    # The large hole should remain
    assert np.any(d1[large_pix] == 0)


def test_dist2holes_hole_min_surf_arcmin2():
    nside = 8
    mask, pix_small = make_small_hole_mask(nside)
    npix = mask.size
    # Compute area per pixel
    area_per_pix = 4 * np.pi / npix * (180 * 60 / np.pi) ** 2
    # Set threshold just above one pixel
    d2 = hp.dist2holes(mask, hole_min_surf_arcmin2=area_per_pix * 1.1)
    # The hole should be filled
    assert d2[pix_small] > 0


def test_dist2holes_no_hole_filter():
    nside = 8
    npix = hp.nside2npix(nside)
    mask = np.zeros(npix, dtype=bool)
    mask[100] = True
    d = hp.dist2holes(mask)
    assert d[100] == 0
    assert np.all(d >= 0)


def test_dist2holes_all_valid():
    nside = 8
    npix = hp.nside2npix(nside)
    mask = np.zeros(npix, dtype=bool)
    d = hp.dist2holes(mask)
    # With no holes the distances saturate at maxdist (default pi)
    assert np.allclose(d, np.pi)


def test_dist2holes_all_invalid():
    nside = 8
    npix = hp.nside2npix(nside)
    mask = np.ones(npix, dtype=bool)
    d = hp.dist2holes(mask)
    # All distances should be zero (all holes)
    assert np.all(d == 0)


def test_fill_small_holes_returns_copy_and_dtype():
    nside = 8
    mask, _, _ = make_mask_with_small_and_large_holes(nside)
    filled = hp._masktools.fill_small_holes(mask, nside, min_size=2)
    assert filled.dtype == mask.dtype
    assert not np.shares_memory(mask, filled)


@pytest.mark.parametrize(
    "kwargs, exc",
    [
        ({"min_size": -1}, ValueError),
        ({"min_area_arcmin2": -0.5}, ValueError),
        ({"min_size": 1.5}, TypeError),
    ],
)
def test_fill_small_holes_rejects_invalid_thresholds(kwargs, exc):
    nside = 8
    mask, *_ = make_mask_with_small_and_large_holes(nside)
    with pytest.raises(exc):
        hp._masktools.fill_small_holes(mask, nside, **kwargs)


def test_fill_small_holes_checks_nside():
    nside = 8
    mask, *_ = make_mask_with_small_and_large_holes(nside)
    with pytest.raises(ValueError):
        hp._masktools.fill_small_holes(mask, nside + 1, min_size=2)


def test_dist2holes_rejects_invalid_thresholds():
    nside = 8
    mask, *_ = make_mask_with_small_and_large_holes(nside)
    with pytest.raises(ValueError):
        hp.dist2holes(mask, hole_min_size=-1)
    with pytest.raises(TypeError):
        hp.dist2holes(mask, hole_min_size=1.5)
    with pytest.raises(ValueError):
        hp.dist2holes(mask, hole_min_surf_arcmin2=-1.0)


def test_dist2holes_combined_filters():
    nside = 8
    mask, pix_small, large_pix = make_mask_with_small_and_large_holes(nside)
    # Thresholds are generous individually but together still fill the small hole
    area_per_pix = 4 * np.pi / mask.size * (180 * 60 / np.pi) ** 2
    d = hp.dist2holes(
        mask,
        hole_min_size=2,
        hole_min_surf_arcmin2=area_per_pix * 2.0,
    )
    assert d[pix_small] > 0
    assert np.any(d[large_pix] == 0)
