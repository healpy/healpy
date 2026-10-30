import numpy as np
import pytest

import healpy as hp
from healpy import superpixels
from healpy.superpixels import (
    _adjust_centers,
    _canonical_valid_pix,
    _dedup_centers,
)


def make_analytic_map(nside):
    """Large-scale gradient plus Gaussian blobs, smooth in nside."""
    npix = hp.nside2npix(nside)
    theta, phi = hp.pix2ang(nside, np.arange(npix))
    vec = np.asarray(hp.pix2vec(nside, np.arange(npix)))
    m = 4.0 * np.cos(theta)
    for amp, th0, ph0, sig in [
        (3.0, 1.0, 0.5, 0.4),
        (-2.0, 2.2, 3.5, 0.3),
        (2.5, 1.4, 5.0, 0.25),
    ]:
        v0 = hp.ang2vec(th0, ph0)
        d = np.arccos(np.clip(vec.T @ v0, -1.0, 1.0))
        m += amp * np.exp(-0.5 * (d / sig) ** 2)
    return m


def ari(a, b):
    """Adjusted Rand index of two integer label arrays."""
    a = np.asarray(a)
    b = np.asarray(b)
    n = a.size
    ua, ia = np.unique(a, return_inverse=True)
    ub, ib = np.unique(b, return_inverse=True)
    C = np.zeros((ua.size, ub.size), dtype=np.int64)
    np.add.at(C, (ia, ib), 1)

    def c2(x):
        x = np.asarray(x, dtype=np.float64)
        return float((x * (x - 1) / 2).sum())

    sij = c2(C)
    sa = c2(C.sum(axis=1))
    sb = c2(C.sum(axis=0))
    tot = n * (n - 1) / 2
    expected = sa * sb / tot
    denom = (sa + sb) / 2 - expected
    if denom == 0:
        return 1.0 if sij == expected else 0.0
    return (sij - expected) / denom


NSIDE = 32


@pytest.fixture(scope="module")
def analytic_map():
    return make_analytic_map(NSIDE)


@pytest.mark.parametrize("init", ["hierarchical", "farthest", "greedy"])
@pytest.mark.parametrize("nest", [False, True])
def test_output_basic(analytic_map, init, nest):
    """Output has the right shape, dtype, labels, and uses all segments."""
    kwargs = dict(seed=1234) if init == "greedy" else {}
    labels = superpixels(analytic_map, 16, nest=nest, init=init, **kwargs)
    npix = hp.nside2npix(NSIDE)
    assert labels.shape == (npix,)
    assert labels.dtype == np.int32
    assert np.all(labels >= 0)
    assert np.unique(labels).size == 16
    assert labels.max() == 15


def test_invalid_pixels_unseen(analytic_map):
    m = analytic_map.copy()
    m[::10] = hp.UNSEEN
    labels = superpixels(m, 16)
    assert np.all(labels[::10] == -1)
    assert np.all(labels[np.mod(np.arange(m.size), 10) == 3] >= 0)


def test_invalid_pixels_nan(analytic_map):
    m = analytic_map.copy()
    m[::10] = np.nan
    labels = superpixels(m, 16)
    assert np.all(labels[::10] == -1)


def test_invalid_pixels_masked_array(analytic_map):
    mask = np.mod(np.arange(analytic_map.size), 7) == 0
    m = np.ma.MaskedArray(analytic_map, mask=mask)
    labels = superpixels(m, 16)
    assert np.all(labels[mask] == -1)


def test_invalid_pixels_any_feature(analytic_map):
    """A pixel is invalid if any feature is bad."""
    f2 = analytic_map.copy()
    bad = np.mod(np.arange(analytic_map.size), 5) == 0
    f2[bad] = hp.UNSEEN
    labels = superpixels(np.stack([analytic_map, f2]), 16)
    assert np.all(labels[bad] == -1)
    good = np.mod(np.arange(analytic_map.size), 5) == 2
    assert np.all(labels[good] >= 0)


@pytest.mark.parametrize("init", ["hierarchical", "farthest", "greedy"])
def test_ring_nested_identity(analytic_map, init):
    """RING and NESTED input give bit-identical partitions, any init."""
    kwargs = dict(seed=99) if init == "greedy" else {}
    allpix = np.arange(analytic_map.size)
    labels_ring = superpixels(analytic_map, 32, init=init, **kwargs)
    m_nest = analytic_map[hp.nest2ring(NSIDE, allpix)]
    labels_nest = superpixels(m_nest, 32, nest=True, init=init, **kwargs)
    # RING->NEST conversion is m[nest2ring(arange)]
    assert np.array_equal(labels_nest, labels_ring[hp.nest2ring(NSIDE, allpix)])


def test_ring_nested_identity_masked_multifeature():
    """Ordering independence also holds with masks and feature stacks."""
    nside = 64
    npix = hp.nside2npix(nside)
    m = make_analytic_map(nside)
    f2 = 1000.0 * np.cos(2 * hp.pix2ang(nside, np.arange(npix))[0])
    bad = np.mod(np.arange(npix), 5) == 0
    f2[bad] = hp.UNSEEN
    stack = np.stack([m, f2])
    allpix = np.arange(npix)
    labels_ring = superpixels(stack, 32, init="farthest")
    stack_nest = stack[:, hp.nest2ring(nside, allpix)]
    labels_nest = superpixels(stack_nest, 32, nest=True, init="farthest")
    assert np.array_equal(labels_nest, labels_ring[hp.nest2ring(nside, allpix)])


def test_determinism(analytic_map):
    l1 = superpixels(analytic_map, 16, init="greedy", seed=42)
    l2 = superpixels(analytic_map, 16, init="greedy", seed=42)
    l3 = superpixels(analytic_map, 16, init="greedy", seed=43)
    assert np.array_equal(l1, l2)
    assert not np.array_equal(l1, l3)
    l1 = superpixels(analytic_map, 16, init="farthest")
    l2 = superpixels(analytic_map, 16, init="farthest")
    assert np.array_equal(l1, l2)


@pytest.mark.parametrize("init", ["hierarchical", "farthest", "greedy"])
def test_single_segment(analytic_map, init):
    kwargs = dict(seed=1) if init == "greedy" else {}
    labels = superpixels(analytic_map, 1, init=init, **kwargs)
    assert np.all(labels == 0)


def test_k_equals_n_valid():
    """With one segment per valid pixel every pixel is its own cluster."""
    nside = 16
    m = make_analytic_map(nside)
    labels = superpixels(m, hp.nside2npix(nside), init="farthest")
    counts = np.bincount(labels[labels >= 0], minlength=m.size)
    assert np.all(counts == 1)


@pytest.mark.parametrize(
    "n_segments", [0, -5, 48.5, np.inf, np.nan, None, 1j, hp.nside2npix(NSIDE) + 1]
)
def test_n_segments_validation(analytic_map, n_segments):
    with pytest.raises(ValueError, match="n_segments"):
        superpixels(analytic_map, n_segments)


def test_n_segments_array_rejected(analytic_map):
    # size-1 arrays only fail int() on recent numpy, so reject explicitly
    with pytest.raises(ValueError, match="n_segments"):
        superpixels(analytic_map, np.array([16]))
    with pytest.raises(ValueError, match="n_segments"):
        superpixels(analytic_map, np.array([16, 17]))


def test_n_segments_integral_float_accepted(analytic_map):
    labels = superpixels(analytic_map, 16.0)
    assert np.all(labels >= 0)


@pytest.mark.parametrize("max_iter", [0, -1, 2.5, np.inf, np.nan])
def test_max_iter_validation(analytic_map, max_iter):
    with pytest.raises(ValueError, match="max_iter"):
        superpixels(analytic_map, 16, max_iter=max_iter)


def test_max_iter_integral_float_accepted(analytic_map):
    labels = superpixels(analytic_map, 16, max_iter=10.0)
    assert np.all(labels >= 0)


@pytest.mark.parametrize("compactness", [0.0, -2.0, np.nan, np.inf, None, "small"])
def test_compactness_validation(analytic_map, compactness):
    with pytest.raises(ValueError, match="compactness"):
        superpixels(analytic_map, 16, compactness=compactness)


def test_compactness_array_rejected(analytic_map):
    with pytest.raises(ValueError, match="compactness"):
        superpixels(analytic_map, 16, compactness=np.array([0.5]))


def test_max_iter_one_runs(analytic_map):
    labels = superpixels(analytic_map, 16, max_iter=1)
    assert np.all(labels >= 0)


def test_init_validation(analytic_map):
    with pytest.raises(ValueError, match="init"):
        superpixels(analytic_map, 16, init="bogus")


def test_bad_map_size():
    with pytest.raises(ValueError, match="HEALPix"):
        superpixels(np.arange(1000.0), 8)


def test_3d_input_rejected():
    with pytest.raises(ValueError, match="shape"):
        superpixels(np.zeros((2, 3, 100)), 8)


def test_transposed_stack_rejected():
    """A (npix, nfeatures) stack is not a valid input."""
    with pytest.raises(ValueError):
        superpixels(np.zeros((100, 2)), 8)


def test_all_invalid_rejected():
    with pytest.raises(ValueError, match="n_segments"):
        superpixels(np.full(hp.nside2npix(NSIDE), np.nan), 8)


def test_return_centers(analytic_map):
    """Center values are the cluster means, in the units of the input map."""
    labels, center_pix, center_values = superpixels(
        analytic_map, 16, return_centers=True
    )
    assert center_pix.shape == (16,)
    assert center_values.shape == (16, 1)
    means = np.array([analytic_map[labels == i].mean() for i in range(16)])
    assert np.allclose(center_values[:, 0], means, rtol=1e-9)


def test_compactness_scale_invariance(analytic_map):
    """Normalisation makes the partition independent of the map units.

    Scaling by a power of two is exact in floating point, so the
    partitions must match bit for bit.
    """
    l1 = superpixels(analytic_map, 16, compactness=0.5, init="farthest")
    l2 = superpixels(1024.0 * analytic_map, 16, compactness=0.5, init="farthest")
    assert np.array_equal(l1, l2)


def test_feature_scale_invariance(analytic_map):
    """Scaling one feature of a stack does not change the partition."""
    stack = np.stack([analytic_map, 1024.0 * analytic_map])
    l1 = superpixels(stack, 16, compactness=0.5, init="farthest")
    stack2 = np.stack([analytic_map, analytic_map])
    l2 = superpixels(stack2, 16, compactness=0.5, init="farthest")
    assert np.array_equal(l1, l2)


def test_multifeature_balances_scales(analytic_map):
    """Two different features both contribute to the partition."""
    theta = hp.pix2ang(NSIDE, np.arange(analytic_map.size))[0]
    feat2 = 1000.0 * np.cos(2 * theta)
    rng = np.random.default_rng(3)
    feat2 = feat2 + rng.standard_normal(feat2.size)
    l_both = superpixels(np.stack([analytic_map, feat2]), 64, init="farthest")
    l_f1 = superpixels(analytic_map, 64, init="farthest")
    l_f2 = superpixels(feat2, 64, init="farthest")
    # neither single-feature partition dominates the joint one
    assert ari(l_both, l_f1) < 0.95
    assert ari(l_both, l_f2) < 0.95


@pytest.mark.parametrize("init", ["hierarchical", "farthest", "greedy"])
def test_constant_map(init):
    """Constant maps (zero feature range) work and use all segments."""
    nside = 64
    m = np.zeros(hp.nside2npix(nside))
    kwargs = dict(seed=7) if init == "greedy" else {}
    labels = superpixels(m, 16, init=init, **kwargs)
    assert np.all(labels >= 0)
    assert np.unique(labels).size == 16


def test_constant_feature_in_stack(analytic_map):
    stack = np.stack([np.full(analytic_map.size, 3.14), analytic_map])
    labels = superpixels(stack, 16)
    assert np.all(labels >= 0)


def test_leftover_assignment_band():
    """Elongated masks with few segments leave no valid pixel behind.

    The nominal spacing assumes a roughly isotropic valid region; a thin
    band with few segments has valid pixels beyond 2S of every centre,
    which are assigned to the globally closest centre.
    """
    nside = 64
    npix = hp.nside2npix(nside)
    theta, phi = hp.pix2ang(nside, np.arange(npix))
    band = np.abs(np.pi / 2 - theta) < np.radians(5)
    m = np.where(band, 1.0, hp.UNSEEN)
    labels = superpixels(m, 4, init="farthest")
    assert np.all(labels[band] >= 0)
    assert np.all(labels[~band] == -1)
    assert np.unique(labels[band]).size == 4


def test_duplicate_centers_dedup():
    """The seed nudge can duplicate centres; dedup keeps them distinct."""
    nside = 64
    npix = hp.nside2npix(nside)
    m = np.zeros((1, npix))
    valid = np.ones(npix, dtype=bool)
    xyz = np.asarray(hp.pix2vec(nside, np.arange(npix))).T
    # two adjacent seeds on a constant map nudge to the same pixel
    seeds = np.array([100, 101])
    adjusted = _adjust_centers(m, valid, xyz, nside, seeds, 2, nest=False)
    assert adjusted[0] == adjusted[1]
    fixed = _dedup_centers(
        _canonical_valid_pix(valid, nside, False), xyz, nside, adjusted, nest=False
    )
    assert np.unique(fixed).size == 2
    # the first centre is kept, the duplicate is re-placed
    assert fixed[0] == adjusted[0]


def test_duplicate_centers_dedup_ordering():
    """Dedup re-places duplicates at the same physical pixels in RING and NESTED."""
    nside = 64
    npix = hp.nside2npix(nside)
    valid = np.ones(npix, dtype=bool)
    allpix = np.arange(npix)
    xyz_ring = np.asarray(hp.pix2vec(nside, allpix, nest=False)).T
    xyz_nest = xyz_ring[hp.nest2ring(nside, allpix)]

    # duplicated centre at ring pixel 100, expressed in each ordering
    dup_ring = np.array([100, 100])
    dup_nest = hp.ring2nest(nside, dup_ring)
    fixed_ring = _dedup_centers(
        _canonical_valid_pix(valid, nside, False), xyz_ring, nside,
        dup_ring, nest=False,
    )
    fixed_nest = _dedup_centers(
        _canonical_valid_pix(valid, nside, True), xyz_nest, nside,
        dup_nest, nest=True,
    )
    assert np.array_equal(fixed_nest, hp.ring2nest(nside, fixed_ring))


@pytest.mark.parametrize("nside", [1, 2])
@pytest.mark.parametrize("init", ["hierarchical", "farthest", "greedy"])
def test_tiny_maps(nside, init):
    m = make_analytic_map(nside)
    kwargs = dict(seed=1) if init == "greedy" else {}
    labels = superpixels(m, 4, init=init, **kwargs)
    assert np.all(labels >= 0)
    assert np.unique(labels).size == 4


def test_input_types(analytic_map):
    l_ref = superpixels(analytic_map, 8)
    assert np.array_equal(superpixels(list(analytic_map), 8), l_ref)
    mi = np.round(100 * analytic_map).astype(np.int64)
    assert np.array_equal(
        superpixels(mi, 8), superpixels(mi.astype(np.float64), 8)
    )
    labels = superpixels(analytic_map.astype(np.float32), 8)
    assert np.all(labels >= 0)


def test_resolution_consistency():
    """Partitions are equally good across nside but not identical.

    SLIC is not resolution-invariant: the same map at different nside
    gives statistically equivalent partitions (same within-superpixel
    scatter) but spatially different ones. This test documents that
    behaviour; see the superpixels docstring.
    """
    k = 64
    labels = {}
    rms = {}
    for nside in [32, 64]:
        m = make_analytic_map(nside)
        labels[nside] = superpixels(m, k, init="farthest")
        vals = m[labels[nside] >= 0]
        labs = labels[nside][labels[nside] >= 0]
        means = np.bincount(labs, weights=vals, minlength=k) / np.maximum(
            np.bincount(labs, minlength=k), 1
        )
        rms[nside] = np.sqrt(np.mean((vals - means[labs]) ** 2))
    # partition quality is the same at both resolutions (within 10%)
    assert rms[64] == pytest.approx(rms[32], rel=0.1)
    # but the partitions themselves are far from identical
    m64 = make_analytic_map(64)
    theta, phi = hp.pix2ang(64, np.arange(m64.size))
    sampled = labels[32][hp.ang2pix(32, theta, phi)]
    agreement = ari(sampled, labels[64])
    assert 0.2 < agreement < 0.95


def test_compactness_saturation(analytic_map):
    """Past the normalised feature range the partition barely changes."""
    l10 = superpixels(analytic_map, 64, compactness=10.0, init="farthest")
    l100 = superpixels(analytic_map, 64, compactness=100.0, init="farthest")
    assert ari(l10, l100) > 0.9
