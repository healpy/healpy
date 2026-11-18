#
#  This file is part of Healpy.
#
#  Healpy is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  Healpy is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Healpy; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#
#  For more information about Healpy, see https://github.com/healpy/healpy
#
"""Superpixel segmentation of HEALPix maps with spherical SLIC.

Original implementation by Nina Öhlckers, proposed in
https://github.com/healpy/healpy/issues/1131, adapting the SLIC
algorithm of Achanta et al. (2012) to the HEALPix grid.
"""

import logging

import numpy as np

from . import pixelfunc
from ._query_disc import query_disc

log = logging.getLogger("healpy")

__all__ = ["superpixels"]

_INIT_STRATEGIES = ("hierarchical", "greedy", "farthest")


def _as_int(value, name):
    """Convert to int, raising ValueError for invalid or non-integral input."""
    if np.ndim(value) > 0:
        # int() of a size-1 array is only an error on recent numpy, and
        # comparing the result against the array is ambiguous: reject
        # arrays of any shape explicitly
        raise ValueError(f"{name} must be an integer, got {value!r}")
    try:
        result = int(value)
    except (TypeError, ValueError, OverflowError):
        raise ValueError(f"{name} must be an integer, got {value!r}") from None
    if result != value:
        raise ValueError(f"{name} must be an integer, got {value!r}")
    return result


def superpixels(
    m,
    n_segments,
    compactness=1.0,
    nest=False,
    init="hierarchical",
    max_iter=100,
    seed=None,
    return_centers=False,
):
    """Segment a HEALPix map into superpixels with spherical SLIC.

    Superpixels oversegment a map into a small number of regions that
    group neighbouring pixels with similar values, so that downstream
    tasks can work with a few hundred coherent regions instead of
    millions of pixels. This function adapts the Simple Linear
    Iterative Clustering (SLIC) algorithm of Achanta et al. (2012) [1]_
    to the sphere: spatial distances are angular distances between pixel
    centres and the clustering works directly on the HEALPix grid.
    Superpixel segmentation on the sphere was also studied for
    panoramic images as SphSLIC (Zhao et al. 2018 [2]_). This function
    is based on an implementation by Nina Öhlckers.

    Parameters
    ----------
    m : array-like, shape (npix,) or (nfeatures, npix)
        Input map, or a stack of feature maps sharing the same
        pixelisation. Features are normalised to ``[0, 1]`` over the
        valid pixels before clustering, so they may have arbitrary units
        and scales. Pixels that are ``UNSEEN``, NaN or masked in any
        feature are not assigned to a superpixel.
    n_segments : int
        Number of superpixels (cluster centres) to create.
    compactness : float, optional
        Scale of the feature distance relative to the spatial distance.
        Features are normalised internally, so ``compactness`` is
        dimensionless: values close to 0 let superpixels follow the map
        values closely, intermediate values give roughly round
        superpixels that snap to strong edges, and values of order 1
        and above give compact, regular superpixels, saturating to an
        almost purely spatial partition (a fancier ``ud_grade``) for
        much larger values. Default: 1.0
    nest : bool, optional
        If True, assume NESTED pixel ordering, otherwise RING ordering.
        For the same random draws (a fixed `seed` with ``init='greedy'``),
        the returned partition is identical for the two orderings.
        Default: False
    init : {'hierarchical', 'greedy', 'farthest'}, optional
        Strategy used to place the initial cluster centres:

        * ``'hierarchical'``: deterministic seeding from the pixel
          centres of the coarser HEALPix level whose pixel count is
          closest to `n_segments`, subsampled or supplemented to
          exactly `n_segments`.
        * ``'greedy'``: greedy k-means++ [3]_ seeding on angular distance
          (random, controlled by `seed`).
        * ``'farthest'``: deterministic farthest-point sampling,
          starting from the lowest valid RING pixel.

        Default: 'hierarchical'
    max_iter : int, optional
        Maximum number of assignment/update iterations. The iteration
        stops early once the cluster centre pixels no longer change.
        Default: 100
    seed : None, int or numpy.random.Generator, optional
        Seed for the random number generator used by ``init='greedy'``.
    return_centers : bool, optional
        If True, also return the centre pixels and centre feature
        values. Default: False

    Returns
    -------
    labels : ndarray of int32, shape (npix,)
        Superpixel index of every pixel, in ``[0, n_segments)``.
        Invalid pixels are labelled ``-1``. In rare cases a cluster can
        end up without pixels, leaving unused indices.
    center_pix : ndarray of int, shape (n_segments,)
        Pixel index of each cluster centre. Only returned if
        `return_centers` is True.
    center_values : ndarray, shape (n_segments, nfeatures)
        Mean feature values of each cluster, in the units of the input
        map. A cluster that ends up without pixels keeps its previous
        centre values. Only returned if `return_centers` is True.

    Raises
    ------
    ValueError
        If the map size is not a valid HEALPix pixel number, if
        `n_segments` is not an integer in ``[1, number of valid
        pixels]``, if `compactness` is not a finite positive number, if
        `max_iter` is not an integer of at least 1, or if `init` is
        unknown.

    Notes
    -----
    Each valid pixel is assigned to the centre minimising

    .. math::

        D = \\sqrt{(d_f / c)^2 + (d_s / S)^2},

    where :math:`d_f` is the Euclidean distance in the normalised
    feature space, :math:`c` is `compactness`, :math:`d_s` is the
    angular distance between pixel and centre in radians, and
    :math:`S = 2 \\sqrt{f / k}` is the nominal cluster scale for
    :math:`k` superpixels covering a valid fraction :math:`f` of the
    sphere (the radius of a disc with the area of one superpixel), so
    masked maps keep the same superpixel density. Only
    pixels within an angular radius :math:`2S` of a centre are compared
    in the iteration (this local search is what makes SLIC fast); valid
    pixels not reached by any search disc are assigned to the globally
    closest centre afterwards, so every valid pixel gets a label.
    Before iterating, every initial centre is moved to the pixel with
    the lowest local gradient within a radius :math:`S / 4`, so that
    centres do not start on edges; centres that end up duplicated by
    this nudge are re-placed deterministically.

    The result is exactly reproducible: ``init='hierarchical'`` and
    ``init='farthest'`` are fully deterministic, and with a fixed `seed`
    so is ``init='greedy'``. Given the same random draws, the partition
    is bit-identical for RING and NESTED input for any `init`.

    Two caveats are inherent to SLIC rather than to this
    implementation:

    * The segmentation is not resolution-invariant: the same map
      evaluated at several ``nside`` with the same `n_segments` and
      `compactness` gives partitions of the same quality (same
      within-superpixel scatter and superpixel size) but spatially
      different, with an adjusted Rand index of roughly 0.5 to 0.7
      between adjacent resolutions. Use superpixels for region
      statistics at one resolution, not to track structures across
      resolutions.
    * Superpixels are not guaranteed to be spatially connected, in
      particular on noisy maps at low `compactness`. See SNIC
      (Achanta & Süsstrunk 2017) for a superpixel algorithm that
      enforces connectivity.

    References
    ----------
    .. [1] R. Achanta, A. Shaji, K. Smith, A. Lucchi, P. Fua
       and S. Süsstrunk, "SLIC Superpixels Compared to
       State-of-the-Art Superpixel Methods", IEEE TPAMI 34(11),
       2274-2282 (2012).
    .. [2] Q. Zhao, D. Yang, H. Lu, T. Huang and Z. Wang, "Spherical
       Superpixel Segmentation", IEEE Trans. Multimedia 20(6),
       1406-1417 (2018).
    .. [3] D. Arthur and S. Vassilvitskii, "k-means++: The Advantages
       of Careful Seeding", ACM-SIAM SODA (2007).

    Examples
    --------
    >>> import numpy as np
    >>> import healpy as hp
    >>> nside = 16
    >>> theta, phi = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)))
    >>> m = np.cos(theta)
    >>> labels = hp.superpixels(m, 48, compactness=0.1)
    >>> labels.shape
    (3072,)

    Masked maps and ``return_centers``:

    >>> m = np.ma.MaskedArray(m, mask=(theta < 0.5))
    >>> labels, center_pix, center_values = hp.superpixels(
    ...     m, 12, compactness=0.5, return_centers=True
    ... )
    >>> center_values.shape
    (12, 1)
    >>> bool(np.all(labels[np.ma.getmaskarray(m)] == -1))
    True
    """
    if init not in _INIT_STRATEGIES:
        raise ValueError(
            f"Unknown init strategy {init!r}, expected one of {_INIT_STRATEGIES}"
        )

    if np.ndim(compactness) > 0:
        raise ValueError(
            f"compactness must be a finite positive number, got {compactness!r}"
        )
    try:
        compactness = float(compactness)
    except (TypeError, ValueError):
        raise ValueError(
            f"compactness must be a finite positive number, got {compactness!r}"
        ) from None
    if not np.isfinite(compactness) or compactness <= 0:
        raise ValueError(
            f"compactness must be a finite positive number, got {compactness!r}"
        )

    invalid = np.ma.getmaskarray(m)
    m_raw = np.atleast_2d(np.asarray(np.ma.getdata(m), dtype=np.float64))
    if m_raw.ndim != 2:
        raise ValueError("m must have shape (npix,) or (nfeatures, npix)")
    invalid = np.atleast_2d(invalid) | ~np.isfinite(m_raw) | pixelfunc.mask_bad(m_raw)
    valid = ~np.any(invalid, axis=0)

    npix = m_raw.shape[1]
    if not pixelfunc.isnpixok(npix):
        raise ValueError(f"Map has {npix} pixels, which is not a valid HEALPix size")
    nside = pixelfunc.npix2nside(npix)

    k = _as_int(n_segments, "n_segments")
    n_valid = int(valid.sum())
    if not 1 <= k <= n_valid:
        raise ValueError(
            f"n_segments must be between 1 and the number of valid pixels "
            f"({n_valid}), got {n_segments}"
        )

    max_iter = _as_int(max_iter, "max_iter")
    if max_iter < 1:
        raise ValueError(f"max_iter must be at least 1, got {max_iter}")

    # normalise each feature to [0, 1] over the valid pixels, so that
    # compactness is dimensionless; constant features are mapped to 0
    # and contribute no feature distance
    m_min = m_raw[:, valid].min(axis=1)
    m_max = m_raw[:, valid].max(axis=1)
    feat_range = np.where(m_max > m_min, m_max - m_min, 1.0)
    m = (m_raw - m_min[:, None]) / feat_range[:, None]

    rng = np.random.default_rng(seed)
    # always compute the pixel vectors in RING ordering and remap for NESTED,
    # so the geometry is bitwise identical regardless of the map ordering
    # (pix2vec itself is not bitwise reproducible across orderings)
    xyz = np.asarray(pixelfunc.pix2vec(nside, np.arange(npix), nest=False)).T
    if nest:
        xyz = xyz[pixelfunc.nest2ring(nside, np.arange(npix))]

    center_pix, center_values = _init_centers(m, valid, xyz, nside, k, init, nest, rng)
    if np.unique(center_pix).size < k:
        # the low-gradient nudge can move two seeds to the same pixel:
        # re-place the duplicates deterministically so that all
        # n_segments clusters are distinct
        center_pix = _dedup_centers(
            _canonical_valid_pix(valid, nside, nest), xyz, nside, center_pix, nest
        )
        center_values = m[:, center_pix].T.copy()

    labels = -np.ones(npix, dtype=np.int32)
    for it in range(max_iter):
        labels = _assign(
            m, valid, xyz, nside, center_pix, center_values, compactness, nest
        )
        new_pix, new_values = _update_centers(
            m, labels, xyz, nside, center_pix, center_values, nest
        )
        converged = np.array_equal(center_pix, new_pix)
        center_pix, center_values = new_pix, new_values
        if converged:
            break
    log.info("SLIC stopped after %d iteration(s)", it + 1)

    # return the centre values in the units of the input map
    center_values = center_values * feat_range[None, :] + m_min[None, :]

    if return_centers:
        return labels, center_pix, center_values
    return labels


def _angular_distance(a, b):
    """Angular distance in radians between unit vectors ``a`` (n, 3) and ``b`` (3,)."""
    return np.arccos(np.clip(a @ b, -1.0, 1.0))


def _feature_distance(m, pix, center_values):
    """Euclidean feature distance between pixels ``pix`` and a centre."""
    diff = m[:, pix] - center_values[:, None]
    return np.sqrt(np.sum(diff**2, axis=0))


def _canonical_valid_pix(valid, nside, nest):
    """Valid pixels sorted by RING index, so that sampling and
    tie-breaking do not depend on the map ordering."""
    vp = np.flatnonzero(valid)
    if nest:
        vp = vp[np.argsort(pixelfunc.nest2ring(nside, vp), kind="stable")]
    return vp


def _seeds_hierarchical(valid, nside, k, nest):
    """Seeds at the centres of the closest coarser HEALPix level.

    Coarse pixel centres fall on fine pixel corners, where ``ang2pix``
    breaks ties differently per ordering, so seeding is always done in
    RING indices and converted afterwards to keep it ordering
    independent.
    """
    npix = valid.size
    if nest:
        valid = valid[pixelfunc.ring2nest(nside, np.arange(npix))]

    seed_nside = min(
        (2**i for i in range(int(np.log2(nside)) + 1)),
        key=lambda ns: abs(pixelfunc.nside2npix(ns) - k),
    )
    theta, phi = pixelfunc.pix2ang(
        seed_nside, np.arange(pixelfunc.nside2npix(seed_nside))
    )
    seeds = pixelfunc.ang2pix(nside, theta, phi)
    seeds = seeds[valid[seeds]]

    if seeds.size > k:
        # uniformly subsample the coarse grid
        seeds = seeds[np.linspace(0, seeds.size - 1, k, dtype=int)]
    elif seeds.size < k:
        # supplement with approximately uniformly spread fine pixels
        available_mask = valid.copy()
        available_mask[seeds] = False
        available = np.arange(npix)[available_mask]
        extra = np.linspace(0, available.size - 1, k - seeds.size, dtype=int)
        seeds = np.concatenate([seeds, available[extra]])

    if nest:
        seeds = pixelfunc.ring2nest(nside, seeds)
    return seeds


def _seeds_farthest(valid_pix, xyz, k):
    """Deterministic farthest-point sampling on squared angular distance."""
    pts = xyz[valid_pix]
    seeds = [0]
    min_dist = _angular_distance(pts, pts[0]) ** 2
    for _ in range(1, k):
        p = np.argmax(min_dist)
        seeds.append(p)
        min_dist = np.minimum(min_dist, _angular_distance(pts, pts[p]) ** 2)
    return valid_pix[np.asarray(seeds)]


def _seeds_kmeanspp(valid_pix, xyz, k, rng, n_local_trials=None):
    """Greedy k-means++ seeding on squared angular distance."""
    if n_local_trials is None:
        n_local_trials = 2 + int(np.log(k))
    pts = xyz[valid_pix]

    first = rng.integers(valid_pix.size)
    seeds = [first]
    min_dist = _angular_distance(pts, pts[first]) ** 2

    for _ in range(1, k):
        probabilities = min_dist / min_dist.sum()
        n_trials = min(n_local_trials, np.count_nonzero(probabilities))
        candidates = rng.choice(
            valid_pix.size, size=n_trials, replace=False, p=probabilities
        )

        best_potential = np.inf
        for candidate in candidates:
            new_min_dist = np.minimum(
                min_dist, _angular_distance(pts, pts[candidate]) ** 2
            )
            potential = new_min_dist.sum()
            if potential < best_potential:
                best_potential = potential
                best_candidate = candidate
                best_dist = new_min_dist

        seeds.append(best_candidate)
        min_dist = best_dist

    return valid_pix[np.asarray(seeds)]


def _pixel_gradient(m, valid, nside, pix, nest):
    """Mean feature distance between each pixel in ``pix`` and its valid neighbours.

    The HEALPix neighbours are roughly equidistant, so the feature
    difference is used directly as the gradient magnitude.
    """
    neighbours = pixelfunc.get_all_neighbours(nside, pix, nest=nest)
    usable = neighbours >= 0
    usable[usable] = valid[neighbours[usable]]
    neighbours = np.where(usable, neighbours, 0)

    diff = m[:, neighbours] - m[:, pix][:, None, :]
    dist = np.sqrt(np.sum(diff**2, axis=0))
    n_usable = usable.sum(axis=0)
    total = np.where(usable, dist, 0.0).sum(axis=0)
    return np.divide(total, n_usable, out=np.zeros_like(total), where=n_usable > 0)


def _spacing(valid, k):
    """Nominal cluster spacing S for k centres spread over the valid area."""
    return 2.0 * np.sqrt(valid.mean() / k)


def _adjust_centers(m, valid, xyz, nside, seeds, k, nest, adjustment_fraction=0.25):
    """Move every seed to the lowest-gradient valid pixel nearby."""
    search_radius = adjustment_fraction * _spacing(valid, k)
    adjusted = np.empty_like(seeds)
    for i, p in enumerate(seeds):
        candidates = query_disc(nside, xyz[p], search_radius, nest=nest)
        candidates = candidates[valid[candidates]]
        if candidates.size == 0:
            adjusted[i] = p
            continue
        gradients = _pixel_gradient(m, valid, nside, candidates, nest)
        # deterministic, ordering-independent tie-break: among candidates
        # with the exact same minimum gradient, pick the lowest RING index
        cand_ring = pixelfunc.nest2ring(nside, candidates) if nest else candidates
        adjusted[i] = candidates[np.lexsort((cand_ring, gradients))[0]]
    return adjusted


def _dedup_centers(valid_pix_canon, xyz, nside, center_pix, nest):
    """Re-place duplicate centres at the valid pixel farthest from all centres.

    The first occurrence of each duplicated pixel is kept; ties are
    broken by lowest RING index so the result is ordering-independent.
    """
    center_pix = np.asarray(center_pix)
    _, first = np.unique(center_pix, return_index=True)
    out = center_pix.copy()
    keep = np.zeros(center_pix.size, dtype=bool)
    keep[first] = True

    used = np.zeros(xyz.shape[0], dtype=bool)
    used[center_pix] = True
    pts = xyz[valid_pix_canon]
    ring = pixelfunc.nest2ring(nside, valid_pix_canon) if nest else valid_pix_canon

    min_d = np.full(valid_pix_canon.size, np.inf)
    for c in center_pix[np.sort(first)]:
        min_d = np.minimum(min_d, _angular_distance(pts, xyz[c]) ** 2)

    for i in np.flatnonzero(~keep):
        cand_d = np.where(used[valid_pix_canon], -np.inf, min_d)
        best = np.lexsort((ring, -cand_d))[0]
        out[i] = valid_pix_canon[best]
        used[out[i]] = True
        min_d = np.minimum(min_d, _angular_distance(pts, xyz[out[i]]) ** 2)
    return out


def _init_centers(m, valid, xyz, nside, k, init, nest, rng):
    """Initial centre pixels and values for the chosen strategy."""
    if init == "hierarchical":
        seeds = _seeds_hierarchical(valid, nside, k, nest)
    else:
        # greedy and farthest sample over the valid pixels: iterate them
        # in canonical RING order so the sampled physical pixels (and
        # the RNG stream for greedy) do not depend on the map ordering
        valid_pix_canon = _canonical_valid_pix(valid, nside, nest)
        if init == "greedy":
            seeds = _seeds_kmeanspp(valid_pix_canon, xyz, k, rng)
        else:
            seeds = _seeds_farthest(valid_pix_canon, xyz, k)

    seeds = _adjust_centers(m, valid, xyz, nside, seeds, k, nest)
    return seeds, m[:, seeds].T.copy()


def _assign_nearest_chunked(
    m, xyz, pix, center_pix, center_values, compactness, S,
    distances, labels,
):
    """Assign each pixel in ``pix`` to the centre minimising D, with no
    radius cut, in chunks sized to bound memory for large k."""
    k = center_pix.size
    nfeatures = m.shape[0]
    # bound the (chunk, k, nfeatures) temporary to roughly 20M elements
    chunk = max(1, int(2e7 / max(k * nfeatures, 1)))
    center_vec = xyz[center_pix]
    for start in range(0, pix.size, chunk):
        p = pix[start : start + chunk]
        diff = m[:, p].T[:, None, :] - center_values[None, :, :]
        d_value = np.sqrt(np.sum(diff**2, axis=-1))
        d_spatial = np.arccos(np.clip(xyz[p] @ center_vec.T, -1.0, 1.0))
        D = np.sqrt((d_value / compactness) ** 2 + (d_spatial / S) ** 2)
        best = np.argmin(D, axis=1)
        best_D = D[np.arange(p.size), best]
        better = best_D < distances[p]
        labels[p[better]] = best[better]
        distances[p[better]] = best_D[better]
    return labels


def _assign(m, valid, xyz, nside, center_pix, center_values, compactness, nest):
    """Assign every valid pixel within 2S of a centre to it.

    Valid pixels not reached by any search disc are assigned to the
    globally closest centre, so that elongated masks with few segments
    do not leave holes (the nominal spacing assumes a roughly isotropic
    valid region and underestimates it there).
    """
    npix = valid.size
    k = center_pix.size
    labels = -np.ones(npix, dtype=np.int32)
    distances = np.full(npix, np.inf)

    S = _spacing(valid, k)
    search_radius = 2.0 * S
    for i in range(k):
        center_vec = xyz[center_pix[i]]
        candidates = query_disc(nside, center_vec, search_radius, nest=nest)
        candidates = candidates[valid[candidates]]

        d_value = _feature_distance(m, candidates, center_values[i])
        d_spatial = _angular_distance(xyz[candidates], center_vec)
        D = np.sqrt((d_value / compactness) ** 2 + (d_spatial / S) ** 2)

        better = D < distances[candidates]
        distances[candidates[better]] = D[better]
        labels[candidates[better]] = i

    remaining = np.flatnonzero(valid & (labels < 0))
    if remaining.size:
        labels = _assign_nearest_chunked(
            m, xyz, remaining, center_pix, center_values, compactness, S,
            distances, labels,
        )
    return labels


def _update_centers(m, labels, xyz, nside, center_pix, center_values, nest):
    """Move each centre to the mean position and mean values of its pixels."""
    k = center_pix.size
    nfeatures = m.shape[0]
    assigned = labels >= 0
    pix = np.flatnonzero(assigned)
    lab_all = labels[pix]

    # accumulate the per-cluster sums in a canonical order (cluster, RING
    # index) so that the floating point sums are bitwise identical
    # regardless of the map ordering
    pix_ring = pixelfunc.nest2ring(nside, pix) if nest else pix
    order = np.lexsort((pix_ring, lab_all))
    lab = lab_all[order]
    sel = m[:, pix][:, order]
    vec_sel = xyz[pix][order]

    count = np.bincount(lab, minlength=k)
    value_sum = np.empty((k, nfeatures))
    for f in range(nfeatures):
        value_sum[:, f] = np.bincount(lab, weights=sel[f], minlength=k)
    vector_sum = np.empty((k, 3))
    for d in range(3):
        vector_sum[:, d] = np.bincount(lab, weights=vec_sel[:, d], minlength=k)

    new_pix = center_pix.copy()
    new_values = center_values.copy()
    for i in np.flatnonzero(count):
        new_values[i] = value_sum[i] / count[i]
        norm = np.linalg.norm(vector_sum[i])
        if norm > 0:
            v = vector_sum[i] / norm
            new_pix[i] = pixelfunc.vec2pix(nside, v[0], v[1], v[2], nest=nest)
    return new_pix, new_values
