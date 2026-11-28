"""The pixel ufuncs split long inputs across OpenMP threads.

Calls on short slices stay on the calling thread, so comparing one long call
against the same input cut into short slices checks the threaded path against
the serial one. These tests pass the same way whether healpy is built with or
without OpenMP, and for any ``OMP_NUM_THREADS``.
"""

import os
import subprocess
import sys
import textwrap

import numpy as np
import pytest

from healpy import _healpy_pixel_lib as pixlib

# Long enough to be split across threads; slices of SHORT run serially
# (the threshold in src/_healpy_pixel_lib.cc is 1000 elements).
LONG = 20_000
SHORT = 500


def _serial(func, *args):
    """Call ``func`` on slices short enough to stay on one thread."""
    args = np.broadcast_arrays(*args)
    pieces = [
        func(*(a[i : i + SHORT] for a in args)) for i in range(0, len(args[0]), SHORT)
    ]
    if isinstance(pieces[0], tuple):
        return tuple(np.concatenate(p) for p in zip(*pieces))
    return np.concatenate(pieces)


def _assert_same(actual, expected):
    if not isinstance(expected, tuple):
        actual, expected = (actual,), (expected,)
    assert len(actual) == len(expected)
    for a, e in zip(actual, expected):
        # Each element is computed by the same code whichever thread runs it,
        # so the results must be identical, not merely close.
        np.testing.assert_array_equal(a, e, strict=True)


def _inputs(nside, rng):
    npix = 12 * nside**2
    theta = np.arccos(rng.uniform(-1, 1, LONG))
    phi = rng.uniform(0, 2 * np.pi, LONG)
    ipix = rng.integers(0, npix, LONG)
    vec = rng.normal(size=(3, LONG))
    x, y = rng.integers(0, nside, (2, LONG))
    face = rng.integers(0, 12, LONG)
    return theta, phi, ipix, vec, x, y, face


def _cases(nside, scheme, theta, phi, ipix, vec, x, y, face):
    f = lambda name: getattr(pixlib, f"_{name}_{scheme}")
    return {
        "ang2pix": (f("ang2pix"), nside, theta, phi),
        "pix2ang": (f("pix2ang"), nside, ipix),
        "xyf2pix": (f("xyf2pix"), nside, x, y, face),
        "pix2xyf": (f("pix2xyf"), nside, ipix),
        "vec2pix": (f("vec2pix"), nside, *vec),
        "pix2vec": (f("pix2vec"), nside, ipix),
        "get_interpol": (f("get_interpol"), nside, theta, phi),
        "get_neighbors": (f("get_neighbors"), nside, ipix),
    }


@pytest.mark.parametrize("scheme", ["ring", "nest"])
@pytest.mark.parametrize("nside", [1, 16, 4096, 2**29])
def test_threaded_matches_serial(nside, scheme):
    rng = np.random.default_rng(nside)
    for name, (func, *args) in _cases(nside, scheme, *_inputs(nside, rng)).items():
        _assert_same(func(*args), _serial(func, *args))


@pytest.mark.parametrize("nside", [3, 100])
def test_threaded_matches_serial_ring_non_power_of_two(nside):
    rng = np.random.default_rng(nside)
    for name, (func, *args) in _cases(nside, "ring", *_inputs(nside, rng)).items():
        _assert_same(func(*args), _serial(func, *args))


@pytest.mark.parametrize("nside", [64, 2**20])
def test_ring2nest_nest2ring(nside):
    ipix = np.random.default_rng(0).integers(0, 12 * nside**2, LONG)
    for func in (pixlib._ring2nest, pixlib._nest2ring):
        _assert_same(func(nside, ipix), _serial(func, nside, ipix))
    np.testing.assert_array_equal(
        pixlib._nest2ring(nside, pixlib._ring2nest(nside, ipix)), ipix
    )


def test_per_element_nside():
    # Each thread keeps its own nside cache; a different nside per element
    # makes every thread reset it many times.
    rng = np.random.default_rng(1)
    nside = 2 ** rng.integers(0, 20, LONG)
    theta = np.arccos(rng.uniform(-1, 1, LONG))
    phi = rng.uniform(0, 2 * np.pi, LONG)
    for scheme in ["ring", "nest"]:
        ang2pix = getattr(pixlib, f"_ang2pix_{scheme}")
        pix2ang = getattr(pixlib, f"_pix2ang_{scheme}")
        neighbors = getattr(pixlib, f"_get_neighbors_{scheme}")
        ipix = ang2pix(nside, theta, phi)
        _assert_same(ipix, _serial(ang2pix, nside, theta, phi))
        _assert_same(pix2ang(nside, ipix), _serial(pix2ang, nside, ipix))
        _assert_same(neighbors(nside, ipix), _serial(neighbors, nside, ipix))


def test_strided_and_broadcast_inputs():
    rng = np.random.default_rng(2)
    theta = np.arccos(rng.uniform(-1, 1, 3 * LONG))[::3]
    phi = rng.uniform(0, 2 * np.pi, (LONG, 2))[:, 1]
    ipix = pixlib._ang2pix_nest(256, theta, phi)
    _assert_same(ipix, pixlib._ang2pix_nest(256, theta.copy(), phi.copy()))
    # 2-d broadcast: numpy calls the loop once per row
    grid = pixlib._ang2pix_ring(256, theta[:, None], phi[None, :8])
    _assert_same(grid[:, 3], pixlib._ang2pix_ring(256, theta, phi[3]))


def test_invalid_theta_gives_minus_one_on_every_thread():
    theta = np.full(LONG, 0.5)
    theta[::7] = np.nan
    phi = np.zeros(LONG)
    ipix = pixlib._ang2pix_ring(64, theta, phi)
    assert (ipix[::7] == -1).all()
    assert (np.delete(ipix, np.s_[::7]) >= 0).all()
    pix, wgt = np.split(np.array(pixlib._get_interpol_ring(64, theta, phi)), 2)
    assert (pix[:, ::7] == -1).all() and np.isnan(wgt[:, ::7]).all()


@pytest.mark.skipif(not hasattr(os, "fork"), reason="needs fork()")
def test_forked_child_does_not_hang():
    # GNU libgomp hangs when a forked child starts a team after the parent
    # used one, so the child must fall back to the serial loop.
    script = textwrap.dedent(
        """
        import multiprocessing, warnings
        import numpy as np
        import healpy as hp

        warnings.simplefilter("ignore", DeprecationWarning)  # fork with threads
        theta = np.linspace(0, np.pi, 200_000)

        def work(_):
            return int(hp.ang2pix(1024, theta, theta).sum())

        expected = work(0)  # the parent runs the loop on threads first
        with multiprocessing.get_context("fork").Pool(2) as pool:
            assert pool.map_async(work, range(4)).get(timeout=60) == [expected] * 4
        """
    )
    env = dict(os.environ, OMP_NUM_THREADS="2")
    subprocess.run([sys.executable, "-c", script], env=env, check=True, timeout=120)


@pytest.mark.parametrize("threads", ["1", "3"])
def test_results_do_not_depend_on_omp_num_threads(threads, tmp_path):
    script = textwrap.dedent(
        """
        import sys
        import numpy as np
        import healpy as hp

        rng = np.random.default_rng(3)
        theta = np.arccos(rng.uniform(-1, 1, 100_000))
        phi = rng.uniform(0, 2 * np.pi, 100_000)
        ipix = hp.ang2pix(2048, theta, phi, nest=True)
        np.savez(
            sys.argv[1],
            ipix=ipix,
            ang=hp.pix2ang(2048, ipix, nest=True),
            vec=hp.pix2vec(2048, ipix, nest=True),
            ring=hp.nest2ring(2048, ipix),
            interp=hp.get_interp_weights(2048, theta, phi)[1],
        )
        """
    )
    out = tmp_path / "out.npz"
    env = dict(os.environ, OMP_NUM_THREADS=threads)
    subprocess.run([sys.executable, "-c", script, str(out)], env=env, check=True)
    got = np.load(out)

    rng = np.random.default_rng(3)
    theta = np.arccos(rng.uniform(-1, 1, 100_000))
    phi = rng.uniform(0, 2 * np.pi, 100_000)
    ipix = _serial(pixlib._ang2pix_nest, 2048, theta, phi)
    np.testing.assert_array_equal(got["ipix"], ipix)
    np.testing.assert_array_equal(got["ang"], _serial(pixlib._pix2ang_nest, 2048, ipix))
    np.testing.assert_array_equal(got["vec"], _serial(pixlib._pix2vec_nest, 2048, ipix))
    np.testing.assert_array_equal(got["ring"], _serial(pixlib._nest2ring, 2048, ipix))
    np.testing.assert_array_equal(
        got["interp"], _serial(pixlib._get_interpol_nest, 2048, theta, phi)[4:]
    )
