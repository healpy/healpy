import copy

import healpy as hp
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from healpy.newvisufunc import projview
from healpy.visufunc import mollview

matplotlib.use("agg")


def _make_map():
    nside = 8
    npix = hp.nside2npix(nside)
    m = np.full(npix, 1.0, dtype=np.double)
    strip = hp.query_strip(nside, np.pi / 2 - np.radians(15), np.pi / 2 + np.radians(15))
    m[strip] = hp.UNSEEN
    return m


def _get_plot_cmap(ax):
    if ax.images:
        return ax.images[0].get_cmap()
    if ax.collections:
        return ax.collections[0].get_cmap()
    raise AssertionError("No colormap object found in plot output")


def _assert_rgb_close(actual, expected, msg):
    np.testing.assert_allclose(actual[:3], expected[:3], atol=1e-2, err_msg=msg)


def test_mollview_cmap_object_applies_bad_bgcolor_args():
    cmap_obj = copy.copy(plt.get_cmap("viridis"))
    cmap_obj.set_bad("red")
    cmap_obj.set_under("black")

    fig = plt.figure()
    try:
        mollview(_make_map(), cmap=cmap_obj, badcolor="gray", bgcolor="white", cbar=False)
        plot_cmap = _get_plot_cmap(plt.gca())
        _assert_rgb_close(plot_cmap._rgba_bad, (0.5, 0.5, 0.5, 1.0), "bad color mismatch")
        _assert_rgb_close(plot_cmap._rgba_under, (1.0, 1.0, 1.0, 1.0), "under color mismatch")
    finally:
        plt.close(fig)


def test_mollview_cmap_string_applies_bad_bgcolor_args():
    fig = plt.figure()
    try:
        mollview(_make_map(), cmap="viridis", badcolor="gray", bgcolor="white", cbar=False)
        plot_cmap = _get_plot_cmap(plt.gca())
        _assert_rgb_close(plot_cmap._rgba_bad, (0.5, 0.5, 0.5, 1.0), "bad color mismatch")
        _assert_rgb_close(plot_cmap._rgba_under, (1.0, 1.0, 1.0, 1.0), "under color mismatch")
    finally:
        plt.close(fig)


def test_projview_cmap_object_applies_bad_bgcolor_args():
    cmap_obj = copy.copy(plt.get_cmap("viridis"))
    cmap_obj.set_bad("red")
    cmap_obj.set_under("black")

    fig = plt.figure()
    try:
        projview(
            _make_map(),
            projection_type="mollweide",
            cmap=cmap_obj,
            badcolor="gray",
            bgcolor="white",
            cbar=False,
        )
        plot_cmap = _get_plot_cmap(plt.gca())
        _assert_rgb_close(plot_cmap._rgba_bad, (0.5, 0.5, 0.5, 1.0), "bad color mismatch")
        _assert_rgb_close(plot_cmap._rgba_under, (1.0, 1.0, 1.0, 1.0), "under color mismatch")
    finally:
        plt.close(fig)


def test_projview_masks_unseen_values():
    fig = plt.figure()
    try:
        projview(_make_map(), projection_type="mollweide", cmap="viridis", cbar=False)
        ax = plt.gca()
        mesh = ax.collections[0]
        arr = mesh.get_array()
        assert np.ma.isMaskedArray(arr)
        assert np.any(arr.mask)
    finally:
        plt.close(fig)
