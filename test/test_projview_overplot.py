"""
Test that projplot, projscatter, and projtext work on axes created by projview.

Regression test for issue #1132: these functions used to return None and draw
nothing when the current axes was created by projview, because they only
supported the SphericalProjAxes used by mollview and friends.
"""
import matplotlib

matplotlib.use("agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest

import healpy as hp

NSIDE = 32


@pytest.fixture(autouse=True)
def close_figures():
    yield
    plt.close("all")


def one_pixel_map(theta, phi):
    """Map with a single hot pixel at the pixel containing (theta, phi)."""
    m = np.zeros(hp.nside2npix(NSIDE))
    pix = hp.ang2pix(NSIDE, theta, phi)
    m[pix] = 1.0
    # mark the pixel center, so the marker drawn at the center of the
    # projected pixel is guaranteed to fall on a hot grid cell
    theta_p, phi_p = hp.pix2ang(NSIDE, pix)
    return m, theta_p, phi_p


def marker_lands_on_hot_pixel(m, theta, phi, **projview_kwargs):
    """Projview a map, projscatter a marker at (theta, phi) and check that the
    marker falls on a hot pixel of the displayed map."""
    lon_d, lat_d, grid = hp.projview(
        m, return_only_data=True, xsize=800, **projview_kwargs
    )
    hp.projview(m, xsize=800, **projview_kwargs)
    offsets = hp.projscatter(theta, phi).get_offsets()[0]
    ix = np.argmin(np.abs(lon_d - offsets[0]))
    iy = np.argmin(np.abs(lat_d - offsets[1]))
    return grid[iy, ix] == 1.0


@pytest.mark.parametrize(
    "projection_type", ["mollweide", "hammer", "aitoff", "cart", "polar", "lambert"]
)
def test_projscatter_on_projview_all_projections(projection_type):
    # spot at lon -30 deg, lat 40 deg, visible in all projections
    m, theta_p, phi_p = one_pixel_map(np.radians(50), np.radians(30))
    assert marker_lands_on_hot_pixel(m, theta_p, phi_p, projection_type=projection_type)


def test_projscatter_on_projview():
    hp.projview(np.zeros(hp.nside2npix(4)))
    ret = hp.projscatter(np.radians(45), np.radians(60), color="red", s=50)
    # before the fix, projscatter returned None and drew nothing
    assert ret is not None
    ax = plt.gca()
    # the pcolormesh of the map is also a collection, check ours was added
    assert ret in ax.collections
    # astro convention (default): longitude is mirrored
    assert np.allclose(ret.get_offsets()[0], [-np.radians(60), np.radians(45)])


def test_projscatter_lonlat():
    hp.projview(np.zeros(hp.nside2npix(4)))
    ret = hp.projscatter([45], [60], lonlat=True)
    assert np.allclose(ret.get_offsets()[0], [-np.radians(45), np.radians(60)])


def test_projscatter_geo_flip():
    hp.projview(np.zeros(hp.nside2npix(4)), flip="geo")
    ret = hp.projscatter(np.radians(45), np.radians(60))
    assert np.allclose(ret.get_offsets()[0], [np.radians(60), np.radians(45)])


def test_projscatter_wraps_longitude():
    hp.projview(np.zeros(hp.nside2npix(4)))
    ret = hp.projscatter(np.radians(45), np.radians(270))
    # -270 deg wrapped into [-180, 180)
    assert np.allclose(ret.get_offsets()[0], [np.radians(90), np.radians(45)])


def test_projplot_on_projview():
    hp.projview(np.zeros(hp.nside2npix(4)))
    ret = hp.projplot(
        [np.radians(45), np.radians(90)], [np.radians(60), np.radians(0)]
    )
    assert ret is not None
    line = ret[0]
    assert np.allclose(np.degrees(line.get_xdata()), [-60.0, 0.0])
    assert np.allclose(np.degrees(line.get_ydata()), [45.0, 0.0])


def test_projplot_thetaphi_array():
    hp.projview(np.zeros(hp.nside2npix(4)))
    ret = hp.projplot(np.array([[np.radians(45)], [np.radians(60)]]), "ro")
    assert len(ret) == 1
    assert np.allclose(np.degrees(ret[0].get_xdata()), [-60.0])
    assert np.allclose(np.degrees(ret[0].get_ydata()), [45.0])


def test_projplot_splits_line_at_boundary():
    hp.projview(np.zeros(hp.nside2npix(4)))
    # a path crossing the 180 deg meridian: phi = 160, 180, 200 deg
    ret = hp.projplot(
        [np.radians(90)] * 3, np.radians([160, 180, 200]), "r-", label="path"
    )
    # the path is split at the boundary, the label only on the first segment
    assert len(ret) == 2
    assert np.allclose(np.degrees(ret[0].get_xdata()), [-160.0, -180.0])
    assert np.allclose(np.degrees(ret[1].get_xdata()), [160.0])
    assert ret[0].get_label() == "path"
    assert ret[1].get_label().startswith("_")


def test_projtext_on_projview():
    hp.projview(np.zeros(hp.nside2npix(4)))
    ret = hp.projtext(np.radians(45), np.radians(60), "x")
    assert ret is not None
    assert ret.get_text() == "x"
    assert np.allclose(np.degrees(ret.get_position()), [-60.0, 45.0])
    assert len(plt.gca().texts) == 1


def test_projscatter_on_rotated_projview():
    # coordinates are interpreted in the map frame: the rot of the view
    # is applied automatically
    m, theta_p, phi_p = one_pixel_map(np.radians(60), np.radians(45))
    assert marker_lands_on_hot_pixel(m, theta_p, phi_p, rot=[45, 30])


def test_projscatter_on_projview_coord_transform():
    m, theta_p, phi_p = one_pixel_map(np.radians(60), np.radians(45))
    assert marker_lands_on_hot_pixel(m, theta_p, phi_p, coord=["G", "C"])


def test_projscatter_rotated_lonlat():
    m = np.zeros(hp.nside2npix(NSIDE))
    m[hp.ang2pix(NSIDE, np.radians(60), np.radians(45))] = 1.0
    hp.projview(m, rot=[45, 30])
    ret_theta_phi = hp.projscatter(np.radians(60), np.radians(45))
    plt.close("all")
    hp.projview(m, rot=[45, 30])
    ret_lonlat = hp.projscatter(45, 30, lonlat=True)
    assert np.allclose(ret_theta_phi.get_offsets(), ret_lonlat.get_offsets())


def test_newprojplot_applies_view_rotation():
    m = np.zeros(hp.nside2npix(NSIDE))
    m[hp.ang2pix(NSIDE, np.radians(60), np.radians(45))] = 1.0
    hp.projview(m, rot=[45, 30])
    lines = hp.newprojplot([np.radians(60)], [np.radians(45)], marker="o")
    plt.close("all")
    hp.projview(m, rot=[45, 30])
    scatter = hp.projscatter(np.radians(60), np.radians(45))
    assert np.allclose(
        np.column_stack([lines[0].get_xdata(), lines[0].get_ydata()]),
        scatter.get_offsets(),
    )


def test_projview_3d_overplot_unsupported():
    hp.projview(np.zeros(hp.nside2npix(4)), projection_type="3d")
    with pytest.raises(ValueError, match="3d"):
        hp.projscatter(np.radians(45), np.radians(60))
    with pytest.raises(ValueError, match="3d"):
        hp.projplot([np.radians(45)], [np.radians(60)])
    with pytest.raises(ValueError, match="3d"):
        hp.projtext(np.radians(45), np.radians(60), "x")


def test_rot_coord_keywords_unsupported_on_projview():
    hp.projview(np.zeros(hp.nside2npix(4)))
    with pytest.raises(ValueError, match="rot and coord"):
        hp.projscatter(np.radians(45), np.radians(60), rot=[45, 30])
    with pytest.raises(ValueError, match="rot and coord"):
        hp.projplot([np.radians(45)], [np.radians(60)], coord="G")
    with pytest.raises(ValueError, match="rot and coord"):
        hp.projtext(np.radians(45), np.radians(60), "x", coord="G")


@pytest.mark.parametrize("func", ["projplot", "projscatter", "projtext"])
def test_warns_on_axes_not_a_healpy_map(func):
    fig = plt.figure()
    fig.add_subplot(111)
    args = {
        "projplot": ([np.radians(45)], [np.radians(60)]),
        "projscatter": (np.radians(45), np.radians(60)),
        "projtext": (np.radians(45), np.radians(60), "x"),
    }[func]
    with pytest.warns(RuntimeWarning, match="not a healpy map"):
        ret = getattr(hp, func)(*args)
    # no axes to draw on: nothing is drawn, still returns None
    assert ret is None
