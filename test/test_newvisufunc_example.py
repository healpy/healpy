import healpy as hp
from healpy.newvisufunc import projview, newprojplot
import matplotlib.pyplot as plt
import numpy as np
import pytest
import os.path

path = os.path.dirname(os.path.realpath(__file__))


@pytest.fixture
def map_data():
    return hp.read_map(
        os.path.join(
            path,
            "data",
            "wmap_band_iqumap_r9_7yr_W_v4_udgraded32_masked_smoothed10deg_fortran.fits",
        )
    )


def test_projview_mollweide(map_data):
    projview(map_data, coord=["G"], flip="astro", projection_type="mollweide")


def test_projview_mollweide_graticule(map_data):
    projview(
        map_data,
        coord=["G"],
        graticule=True,
        graticule_labels=True,
        projection_type="mollweide",
    )


def test_projview_mollweide_graticule_labels(map_data):
    projview(
        map_data,
        coord=["G"],
        graticule=True,
        graticule_labels=True,
        unit="cbar label",
        xlabel="longitude",
        ylabel="latitude",
        cb_orientation="horizontal",
        projection_type="mollweide",
    )


def test_projview_mollweide_vertical_cbar(map_data):
    projview(
        map_data,
        coord=["G"],
        graticule=True,
        graticule_labels=True,
        unit="cbar label",
        xlabel="longitude",
        ylabel="latitude",
        cb_orientation="vertical",
        projection_type="mollweide",
    )


def test_projview_mollweide_extended_cbar(map_data):
    projview(
        map_data,
        coord=["G"],
        graticule=True,
        graticule_labels=True,
        unit="cbar label",
        xlabel="longitude",
        ylabel="latitude",
        cb_orientation="vertical",
        min=-0.05,
        max=0.05,
        latitude_grid_spacing=45,
        projection_type="mollweide",
        title="Mollweide projection, astro convention (default)",
    )
    newprojplot(
        theta=np.radians(50), phi=np.radians(60), marker="o", color="r", markersize=10
    )


def test_projview_mollweide_geo_convention(map_data):
    projview(
        map_data,
        coord=["G"],
        graticule=True,
        graticule_labels=True,
        unit="cbar label",
        xlabel="longitude",
        ylabel="latitude",
        cb_orientation="vertical",
        min=-0.05,
        max=0.05,
        latitude_grid_spacing=45,
        projection_type="mollweide",
        title="Mollweide projection, geo convention",
        flip="geo",
        phi_convention="clockwise",
    )
    newprojplot(
        theta=np.radians(50), phi=np.radians(60), marker="o", color="r", markersize=10
    )


def test_projview_hammer(map_data):
    projview(
        map_data,
        coord=["G"],
        graticule=True,
        graticule_labels=True,
        unit="cbar label",
        xlabel="longitude",
        ylabel="latitude",
        cb_orientation="vertical",
        min=-0.05,
        max=0.05,
        latitude_grid_spacing=30,
        projection_type="hammer",
        title="Hammer projection",
    )


def test_projview_aitoff(map_data):
    projview(
        map_data,
        coord=["G"],
        graticule=True,
        graticule_labels=True,
        unit="cbar label",
        xlabel="longitude",
        ylabel="latitude",
        cb_orientation="vertical",
        min=-0.05,
        max=0.05,
        latitude_grid_spacing=30,
        projection_type="aitoff",
        title="Aitoff projection",
    )


def test_projview_cart(map_data):
    projview(map_data, coord=["G"], projection_type="cart")


def test_projview_cart_labels(map_data):
    projview(
        map_data, coord=["G"], projection_type="cart", xlabel="xlabel", ylabel="ylabel"
    )


def test_projview_cart_graticule(map_data):
    projview(
        map_data,
        coord=["G"],
        graticule=True,
        graticule_labels=True,
        unit="cbar label",
        xlabel="longitude",
        ylabel="latitude",
        cb_orientation="horizontal",
        projection_type="cart",
        title="Cart projection",
    )


def test_projview_cart_graticule_vertical_cbar(map_data):
    projview(
        map_data,
        coord=["G"],
        graticule=True,
        graticule_labels=True,
        unit="cbar label",
        xlabel="longitude",
        ylabel="latitude",
        cb_orientation="vertical",
        projection_type="cart",
    )


def test_projview_3d(map_data):
    projview(
        map_data,
        coord=["G"],
        hold=False,
        graticule=True,
        graticule_labels=True,
        unit="cbar label",
        xlabel="longitude",
        ylabel="latitude",
        cb_orientation="horizontal",
        projection_type="3d",
        title="3D projection",
    )


def test_projview_3d_vertical_cbar(map_data):
    projview(
        map_data,
        coord=["G"],
        hold=False,
        graticule=True,
        graticule_labels=True,
        projection_type="3d",
        unit="cbar label",
        xlabel="longitude",
        ylabel="latitude",
        cb_orientation="vertical",
        cmap="viridis",
    )


def test_projview_polar(map_data):
    projview(
        map_data,
        coord=["G"],
        graticule=True,
        graticule_labels=True,
        unit="cbar label",
        cb_orientation="horizontal",
        projection_type="polar",
        title="Polar projection",
    )


def test_projview_polar_vertical_cbar(map_data):
    projview(
        map_data,
        coord=["G"],
        graticule=True,
        graticule_labels=True,
        unit="cbar label",
        cb_orientation="vertical",
        projection_type="polar",
    )


def test_projview_polar_override(map_data):
    projview(
        map_data,
        coord=["G"],
        hold=False,
        graticule=True,
        graticule_labels=True,
        flip="astro",
        projection_type="polar",
        unit="cbar label",
        cb_orientation="horizontal",
        override_plot_properties={
            "cbar_shrink": 0.5,
            "cbar_pad": 0.02,
            "cbar_label_pad": -35,
            "figure_width": 16,
            "figure_size_ratio": 0.63,
        },
    )


def test_projview_hammer_fontsize(map_data):
    projview(
        map_data,
        coord=["G"],
        graticule=True,
        graticule_labels=True,
        unit=r"cbar label $\alpha$",
        xlabel="longitude",
        ylabel="latitude",
        cb_orientation="vertical",
        min=-0.05,
        max=0.05,
        latitude_grid_spacing=30,
        projection_type="hammer",
        title="Hammer projection",
        fontsize={
            "xlabel": 20,
            "ylabel": 20,
            "xtick_label": 20,
            "ytick_label": 20,
            "title": 20,
            "cbar_label": 20,
            "cbar_tick_label": 20,
        },
        xtick_label_color="r",
        ytick_label_color="g",
        graticule_color="black",
    )


def test_projview_hammer_no_phi_tick_label_shift(map_data):
    projview(
        map_data,
        coord=["G"],
        graticule=True,
        graticule_labels=True,
        unit=r"cbar label",
        xlabel="longitude",
        ylabel="latitude",
        cb_orientation="vertical",
        min=-0.05,
        max=0.05,
        latitude_grid_spacing=30,
        projection_type="hammer",
        title="Hammer projection",
        phi_convention="symmetrical",
    )
    s = 500
    plt.scatter(np.deg2rad(0), np.deg2rad(0), color="r", marker="x", linewidth=10, s=s)
    plt.scatter(
        np.deg2rad(120), np.deg2rad(0), color="r", marker="x", linewidth=10, s=s
    )
    plt.scatter(
        np.deg2rad(-120), np.deg2rad(0), color="r", marker="x", linewidth=10, s=s
    )
    plt.scatter(np.deg2rad(0), np.deg2rad(60), color="r", marker="x", linewidth=10, s=s)
    plt.scatter(
        np.deg2rad(0), np.deg2rad(-60), color="r", marker="x", linewidth=10, s=s
    )


def test_projview_hammer_override_axis_tick_labels(map_data):
    projview(
        map_data,
        coord=["G"],
        graticule=True,
        graticule_labels=True,
        unit=r"cbar label",
        xlabel="longitude",
        ylabel="latitude",
        cb_orientation="vertical",
        min=-0.05,
        max=0.05,
        latitude_grid_spacing=30,
        projection_type="hammer",
        title="Hammer projection",
        custom_xtick_labels=["A", "B", "C", "D", "E"],
        custom_ytick_labels=["F", "G", "H", "I", "J"],
    )


def test_projview_hammer_equatorial(map_data):
    projview(
        map_data,
        coord=["G", "C"],
        graticule=True,
        graticule_labels=True,
        unit=r"cbar label",
        xlabel="RA",
        ylabel="DEC",
        cb_orientation="vertical",
        min=-0.05,
        max=0.05,
        latitude_grid_spacing=30,
        projection_type="hammer",
        title="Hammer projection",
    )


def test_projview_cart_local_azimuth_counterclockwise(map_data):
    local_sidereal_time = 18
    altitude = -35.206667
    rotAngles = [(180 + (local_sidereal_time * 15)) % 360, -(altitude - 90)]
    projview(
        map_data,
        coord=["G", "C"],
        graticule=True,
        graticule_labels=True,
        unit=r"cbar label",
        xlabel="azimuth",
        ylabel="zenith",
        cb_orientation="vertical",
        min=-0.05,
        max=0.05,
        latitude_grid_spacing=30,
        projection_type="cart",
        title="Cartesian projection",
        rot=rotAngles,
        fontsize={"xtick_label": 20},
        phi_convention="counterclockwise",
    )


def test_projview_cart_local_azimuth_clockwise(map_data):
    local_sidereal_time = 18
    altitude = -35.206667
    rotAngles = [(180 + (local_sidereal_time * 15)) % 360, -(altitude - 90)]
    projview(
        map_data,
        coord=["G", "C"],
        graticule=True,
        graticule_labels=True,
        unit=r"cbar label",
        xlabel="azimuth",
        ylabel="zenith",
        cb_orientation="vertical",
        min=-0.05,
        max=0.05,
        latitude_grid_spacing=30,
        projection_type="cart",
        title="Cartesian projection",
        rot=rotAngles,
        fontsize={"xtick_label": 20},
        phi_convention="clockwise",
    )


def test_projview_cart_local_azimuth_symmetrical(map_data):
    local_sidereal_time = 18
    altitude = -35.206667
    rotAngles = [(180 + (local_sidereal_time * 15)) % 360, -(altitude - 90)]
    projview(
        map_data,
        coord=["G", "C"],
        graticule=True,
        graticule_labels=True,
        unit=r"cbar label",
        xlabel="azimuth",
        ylabel="zenith",
        cb_orientation="vertical",
        min=-0.05,
        max=0.05,
        latitude_grid_spacing=30,
        projection_type="cart",
        title="Cartesian projection",
        rot=rotAngles,
        fontsize={"xtick_label": 20},
        phi_convention="symmetrical",
    )


def test_projview_return_only_data(map_data):
    longitude, latitude, grid_map = projview(
        map_data, coord=["G"], return_only_data=True
    )


def test_projview_planck_colormap(map_data):
    m_scaled = map_data * 3000
    projview(
        m_scaled,
        title="Planck colormap",
        cmap="planck",
        rlabel=r"A$_{\mathsf{ex. 1}}$",
        llabel=r"$Q$",
        unit=r"$\mu$K",
        fontname="serif",
        width=10,
        show_tickmarkers=True,
        cbar_ticks=[-300, 0, 300],
        cb_orientation="vertical",
        sub=121,
        override_plot_properties={"cbar_tick_direction": "in"},
    )
    projview(
        m_scaled,
        title="Planck logarithmic colormap",
        cmap="planck_log",
        norm="symlog2",
        rlabel=r"A$_{\mathrm{ex. 3}}$",
        llabel=r"$I$",
        unit=r"$\mu$K",
        min=-1e3,
        max=1e7,
        cb_orientation="vertical",
        sub=122,
    )
    plt.tight_layout()


def test_projview_symlog_normalization(map_data):
    m_scaled = map_data * 3000
    projview(
        m_scaled,
        title="symlog normalization",
        cmap="planck",
        norm="symlog",
        rlabel=r"A$_{\mathsf{ex. 2}}$",
        llabel=r"$Q$",
        unit=r"$\mu$K",
        cbar_ticks=[-3000, -30, 0, 30, 3000],
        remove_mono=True,
        show_tickmarkers=True,
        sub=121,
        override_plot_properties={"cbar_tick_direction": "in"},
        norm_dict={"linscale": 0.5},
    )
    projview(
        m_scaled,
        title="WMAP colormap",
        cmap="wmap",
        rlabel=r"A$_{\mathrm{ex. 4}}$",
        llabel=r"$I$",
        unit=r"$\mu$K",
        fontname="serif",
        min=-100,
        max=100,
        sub=122,
    )
    plt.tight_layout()

    def test_projview_mollweide_badcolor(map_data):
        projview(
            map_data,
            coord=["G"],
            graticule=True,
            graticule_labels=True,
            unit="cbar label",
            xlabel="longitude",
            ylabel="latitude",
            cb_orientation="vertical",
            min=-0.05,
            max=0.05,
            latitude_grid_spacing=45,
            projection_type="mollweide",
            title="Mollweide projection with badcolor",
            badcolor="red",
        )

    def test_projview_cart_bgcolor(map_data):
        projview(
            map_data,
            coord=["G"],
            graticule=True,
            graticule_labels=True,
            unit="cbar label",
            xlabel="longitude",
            ylabel="latitude",
            cb_orientation="horizontal",
            projection_type="cart",
            title="Cart projection with bgcolor",
            bgcolor="lightblue",
        )
