from __future__ import division

import os.path

import pytest

import numpy as np

import healpy as hp
from healpy import Rotator
from healpy.rotator import euler, euler_matrix_new

path = os.path.dirname(os.path.realpath(__file__))

# A pytest fixture with autouse=True is run before each of the other tests
@pytest.fixture(autouse=True)
def set_random_seed():
    """Set same seed so tests are reproducible"""
    seed = 12345
    np.random.seed(seed)


def test_rotate_map_polarization():
    """Compare to a rotation from Galactic to Ecliptic of
    a map of pure Q polarization, the expected value was computed with HEALPix IDL:
    https://gist.github.com/zonca/401069e1c520e02eaff8cd86149d5900
    """
    nside = 32
    npix = hp.nside2npix(nside)
    QU_gal = np.zeros((2, npix), dtype=np.double)
    QU_gal[0, : npix // 2] = 1
    gal2ecl = Rotator(coord=["G", "E"])
    QU_ecl = gal2ecl.rotate_map_pixel(QU_gal)

    expected = hp.ma(
        hp.read_map(os.path.join(path, "data", "justq_gal2ecl.fits.gz"), [0, 1])
    )

    expected.mask = expected == 0

    for i_pol, pol in enumerate("QU"):
        assert (
            (np.abs(expected[i_pol] - QU_ecl[i_pol]) < 0.05).sum()
            / np.logical_not(expected[i_pol].mask).sum()
        ) > 0.9, (pol + " comparison failed in rotate_map")


def test_rotate_map_polarization_alms():
    lmax = 64
    path = os.path.dirname(os.path.realpath(__file__))
    map1 = hp.read_map(
        os.path.join(path, "data", "wmap_band_iqumap_r9_7yr_W_v4_udgraded32.fits"),
        (0, 1, 2),
    )

    # do the rotation with hp.rotate_alm
    angles = hp.rotator.coordsys2euler_zyz(coord=["G", "E"])
    alm = hp.map2alm(map1, lmax=lmax, use_pixel_weights=True)
    hp.rotate_alm(alm, *angles)
    rotated_map1 = hp.alm2map(alm, nside=hp.get_nside(map1), lmax=lmax)

    # do the rotation with hp.Rotator
    gal2ecl = hp.Rotator(coord=["G", "E"])
    rotate_map1_rotate_map_alms = gal2ecl.rotate_map_alms(map1, lmax=lmax)

    np.testing.assert_allclose(rotate_map1_rotate_map_alms, rotated_map1, rtol=1e-5)


def test_rotate_map_polarization_with_spectrum():
    """Rotation of reference frame should not change the angular power spectrum.
    In this test we create a map from a spectrum with a pure EE signal and check
    that the spectrum of this map and the spectrum of the same map rotated
    from Galactic to Ecliptic agrees.
    This test checks if the QU rotation is correct"""
    nside = 32
    cl = np.zeros((6, 96), dtype=np.double)
    # Set ell=1 for EE to 1
    cl[1][2] = 1
    gal2ecl = Rotator(coord=["G", "E"])
    m_original = hp.synfast(cl, nside=nside, new=True)
    cl_from_m_original = hp.anafast(m_original)
    m_rotated = gal2ecl.rotate_map_pixel(m_original)
    cl_from_m_rotated = hp.anafast(m_rotated)

    assert np.abs(cl_from_m_rotated - cl_from_m_original).sum() < 1e-2


def test_rotate_dipole_and_back():
    """Rotate a smooth signal (dipole) from Galactic to Ecliptic and back"""
    nside = 64
    npix = hp.nside2npix(nside)
    pix = np.arange(npix)
    vec = np.array(hp.pix2vec(nside, pix))
    # dipole max is at North Pole
    dip_dir = np.array([0, 0, 1])
    m_gal = np.dot(vec.T, dip_dir)
    gal2ecl = Rotator(coord=["C", "E"])
    ecl2gal = Rotator(coord=["E", "C"])
    m_ecl = gal2ecl.rotate_map_pixel(m_gal)
    # Remove 10 deg along equator because dipole signal is so low that relative error is
    # too large
    cut_equator_deg = 5
    no_equator = hp.query_strip(
        nside, np.radians(90 + cut_equator_deg), np.radians(90 - cut_equator_deg)
    )
    np.testing.assert_allclose(
        m_gal[no_equator], ecl2gal.rotate_map_pixel(m_ecl)[no_equator], rtol=1e-3
    )


def test_rotator_input_lengths():
    with pytest.raises(ValueError):
        Rotator(coord=[("C", "E"), ("E", "G")], rot=[(0, 0, 90)])


def test_rotator_input_type():
    with pytest.raises(ValueError):
        Rotator(coord="CE", rot=[(0, 0, 90)])


def test_rotator_input_lengths_inv():
    with pytest.raises(ValueError):
        Rotator(
            coord=[("C", "E"), ("E", "G")], rot=[(0, 0, 90), (0, 90, 0)], inv=[True]
        )


def test_rotator_eq():
    rot_1 = Rotator(coord=("G", "E"))
    assert rot_1 == rot_1.get_inverse().get_inverse()


def test_rotate_vector():
    gal2ecl = Rotator(coord=("G", "E"))
    gal_vec = np.array([0.1, -0.1, 1])
    gal_vec_in_ecl = gal2ecl(*gal_vec)
    np.testing.assert_allclose(gal_vec, gal2ecl.I(gal_vec_in_ecl))


@pytest.mark.parametrize("select", list(range(1, 6 + 1)))
@pytest.mark.parametrize("FK4", [0, 1])
def test_euler(select, FK4):
    out = euler(30, 20, select=select, FK4=FK4)
    np.testing.assert_array_equal(np.isnan(out), 0)


@pytest.mark.parametrize("X,Y,ZYX", [(1, 0, 0), (0, 1, 0), (0, 0, 1)])
def test_euler_matrix_new(X, Y, ZYX):
    out = euler_matrix_new(10, 10, 10, X=X, Y=Y, ZYX=ZYX, deg=True)
    np.testing.assert_array_equal(np.isnan(out), 0)
