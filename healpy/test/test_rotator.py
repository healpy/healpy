import os.path

import numpy as np

import healpy as hp
from healpy import Rotator

path = os.path.dirname(os.path.realpath(__file__))


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
    QU_ecl = gal2ecl.rotate_map(QU_gal)

    expected = hp.ma(
        hp.read_map(os.path.join(path, "data", "justq_gal2ecl.fits.gz"), [0, 1])
    )

    expected.mask = expected == 0

    for i_pol, pol in enumerate("QU"):
        assert (
            (np.abs(expected[i_pol] - QU_ecl[i_pol]) < .05).sum()
            / (~expected[i_pol].mask).sum()
        ) > .9, (i_pol + " comparison failed in rotate_map")


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
    m_ecl = gal2ecl.rotate_map(m_gal)
    # Remove 10 deg along equator because dipole signal is so low that relative error is
    # too large
    cut_equator_deg = 5
    no_equator = hp.query_strip(
        nside, np.radians(90 + cut_equator_deg), np.radians(90 - cut_equator_deg)
    )
    np.testing.assert_allclose(
        m_gal[no_equator], ecl2gal.rotate_map(m_ecl)[no_equator], rtol=1e-3
    )
