import os.path

import numpy as np

import healpy as hp
from healpy import Rotator

path = os.path.dirname( os.path.realpath( __file__ ) )

def test_rotate_map_polarization():
    nside = 16
    npix = hp.nside2npix(nside)
    QU_gal = [np.ones(npix, dtype=np.double), np.zeros(npix, dtype=np.double)]
    gal2ecl = Rotator(coord=["G","E"])
    QU_ecl = gal2ecl.rotate_map(QU_gal)

    expected = hp.read_map(os.path.join(path, "data", "justq_gal2ecl.fits.gz"), [0,1])

    np.testing.assert_allclose(QU_ecl[0], expected[0], rtol=1e-4)
    np.testing.assert_allclose(QU_ecl[1], expected[1], rtol=1e-4)
