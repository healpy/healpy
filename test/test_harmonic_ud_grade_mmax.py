import numpy as np
import pytest

import healpy as hp


@pytest.mark.parametrize("input_type", ["map", "alm"])
@pytest.mark.parametrize("mmax", [0, 1, 3, 5])
@pytest.mark.parametrize("filtered", [False, True])
@pytest.mark.parametrize("polarized", [False, True])
def test_harmonic_ud_grade_reduced_mmax(input_type, mmax, filtered, polarized):
    lmax = 5
    alm = np.arange(1, hp.Alm.getsize(lmax) + 1, dtype=np.complex128)
    if polarized:
        alm = np.array([alm, 2 * alm, 3 * alm])
    maps = hp.alm2map(alm, 4, lmax=lmax, pol=polarized)
    source = maps if input_type == "map" else alm
    original = source.copy()
    if input_type == "map":
        limited = hp.map2alm(maps, lmax=lmax, mmax=mmax, iter=0, pol=polarized)
    else:
        limited = hp.resize_alm(alm, lmax, lmax, lmax, mmax)
    fwhm = 0.1 if filtered else 0
    if filtered:
        beam = hp.gauss_beam(fwhm, lmax=lmax, pol=polarized)
        if polarized:
            limited = [hp.almxfl(a, beam[:, min(i, 1)], mmax=mmax) for i, a in enumerate(limited)]
        else:
            limited = hp.almxfl(limited, beam, mmax=mmax)
    expected = hp.alm2map(limited, 2, lmax=lmax, mmax=mmax, pol=polarized)
    kwargs = {"nside_in": 4} if input_type == "alm" else {"iter": 0, "use_pixel_weights": False}
    actual = hp.harmonic_ud_grade(
        source, 2, lmax=lmax, mmax=mmax, input_type=input_type,
        pixwin=False, fwhm_out=fwhm, pol=polarized, **kwargs
    )
    np.testing.assert_allclose(actual, expected, rtol=1e-13, atol=1e-13)
    np.testing.assert_array_equal(source, original)
