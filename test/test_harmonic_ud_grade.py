import numpy as np

import healpy as hp


def _single_mode_map(nside, ell, emm=0):
    alm = np.zeros(hp.Alm.getsize(ell), dtype=np.complex128)
    alm[hp.Alm.getidx(ell, ell, emm)] = 1.0
    return hp.alm2map(alm, nside=nside, pixwin=False)


def test_harmonic_ud_grade_preserves_low_ell_mode():
    nside_in = 128
    nside_out = 32
    input_map = _single_mode_map(nside_in, ell=20)
    expected = _single_mode_map(nside_out, ell=20)

    output = hp.harmonic_ud_grade(input_map, nside_out=nside_out)

    np.testing.assert_allclose(output, expected, rtol=1e-3, atol=1e-6)


def test_harmonic_ud_grade_suppresses_aliasing_vs_ud_grade():
    nside_in = 128
    nside_out = 32
    input_map = _single_mode_map(nside_in, ell=120)

    output_harmonic = hp.harmonic_ud_grade(input_map, nside_out=nside_out)
    output_ud_grade = hp.ud_grade(input_map, nside_out=nside_out)

    std_harmonic = np.std(output_harmonic)
    std_ud_grade = np.std(output_ud_grade)

    assert std_harmonic < 0.2 * std_ud_grade


def test_harmonic_ud_grade_multimap_pol_false_and_dtype():
    nside_in = 64
    nside_out = 16
    npix = hp.nside2npix(nside_in)
    rng = np.random.default_rng(123)
    input_maps = rng.normal(size=(2, npix))

    output = hp.harmonic_ud_grade(
        input_maps,
        nside_out=nside_out,
        pol=False,
        dtype=np.float32,
    )

    assert output.shape == (2, hp.nside2npix(nside_out))
    assert output.dtype == np.float32
