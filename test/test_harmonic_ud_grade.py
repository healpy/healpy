import numpy as np
import pytest

import healpy as hp


def _single_mode_map(nside, ell, emm=0):
    alm = np.zeros(hp.Alm.getsize(ell), dtype=np.complex128)
    alm[hp.Alm.getidx(ell, ell, emm)] = 1.0
    return hp.alm2map(alm, nside=nside, pixwin=False)


# ------------------------------------------------------------------
# Backward-compatible tests (plain bandlimit clipping, no correction)
# ------------------------------------------------------------------


def test_harmonic_ud_grade_preserves_low_ell_mode():
    """Plain clipping preserves a low-ell mode below the output bandlimit."""
    nside_in = 128
    nside_out = 32
    input_map = _single_mode_map(nside_in, ell=20)
    expected = _single_mode_map(nside_out, ell=20)

    output = hp.harmonic_ud_grade(
        input_map, nside_out=nside_out, use_pixel_weights=False,
        pixwin=False, fwhm_out=0,
    )

    np.testing.assert_allclose(output, expected, rtol=1e-3, atol=1e-6)


def test_harmonic_ud_grade_suppresses_aliasing_vs_ud_grade():
    nside_in = 128
    nside_out = 32
    input_map = _single_mode_map(nside_in, ell=120)

    output_harmonic = hp.harmonic_ud_grade(
        input_map, nside_out=nside_out, use_pixel_weights=False,
        pixwin=False, fwhm_out=0,
    )
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
        use_pixel_weights=False,
        pixwin=False,
        fwhm_out=0,
        dtype=np.float32,
    )

    assert output.shape == (2, hp.nside2npix(nside_out))
    assert output.dtype == np.float32


def test_harmonic_ud_grade_missing_pixel_weights_file_raises(tmp_path):
    nside_in = 32
    nside_out = 16
    input_map = np.ones(hp.nside2npix(nside_in), dtype=np.float64)
    expected = tmp_path / "full_weights" / "healpix_full_weights_nside_0032.fits"

    with pytest.raises(RuntimeError, match="Pixel weights are required by default"):
        hp.harmonic_ud_grade(
            input_map,
            nside_out=nside_out,
            use_pixel_weights=True,
            datapath=str(tmp_path),
        )

    with pytest.raises(RuntimeError, match=str(expected)):
        hp.harmonic_ud_grade(
            input_map,
            nside_out=nside_out,
            use_pixel_weights=True,
            datapath=str(tmp_path),
        )

    with pytest.raises(RuntimeError, match="use_pixel_weights=False"):
        hp.harmonic_ud_grade(
            input_map,
            nside_out=nside_out,
            use_pixel_weights=True,
            datapath=str(tmp_path),
        )


# ------------------------------------------------------------------
# New tests for pixwin and fwhm_out parameters
# ------------------------------------------------------------------


def test_pixwin_changes_result():
    """pixwin=True produces a different result from pixwin=False."""
    nside_in = 128
    nside_out = 32
    input_map = _single_mode_map(nside_in, ell=20)

    out_no_pw = hp.harmonic_ud_grade(
        input_map, nside_out=nside_out, use_pixel_weights=False,
        pixwin=False, fwhm_out=0,
    )
    out_pw = hp.harmonic_ud_grade(
        input_map, nside_out=nside_out, use_pixel_weights=False,
        pixwin=True, fwhm_out=0,
    )

    assert not np.allclose(out_no_pw, out_pw, rtol=1e-4), \
        "pixwin correction should change the result"


def test_fwhm_out_smooths_output():
    """fwhm_out > 0 should suppress high-ell power compared to fwhm_out=0."""
    nside_in = 128
    nside_out = 32
    lmax_out = 3 * nside_out - 1
    rng = np.random.default_rng(42)
    cl = np.zeros(3 * nside_in)
    cl[2:] = np.arange(2, 3 * nside_in)[..., ::-1].astype(float)
    input_map = hp.synfast(cl, nside_in, new=True)

    out_no_beam = hp.harmonic_ud_grade(
        input_map, nside_out=nside_out, use_pixel_weights=False,
        pixwin=False, fwhm_out=0,
    )
    out_beam = hp.harmonic_ud_grade(
        input_map, nside_out=nside_out, use_pixel_weights=False,
        pixwin=False, fwhm_out=3 * hp.nside2resol(nside_out),
    )

    cl_no_beam = hp.anafast(out_no_beam, lmax=lmax_out)
    cl_beam = hp.anafast(out_beam, lmax=lmax_out)

    # High-ell power should be suppressed by the beam
    high_ell = slice(lmax_out // 2, None)
    assert np.mean(cl_beam[high_ell]) < 0.5 * np.mean(cl_no_beam[high_ell]), \
        "Beam should suppress high-ell power"


def test_fwhm_out_default_applies_beam():
    """Default fwhm_out=None auto-computes 3*nside2resol and differs from 0."""
    nside_in = 128
    nside_out = 32
    input_map = _single_mode_map(nside_in, ell=20)

    out_default = hp.harmonic_ud_grade(
        input_map, nside_out=nside_out, use_pixel_weights=False,
        pixwin=False,  # fwhm_out defaults to None = auto
    )
    out_no_beam = hp.harmonic_ud_grade(
        input_map, nside_out=nside_out, use_pixel_weights=False,
        pixwin=False, fwhm_out=0,
    )

    assert not np.allclose(out_default, out_no_beam, rtol=1e-4), \
        "Default fwhm_out should apply a beam"


def test_fwhm_in_deconvolves_input_beam():
    """fwhm_in should deconvolve, partially undoing prior smoothing."""
    nside_in = 128
    nside_out = 32
    fwhm = 0.05  # radians

    rng = np.random.default_rng(7)
    cl = np.zeros(3 * nside_in)
    cl[2:] = 1.0
    input_map = hp.synfast(cl, nside_in, new=True)

    # Smooth the input to simulate a beamed map
    beamed_map = hp.smoothing(input_map, fwhm=fwhm)

    # Downgrade with deconvolution of the input beam
    out_deconv = hp.harmonic_ud_grade(
        beamed_map, nside_out=nside_out, use_pixel_weights=False,
        pixwin=False, fwhm_in=fwhm, fwhm_out=0,
    )
    # Downgrade without deconvolution
    out_no_deconv = hp.harmonic_ud_grade(
        beamed_map, nside_out=nside_out, use_pixel_weights=False,
        pixwin=False, fwhm_in=0, fwhm_out=0,
    )

    # Deconvolved output should have more high-ell power
    lmax_out = 3 * nside_out - 1
    cl_deconv = hp.anafast(out_deconv, lmax=lmax_out)
    cl_no_deconv = hp.anafast(out_no_deconv, lmax=lmax_out)

    high_ell = slice(lmax_out // 2, None)
    assert np.mean(cl_deconv[high_ell]) > 1.5 * np.mean(cl_no_deconv[high_ell]), \
        "Deconvolving input beam should restore high-ell power"

def test_harmonic_ud_grade_upgrading_caps_lmax():
    """When upgrading, lmax should be capped to the input resolution limit."""
    nside_in = 16
    nside_out = 32
    # If lmax is not capped, it will try to extract up to 3*32-1=95
    # which map2alm will complain about or return noisy modes.
    # The fix ensures lmax is min(3*32-1, 3*16-1) = 47.
    input_map = _single_mode_map(nside_in, ell=10)
    
    # Should not raise an error or complain about lmax > lmax_in
    output = hp.harmonic_ud_grade(
        input_map, nside_out=nside_out, use_pixel_weights=False,
        pixwin=False, fwhm_out=0,
    )
    assert hp.get_nside(output) == nside_out


def test_harmonic_ud_grade_pol_beam():
    """Polarized beam deconvolution should use the polarized beam transfer function."""
    nside_in = 32
    nside_out = 32
    npix = hp.nside2npix(nside_in)
    rng = np.random.default_rng(123)
    input_maps = rng.normal(size=(3, npix))
    
    # Run with a large beam to exaggerate differences
    out_pol = hp.harmonic_ud_grade(
        input_maps, nside_out=nside_out, pol=True,
        fwhm_in=np.radians(10.0), fwhm_out=np.radians(10.0),
        use_pixel_weights=False, pixwin=False
    )
    out_unpol = hp.harmonic_ud_grade(
        input_maps, nside_out=nside_out, pol=False,
        fwhm_in=np.radians(10.0), fwhm_out=np.radians(10.0),
        use_pixel_weights=False, pixwin=False
    )
    
    # T component should be identical
    np.testing.assert_allclose(out_pol[0], out_unpol[0], rtol=1e-5)
    
    # Q and U components should be slightly different because out_pol used the E/B beam
    assert not np.allclose(out_pol[1], out_unpol[1], rtol=1e-5)
    assert not np.allclose(out_pol[2], out_unpol[2], rtol=1e-5)

