import re

import numpy as np
import pytest

import healpy as hp


def _single_mode_map(nside, ell, emm=0):
    """Create a map containing a single spherical-harmonic mode.

    Generates a HEALPix map at *nside* with a single (ell, emm) mode
    set to unit amplitude.  Useful for controlled tests where the
    expected output can be computed analytically.
    """
    alm = np.zeros(hp.Alm.getsize(ell), dtype=np.complex128)
    alm[hp.Alm.getidx(ell, ell, emm)] = 1.0
    return hp.alm2map(alm, nside=nside, pixwin=False)


# ------------------------------------------------------------------
# Basic functionality tests
# ------------------------------------------------------------------


def test_harmonic_ud_grade_preserves_low_ell_mode():
    """Plain bandlimit truncation preserves a low-ell mode below the output Nyquist.

    A single mode at ell=20 (well below the output lmax of 3*32-1=95)
    should survive the downgrade from NSIDE 128 to 32 with sub-percent
    relative error, since no aliasing or beam correction is involved.
    """
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
    """harmonic_ud_grade suppresses aliasing that ud_grade introduces for high-ell modes.

    A single mode at ell=120 (above the output NSIDE 32 bandlimit of
    ~95) will alias strongly under pixel averaging (ud_grade) but should
    be suppressed to <20% of the ud_grade amplitude by harmonic
    bandlimit truncation.
    """
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
    """Multi-map input with pol=False processes each map independently and respects dtype.

    Two random maps downgraded from NSIDE 64 to 16 should produce a
    (2, Npix_out) output in float32 when dtype=np.float32 is specified.
    With pol=False, each map is treated as spin-0 independently (not
    as TQU components).
    """
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
    """Missing pixel-weights file with use_pixel_weights=True raises RuntimeError.

    When pixel weights are required but the weight file does not exist
    in the specified datapath, the function must raise RuntimeError with
    a helpful message that mentions the expected file path and suggests
    passing use_pixel_weights=False.
    """
    nside_in = 32
    nside_out = 16
    input_map = np.ones(hp.nside2npix(nside_in), dtype=np.float64)
    expected_path = tmp_path / "full_weights" / "healpix_full_weights_nside_0032.fits"

    with pytest.raises(
        RuntimeError,
        match=r"Pixel weights are required by default.*"
        + re.escape(str(expected_path))
        + r".*use_pixel_weights=False",
    ):
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
    """pixwin=True produces a different result from pixwin=False.

    The pixel-window correction deconvolves the input pixel window and
    applies the output pixel window.  For a single-mode map, this
    correction should visibly change the output because the pixel
    windows at different NSIDEs have different shapes.
    """
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
    """fwhm_out > 0 should suppress high-ell power compared to fwhm_out=0.

    Applying a Gaussian output beam (fwhm_out = 3 * nside2resol) should
    reduce power at high multipoles relative to plain bandlimit
    truncation (fwhm_out=0).  This verifies that the output beam is
    actually being applied to the harmonic coefficients.
    """
    nside_in = 128
    nside_out = 32
    lmax_out = 3 * nside_out - 1
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
    """Passing fwhm_out=None auto-computes the effective resolution beam and differs from 0.

    The effective resolution beam (fwhm_out=None) computes
    effective_resolution_fwhm(nside_out), which should produce a visibly
    different (smoother) result than fwhm_out=0 (plain truncation).
    """
    nside_in = 128
    nside_out = 32
    input_map = _single_mode_map(nside_in, ell=20)

    out_auto = hp.harmonic_ud_grade(
        input_map, nside_out=nside_out, use_pixel_weights=False,
        pixwin=False, fwhm_out=None,  # effective resolution beam
    )
    out_no_beam = hp.harmonic_ud_grade(
        input_map, nside_out=nside_out, use_pixel_weights=False,
        pixwin=False, fwhm_out=0,  # no output beam (default)
    )

    assert not np.allclose(out_auto, out_no_beam, rtol=1e-4), \
        "fwhm_out=None should apply a Planck-scaled beam"


def test_fwhm_in_deconvolves_input_beam():
    """fwhm_in should deconvolve the input beam, partially restoring high-ell power.

    A map that has been smoothed with a Gaussian beam will have
    suppressed high-ell power.  Downgrading with fwhm_in set to that
    beam's FWHM should partially undo the smoothing, restoring more
    high-ell power than downgrading without deconvolution.
    """
    nside_in = 128
    nside_out = 32
    fwhm = 0.05  # radians

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
    """When upgrading NSIDE, lmax is capped to the input resolution limit.

    Upgrading from NSIDE 16 to 32 should not try to extract modes above
    3*16-1=47.  Without capping, map2alm would attempt to recover
    modes up to 3*32-1=95 which do not exist in the input, causing
    errors or noisy results. The output must also preserve a low-ell
    mode rather than just having the correct nside.
    """
    nside_in = 16
    nside_out = 32
    ell = 10
    input_map = _single_mode_map(nside_in, ell=ell)
    expected = _single_mode_map(nside_out, ell=ell)

    output = hp.harmonic_ud_grade(
        input_map, nside_out=nside_out, use_pixel_weights=False,
        pixwin=False, fwhm_out=0,
    )
    assert hp.get_nside(output) == nside_out
    # The low-ell mode is below the input bandlimit and must survive
    # upgrading. The map2alm round-trip at coarse nside without pixel
    # weights leaves a ~1% pixelization residual, so use a 2% tolerance.
    np.testing.assert_allclose(output, expected, rtol=2e-2, atol=2e-2)


def test_harmonic_ud_grade_pol_beam():
    """Polarized beam deconvolution uses the E/B beam column for Q and U.

    For a 3-component (TQU) input with pol=True, the temperature
    component uses the T beam column while Q and U use the E/B column
    (which differs for large beams).  With pol=False, all three maps
    use the same T beam.  The T outputs should match; Q and U should
    differ between the two calls.
    """
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
    
    # Q and U components should differ because out_pol used the E/B beam
    assert not np.allclose(out_pol[1], out_unpol[1], rtol=1e-5)
    assert not np.allclose(out_pol[2], out_unpol[2], rtol=1e-5)


# ------------------------------------------------------------------
# Tests for beam_window_in / beam_window_out parameters
# ------------------------------------------------------------------


def test_beam_window_in_overrides_fwhm_in():
    """beam_window_in with a Gaussian array produces the same result as the matching fwhm_in.

    Passing beam_window_in=gauss_beam(fwhm, lmax) should be equivalent
    to passing fwhm_in=fwhm, verifying that the custom beam array path
    correctly overrides the FWHM path.
    """
    nside_in = 64
    nside_out = 16
    fwhm = 0.05  # radians
    lmax = 3 * nside_out - 1

    input_map = hp.synfast(np.ones(3 * nside_in), nside_in, new=True)

    beam_in = hp.gauss_beam(fwhm, lmax=lmax)

    out_fwhm = hp.harmonic_ud_grade(
        input_map, nside_out=nside_out, use_pixel_weights=False,
        pixwin=False, fwhm_in=fwhm, fwhm_out=0,
    )
    out_beam = hp.harmonic_ud_grade(
        input_map, nside_out=nside_out, use_pixel_weights=False,
        pixwin=False, fwhm_in=0, fwhm_out=0, beam_window_in=beam_in,
    )

    np.testing.assert_allclose(out_fwhm, out_beam, rtol=1e-10)


def test_beam_window_out_overrides_fwhm_out():
    """beam_window_out with a Gaussian array produces the same result as the matching fwhm_out.

    Passing beam_window_out=gauss_beam(fwhm, lmax) should be equivalent
    to passing fwhm_out=fwhm, verifying that the custom output beam
    array correctly overrides the FWHM parameter.
    """
    nside_in = 64
    nside_out = 16
    fwhm = 0.05  # radians
    lmax = 3 * nside_out - 1

    input_map = hp.synfast(np.ones(3 * nside_in), nside_in, new=True)

    beam_out = hp.gauss_beam(fwhm, lmax=lmax)

    out_fwhm = hp.harmonic_ud_grade(
        input_map, nside_out=nside_out, use_pixel_weights=False,
        pixwin=False, fwhm_out=fwhm,
    )
    out_beam = hp.harmonic_ud_grade(
        input_map, nside_out=nside_out, use_pixel_weights=False,
        pixwin=False, beam_window_out=beam_out,
    )

    np.testing.assert_allclose(out_fwhm, out_beam, rtol=1e-10)


def test_beam_window_custom_shape():
    """A non-Gaussian beam_window produces different results from any Gaussian.

    A top-hat beam (1 up to lmax/2, then 0) should differ from a
    Gaussian beam with similar width, verifying that arbitrary beam
    shapes are actually being applied rather than ignored.
    """
    nside_in = 64
    nside_out = 16
    lmax = 3 * nside_out - 1

    input_map = hp.synfast(np.ones(3 * nside_in), nside_in, new=True)

    # Top-hat beam: 1 up to lmax/2, then 0
    custom_beam = np.ones(lmax + 1)
    custom_beam[lmax // 2:] = 0.0

    out_custom = hp.harmonic_ud_grade(
        input_map, nside_out=nside_out, use_pixel_weights=False,
        pixwin=False, fwhm_out=0, beam_window_out=custom_beam,
    )
    # Gaussian beam with similar width
    out_gauss = hp.harmonic_ud_grade(
        input_map, nside_out=nside_out, use_pixel_weights=False,
        pixwin=False, fwhm_out=hp.nside2resol(nside_out),
    )

    assert not np.allclose(out_custom, out_gauss, rtol=1e-2), \
        "Custom beam should differ from Gaussian"


def test_beam_window_polarized_2d():
    """2D beam_window with gauss_beam(pol=True) format works for TQU maps.

    A 2D beam array from gauss_beam(fwhm, pol=True) — shape
    (lmax+1, 4) with [T, E/B, T→E, T→B] columns — should produce the
    same result as passing fwhm_in/fwhm_out for a polarized 3-component
    input with pol=True.
    """
    nside_in = 32
    nside_out = 32
    fwhm = np.radians(10.0)
    lmax = 3 * nside_out - 1

    rng = np.random.default_rng(77)
    input_maps = rng.normal(size=(3, hp.nside2npix(nside_in)))

    beam = hp.gauss_beam(fwhm, lmax=lmax, pol=True)

    out_fwhm = hp.harmonic_ud_grade(
        input_maps, nside_out=nside_out, pol=True,
        fwhm_in=fwhm, fwhm_out=fwhm,
        use_pixel_weights=False, pixwin=False,
    )
    out_beam = hp.harmonic_ud_grade(
        input_maps, nside_out=nside_out, pol=True,
        beam_window_in=beam, beam_window_out=beam,
        use_pixel_weights=False, pixwin=False,
    )

    np.testing.assert_allclose(out_fwhm, out_beam, rtol=1e-10)


def test_beam_window_1d_raises_for_polarized():
    """1D beam_window for a 3-component polarized map raises ValueError.

    A 1D beam array does not have separate T and E/B columns, so it
    cannot be used for polarized transforms.  The function must reject
    this with a clear error message.
    """
    nside_in = 32
    nside_out = 32
    lmax = 3 * nside_out - 1

    rng = np.random.default_rng(55)
    input_maps = rng.normal(size=(3, hp.nside2npix(nside_in)))

    beam_1d = hp.gauss_beam(np.radians(10.0), lmax=lmax, pol=False)

    with pytest.raises(ValueError, match="1D.*polarization"):
        hp.harmonic_ud_grade(
            input_maps, nside_out=nside_out, pol=True,
            beam_window_in=beam_1d,
            use_pixel_weights=False, pixwin=False,
        )


def test_beam_window_too_short_raises():
    """beam_window shorter than lmax+1 raises ValueError.

    The beam transfer function must cover all multipoles up to lmax.
    Passing an array that is too short should raise immediately rather
    than producing silently wrong results from out-of-bounds access.
    """
    nside_in = 32
    nside_out = 16

    input_map = np.ones(hp.nside2npix(nside_in))

    with pytest.raises(ValueError, match="beam_window_in.*length"):
        hp.harmonic_ud_grade(
            input_map, nside_out=nside_out,
            beam_window_in=np.ones(5),
            use_pixel_weights=False, pixwin=False,
        )

    with pytest.raises(ValueError, match="beam_window_out.*length"):
        hp.harmonic_ud_grade(
            input_map, nside_out=nside_out,
            beam_window_out=np.ones(5),
            use_pixel_weights=False, pixwin=False,
        )


# ------------------------------------------------------------------
# Tests for reconvolution (nside_out == nside_in)
# ------------------------------------------------------------------


def test_reconvolution_same_nside():
    """harmonic_ud_grade with nside_out==nside_in reconvolves to a new beam.

    A map smoothed with a 5-arcmin beam, when reconvolved to a
    15-arcmin beam at the same NSIDE, should have suppressed high-ell
    power.  This tests the single-step transfer ratio b_out/b_in at
    fixed pixelisation.
    """
    nside = 32
    fwhm_in = np.radians(5.0)
    fwhm_out = np.radians(15.0)

    cl = np.zeros(3 * nside)
    cl[2:] = 1.0
    input_map = hp.synfast(cl, nside, new=True)

    # Smooth with input beam
    beamed_map = hp.smoothing(input_map, fwhm=fwhm_in)

    # Reconvolve: deconvolve input beam, apply wider output beam
    output = hp.harmonic_ud_grade(
        beamed_map, nside_out=nside,
        fwhm_in=fwhm_in, fwhm_out=fwhm_out,
        use_pixel_weights=False, pixwin=False,
    )

    # Output should be smoother (wider beam) than the input beamed map
    lmax = 3 * nside - 1
    cl_input = hp.anafast(beamed_map, lmax=lmax)
    cl_output = hp.anafast(output, lmax=lmax)

    # High-ell power should be suppressed by the wider output beam
    high_ell = slice(lmax // 2, None)
    assert np.mean(cl_output[high_ell]) < 0.5 * np.mean(cl_input[high_ell]), \
        "Reconvolution with wider beam should suppress high-ell power"


def test_reconvolution_with_beam_window():
    """Reconvolution with beam_window arrays matches FWHM-based result.

    Passing beam_window_in and beam_window_out arrays (from gauss_beam)
    should produce the same result as passing the equivalent fwhm_in and
    fwhm_out, verifying consistency between the two beam-specification
    APIs for same-NSIDE reconvolution.
    """
    nside = 32
    fwhm_in = np.radians(5.0)
    fwhm_out = np.radians(15.0)
    lmax = 3 * nside - 1

    cl = np.zeros(3 * nside)
    cl[2:] = 1.0
    input_map = hp.synfast(cl, nside, new=True)

    beamed_map = hp.smoothing(input_map, fwhm=fwhm_in)

    beam_in = hp.gauss_beam(fwhm_in, lmax=lmax)
    beam_out = hp.gauss_beam(fwhm_out, lmax=lmax)

    out_fwhm = hp.harmonic_ud_grade(
        beamed_map, nside_out=nside,
        fwhm_in=fwhm_in, fwhm_out=fwhm_out,
        use_pixel_weights=False, pixwin=False,
    )
    out_beam = hp.harmonic_ud_grade(
        beamed_map, nside_out=nside,
        beam_window_in=beam_in, beam_window_out=beam_out,
        use_pixel_weights=False, pixwin=False,
    )

    np.testing.assert_allclose(out_fwhm, out_beam, rtol=1e-10)


# ------------------------------------------------------------------
# input_type='alm' tests
# ------------------------------------------------------------------


def test_input_type_alm_requires_nside_in():
    """input_type='alm' without nside_in raises ValueError.

    When passing a_lm coefficients directly, nside_in cannot be inferred
    from the input (unlike maps where it comes from Npix), so it must
    be provided explicitly for the pixel-window correction to work.
    """
    alm = hp.synalm(np.ones(100), lmax=99)
    with pytest.raises(ValueError, match="nside_in must be provided"):
        hp.harmonic_ud_grade(alm, nside_out=16, input_type="alm")


def test_input_type_invalid():
    """Invalid input_type raises ValueError with a clear message."""
    m = np.zeros(hp.nside2npix(32))
    with pytest.raises(ValueError, match="input_type must be 'map' or 'alm'"):
        hp.harmonic_ud_grade(m, nside_out=16, input_type="bad")


def test_negative_fwhm_raises():
    """Negative fwhm_in and fwhm_out raise ValueError.

    Beam FWHM must be non-negative (0 means no beam).  Negative values
    are not physically meaningful and must be caught early with a clear
    error message.
    """
    m = np.zeros(hp.nside2npix(32))
    with pytest.raises(ValueError, match="fwhm_in must be >= 0"):
        hp.harmonic_ud_grade(m, nside_out=16, fwhm_in=-1.0,
                             use_pixel_weights=False)
    with pytest.raises(ValueError, match="fwhm_out must be >= 0 or None"):
        hp.harmonic_ud_grade(m, nside_out=16, fwhm_out=-1.0,
                             use_pixel_weights=False)


def test_input_type_alm_matches_map_input():
    """input_type='alm' produces the same output as map input for a single-mode signal.

    For a single-mode a_lm at ell=10 (no iterative SHT errors), the
    two code paths should give identical results: the map path does
    map→map2alm→transfer→alm2map, while the alm path skips map2alm and
    applies the transfer directly.
    """
    nside_in = 64
    nside_out = 16
    ell = 10
    lmax = 3 * nside_out - 1

    # Single-mode alm is exact — no iterative SHT errors
    alm = np.zeros(hp.Alm.getsize(lmax), dtype=np.complex128)
    alm[hp.Alm.getidx(lmax, ell, 0)] = 1.0
    m_in = hp.alm2map(alm, nside=nside_in, lmax=lmax, pixwin=False)

    # Map input path (plain truncation, no pixwin, no beam)
    m_out_map = hp.harmonic_ud_grade(
        m_in, nside_out=nside_out, lmax=lmax,
        pixwin=False, fwhm_out=0, use_pixel_weights=False,
    )

    # Alm input path (skips map2alm)
    m_out_alm = hp.harmonic_ud_grade(
        alm, nside_out=nside_out, nside_in=nside_in, lmax=lmax,
        input_type="alm",
        pixwin=False, fwhm_out=0,
    )

    np.testing.assert_allclose(m_out_map, m_out_alm, atol=1e-12)


def test_input_type_alm_with_pixwin():
    """input_type='alm' with pixwin=True applies pixel-window transfer.

    The pixel-window correction should change the output compared to
    pixwin=False even when the input is a_lm (not a map), because
    the transfer function includes the p_out/p_in ratio.
    """
    nside_in = 64
    nside_out = 16
    lmax = 3 * nside_out - 1

    cl = np.zeros(lmax + 1)
    cl[2:] = 1.0 / np.arange(2, lmax + 1) ** 2
    alm = hp.synalm(cl, lmax=lmax)

    # With pixwin
    m_out_pw = hp.harmonic_ud_grade(
        alm, nside_out=nside_out, nside_in=nside_in,
        input_type="alm",
        pixwin=True, fwhm_out=0,
    )

    # Without pixwin
    m_out_no_pw = hp.harmonic_ud_grade(
        alm, nside_out=nside_out, nside_in=nside_in,
        input_type="alm",
        pixwin=False, fwhm_out=0,
    )

    # Pixel window correction should change the result
    assert not np.allclose(m_out_pw, m_out_no_pw, atol=1e-10)


def test_input_type_alm_polarized():
    """input_type='alm' with polarized 3-component TEB input produces TQU output.

    A 3-component a_lm list [alm_T, alm_E, alm_B] with pol=True and
    pixwin=True should produce a (3, Npix) TQU output where the T and
    E/B components differ because the polarization pixel window differs
    from the temperature one.
    """
    nside_in = 32
    nside_out = 16
    lmax = 3 * nside_out - 1

    alm = np.zeros(hp.Alm.getsize(lmax), dtype=np.complex128)
    alm[hp.Alm.getidx(lmax, 5, 0)] = 1.0
    alm_teb = [alm, alm * 0.5, alm * 0.3]

    output = hp.harmonic_ud_grade(
        alm_teb, nside_out=nside_out, nside_in=nside_in,
        input_type="alm", pol=True,
        pixwin=True, fwhm_out=0,
    )

    assert output.shape == (3, hp.nside2npix(nside_out))
    # T, E, B should differ because polarization pixel window differs from T
    assert not np.allclose(output[0], output[1])


def test_input_type_alm_does_not_mutate_input():
    """input_type='alm' must not modify the caller's input arrays.

    Regression test: an earlier version of _apply_harmonic_transfer
    assigned to alm[0], alm[1], alm[2] in place, mutating the user's
    TEB list through aliasing. Now the function must return a new list
    and leave the inputs untouched.
    """
    nside_in = 32
    nside_out = 16
    lmax = 3 * nside_out - 1

    base = np.zeros(hp.Alm.getsize(lmax), dtype=np.complex128)
    base[hp.Alm.getidx(lmax, 5, 0)] = 1.0
    alm_teb = [base.copy(), (base * 0.5).copy(), (base * 0.3).copy()]
    expected = [a.copy() for a in alm_teb]

    hp.harmonic_ud_grade(
        alm_teb, nside_out=nside_out, nside_in=nside_in,
        input_type="alm", pol=True,
        fwhm_in=np.radians(5.0), fwhm_out=np.radians(10.0),
        pixwin=True,
    )

    for got, want in zip(alm_teb, expected):
        np.testing.assert_array_equal(got, want)

    # Same check for the 1D / spin-0 path
    alm_1d = base.copy()
    snapshot = alm_1d.copy()
    hp.harmonic_ud_grade(
        alm_1d, nside_out=nside_out, nside_in=nside_in,
        input_type="alm", pol=False,
        fwhm_in=np.radians(5.0), fwhm_out=np.radians(10.0),
        pixwin=True,
    )
    np.testing.assert_array_equal(alm_1d, snapshot)


def test_input_type_alm_rejects_unsupported_multi_alm():
    """input_type='alm' rejects 2 or 4+ alm arrays with a clear ValueError.

    A list of 2 (T+E) or 4+ alm arrays has no map-path analogue (the map
    path supports either 1 spin-0 map or a TQU triplet), so the function
    must reject it explicitly rather than silently demoting to spin-0.
    """
    nside_in = 32
    nside_out = 16
    lmax = 3 * nside_out - 1
    base = np.zeros(hp.Alm.getsize(lmax), dtype=np.complex128)

    with pytest.raises(ValueError, match="supports a single 1D alm array or a TEB triplet"):
        hp.harmonic_ud_grade(
            [base, base], nside_out=nside_out, nside_in=nside_in,
            input_type="alm",
        )

    with pytest.raises(ValueError, match="supports a single 1D alm array or a TEB triplet"):
        hp.harmonic_ud_grade(
            [base, base, base, base], nside_out=nside_out, nside_in=nside_in,
            input_type="alm",
        )


def test_input_type_alm_teb_mismatched_sizes_raises():
    """input_type='alm' rejects a TEB triplet whose components have different sizes.

    The function assumes all three TEB alms share the same lmax (true for
    everything healpy itself produces). A hostile caller could pass
    mismatched sizes; the function must catch this early with a clear
    ValueError rather than crashing later with a confusing shape error.
    """
    nside_in = 32
    nside_out = 16
    lmax_a = 3 * nside_out - 1
    lmax_b = 3 * nside_out  # different lmax
    alm_a = np.zeros(hp.Alm.getsize(lmax_a), dtype=np.complex128)
    alm_b = np.zeros(hp.Alm.getsize(lmax_b), dtype=np.complex128)

    with pytest.raises(ValueError, match="TEB.*same size|same.*size"):
        hp.harmonic_ud_grade(
            [alm_a, alm_b, alm_a], nside_out=nside_out, nside_in=nside_in,
            input_type="alm",
        )


def test_input_type_alm_warns_on_ignored_args():
    """input_type='alm' warns when args that only apply to the map path are set.

    The ``iter``, ``use_weights``, ``use_pixel_weights``, and ``datapath``
    arguments are only used by the map2alm step. When ``input_type='alm'``
    skips that step, setting any of these to a non-default value is a
    debugging pitfall and must emit a UserWarning.
    """
    nside_in = 32
    nside_out = 16
    lmax = 3 * nside_out - 1
    alm = np.zeros(hp.Alm.getsize(lmax), dtype=np.complex128)
    alm[hp.Alm.getidx(lmax, 5, 0)] = 1.0

    common = dict(
        nside_out=nside_out, nside_in=nside_in, input_type="alm",
        pixwin=False, fwhm_out=0,
    )

    # iter (default is None for harmonic_ud_grade)
    with pytest.warns(UserWarning, match="iter.*ignored.*input_type='alm'"):
        hp.harmonic_ud_grade(alm, iter=3, **common)

    # use_weights (default False)
    with pytest.warns(UserWarning, match="use_weights.*ignored.*input_type='alm'"):
        hp.harmonic_ud_grade(alm, use_weights=True, **common)

    # use_pixel_weights (default True) — explicit True is fine (matches default)
    # but explicit False should warn since the user is asking for non-default behavior
    with pytest.warns(UserWarning, match="use_pixel_weights.*ignored.*input_type='alm'"):
        hp.harmonic_ud_grade(alm, use_pixel_weights=False, **common)

    # datapath (default None)
    with pytest.warns(UserWarning, match="datapath.*ignored.*input_type='alm'"):
        hp.harmonic_ud_grade(alm, datapath="/tmp", **common)


def test_input_type_alm_no_warning_when_defaults():
    """input_type='alm' is silent when only default values are passed for ignored args.

    Default values for ``iter``, ``use_weights``, ``use_pixel_weights``,
    and ``datapath`` must not trigger the ignored-args warning — otherwise
    every alm-path call would be noisy.
    """
    nside_in = 32
    nside_out = 16
    lmax = 3 * nside_out - 1
    alm = np.zeros(hp.Alm.getsize(lmax), dtype=np.complex128)
    alm[hp.Alm.getidx(lmax, 5, 0)] = 1.0

    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # any warning becomes an error
        hp.harmonic_ud_grade(
            alm, nside_out=nside_out, nside_in=nside_in,
            input_type="alm", pixwin=False, fwhm_out=0,
        )


def test_input_type_alm_with_beam():
    """input_type='alm' correctly applies beam transfer functions.

    The a_lm path (no map2alm) and the map path (with map2alm) should
    agree to within the map2alm round-trip tolerance when both are given
    the same fwhm_in and fwhm_out.  This verifies that beam deconvolution
    and reconvolution work correctly when the forward SHT is skipped.
    """
    nside_in = 64
    nside_out = 32
    lmax = 3 * nside_out - 1

    cl = np.zeros(lmax + 1)
    cl[2:] = 1.0 / np.arange(2, lmax + 1) ** 2
    alm = hp.synalm(cl, lmax=lmax)

    fwhm_in = np.radians(30.0)
    fwhm_out = np.radians(60.0)

    # Map path: map -> map2alm -> transfer -> alm2map
    m_in = hp.alm2map(alm, nside=nside_in, lmax=lmax, pixwin=False)
    m_out_map = hp.harmonic_ud_grade(
        m_in, nside_out=nside_out,
        fwhm_in=fwhm_in, fwhm_out=fwhm_out,
        pixwin=False, use_pixel_weights=False,
    )

    # Alm path: alm -> transfer -> alm2map (no map2alm rounding)
    m_out_alm = hp.harmonic_ud_grade(
        alm, nside_out=nside_out, nside_in=nside_in,
        input_type="alm",
        fwhm_in=fwhm_in, fwhm_out=fwhm_out,
        pixwin=False,
    )

    # Should be close; difference only from map2alm round-trip in map path
    np.testing.assert_allclose(m_out_map, m_out_alm, atol=1e-4)


def test_input_type_alm_same_nside_reconvolution():
    """input_type='alm' with same NSIDE performs beam reconvolution.

    Reconvolving from a 10-arcmin to 30-arcmin beam via a_lm input
    should suppress high-ell power relative to the original unsmoothed
    signal, confirming that the transfer ratio b_out/b_in is applied
    correctly even when the map2alm step is skipped.
    """
    nside = 32
    lmax = 3 * nside - 1

    cl = np.zeros(lmax + 1)
    cl[2:] = 1.0 / np.arange(2, lmax + 1) ** 2
    alm = hp.synalm(cl, lmax=lmax)

    # Use small beams so the division doesn't overflow at high ℓ
    fwhm_in = np.radians(10.0)
    fwhm_out = np.radians(30.0)

    m_out = hp.harmonic_ud_grade(
        alm, nside_out=nside, nside_in=nside,
        input_type="alm",
        fwhm_in=fwhm_in, fwhm_out=fwhm_out,
        pixwin=False,
    )

    # Output should be at the same nside
    assert m_out.shape == (hp.nside2npix(nside),)

    # The output beam should be broader than the input
    # Compare power: smoothing with a wider beam reduces high-ell power
    m_original = hp.alm2map(alm, nside=nside, lmax=lmax, pixwin=False)
    cl_original = hp.anafast(m_original, lmax=lmax)
    cl_out = hp.anafast(m_out, lmax=lmax)
    # At high ell, reconvolved power should be suppressed
    assert cl_out[lmax] < cl_original[lmax]


def test_input_type_alm_truncates_higher_lmax():
    """input_type='alm' with alm lmax > output lmax truncates via resize_alm.

    When the input a_lm have lmax higher than the output bandlimit
    (3*nside_out-1), the function must truncate the a_lm using
    resize_alm before calling alm2map.  This test verifies the
    truncation by comparing against a manual resize_alm + alm2map
    reference.
    """
    nside_in = 256
    nside_out = 64
    lmax_out = 3 * nside_out - 1

    cl = np.zeros(3 * nside_in)
    cl[2:] = 1.0 / np.arange(2, 3 * nside_in) ** 2
    # synalm uses the global numpy RNG; save and restore around the call
    # so this test doesn't perturb state seen by other tests.
    saved_state = np.random.get_state()
    try:
        np.random.seed(42)
        # Create alm at lmax_in = 3*nside_in - 1, much higher than lmax_out
        alm = hp.synalm(cl, lmax=3 * nside_in - 1)
    finally:
        np.random.set_state(saved_state)

    assert hp.Alm.getlmax(len(alm)) > lmax_out  # alm has higher lmax

    # Should work: truncates alm from high lmax to lmax_out
    m_out = hp.harmonic_ud_grade(
        alm, nside_out=nside_out, nside_in=nside_in,
        input_type="alm", pol=False, pixwin=False, fwhm_out=0,
    )

    # Verify by manually truncating and calling alm2map
    alm_truncated = hp.resize_alm(alm, 3 * nside_in - 1, 3 * nside_in - 1, lmax_out, lmax_out)
    m_ref = hp.alm2map(alm_truncated, nside=nside_out, lmax=lmax_out, pixwin=False)

    np.testing.assert_allclose(m_out, m_ref, atol=1e-12)


def test_input_type_alm_uses_lower_input_lmax_by_default():
    """input_type='alm' accepts alm arrays below the default output bandlimit."""
    nside_in = 64
    nside_out = 32
    lmax = 8

    alm = np.zeros(hp.Alm.getsize(lmax), dtype=np.complex128)
    alm[hp.Alm.getidx(lmax, 4, 0)] = 1.0

    m_out = hp.harmonic_ud_grade(
        alm, nside_out=nside_out, nside_in=nside_in,
        input_type="alm", pol=False, pixwin=False, fwhm_out=0,
    )
    m_ref = hp.alm2map(alm, nside=nside_out, lmax=lmax, pixwin=False)

    np.testing.assert_allclose(m_out, m_ref, atol=1e-12)


def test_reconvolution_suppresses_high_ell_power():
    """Reconvolution to a wider beam suppresses high-ell power, not amplifies it.

    The transfer function is applied as a single-step ratio
    fl = b_out / b_in, so going from a narrow to a wide beam should
    reduce power at high ell.  This also verifies that pixel-weight
    vs ring-weight SHT differences are small relative to the signal
    (the residual is not amplified by the transfer).
    """
    nside = 64
    lmax = 3 * nside - 1

    fwhm_narrow = np.radians(10.0 / 60)  # 10 arcmin
    fwhm_wide = np.radians(30.0 / 60)    # 30 arcmin

    cl = np.zeros(lmax + 1)
    cl[2:] = 1.0 / np.arange(2, lmax + 1) ** 2
    saved_state = np.random.get_state()
    try:
        np.random.seed(42)
        alm = hp.synalm(cl, lmax=lmax)
    finally:
        np.random.set_state(saved_state)
    map_narrow = hp.alm2map(alm, nside=nside, fwhm=fwhm_narrow)

    # Reconvolve from narrow to wide at same NSIDE
    map_reconv = hp.harmonic_ud_grade(
        map_narrow, nside_out=nside, pol=False, pixwin=False,
        fwhm_in=fwhm_narrow, fwhm_out=fwhm_wide,
        use_pixel_weights=False,
    )

    # Reference: directly smooth from alm at wide beam
    map_wide_direct = hp.alm2map(alm, nside=nside, fwhm=fwhm_wide)

    # 1. Reconvolved map should have LESS power than input at high ell
    cl_narrow = hp.anafast(map_narrow, lmax=lmax)
    cl_reconv = hp.anafast(map_reconv, lmax=lmax)

    # At high ell (ell > 200), reconvolved power should be suppressed
    high_ell = lmax // 2
    assert cl_reconv[high_ell] < cl_narrow[high_ell], (
        f"Reconvolution should suppress high-ell power: "
        f"cl_reconv[{high_ell}]={cl_reconv[high_ell]:.2e} >= "
        f"cl_narrow[{high_ell}]={cl_narrow[high_ell]:.2e}"
    )

    # 2. Reconvolved spectrum should match direct smoothing closely
    #    (within the map2alm round-trip accuracy)
    cl_wide = hp.anafast(map_wide_direct, lmax=lmax)
    np.testing.assert_allclose(
        cl_reconv[10:high_ell], cl_wide[10:high_ell],
        rtol=0.05,
        err_msg="Reconvolved power spectrum should match direct smoothing",
    )

    # 3. Pixel-level residual should be small relative to signal
    #    The residual comes from the map2alm round-trip (estimating alm
    #    from a pixelized map, then resynthesizing). This is typically
    #    a few percent for a beam-smoothed map at NSIDE=64.
    residual = map_reconv - map_wide_direct
    relative_std = residual.std() / map_wide_direct.std()
    assert relative_std < 0.05, (
        f"Pixel-level residual too large: relative std = {relative_std:.4f}"
    )


# ── Corner-case tests for internal helpers ──────────────────────────


class TestResolveBeam:
    """Tests for the _resolve_beam helper."""

    def test_beam_window_takes_precedence(self):
        """beam_window overrides gauss_bl even when both are provided."""
        from healpy.sphtfunc import _resolve_beam

        lmax = 10
        bw = np.ones((lmax + 1, 1)) * 0.5  # custom beam window
        gb = np.ones(lmax + 1)  # gauss beam (would return 1.0)
        result = _resolve_beam(bw, gb, polar_component=False)
        np.testing.assert_array_equal(result, np.full(lmax + 1, 0.5))

    def test_none_both_returns_none(self):
        """None beam_window and None gauss_bl → None."""
        from healpy.sphtfunc import _resolve_beam

        assert _resolve_beam(None, None, polar_component=False) is None

    def test_1d_gauss_bl_temperature(self):
        """1D gauss_bl with polar_component=False returns the array as-is."""
        from healpy.sphtfunc import _resolve_beam

        gb = np.array([1.0, 0.9, 0.8])
        result = _resolve_beam(None, gb, polar_component=False)
        np.testing.assert_array_equal(result, gb)

    def test_1d_gauss_bl_polarization_raises(self):
        """1D gauss_bl with polar_component=True raises ValueError.

        The polarization column does not exist when gauss_beam was called
        with pol=False, so attempting to extract it is an error.
        """
        from healpy.sphtfunc import _resolve_beam

        gb = np.array([1.0, 0.9, 0.8])
        with pytest.raises(ValueError, match="polar_component=True but gauss_bl is 1D"):
            _resolve_beam(None, gb, polar_component=True)

    def test_2d_gauss_bl_temperature_column(self):
        """2D gauss_bl with polar_component=False returns column 0."""
        from healpy.sphtfunc import _resolve_beam

        gb = np.array([[1.0, 0.5], [0.9, 0.4], [0.8, 0.3]])
        result = _resolve_beam(None, gb, polar_component=False)
        np.testing.assert_array_equal(result, np.array([1.0, 0.9, 0.8]))

    def test_2d_gauss_bl_polarization_column(self):
        """2D gauss_bl with polar_component=True returns column 1."""
        from healpy.sphtfunc import _resolve_beam

        gb = np.array([[1.0, 0.5], [0.9, 0.4], [0.8, 0.3]])
        result = _resolve_beam(None, gb, polar_component=True)
        np.testing.assert_array_equal(result, np.array([0.5, 0.4, 0.3]))


class TestApplyHarmonicTransfer:
    """Tests for the _apply_harmonic_transfer helper."""

    @pytest.fixture
    def lmax(self):
        return 10

    @pytest.fixture
    def fl_T(self, lmax):
        return np.ones(lmax + 1)

    @pytest.fixture
    def fl_P(self, lmax):
        return np.ones(lmax + 1) * 0.5

    @pytest.fixture
    def alm_1d(self, lmax):
        """A single 1D alm array."""
        return np.ones(hp.Alm.getsize(lmax), dtype=np.complex128)

    @pytest.fixture
    def alm_teb_list(self, lmax):
        """A TEB triplet as a list of 3 arrays."""
        base = np.ones(hp.Alm.getsize(lmax), dtype=np.complex128)
        return [base.copy(), base.copy() * 0.5, base.copy() * 0.3]

    @pytest.fixture
    def alm_teb_2d(self, lmax):
        """A TEB triplet as a 2D ndarray with 3 rows."""
        base = np.ones(hp.Alm.getsize(lmax), dtype=np.complex128)
        return np.array([base, base * 0.5, base * 0.3])

    def test_1d_alm_returns_1d(self, alm_1d, fl_T):
        """1D alm input returns a 1D array."""
        from healpy.sphtfunc import _apply_harmonic_transfer

        result = _apply_harmonic_transfer(alm_1d, fl_T)
        assert isinstance(result, np.ndarray)
        assert result.ndim == 1

    def test_1d_alm_does_not_mutate(self, alm_1d, fl_T):
        """1D alm input is not modified in place."""
        from healpy.sphtfunc import _apply_harmonic_transfer

        snapshot = alm_1d.copy()
        fl_T_scaled = np.ones_like(fl_T) * 0.8
        _apply_harmonic_transfer(alm_1d, fl_T_scaled)
        np.testing.assert_array_equal(alm_1d, snapshot)

    def test_teb_list_with_fl_P(self, alm_teb_list, fl_T, fl_P):
        """TEB list with fl_P: T gets fl_T, E/B get fl_P."""
        from healpy.sphtfunc import _apply_harmonic_transfer

        result = _apply_harmonic_transfer(alm_teb_list, fl_T, fl_P)
        assert isinstance(result, list)
        assert len(result) == 3
        # T component: fl_T (all ones) * original
        np.testing.assert_allclose(result[0], alm_teb_list[0])
        # E component: fl_P (0.5) * original
        np.testing.assert_allclose(result[1], alm_teb_list[1] * 0.5)
        # B component: fl_P (0.5) * original
        np.testing.assert_allclose(result[2], alm_teb_list[2] * 0.5)

    def test_teb_list_without_fl_P(self, alm_teb_list, fl_T):
        """TEB list without fl_P: all components get fl_T."""
        from healpy.sphtfunc import _apply_harmonic_transfer

        fl_T_scaled = np.ones_like(fl_T) * 0.7
        result = _apply_harmonic_transfer(alm_teb_list, fl_T_scaled, fl_P=None)
        assert isinstance(result, list)
        assert len(result) == 3
        # All three get the same fl_T
        np.testing.assert_allclose(result[0], alm_teb_list[0] * 0.7)
        np.testing.assert_allclose(result[1], alm_teb_list[1] * 0.7)
        np.testing.assert_allclose(result[2], alm_teb_list[2] * 0.7)

    def test_teb_2d_with_fl_P(self, alm_teb_2d, fl_T, fl_P):
        """2D ndarray TEB with fl_P: T gets fl_T, E/B get fl_P."""
        from healpy.sphtfunc import _apply_harmonic_transfer

        result = _apply_harmonic_transfer(alm_teb_2d, fl_T, fl_P)
        assert isinstance(result, list)
        assert len(result) == 3

    def test_teb_list_does_not_mutate(self, alm_teb_list, fl_T, fl_P):
        """TEB list input arrays are not modified in place."""
        from healpy.sphtfunc import _apply_harmonic_transfer

        snapshots = [a.copy() for a in alm_teb_list]
        _apply_harmonic_transfer(alm_teb_list, fl_T, fl_P)
        for got, want in zip(alm_teb_list, snapshots):
            np.testing.assert_array_equal(got, want)

    def test_teb_2d_does_not_mutate(self, alm_teb_2d, fl_T, fl_P):
        """2D ndarray TEB input is not modified in place."""
        from healpy.sphtfunc import _apply_harmonic_transfer

        snapshot = alm_teb_2d.copy()
        _apply_harmonic_transfer(alm_teb_2d, fl_T, fl_P)
        np.testing.assert_array_equal(alm_teb_2d, snapshot)

    def test_single_element_list_unwraps(self, lmax, fl_T):
        """List of 1 alm array is treated as spin-0 (not TEB)."""
        from healpy.sphtfunc import _apply_harmonic_transfer

        alm = [np.ones(hp.Alm.getsize(lmax), dtype=np.complex128)]
        result = _apply_harmonic_transfer(alm, fl_T)
        # Should return a list with 1 element (spin-0 fallback)
        assert isinstance(result, list)
        assert len(result) == 1


class TestHarmonicUdGradeAlmInputShapes:
    """Tests for harmonic_ud_grade input_type='alm' shape validation."""

    @pytest.fixture
    def nside_in(self):
        return 32

    @pytest.fixture
    def nside_out(self):
        return 16

    @pytest.fixture
    def lmax(self, nside_out):
        return 3 * nside_out - 1

    @pytest.fixture
    def base_alm(self, lmax):
        """A base alm array with a single mode at ell=5."""
        alm = np.zeros(hp.Alm.getsize(lmax), dtype=np.complex128)
        alm[hp.Alm.getidx(lmax, 5, 0)] = 1.0
        return alm

    def test_2d_three_rows_polarized(self, base_alm, nside_in, nside_out, lmax):
        """2D ndarray with 3 rows: polarized input_type='alm' works."""
        alm_2d = np.array([base_alm, base_alm * 0.5, base_alm * 0.3])
        output = hp.harmonic_ud_grade(
            alm_2d, nside_out=nside_out, nside_in=nside_in,
            input_type="alm", pol=True, pixwin=False, fwhm_out=0,
        )
        assert len(output) == 3

    def test_2d_three_rows_unpolarized(self, base_alm, nside_in, nside_out, lmax):
        """2D ndarray with 3 rows, pol=False: treated as unpolarized TEB."""
        alm_2d = np.array([base_alm, base_alm * 0.5, base_alm * 0.3])
        output = hp.harmonic_ud_grade(
            alm_2d, nside_out=nside_out, nside_in=nside_in,
            input_type="alm", pol=False, pixwin=False, fwhm_out=0,
        )
        # pol=False with 3 alm arrays: is_polarized=False, so
        # all three get fl_T (same transfer). Output is still 3 maps.
        assert len(output) == 3

    def test_2d_one_row_spin0(self, base_alm, nside_in, nside_out, lmax):
        """2D ndarray with 1 row: unwrapped to spin-0."""
        alm_2d = base_alm.reshape(1, -1)
        output = hp.harmonic_ud_grade(
            alm_2d, nside_out=nside_out, nside_in=nside_in,
            input_type="alm", pol=False, pixwin=False, fwhm_out=0,
        )
        assert hp.get_nside(output) == nside_out

    def test_2d_two_rows_raises(self, base_alm, nside_in, nside_out, lmax):
        """2D ndarray with 2 rows: rejected (no map-path analogue)."""
        alm_2d = np.array([base_alm, base_alm])
        with pytest.raises(ValueError, match="2D ndarray with 1 or 3 rows"):
            hp.harmonic_ud_grade(
                alm_2d, nside_out=nside_out, nside_in=nside_in,
                input_type="alm",
            )

    def test_2d_four_rows_raises(self, base_alm, nside_in, nside_out, lmax):
        """2D ndarray with 4 rows: rejected."""
        alm_2d = np.array([base_alm, base_alm, base_alm, base_alm])
        with pytest.raises(ValueError, match="2D ndarray with 1 or 3 rows"):
            hp.harmonic_ud_grade(
                alm_2d, nside_out=nside_out, nside_in=nside_in,
                input_type="alm",
            )

    def test_list_one_element_spin0(self, base_alm, nside_in, nside_out, lmax):
        """List of 1 alm array: unwrapped to spin-0."""
        output = hp.harmonic_ud_grade(
            [base_alm], nside_out=nside_out, nside_in=nside_in,
            input_type="alm", pixwin=False, fwhm_out=0,
        )
        assert hp.get_nside(output) == nside_out

    def test_list_three_elements_polarized(self, base_alm, nside_in, nside_out, lmax):
        """List of 3 alm arrays with pol=True: polarized path."""
        output = hp.harmonic_ud_grade(
            [base_alm, base_alm * 0.5, base_alm * 0.3],
            nside_out=nside_out, nside_in=nside_in,
            input_type="alm", pol=True, pixwin=False, fwhm_out=0,
        )
        assert len(output) == 3
