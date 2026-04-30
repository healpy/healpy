# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.1
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# # `harmonic_ud_grade` vs `ud_grade`: A Focused Comparison
#
# This notebook demonstrates why `healpy.harmonic_ud_grade` is the
# preferred method for changing HEALPix resolution, compared to the
# pixel-space `ud_grade`.
#
# The key improvement is **anti-aliasing**: `ud_grade` simply averages
# sub-pixels, so high-multipole power from the input grid folds into
# the output. `harmonic_ud_grade` works in spherical-harmonic space,
# applying the Planck 2015 XVI Eq. 1 prescription that corrects for
# pixel windows and optionally downgrades the beam.

# %%
import time
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt

plt.rcParams.update({"figure.dpi": 120, "font.size": 11})

# Planck 2013 XXIII Table 1: exact FWHM-to-pixel ratio
PLANCK_K = 160.0 / (np.degrees(hp.nside2resol(64)) * 60)
print(f"healpy {hp.__version__}")
print(f"Planck FWHM/pixel ratio: {PLANCK_K:.4f}")

# %% [markdown]
# ## 1. Aliasing Stress Test: A Single High-$ℓ$ Mode
#
# We create a map at `nside_in = 128` containing **only** one mode at
# $ℓ = 120$ — well above the output bandlimit of
# `3 * 32 - 1 = 95`.  Any non-zero output is aliasing.

# %%
nside_in = 128
nside_out = 32
lmax_out = 3 * nside_out - 1

# Single mode at ell=120
alm_single = np.zeros(hp.Alm.getsize(120), dtype=np.complex128)
alm_single[hp.Alm.getidx(120, 120, 0)] = 1.0
m_in = hp.alm2map(alm_single, nside=nside_in, pixwin=False)

# Downgrade methods
m_ud = hp.ud_grade(m_in, nside_out=nside_out)

# harmonic_ud_grade with new API: fwhm_in is now required
# For this synthetic test there is no input beam → pass 0
m_harm = hp.harmonic_ud_grade(
    m_in,
    nside_out=nside_out,
    fwhm_in=0,
    use_pixel_weights=False,
    pixwin=False,
    fwhm_out=0,
)

# %%
fig, axes = plt.subplots(1, 3, figsize=(14, 3.5))
for ax, m, title in zip(
    axes,
    [m_in, m_ud, m_harm],
    [f"Input (Nside={nside_in})", f"ud_grade (Nside={nside_out})", f"harmonic_ud_grade (Nside={nside_out})"],
):
    hp.mollview(m, title=title, hold=True, cbar=True, cmap="RdBu_r")
plt.tight_layout()
plt.show()

print(f"Input map RMS:              {np.std(m_in):.6e}")
print(f"ud_grade output RMS:        {np.std(m_ud):.6e}  <-- aliasing")
print(f"harmonic_ud_grade RMS:      {np.std(m_harm):.6e}  <-- suppressed")
print(f"Ratio ud/harm:              {np.std(m_ud) / np.std(m_harm):.1f}x")

# %% [markdown]
# ## 2. Power Spectrum Recovery — Why It Matters
#
# Real maps have broad power spectra. `ud_grade` folds high-$\ell$
# power into all multipoles, corrupting the spectrum. Here we
# synthesise a smooth $C_\ell \propto \ell^{-2}$ map at
# `Nside_in = 256`, downgrade to `Nside_out = 64`, and recover the
# spectrum.

# %%
nside_in = 256
nside_out = 64
lmax_in = 3 * nside_in - 1
lmax_out = 3 * nside_out - 1

np.random.seed(42)
ell = np.arange(lmax_in + 1, dtype=float)
cl_in = np.zeros(lmax_in + 1)
cl_in[2:] = ell[2:] ** (-2.0)
elm = hp.synalm(cl_in, lmax=lmax_in)
m_in = hp.alm2map(elm, nside=nside_in, lmax=lmax_in, pixwin=False)

# Ground-truth reference: directly synthesise at output resolution
elm_ref = hp.resize_alm(elm, lmax_in, lmax_in, lmax_out, lmax_out)
m_ref = hp.alm2map(elm_ref, nside=nside_out, lmax=lmax_out, pixwin=False)

m_ud = hp.ud_grade(m_in, nside_out=nside_out)

# New API: fwhm_in=0 for this synthetic simulation with no beam
m_harm = hp.harmonic_ud_grade(
    m_in,
    nside_out=nside_out,
    fwhm_in=0,
    use_pixel_weights=False,
    pixwin=False,
    fwhm_out=0,
)

# Measure spectra
cl_ref = hp.anafast(m_ref, lmax=lmax_out)
cl_ud = hp.anafast(m_ud, lmax=lmax_out)
cl_harm = hp.anafast(m_harm, lmax=lmax_out)
ell_out = np.arange(lmax_out + 1)

# %%
fig, axes = plt.subplots(1, 2, figsize=(12, 4))

axes[0].loglog(ell_out[2:], cl_in[:lmax_out + 1][2:], "k--", lw=2, label="Input $C_\ell$ (truth)")
axes[0].loglog(ell_out[2:], cl_ud[2:], alpha=0.8, label="ud_grade")
axes[0].loglog(ell_out[2:], cl_harm[2:], alpha=0.8, label="harmonic_ud_grade")
axes[0].axvline(lmax_out, color="red", ls=":", label=r"$\ell_{\max}$")
axes[0].set_xlabel(r"$\ell$")
axes[0].set_ylabel(r"$C_\ell$")
axes[0].set_title("Power spectra after downgrade")
axes[0].legend()
axes[0].grid(alpha=0.3)

# Fractional error
frac_ud = (cl_ud[2:] - cl_ref[2:]) / cl_ref[2:]
frac_harm = (cl_harm[2:] - cl_ref[2:]) / cl_ref[2:]

axes[1].plot(ell_out[2:], frac_ud * 100, alpha=0.8, label="ud_grade")
axes[1].plot(ell_out[2:], frac_harm * 100, alpha=0.8, label="harmonic_ud_grade")
axes[1].axhline(0, color="k", lw=0.5)
axes[1].axvline(lmax_out, color="red", ls=":")
axes[1].set_xlabel(r"$\ell$")
axes[1].set_ylabel("Fractional error [%]")
axes[1].set_title("Spectral error vs truth")
axes[1].legend()
axes[1].grid(alpha=0.3)

plt.tight_layout()
plt.show()

# RMS fractional error
safe = slice(2, 2 * nside_out + 1)
print(
    f"RMS frac. error over safe band ℓ ∈ [2, {2*nside_out}]:\n"
    f"  ud_grade:          {np.sqrt(np.mean(frac_ud[safe]**2))*100:.2f}%\n"
    f"  harmonic_ud_grade: {np.sqrt(np.mean(frac_harm[safe]**2))*100:.2f}%"
)

# %% [markdown]
# ## 3. Noise Aliasing — The Real-World Failure
#
# So far the input was smooth and band-limited. Real data contain
# **white noise** at the pixel level. `ud_grade` has no anti-aliasing
# filter, so this noise leaks into every multipole of the output.
# Here we add flat-spectrum noise to the previous test and show how
# each method handles it.

# %%
nside_in = 256
nside_out = 64
lmax_out = 3 * nside_out - 1

np.random.seed(42)
ell_in = np.arange(3 * nside_in + 1, dtype=float)
cl_signal = np.zeros(3 * nside_in + 1)
cl_signal[2:] = ell_in[2:] ** (-2.0)
map_signal = hp.synfast(cl_signal, nside_in, lmax=3 * nside_in, new=True)

# Flat noise: amplitude equal to signal at ℓ = 50 → noise dominates above
ell_knee = 50
cl_noise = np.zeros(3 * nside_in + 1)
cl_noise[2:] = cl_signal[ell_knee]
map_noise = hp.synfast(cl_noise, nside_in, lmax=3 * nside_in, new=True)
map_noisy = map_signal + map_noise

# Reference: downgrade signal-only with harmonic
m_ref_signal = hp.harmonic_ud_grade(
    map_signal,
    nside_out=nside_out,
    fwhm_in=0,
    use_pixel_weights=False,
    pixwin=False,
    fwhm_out=0,
)

# Downgrade methods
m_ud = hp.ud_grade(map_noisy, nside_out=nside_out)
m_harm = hp.harmonic_ud_grade(
    map_noisy,
    nside_out=nside_out,
    fwhm_in=0,
    use_pixel_weights=False,
    pixwin=False,
    fwhm_out=0,
)

# Isolate noise contribution
cl_noise_ud = hp.anafast(m_ud - hp.ud_grade(map_signal, nside_out=nside_out), lmax=lmax_out)
cl_noise_harm = hp.anafast(m_harm - m_ref_signal, lmax=lmax_out)

# %%
fig, axes = plt.subplots(1, 2, figsize=(12, 4))

axes[0].loglog(ell_out[2:], cl_signal[:lmax_out + 1][2:], "k--", lw=2, label="Signal (truth)")
axes[0].loglog(ell_out[2:], cl_noise_ud[2:], alpha=0.8, label="Noise in ud_grade")
axes[0].loglog(ell_out[2:], cl_noise_harm[2:], alpha=0.8, label="Noise in harmonic_ud_grade")
axes[0].set_xlabel(r"$\ell$")
axes[0].set_ylabel(r"$C_\ell$ (noise)")
axes[0].set_title("Noise contamination after downgrade")
axes[0].legend()
axes[0].grid(alpha=0.3)

# Noise-to-signal ratio
nsr_ud = cl_noise_ud[2:] / (cl_signal[:lmax_out + 1][2:] + 1e-30)
nsr_harm = cl_noise_harm[2:] / (cl_signal[:lmax_out + 1][2:] + 1e-30)

axes[1].semilogy(ell_out[2:], nsr_ud, alpha=0.8, label="ud_grade")
axes[1].semilogy(ell_out[2:], nsr_harm, alpha=0.8, label="harmonic_ud_grade")
axes[1].axhline(1.0, color="k", ls=":", label="noise = signal")
axes[1].axvline(lmax_out, color="red", ls=":", alpha=0.5)
axes[1].set_xlabel(r"$\ell$")
axes[1].set_ylabel("Noise / Signal")
axes[1].set_title("Noise-to-signal ratio")
axes[1].legend()
axes[1].grid(alpha=0.3)

plt.tight_layout()
plt.show()

frac_corrupted_ud = np.mean(nsr_ud[lmax_out // 2:] > 1.0) * 100
frac_corrupted_harm = np.mean(nsr_harm[lmax_out // 2:] > 1.0) * 100
print(
    f"Fraction of high-ℓ bins where noise exceeds signal:\n"
    f"  ud_grade:          {frac_corrupted_ud:.0f}%\n"
    f"  harmonic_ud_grade: {frac_corrupted_harm:.0f}%"
)

# %% [markdown]
# ## 4. The New API: Required `fwhm_in` and `input_type="alm"`
#
# `harmonic_ud_grade` now requires `fwhm_in` so the resolution scaling
# is always explicit.  The default `fwhm_out` preserves the Planck
# FWHM-to-pixel ratio across resolutions.
#
# ### 4.1 Default beam scaling

# %%
nside_in = 2048
nside_out = 64

# Planck beam at Nside=2048 is 5 arcmin
fwhm_in = np.radians(5 / 60)

m_input = np.zeros(hp.nside2npix(nside_in))

# With the new API, fwhm_out is auto-computed to preserve the ratio
# fwhm_out = fwhm_in * (resol_out / resol_in)
expected_fwhm_out = fwhm_in * (
    hp.nside2resol(nside_out) / hp.nside2resol(nside_in)
)

# This should equal 160 arcmin (Planck Nside=64 beam)
print(f"Nside_in = {nside_in}, Nside_out = {nside_out}")
print(f"  fwhm_in  = {np.degrees(fwhm_in)*60:.1f} arcmin")
print(f"  expected fwhm_out = {np.degrees(expected_fwhm_out)*60:.1f} arcmin ← matches Planck Nside=64")

# Verify by running the function (it will work even on zeros)
m_out = hp.harmonic_ud_grade(
    m_input,
    nside_out=nside_out,
    fwhm_in=fwhm_in,
    use_pixel_weights=False,
    pixwin=False,
)
print(f"  ✓ harmonic_ud_grade accepted the call with default fwhm_out")

# %% [markdown]
# ### 4.2 `input_type="alm"` — skip map2alm
#
# If you already have $a_{\ell m}$ coefficients, you can feed them
# directly.

# %%
nside = 64
lmax = 3 * nside - 1

np.random.seed(7)
cl = np.zeros(3 * nside)
cl[2:] = 1.0
m = hp.synfast(cl, nside, new=True)

# Get alm directly
alm = hp.map2alm(m, lmax=lmax, use_pixel_weights=False)

# Pass alm instead of map
m_from_alm = hp.harmonic_ud_grade(
    alm,
    nside_out=nside,
    fwhm_in=0,
    input_type="alm",
    nside_in=nside,
    use_pixel_weights=False,
    pixwin=False,
    fwhm_out=0,
)

# Should match going through the map path
m_from_map = hp.harmonic_ud_grade(
    m,
    nside_out=nside,
    fwhm_in=0,
    use_pixel_weights=False,
    pixwin=False,
    fwhm_out=0,
)

np.testing.assert_allclose(m_from_alm, m_from_map, rtol=1e-5)
print("✓ alm input produced identical result to map input")

# %% [markdown]
# ## 5. Performance

# %%
nside_in = 512
nside_out = 128
np.random.seed(42)
m = hp.synfast(np.ones(3 * nside_in), nside_in, new=True)

# ud_grade
t0 = time.perf_counter()
for _ in range(10):
    _ = hp.ud_grade(m, nside_out=nside_out)
t_ud = (time.perf_counter() - t0) / 10

# harmonic_ud_grade
t0 = time.perf_counter()
for _ in range(10):
    _ = hp.harmonic_ud_grade(
        m,
        nside_out=nside_out,
        fwhm_in=0,
        use_pixel_weights=False,
        pixwin=False,
        fwhm_out=0,
    )
t_harm = (time.perf_counter() - t0) / 10

print(f"ud_grade:           {t_ud*1000:.1f} ms")
print(f"harmonic_ud_grade:  {t_harm*1000:.1f} ms")
print(f"slowdown:           {t_harm/t_ud:.1f}x")

# %% [markdown]
# ## Summary
#
# | Feature | `ud_grade` | `harmonic_ud_grade` |
# |--------|------------|---------------------|
# | Aliasing | Heavy leakage | Suppressed to noise floor |
# | Spectrum fidelity | Corrupted at all $\ell$ | Correct within band |
# | Pixel window | Ignored | Corrected (Eq. 1) |
# | Beam scaling | Not handled | Preserves FWHM/pixel ratio |
# | Alm input | Not supported | Supported |
# | Speed | Fast | ~5-15x slower (SHT cost) |
# |
# **Recommendation**: Use `harmonic_ud_grade` whenever you care about
# the power spectrum, map fidelity, or when downgrading maps that
# contain noise or unresolved sub-pixel structure.  Use `ud_grade` only
# when speed is critical and scientific accuracy is not required.
