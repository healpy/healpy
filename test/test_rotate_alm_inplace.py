"""In-place rotate_alm I/O for non-contiguous complex128 (healpy#702)."""

import numpy as np
import healpy as hp


def test_rotate_alm_strided_complex128_inplace():
    psi, theta, phi = 1.0, 2.0, 3.0
    base = hp.map2alm(np.arange(12, dtype=float))
    ref = base.copy()
    hp.rotate_alm(ref, psi, theta, phi)
    buf = np.empty(base.size * 2, dtype=np.complex128)
    alm = buf[::2]
    alm[:] = base
    hp.rotate_alm(alm, psi, theta, phi)
    np.testing.assert_allclose(alm, ref)


def test_rotate_alm_contiguous_complex128_still_inplace():
    psi, theta, phi = 1.0, 2.0, 3.0
    alm = hp.map2alm(np.arange(12, dtype=float))
    ref = alm.copy()
    hp.rotate_alm(ref, psi, theta, phi)
    hp.rotate_alm(alm, psi, theta, phi)
    np.testing.assert_allclose(alm, ref)


def test_rotate_alm_strided_spin_inplace():
    lmax = 32
    nalm = hp.Alm.getsize(lmax)
    base = np.zeros((3, nalm), dtype=np.complex128)
    base[0, 1] = 1
    base[1, 2] = 1
    ref = base.copy()
    hp.rotate_alm(ref, 0.1, 0.2, 0.3)
    buf = np.empty((3, nalm * 2), dtype=np.complex128)
    alm = buf[:, ::2]
    alm[:] = base
    hp.rotate_alm(alm, 0.1, 0.2, 0.3)
    np.testing.assert_allclose(alm, ref)
