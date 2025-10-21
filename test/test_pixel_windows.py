import numpy as np
import pytest
from astropy.io import fits
from pathlib import Path
from urllib.error import URLError

import healpy as hp


def _write_pixel_window_file(path: Path, nside: int):
    """Create a minimal pixel window FITS file and return the expected arrays."""
    lmax = 3 * nside - 1
    ell = np.arange(lmax + 1, dtype=np.float64)
    temperature = 1.0 / (ell + 1.0)
    polarization = 0.5 / (ell + 1.0)

    cols = fits.ColDefs(
        [
            fits.Column(name="TEMPERATURE", format="D", array=temperature),
            fits.Column(name="POLARIZATION", format="D", array=polarization),
        ]
    )
    table_hdu = fits.BinTableHDU.from_columns(cols, name="PIXEL WINDOW")
    fits.HDUList([fits.PrimaryHDU(), table_hdu]).writeto(path, overwrite=True)
    return temperature, polarization


def test_pixwin_local_datapath_root(tmp_path):
    nside = 32
    data_dir = tmp_path / "pixel_window_functions"
    data_dir.mkdir()
    expected_temp, _ = _write_pixel_window_file(
        data_dir / f"pixel_window_n{nside:04d}.fits", nside
    )

    pixwin = hp.pixwin(nside, datapath=tmp_path)

    np.testing.assert_allclose(pixwin, expected_temp)


def test_pixwin_local_nested_directory(tmp_path):
    nside = 16
    data_dir = tmp_path / "pixel_window_functions"
    data_dir.mkdir()
    expected_temp, expected_pol = _write_pixel_window_file(
        data_dir / f"pixel_window_n{nside:04d}.fits", nside
    )

    pixwin_temp, pixwin_pol = hp.pixwin(nside, pol=True, datapath=data_dir)

    np.testing.assert_allclose(pixwin_temp, expected_temp)
    np.testing.assert_allclose(pixwin_pol, expected_pol)


def test_pixwin_local_file_path(tmp_path):
    nside = 64
    file_path = tmp_path / f"pixel_window_n{nside:04d}.fits"
    expected_temp, _ = _write_pixel_window_file(file_path, nside)

    pixwin = hp.pixwin(nside, datapath=file_path)

    np.testing.assert_allclose(pixwin, expected_temp)


def test_pixwin_custom_lmax(tmp_path):
    nside = 8
    data_dir = tmp_path / "pixel_window_functions"
    data_dir.mkdir()
    expected_temp, _ = _write_pixel_window_file(
        data_dir / f"pixel_window_n{nside:04d}.fits", nside
    )
    lmax = 6

    pixwin = hp.pixwin(nside, lmax=lmax, datapath=tmp_path)

    assert len(pixwin) == lmax + 1
    np.testing.assert_allclose(pixwin, expected_temp[: lmax + 1])


def test_pixwin_missing_local_file(tmp_path):
    with pytest.raises(ValueError, match="Pixel window file not found"):
        hp.pixwin(32, datapath=tmp_path)


def test_pixwin_download_via_astropy(monkeypatch, tmp_path):
    nside = 32
    data_dir = tmp_path / "pixel_window_functions"
    data_dir.mkdir()
    expected_temp, _ = _write_pixel_window_file(
        data_dir / f"pixel_window_n{nside:04d}.fits", nside
    )

    def fake_get_pkg_data_filename(filename, package):
        assert filename.endswith(f"pixel_window_n{nside:04d}.fits")
        return str(data_dir / f"pixel_window_n{nside:04d}.fits")

    monkeypatch.setattr(
        "astropy.utils.data.get_pkg_data_filename", fake_get_pkg_data_filename
    )

    pixwin = hp.pixwin(nside)

    np.testing.assert_allclose(pixwin, expected_temp)


def test_pixwin_download_failure(monkeypatch):
    def raise_urlerror(*args, **kwargs):
        raise URLError("Simulated failure")

    monkeypatch.setattr(
        "astropy.utils.data.get_pkg_data_filename",
        raise_urlerror,
    )

    with pytest.raises(URLError, match="nside 32"):
        hp.pixwin(32)


def test_pixwin_invalid_nside(tmp_path):
    with pytest.raises(ValueError):
        hp.pixwin(15, datapath=tmp_path)
