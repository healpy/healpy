"""Tests for pixel window function loading and caching.

This module tests the pixwin() function's ability to load pixel window
functions from local files and use astropy's caching mechanism for
remote downloads.
"""
from pathlib import Path

import pytest
from astropy.utils import data as astropy_data

TEST_DATA = Path(__file__).resolve().parent / "data"
PIXWIN_NSIDE = 4
PIXWIN_FIXTURE = TEST_DATA / f"pixel_window_n{PIXWIN_NSIDE:04d}.fits"


def test_pixwin_download(monkeypatch):
    """Test that pixwin downloads and caches the file using astropy if not present locally."""
    import healpy as hp

    nside = PIXWIN_NSIDE

    assert PIXWIN_FIXTURE.exists(), "bundled pixel window file missing"

    # Mock astropy's download mechanism to use our bundled fixture
    monkeypatch.setattr(
        astropy_data,
        "get_pkg_data_filename",
        lambda *_, **__: str(PIXWIN_FIXTURE),
    )
    pw = hp.pixwin(nside)
    assert pw is not None
    assert len(pw) == 3 * nside


def test_pixwin_download_polarization(monkeypatch):
    """Test that pixwin returns both temperature and polarization windows when pol=True."""
    import healpy as hp

    nside = PIXWIN_NSIDE

    assert PIXWIN_FIXTURE.exists(), "bundled pixel window file missing"

    # Mock astropy's download mechanism to use our bundled fixture
    monkeypatch.setattr(
        astropy_data,
        "get_pkg_data_filename",
        lambda *_, **__: str(PIXWIN_FIXTURE),
    )
    pw_temp, pw_pol = hp.pixwin(nside, pol=True)
    assert pw_temp is not None
    assert pw_pol is not None
    assert len(pw_temp) == 3 * nside
    assert len(pw_pol) == 3 * nside


def test_pixwin_local_datapath(tmp_path):
    """Test that pixwin loads the file from a local datapath if provided."""
    import healpy as hp

    nside = PIXWIN_NSIDE
    datapath = tmp_path / "pixel_window_functions"
    datapath.mkdir(parents=True)
    assert PIXWIN_FIXTURE.exists(), "bundled pixel window file missing"
    local_file = datapath / f"pixel_window_n{nside:04d}.fits"
    local_file.write_bytes(PIXWIN_FIXTURE.read_bytes())
    pw = hp.pixwin(nside, datapath=tmp_path)
    assert pw is not None
    assert len(pw) == 3 * nside


def test_pixwin_local_datapath_missing_file(tmp_path):
    """Test that pixwin raises an error when the file is not found in the specified datapath."""
    import healpy as hp

    nside = PIXWIN_NSIDE
    # Create an empty datapath without the required file
    datapath = tmp_path / "empty_datapath"
    datapath.mkdir(parents=True)
    
    with pytest.raises(ValueError, match="Pixel window file not found"):
        hp.pixwin(nside, datapath=datapath)
