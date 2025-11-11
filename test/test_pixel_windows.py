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

    # Remove any local file by monkeypatching os.path.isfile to always return False
    monkeypatch.setattr("os.path.isfile", lambda path: False)
    monkeypatch.setattr(
        astropy_data,
        "get_pkg_data_filename",
        lambda *_, **__: str(PIXWIN_FIXTURE),
    )
    pw = hp.pixwin(nside)
    assert pw is not None
    assert len(pw) == 3 * nside - 1 + 1


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
    assert len(pw) == 3 * nside - 1 + 1
