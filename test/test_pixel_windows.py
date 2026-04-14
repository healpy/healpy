from pathlib import Path

import pytest
from astropy.utils import data as astro_data
from urllib.error import URLError


def test_pixwin_download(monkeypatch):
    """Test that pixwin downloads and caches the file using astropy if not present locally."""
    import healpy as hp

    nside = 32
    # Remove any local file by monkeypatching os.path.isfile to always return False
    monkeypatch.setattr("os.path.isfile", lambda path: False)
    pw = hp.pixwin(nside)
    assert pw is not None
    assert len(pw) == 3 * nside - 1 + 1


def test_pixwin_local_datapath(tmp_path):
    """Test that pixwin loads the file from a local datapath if provided."""
    import healpy as hp

    nside = 32
    datapath = tmp_path / "pixel_window_functions"
    datapath.mkdir(parents=True)
    filename = f"pixel_window_functions/pixel_window_n{nside:04d}.fits"
    try:
        with astro_data.conf.set_temp("remote_timeout", 30):
            cached_path = astro_data.get_pkg_data_filename(
                filename, package="healpy"
            )
    except (OSError, URLError) as exc:
        pytest.skip(f"pixel window download unavailable: {exc}")
    local_file = datapath / f"pixel_window_n{nside:04d}.fits"
    local_file.write_bytes(Path(cached_path).read_bytes())
    pw = hp.pixwin(nside, datapath=tmp_path)
    assert pw is not None
    assert len(pw) == 3 * nside - 1 + 1
