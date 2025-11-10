import shutil
from pathlib import Path

import pytest


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
    repo_root = Path(__file__).resolve().parents[1]
    source_file = (
        repo_root / "cextern" / "healpix" / "data" / f"pixel_window_n{nside:04d}.fits"
    )
    if not source_file.is_file():
        pytest.skip("Vendored pixel window file is not available in this checkout.")
    local_file = datapath / source_file.name
    shutil.copyfile(source_file, local_file)
    pw = hp.pixwin(nside, datapath=tmp_path)
    assert pw is not None
    assert len(pw) == 3 * nside - 1 + 1
