import importlib.util
from pathlib import Path
from unittest import mock

import pytest

setuptools = pytest.importorskip("setuptools")


def _load_setup_module():
    """Load setup.py without triggering a real setuptools.setup call."""
    setup_path = Path(__file__).resolve().parents[1] / "setup.py"
    spec = importlib.util.spec_from_file_location("healpy_setup_test", setup_path)
    module = importlib.util.module_from_spec(spec)
    with mock.patch("setuptools.setup"):
        spec.loader.exec_module(module)
    return module


def test_build_external_clib_missing_library(monkeypatch, tmp_path_factory):
    setup_mod = _load_setup_module()

    # Instantiate the command with writable build directories.
    cmd = setup_mod.build_external_clib(mock.Mock())
    cmd.build_temp = str(tmp_path_factory.mktemp("build_temp"))
    cmd.build_clib = str(tmp_path_factory.mktemp("build_clib"))

    def fake_check_output(cmd_args, env=None):
        cmd_list = list(cmd_args)
        if cmd_list[:2] == ["pkg-config", "--version"]:
            return b"1.8.0"
        raise setup_mod.CalledProcessError(1, cmd_args)

    monkeypatch.setattr(setup_mod, "check_output", fake_check_output)

    with pytest.raises(setup_mod.ExecError) as excinfo:
        cmd.build_library(
            "fake_lib",
            pkg_config_name="fake_lib",
            local_source=None,
        )

    assert str(excinfo.value) == "library 'fake_lib' is not installed"
