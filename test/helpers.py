"""Shared helpers for healpy tests."""

import sys
from pathlib import Path

import pytest

DATAPATH = Path(__file__).parent / "data"

# Pyodide has no network access from Python (sockets are unavailable), so
# tests that download pixel windows or pixel weights from healpy-data at
# run time cannot pass there.
requires_network = pytest.mark.skipif(
    sys.platform == "emscripten",
    reason="no network access under Pyodide/Emscripten",
)
