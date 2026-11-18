import numpy as np
import pytest

from astropy.io import fits

import healpy as hp


@pytest.mark.parametrize("column_names", [None, "SIGNAL", "T", ["SIGNAL"]])
@pytest.mark.parametrize("partial", [False, True])
def test_write_map_single_column_name(tmp_path, column_names, partial):
    signal = np.arange(192, dtype=np.float64)
    if partial:
        signal[:4] = hp.UNSEEN
    filename = tmp_path / "map.fits"
    hp.write_map(filename, signal, column_names=column_names, partial=partial)
    expected_name = "TEMPERATURE" if column_names is None else (
        column_names if isinstance(column_names, str) else column_names[0]
    )
    with fits.open(filename) as hdus:
        assert hdus[1].columns.names == (["PIXEL"] if partial else []) + [expected_name]
        np.testing.assert_array_equal(hdus[1].data[expected_name], signal[4:] if partial else signal)
    np.testing.assert_array_equal(hp.read_map(filename), signal)


def test_write_map_string_name_does_not_split_into_multiple_columns(tmp_path):
    signal = np.arange(192, dtype=np.float64)
    with pytest.raises(AssertionError, match="Length column_names"):
        hp.write_map(tmp_path / "map.fits", [signal, signal], column_names="QU")
