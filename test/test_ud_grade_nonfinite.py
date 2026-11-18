import numpy as np
import pytest

import healpy as hp


@pytest.mark.parametrize("invalid", [np.nan, np.inf, -np.inf, hp.UNSEEN])
@pytest.mark.parametrize("pess", [False, True])
@pytest.mark.parametrize("power", [None, -2, 2])
@pytest.mark.parametrize("order", ["RING", "NESTED"])
def test_ud_grade_ignores_nonfinite_children(invalid, pess, power, order):
    signal = np.full(48, 2.0)
    signal[0] = invalid
    signal[4:8] = invalid
    if order == "RING":
        signal = hp.reorder(signal, n2r=True)
    actual = hp.ud_grade(signal, 1, pess=pess, power=power, order_in=order)
    expected = np.full(12, 2.0 * (0.5 ** power if power else 1.0))
    expected[1] = hp.UNSEEN
    if pess:
        expected[0] = hp.UNSEEN
    np.testing.assert_array_equal(actual, expected)
