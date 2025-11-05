import numpy as np
import pytest

from healpy import pixelfunc as pf


def test_check_theta_valid_accepts_bounds():
    pf.check_theta_valid(0.0)
    pf.check_theta_valid(np.pi)
    pf.check_theta_valid(np.array([0.0, np.pi / 2]))


@pytest.mark.parametrize("theta", [-1e-6, np.pi + 2e-5])
def test_check_theta_valid_rejects_out_of_range(theta):
    with pytest.raises(ValueError):
        pf.check_theta_valid(theta)


def test_check_theta_valid_rejects_out_of_range_arrays():
    with pytest.raises(ValueError):
        pf.check_theta_valid(np.array([0.0, -1e-6]))
    with pytest.raises(ValueError):
        pf.check_theta_valid(np.array([np.pi, np.pi + 2e-5]))


def test_lonlat_thetaphi_roundtrip():
    lon = np.array([0.0, 45.0, 120.0])
    lat = np.array([-30.0, 0.0, 60.0])

    theta, phi = pf.lonlat2thetaphi(lon, lat)
    lon_back, lat_back = pf.thetaphi2lonlat(theta, phi)

    assert np.allclose(lon_back, lon)
    assert np.allclose(lat_back, lat)


def test_maptype_single_map():
    npix = pf.nside2npix(1)
    single_map = np.arange(npix, dtype=float)

    assert pf.maptype(single_map) == 0


def test_maptype_sequence_of_maps():
    npix = pf.nside2npix(1)
    stack = [np.arange(npix, dtype=float), np.zeros(npix)]

    assert pf.maptype(stack) == len(stack)


def test_maptype_errors_on_mismatch():
    with pytest.raises(TypeError):
        pf.maptype([np.arange(12), np.zeros(8)])


def test_maptype_errors_on_scalar():
    with pytest.raises(TypeError):
        pf.maptype(42)


def test_get_map_size_array_like():
    npix = pf.nside2npix(2)
    map_vector = np.arange(npix, dtype=float)

    assert pf.get_map_size(map_vector) == npix


def test_get_map_size_dict_like_with_nside():
    npix = pf.nside2npix(2)
    explicit = {0: 1.0, npix - 1: 2.0, "nside": 2}

    assert pf.get_map_size(explicit) == npix
