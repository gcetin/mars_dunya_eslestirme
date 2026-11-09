"""Terrain maths checked against cases with known analytic answers."""

import math

import numpy as np
import pytest

from analog_finder.bodies import EARTH, MARS, Body
from analog_finder.dem import Patch
from analog_finder.terrain import (
    aeqd_to_latlon,
    extract_features,
    radial_profile,
    resample_to_scale,
    scale_features,
    slope_deg,
    window_features,
)


def synthetic_patch(body: Body, lat0: float, lon0: float, half_m: float, pixel_m: float, func) -> Patch:
    """Sample ``func(x_east_m, y_north_m)`` on a lat/lon grid around a point."""
    dlat = math.degrees(pixel_m / body.radius_m)
    dlon = dlat / math.cos(math.radians(lat0))
    n_lat = int(half_m * 1.3 / pixel_m)
    n_lon = n_lat
    lat = lat0 + (n_lat - np.arange(2 * n_lat + 1)) * dlat
    lon = lon0 + (np.arange(2 * n_lon + 1) - n_lon) * dlon
    lon2, lat2 = np.meshgrid(lon, lat)
    x = np.radians(lon2 - lon0) * body.radius_m * math.cos(math.radians(lat0))
    y = np.radians(lat2 - lat0) * body.radius_m
    return Patch(func(x, y), lat, lon, body, "synthetic")


def test_slope_of_plane_is_exact():
    y, x = np.mgrid[:50, :50] * 100.0
    z = x * math.tan(math.radians(10))
    assert slope_deg(z, 100.0) == pytest.approx(10.0)


def test_scale_features_on_plane():
    y, x = np.mgrid[:50, :50] * 1000.0
    features = scale_features(x * math.tan(math.radians(5)), 1000.0)
    assert features["slope_mean"] == pytest.approx(5.0)
    assert features["slope_std"] == pytest.approx(0.0, abs=1e-9)
    assert features["rough"] == pytest.approx(0.0, abs=1e-6)


def test_aeqd_centre_and_north():
    lat, lon = aeqd_to_latlon(np.array(0.0), np.array(0.0), 18.0, 77.0, MARS.radius_m)
    assert lat == pytest.approx(18.0) and lon == pytest.approx(77.0)
    y = math.radians(1.0) * MARS.radius_m
    lat, lon = aeqd_to_latlon(np.array(0.0), np.array(y), 18.0, 77.0, MARS.radius_m)
    assert lat == pytest.approx(19.0) and lon == pytest.approx(77.0)


@pytest.mark.parametrize("body", [MARS, EARTH])
@pytest.mark.parametrize("lat0", [0.0, 18.4, 60.0])
def test_slope_is_recovered_in_metres_on_each_body(body, lat0):
    """A planar 3° ramp must read 3° regardless of radius or latitude."""
    tan = math.tan(math.radians(3))
    patch = synthetic_patch(body, lat0, 30.0, 30_000, 250.0, lambda x, y: x * tan)
    grid = resample_to_scale(patch, lat0, 30.0, 60_000, 1000.0)
    assert np.isfinite(grid).all()
    assert np.median(slope_deg(grid, 1000.0)) == pytest.approx(3.0, abs=0.05)


def test_fine_and_coarse_sources_converge_after_support_matching():
    """Rough terrain seen at 30 m and at 463 m must agree at the 2 km scale."""

    def terrain(x, y):
        return (
            300 * np.sin(2 * np.pi * x / 25_000) * np.cos(2 * np.pi * y / 18_000)
            + 90 * np.sin(2 * np.pi * (x + y) / 7_000)
        )

    fine = synthetic_patch(MARS, 18.0, 77.0, 30_000, 60.0, terrain)
    coarse = synthetic_patch(MARS, 18.0, 77.0, 30_000, 463.0, terrain)
    a = extract_features(fine, 18.0, 77.0, 60_000, (2_000.0,)).values["slope_mean_2km"]
    b = extract_features(coarse, 18.0, 77.0, 60_000, (2_000.0,)).values["slope_mean_2km"]
    assert a == pytest.approx(b, rel=0.03)


def test_scale_finer_than_source_is_refused():
    patch = synthetic_patch(MARS, 0.0, 0.0, 30_000, 463.0, lambda x, y: x * 0.0)
    with pytest.raises(ValueError, match="too fine"):
        resample_to_scale(patch, 0.0, 0.0, 60_000, 500.0)


def _bowl(depth: float, radius_m: float, spacing: float, n: int = 61) -> np.ndarray:
    axis = (np.arange(n) - (n - 1) / 2) * spacing
    x, y = np.meshgrid(axis, axis)
    r = np.hypot(x, y)
    return -depth * np.clip(1 - (r / radius_m) ** 2, 0, None)


def test_basin_features_identify_a_crater():
    z = _bowl(500.0, 22_000.0, 1000.0)
    features = window_features(z, 1000.0)
    assert features["basin_depth"] > 250
    assert features["radial_r2"] > 0.98
    assert window_features(-z, 1000.0)["basin_depth"] < -250


def test_radial_r2_is_low_for_a_linear_ridge():
    n = 61
    axis = (np.arange(n) - (n - 1) / 2) * 1000.0
    x, _ = np.meshgrid(axis, axis)
    ridge = 400 * np.exp(-(x / 6000.0) ** 2)
    assert window_features(ridge, 1000.0)["radial_r2"] < 0.3


def test_regional_tilt_does_not_change_shape_features():
    z = _bowl(500.0, 22_000.0, 1000.0)
    y, x = np.mgrid[:61, :61] * 1000.0
    tilted = z + 0.02 * x + 0.01 * y
    _, profile_a, r2_a = radial_profile(z, 1000.0)
    _, profile_b, r2_b = radial_profile(tilted, 1000.0)
    assert r2_a == pytest.approx(r2_b, abs=1e-9)
    assert profile_a == pytest.approx(profile_b, abs=1e-6)


def test_absolute_elevation_is_ignored():
    z = _bowl(500.0, 22_000.0, 1000.0)
    a, b = window_features(z, 1000.0), window_features(z - 3000.0, 1000.0)
    assert a == pytest.approx(b)


def test_nan_gap_is_tolerated_but_large_holes_invalidate():
    func = lambda x, y: 0.03 * x
    patch = synthetic_patch(MARS, 10.0, 40.0, 30_000, 463.0, func)
    patch.z[100:104, 100:104] = np.nan
    assert np.isfinite(extract_features(patch, 10.0, 40.0).values["slope_mean_1km"])
    patch.z[: patch.z.shape[0] // 2] = np.nan
    result = extract_features(patch, 10.0, 40.0)
    assert result.valid_frac < 0.95
    assert math.isnan(result.values["slope_mean_1km"])
