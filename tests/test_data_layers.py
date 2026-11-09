"""Tile stitching, land-cover sampling and MOLA coverage bookkeeping."""

import numpy as np
import pytest

from analog_finder.bodies import MARS, MOON
from analog_finder.dem import Patch, read_earth_patch, read_global_patch, url_is_absent
from analog_finder.landcover import window_fractions
from analog_finder.mola_support import _tile_for
from analog_finder.validation import resample_mask, with_track_gaps


def elevation(lat, lon):
    return 1000.0 * lat + lon


def fake_tile(lat_floor: int, lon_floor: int, n: int = 240) -> np.ndarray:
    lats = lat_floor + 1 - np.arange(n) / n
    lons = lon_floor + np.arange(n) / n
    lon2, lat2 = np.meshgrid(lons, lats)
    return elevation(lat2, lon2).astype("float32")


def test_earth_tiles_join_without_offset_or_flip():
    # Centred on a corner shared by four tiles.
    patch = read_earth_patch(38.0, 30.0, 40_000, arcsec=15, tile_loader=fake_tile)
    lon2, lat2 = np.meshgrid(patch.lon, patch.lat)
    assert patch.z == pytest.approx(elevation(lat2, lon2), abs=0.02)
    assert patch.lat[0] > patch.lat[-1] and patch.lon[0] < patch.lon[-1]
    assert patch.lat.min() < 38.0 < patch.lat.max() and patch.lon.min() < 30.0 < patch.lon.max()


def test_missing_sea_tile_becomes_nan():
    loader = lambda la, lo: None if lo == 30 else fake_tile(la, lo)
    patch = read_earth_patch(38.5, 30.0, 40_000, arcsec=15, tile_loader=loader)
    assert np.isnan(patch.z[:, patch.lon >= 30.0]).all()
    assert np.isfinite(patch.z[:, patch.lon < 30.0]).all()


def test_antimeridian_is_refused():
    with pytest.raises(ValueError, match="antimeridian"):
        read_earth_patch(10.0, 179.9, 40_000, arcsec=15, tile_loader=fake_tile)


def test_polar_window_needs_a_polar_product():
    with pytest.raises(ValueError, match="polar"):
        read_global_patch(MARS, 85.0, 0.0, 40_000)


def test_landcover_fractions_on_a_known_tile():
    tile = np.full((600, 600), 10, "uint8")  # tree
    tile[:300] = 60  # the northern half is bare
    fractions = window_fractions(37.5, 28.5, 60_000, loader=lambda la, lo: tile if (la, lo) == (36, 27) else None)
    assert fractions["lc_bare"] == pytest.approx(0.5, abs=0.02)
    assert fractions["lc_tree"] == pytest.approx(0.5, abs=0.02)
    assert fractions["lc_valid"] == pytest.approx(1.0)


def test_landcover_with_no_tiles_is_invalid():
    assert window_fractions(0.0, 0.0, 60_000, loader=lambda la, lo: None)["lc_valid"] == 0.0


@pytest.mark.parametrize(
    "lat,lon,name",
    [
        (18.38, 77.58, "megc44n000hb.img"),
        (-5.37, 137.81, "megc00n090hb.img"),
        (-50.0, 200.0, "megc44s180hb.img"),
        (60.0, 10.0, "megc88n000hb.img"),
    ],
)
def test_mola_count_tile_names(lat, lon, name):
    assert _tile_for(lat, lon)[0] == name


def test_track_gaps_rebuild_a_smooth_surface():
    lat = 10 - np.arange(200) * 0.004
    lon = 40 + np.arange(200) * 0.004
    lon2, lat2 = np.meshgrid(lon, lat)
    patch = Patch(500 * lat2 + 300 * lon2, lat, lon, MARS, "plane")
    rng = np.random.default_rng(0)
    gappy = with_track_gaps(patch, rng.random((150, 170)) < 0.33)
    # Linear interpolation is exact on a plane inside the measured hull. The
    # outermost cells may lie beyond it and take the nearest measured value.
    assert gappy.z[6:-6, 6:-6] == pytest.approx(patch.z[6:-6, 6:-6], abs=1e-6)
    assert np.abs(gappy.z - patch.z).max() < 15.0


def test_resample_mask_shape_and_share():
    rng = np.random.default_rng(1)
    mask = rng.random((157, 165)) < 0.33
    resized = resample_mask(mask, (201, 159))
    assert resized.shape == (201, 159)
    assert resized.mean() == pytest.approx(mask.mean(), abs=0.03)


@pytest.mark.network
def test_lola_scale_factor_is_applied():
    """Mare Tranquillitatis lies about 1.5-2 km below datum, not 3-4 km."""
    patch = read_global_patch(MOON, 0.67, 23.47, 15_000)
    assert -2300 < float(np.nanmedian(patch.z)) < -1200


@pytest.mark.network
def test_mola_floor_of_jezero():
    patch = read_global_patch(MARS, 18.38, 77.58, 5_000)
    assert -2800 < float(np.nanmedian(patch.z)) < -2200


@pytest.mark.network
def test_only_404_counts_as_absent():
    assert url_is_absent("https://copernicus-dem-90m.s3.amazonaws.com/does-not-exist.tif")
