"""ESA WorldCover land-cover fractions for Earth windows."""

from __future__ import annotations

import math

import numpy as np
import rasterio

from .bodies import EARTH
from .dem import CACHE_DIR, GDAL_ENV, url_is_absent

URL = "https://esa-worldcover.s3.eu-central-1.amazonaws.com/v200/2021/map/ESA_WorldCover_10m_2021_v200_{name}_Map.tif"
PRODUCT = "ESA WorldCover 10 m 2021 v200"
TILE_DEG = 3

# Reading the 1/64 overview gives ~640 m samples, about 10,000 per 60 km window.
OVERVIEW_INDEX = 5

CLASSES = {
    10: "tree",
    20: "shrub",
    30: "grass",
    40: "crop",
    50: "built",
    60: "bare",
    70: "snow",
    80: "water",
    90: "wetland",
    95: "mangrove",
    100: "moss",
}
COLUMNS = [f"lc_{name}" for name in CLASSES.values()]


def tile_floor(value: float) -> int:
    return math.floor(value / TILE_DEG) * TILE_DEG


def _tile_name(lat: int, lon: int) -> str:
    return f"{'N' if lat >= 0 else 'S'}{abs(lat):02d}{'E' if lon >= 0 else 'W'}{abs(lon):03d}"


def tile_path(lat: int, lon: int):
    return CACHE_DIR / "worldcover" / f"{_tile_name(lat, lon)}.npz"


def fetch_tile(lat: int, lon: int) -> np.ndarray | None:
    """Return the decimated class raster of one 3° tile, or None if there is none."""
    path = tile_path(lat, lon)
    if path.exists():
        with np.load(path) as data:
            return None if int(data["empty"]) else data["z"]

    array = None
    url = URL.format(name=_tile_name(lat, lon))
    try:
        with rasterio.Env(**GDAL_ENV), rasterio.open(f"/vsicurl/{url}", overview_level=OVERVIEW_INDEX) as overview:
            array = overview.read(1)
    except rasterio.errors.RasterioIOError:
        if not url_is_absent(url):
            raise

    path.parent.mkdir(parents=True, exist_ok=True)
    if array is None:
        np.savez_compressed(path, z=np.zeros((0, 0), "uint8"), empty=1)
        return None
    np.savez_compressed(path, z=array, empty=0)
    return array


def window_fractions(lat: float, lon: float, window_m: float, loader=fetch_tile) -> dict[str, float]:
    """Land-cover fractions of a window, as shares of its classified samples.

    ``lc_valid`` is the share of the samples read that have a class at all.
    """
    half_lat = math.degrees(window_m / 2 / EARTH.radius_m)
    half_lon = half_lat / math.cos(math.radians(lat))
    south, north, west, east = lat - half_lat, lat + half_lat, lon - half_lon, lon + half_lon

    samples = []
    for tile_lat in range(tile_floor(south), tile_floor(north) + TILE_DEG, TILE_DEG):
        for tile_lon in range(tile_floor(west), tile_floor(east) + TILE_DEG, TILE_DEG):
            tile = loader(tile_lat, tile_lon)
            if tile is None:
                continue
            height, width = tile.shape
            top = tile_lat + TILE_DEG
            row0 = max(int((top - north) / TILE_DEG * height), 0)
            row1 = min(int(math.ceil((top - south) / TILE_DEG * height)), height)
            col0 = max(int((west - tile_lon) / TILE_DEG * width), 0)
            col1 = min(int(math.ceil((east - tile_lon) / TILE_DEG * width)), width)
            if row1 > row0 and col1 > col0:
                samples.append(tile[row0:row1, col0:col1].ravel())

    result = dict.fromkeys(COLUMNS, 0.0)
    if not samples:
        return result | {"lc_valid": 0.0}
    values = np.concatenate(samples)
    classified = values[values > 0]
    if classified.size == 0:
        return result | {"lc_valid": 0.0}
    for code, name in CLASSES.items():
        result[f"lc_{name}"] = float(np.count_nonzero(classified == code) / classified.size)
    return result | {"lc_valid": float(classified.size / values.size)}
