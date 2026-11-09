"""Read DEM patches for Mars, the Moon and Earth into a common structure.

All remote access is by HTTP range request, so only the bytes covering the
requested area are transferred. Results are cached under ``data/cache``.
"""

from __future__ import annotations

import hashlib
import math
import os
import time
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import numpy as np
import rasterio
import requests
from rasterio.enums import Resampling
from rasterio.windows import Window

from .bodies import EARTH, Body

CACHE_DIR = Path(os.environ.get("ANALOG_CACHE", Path(__file__).resolve().parents[2] / "data" / "cache"))

GDAL_ENV = {
    "GDAL_DISABLE_READDIR_ON_OPEN": "EMPTY_DIR",
    "CPL_VSIL_CURL_ALLOWED_EXTENSIONS": ".tif",
    "GDAL_HTTP_MULTIRANGE": "YES",
    "GDAL_HTTP_MERGE_CONSECUTIVE_RANGES": "YES",
    "GDAL_HTTP_MAX_RETRY": "5",
    "GDAL_HTTP_RETRY_DELAY": "2",
    "VSI_CACHE": "TRUE",
}

# Beyond this latitude a cylindrical global mosaic needs a polar product.
MAX_ABS_LAT = 80.0

NASADEM_URL = "https://nasademeuwest.blob.core.windows.net/nasadem-cog/v001/NASADEM_HGT_{name}.tif"
COPDEM_URLS = {
    30: "https://copernicus-dem-30m.s3.amazonaws.com/Copernicus_DSM_COG_10_{name}_DEM/Copernicus_DSM_COG_10_{name}_DEM.tif",
    90: "https://copernicus-dem-90m.s3.amazonaws.com/Copernicus_DSM_COG_30_{name}_DEM/Copernicus_DSM_COG_30_{name}_DEM.tif",
}
_INT16_NODATA = -32768


@dataclass
class Patch:
    """Elevations on a regular latitude/longitude grid.

    ``z`` is in metres with NaN for missing data. Row 0 is the northern edge.
    ``lat`` and ``lon`` hold pixel-centre coordinates in degrees.
    """

    z: np.ndarray
    lat: np.ndarray
    lon: np.ndarray
    body: Body
    source: str

    @property
    def dlat(self) -> float:
        return float(abs(self.lat[1] - self.lat[0]))

    @property
    def dlon(self) -> float:
        return float(abs(self.lon[1] - self.lon[0]))

    @property
    def native_m(self) -> float:
        """North-south pixel spacing in metres."""
        return math.radians(self.dlat) * self.body.radius_m


def _extent_deg(body: Body, lat: float, half_m: float) -> tuple[float, float]:
    """Half-width in degrees of latitude and longitude covering ``half_m``."""
    half_lat = math.degrees(half_m / body.radius_m) * 1.02
    max_abs = min(abs(lat) + half_lat, 89.0)
    half_lon = half_lat / math.cos(math.radians(max_abs))
    return half_lat, half_lon


def read_global_patch(body: Body, lat: float, lon: float, half_m: float, *, use_cache: bool = True) -> Patch:
    """Read a window from the body's global simple-cylindrical mosaic."""
    if body.dem is None:
        raise ValueError(f"{body.name} has no global mosaic; use read_earth_patch.")
    if abs(lat) > MAX_ABS_LAT:
        raise ValueError(
            f"|lat| > {MAX_ABS_LAT:.0f}° requires a polar stereographic product, "
            "which is not integrated yet."
        )

    key = hashlib.sha1(f"{body.dem.url}|{lat:.4f}|{lon:.4f}|{half_m:.0f}".encode()).hexdigest()[:16]
    cache_file = CACHE_DIR / "patches" / f"{body.key}_{key}.npz"
    if use_cache and cache_file.exists():
        with np.load(cache_file) as data:
            return Patch(data["z"].astype("float64"), data["lat"], data["lon"], body, body.dem.id)

    half_lat, half_lon = _extent_deg(body, lat, half_m)
    if abs(lon) + half_lon > 180.0:
        raise ValueError("Windows crossing the map edge at ±180° are not supported yet.")

    with rasterio.Env(**GDAL_ENV), rasterio.open(f"/vsicurl/{body.dem.url}") as ds:
        pixel_deg = math.degrees(abs(ds.transform.a) / body.radius_m)
        origin_lon = math.degrees(ds.transform.c / body.radius_m)
        origin_lat = math.degrees(ds.transform.f / body.radius_m)
        col0 = max(int(math.floor((lon - half_lon - origin_lon) / pixel_deg)), 0)
        col1 = min(int(math.ceil((lon + half_lon - origin_lon) / pixel_deg)), ds.width)
        row0 = max(int(math.floor((origin_lat - lat - half_lat) / pixel_deg)), 0)
        row1 = min(int(math.ceil((origin_lat - lat + half_lat) / pixel_deg)), ds.height)
        raw = ds.read(1, window=Window(col0, row0, col1 - col0, row1 - row0))
        # Mosaics may store elevation with a scale factor (LOLA uses 0.5 m/DN).
        z = raw.astype("float64") * ds.scales[0] + ds.offsets[0]
        if ds.nodata is not None:
            z[raw == ds.nodata] = np.nan

    lats = origin_lat - (np.arange(row0, row1) + 0.5) * pixel_deg
    lons = origin_lon + (np.arange(col0, col1) + 0.5) * pixel_deg

    cache_file.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(cache_file, z=z.astype("float32"), lat=lats, lon=lons)
    return Patch(z, lats, lons, body, body.dem.id)


# --- Earth ------------------------------------------------------------------

_pc_token: tuple[str, float] | None = None


def _planetary_computer_token() -> str:
    """Anonymous, short-lived read token for the NASADEM COG container."""
    global _pc_token
    if _pc_token is None or time.time() > _pc_token[1]:
        response = requests.get(
            "https://planetarycomputer.microsoft.com/api/sas/v1/token/nasadem", timeout=60
        )
        response.raise_for_status()
        payload = response.json()
        # The service may hand out a token that is already close to expiry.
        expiry = datetime.fromisoformat(payload["msft:expiry"].replace("Z", "+00:00")).timestamp()
        _pc_token = (payload["token"], expiry - 120)
    return _pc_token[0]


def _tile_name_nasadem(lat: int, lon: int) -> str:
    return f"{'n' if lat >= 0 else 's'}{abs(lat):02d}{'e' if lon >= 0 else 'w'}{abs(lon):03d}"


def _tile_name_copdem(lat: int, lon: int) -> str:
    return f"{'N' if lat >= 0 else 'S'}{abs(lat):02d}_00_{'E' if lon >= 0 else 'W'}{abs(lon):03d}_00"


def url_is_absent(url: str) -> bool:
    """True only if the server states the object does not exist.

    GDAL reports throttling and truncated responses with the same exception as
    a missing file. Caching those as "open sea" would silently delete land.
    """
    response = requests.head(url, timeout=60, allow_redirects=True)
    return response.status_code == 404


def default_earth_source(lat_floor: int) -> str:
    """NASADEM where it exists, Copernicus DEM for high latitudes."""
    return "nasadem" if -56 <= lat_floor < 60 else "copdem"


def earth_tile_path(lat: int, lon: int, arcsec: int, source: str) -> Path:
    return CACHE_DIR / "earth" / f"{source}_{arcsec}s" / f"{_tile_name_nasadem(lat, lon)}.npz"


def fetch_earth_tile(lat: int, lon: int, arcsec: int, source: str | None = None) -> np.ndarray | None:
    """Return a 1°×1° tile sampled every ``arcsec`` seconds, or None over open sea.

    Row 0 lies at ``lat + 1`` and column 0 at ``lon``. Tiles of different
    products and latitudes all end up on this same grid.
    """
    source = source or default_earth_source(lat)
    path = earth_tile_path(lat, lon, arcsec, source)
    if path.exists():
        with np.load(path) as data:
            if int(data["empty"]):
                return None
            return _decode(data["z"])

    n = 3600 // arcsec
    if source == "nasadem":
        name = _tile_name_nasadem(lat, lon)
        url = NASADEM_URL.format(name=name) + "?" + _planetary_computer_token()
    elif source == "copdem":
        # GLO-90 is sufficient, and much lighter, for any grid of 3" or coarser.
        url = COPDEM_URLS[30 if arcsec < 3 else 90].format(name=_tile_name_copdem(lat, lon))
    else:
        raise ValueError(f"Unknown Earth DEM source: {source}")

    z = None
    try:
        with rasterio.Env(**GDAL_ENV), rasterio.open(f"/vsicurl/{url}") as ds:
            # NASADEM tiles carry a duplicated last row and column; drop them so
            # that adjacent tiles butt together.
            width = ds.width - 1 if source == "nasadem" else ds.width
            height = ds.height - 1 if source == "nasadem" else ds.height
            raw = ds.read(
                1,
                window=Window(0, 0, width, height),
                out_shape=(n, n),
                resampling=Resampling.average if n < min(width, height) else Resampling.bilinear,
                masked=True,
            )
            z = raw.astype("float64").filled(np.nan)
    except rasterio.errors.RasterioIOError:
        if not url_is_absent(url):
            raise

    path.parent.mkdir(parents=True, exist_ok=True)
    if z is None:
        np.savez_compressed(path, z=np.zeros((0, 0), "int16"), empty=1)
        return None
    np.savez_compressed(path, z=_encode(z), empty=0)
    return z.astype("float32")


def _encode(z: np.ndarray) -> np.ndarray:
    return np.where(np.isfinite(z), np.round(z), _INT16_NODATA).astype("int16")


def _decode(z: np.ndarray) -> np.ndarray:
    return np.where(z == _INT16_NODATA, np.nan, z.astype("float32")).astype("float32")


def read_earth_patch(
    lat: float,
    lon: float,
    half_m: float,
    *,
    arcsec: int = 1,
    source: str | None = None,
    tile_loader=None,
) -> Patch:
    """Mosaic Earth tiles around a point on a regular ``arcsec`` grid."""
    half_lat, half_lon = _extent_deg(EARTH, lat, half_m)
    if abs(lon) + half_lon > 180.0:
        raise ValueError("Windows crossing the antimeridian are not supported yet.")

    loader = tile_loader or (lambda la, lo: fetch_earth_tile(la, lo, arcsec, source))
    lat_floors = range(math.floor(lat - half_lat), math.floor(lat + half_lat) + 1)
    lon_floors = range(math.floor(lon - half_lon), math.floor(lon + half_lon) + 1)
    n = 3600 // arcsec
    mosaic = np.full((len(lat_floors) * n, len(lon_floors) * n), np.nan, dtype="float32")
    used = set()
    for j, tile_lat in enumerate(reversed(lat_floors)):
        for i, tile_lon in enumerate(lon_floors):
            tile = loader(tile_lat, tile_lon)
            if tile is not None:
                mosaic[j * n : (j + 1) * n, i * n : (i + 1) * n] = tile
                used.add(source or default_earth_source(tile_lat))

    step = arcsec / 3600.0
    top = lat_floors[-1] + 1
    west = lon_floors[0]
    lats = top - np.arange(mosaic.shape[0]) * step
    lons = west + np.arange(mosaic.shape[1]) * step
    rows = np.nonzero((lats <= lat + half_lat) & (lats >= lat - half_lat))[0]
    cols = np.nonzero((lons >= lon - half_lon) & (lons <= lon + half_lon))[0]
    z = mosaic[rows[0] : rows[-1] + 1, cols[0] : cols[-1] + 1].astype("float64")
    label = "+".join(sorted(used)) or "no data"
    return Patch(z, lats[rows[0] : rows[-1] + 1], lons[cols[0] : cols[-1] + 1], EARTH, f"{label} @ {arcsec} arc-sec")
