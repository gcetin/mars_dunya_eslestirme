"""Sentinel-2 L2A surface summary for shortlisted Earth candidates.

Only shortlisted windows are read, on demand. Scenes come from the public
Earth Search STAC index and the ``sentinel-cogs`` bucket, which need no account.
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import asdict, dataclass

import numpy as np
import rasterio
import requests
from rasterio.enums import Resampling
from rasterio.warp import transform_bounds
from rasterio.windows import from_bounds

from .bodies import EARTH
from .dem import CACHE_DIR, GDAL_ENV

STAC_SEARCH = "https://earth-search.aws.element84.com/v1/search"
PRODUCT = "Sentinel-2 L2A (Earth Search / AWS Open Data)"
OUT_SIZE = 512
MIN_COVERAGE = 0.85

# SCL classes treated as unusable: no data, saturated, shadow, cloud, cirrus.
_BAD_SCL = {0, 1, 3, 8, 9, 10}


@dataclass
class SurfaceSummary:
    scene_id: str
    date: str
    cloud_pct: float
    coverage: float
    ndvi_median: float
    ndvi_p90: float
    green_frac: float
    ndmi_median: float
    swir_ratio_median: float


def _bbox(lat: float, lon: float, window_m: float) -> tuple[float, float, float, float]:
    half_lat = math.degrees(window_m / 2 / EARTH.radius_m)
    half_lon = half_lat / math.cos(math.radians(lat))
    return lon - half_lon, lat - half_lat, lon + half_lon, lat + half_lat


def _search(bbox, date_range: str) -> list[dict]:
    west, south, east, north = bbox
    # A small box at the centre keeps scenes that contain the window's middle.
    cx, cy = (west + east) / 2, (south + north) / 2
    body = {
        "collections": ["sentinel-2-l2a"],
        "bbox": [cx - 0.01, cy - 0.01, cx + 0.01, cy + 0.01],
        "datetime": date_range,
        "query": {"eo:cloud_cover": {"lt": 15}, "s2:nodata_pixel_percentage": {"lt": 25}},
        "sortby": [{"field": "properties.eo:cloud_cover", "direction": "asc"}],
        "limit": 5,
    }
    response = requests.post(STAC_SEARCH, json=body, timeout=90)
    response.raise_for_status()
    return response.json()["features"]


def _read(href: str, bbox, resampling: Resampling) -> np.ndarray:
    with rasterio.open(f"/vsicurl/{href}") as ds:
        bounds = transform_bounds("EPSG:4326", ds.crs, *bbox)
        window = from_bounds(*bounds, ds.transform)
        return ds.read(
            window=window,
            out_shape=(ds.count, OUT_SIZE, OUT_SIZE),
            boundless=True,
            fill_value=0,
            resampling=resampling,
        )


def _reflectance(item: dict, name: str, bbox) -> np.ndarray:
    """Surface reflectance of one band, with zero (no data) mapped to NaN."""
    asset = item["assets"][name]
    band = asset["raster:bands"][0]
    raw = _read(asset["href"], bbox, Resampling.average)[0].astype("float64")
    # When this flag is set the pixels are already corrected, although the
    # metadata still lists the -0.1 offset. Water near DN 60 confirms it;
    # applying the offset again produces NDVI values above 1.
    already_applied = item["properties"].get("earthsearch:boa_offset_applied", False)
    offset = 0.0 if already_applied else band.get("offset", 0.0)
    out = raw * band.get("scale", 1e-4) + offset
    out[raw == 0] = np.nan
    return out


def _ratio(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    with np.errstate(invalid="ignore", divide="ignore"):
        return (a - b) / (a + b)


def summarize(item: dict, bbox) -> tuple[SurfaceSummary, np.ndarray, np.ndarray]:
    """Return scene statistics, the RGB quicklook and the NDVI raster."""
    with rasterio.Env(**GDAL_ENV):
        red, nir = _reflectance(item, "red", bbox), _reflectance(item, "nir", bbox)
        swir1, swir2 = _reflectance(item, "swir16", bbox), _reflectance(item, "swir22", bbox)
        scl = _read(item["assets"]["scl"]["href"], bbox, Resampling.nearest)[0]
        rgb = _read(item["assets"]["visual"]["href"], bbox, Resampling.average)[:3]

    usable = ~np.isin(scl, list(_BAD_SCL)) & np.isfinite(red) & np.isfinite(nir)
    ndvi = _ratio(nir, red)
    ndmi = _ratio(nir, swir1)
    ndvi[~usable] = np.nan

    def median(values: np.ndarray) -> float:
        selected = values[usable & np.isfinite(values)]
        return float(np.median(selected)) if selected.size else float("nan")

    with np.errstate(invalid="ignore", divide="ignore"):
        swir_ratio = swir1 / swir2
    selected_ndvi = ndvi[usable]
    summary = SurfaceSummary(
        scene_id=item["id"],
        date=item["properties"]["datetime"][:10],
        cloud_pct=float(item["properties"]["eo:cloud_cover"]),
        coverage=float(usable.mean()),
        ndvi_median=median(ndvi),
        ndvi_p90=float(np.percentile(selected_ndvi, 90)) if selected_ndvi.size else float("nan"),
        green_frac=float(np.mean(selected_ndvi > 0.4)) if selected_ndvi.size else float("nan"),
        ndmi_median=median(ndmi),
        swir_ratio_median=median(swir_ratio),
    )
    return summary, np.moveaxis(rgb, 0, -1).astype("uint8"), ndvi


def surface_for_window(
    lat: float,
    lon: float,
    window_m: float = 60_000.0,
    date_range: str = "2024-06-01T00:00:00Z/2026-09-30T23:59:59Z",
):
    """Best available scene for a window, or None if no usable scene exists."""
    key = hashlib.sha1(f"{lat:.4f}|{lon:.4f}|{window_m:.0f}|{date_range}".encode()).hexdigest()[:16]
    cache_dir = CACHE_DIR / "sentinel2"
    json_path, array_path = cache_dir / f"{key}.json", cache_dir / f"{key}.npz"
    if json_path.exists():
        stored = json.loads(json_path.read_text())
        if stored is None:
            return None
        with np.load(array_path) as data:
            return SurfaceSummary(**stored), data["rgb"], data["ndvi"].astype("float64")

    bbox = _bbox(lat, lon, window_m)
    best = None
    for item in _search(bbox, date_range):
        if "raster:bands" not in item["assets"]["red"]:
            continue
        result = summarize(item, bbox)
        if best is None or result[0].coverage > best[0].coverage:
            best = result
        if best[0].coverage >= MIN_COVERAGE:
            break

    cache_dir.mkdir(parents=True, exist_ok=True)
    if best is None:
        json_path.write_text("null")
        return None
    summary, rgb, ndvi = best
    np.savez_compressed(array_path, rgb=rgb, ndvi=ndvi.astype("float32"))
    json_path.write_text(json.dumps(asdict(summary)))
    return best
