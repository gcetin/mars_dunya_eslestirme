"""Build a searchable catalog of Earth terrain windows.

Features are computed once for a grid of windows, so any Mars or Moon query
afterwards is an in-memory comparison.
"""

from __future__ import annotations

import json
import math
import time
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
from rasterio.errors import RasterioIOError

from . import landcover
from .bodies import EARTH
from .dem import _extent_deg, default_earth_source, earth_tile_path, fetch_earth_tile, read_earth_patch
from .provenance import METHOD, data_sources, software_versions, utc_now
from .terrain import DEFAULT_SCALES_M, DEFAULT_WINDOW_M, extract_features, feature_names, margin_m

ARCSEC = 15
MAX_SEA_FRAC = 0.3
CATALOG_DIR = Path(__file__).resolve().parents[2] / "data" / "catalog"

GLOBAL_REGIONS = [
    (-56.0, 83.0, -179.0, 179.0),
    (-79.0, -76.0, 160.0, 165.0),  # McMurdo Dry Valleys
]


def window_centers(regions, window_m: float, stride_m: float) -> list[tuple[float, float]]:
    """Window centres spaced ``stride_m`` apart on the ground."""
    centers = []
    lat_step = math.degrees(stride_m / EARTH.radius_m)
    half_m = window_m / 2 + margin_m()
    for south, north, west, east in regions:
        for lat in np.arange(south, north + 1e-9, lat_step):
            lon_step = lat_step / max(math.cos(math.radians(lat)), 0.02)
            for lon in np.arange(west, east + 1e-9, lon_step):
                _, half_lon = _extent_deg(EARTH, float(lat), half_m)
                if abs(lon) + half_lon < 180.0:
                    centers.append((round(float(lat), 4), round(float(lon), 4)))
    return sorted(set(centers))


def _tiles_for(centers, window_m: float) -> tuple[set, set]:
    """1° DEM tiles and 3° land-cover tiles needed by a list of centres."""
    half_m = window_m / 2 + margin_m()
    dem, cover = set(), set()
    for lat, lon in centers:
        half_lat, half_lon = _extent_deg(EARTH, lat, half_m)
        for la in range(math.floor(lat - half_lat), math.floor(lat + half_lat) + 1):
            for lo in range(math.floor(lon - half_lon), math.floor(lon + half_lon) + 1):
                dem.add((la, lo))
        for la in range(landcover.tile_floor(lat - half_lat), landcover.tile_floor(lat + half_lat) + 3, 3):
            for lo in range(landcover.tile_floor(lon - half_lon), landcover.tile_floor(lon + half_lon) + 3, 3):
                cover.add((la, lo))
    return dem, cover


def _with_retries(function, *args, attempts: int = 4):
    for attempt in range(attempts):
        try:
            return function(*args)
        except (RasterioIOError, OSError):
            if attempt == attempts - 1:
                raise
            time.sleep(2 * (attempt + 1))


def prefetch(tiles: set, kind: str, threads: int = 12, log=print) -> list:
    """Download tiles concurrently; return the ones that kept failing."""
    if kind == "dem":
        pending = [t for t in tiles if not earth_tile_path(*t, ARCSEC, default_earth_source(t[0])).exists()]
        job = lambda t: _with_retries(fetch_earth_tile, t[0], t[1], ARCSEC)
    else:
        pending = [t for t in tiles if not landcover.tile_path(*t).exists()]
        job = lambda t: _with_retries(landcover.fetch_tile, *t)

    log(f"[{kind}] {len(tiles)} tiles, {len(pending)} to download")
    failed = []
    started = time.time()
    with ThreadPoolExecutor(max_workers=threads) as pool:
        futures = {pool.submit(job, t): t for t in pending}
        for i, future in enumerate(as_completed(futures), 1):
            try:
                future.result()
            except Exception as exc:  # recorded and reported in the catalog metadata
                failed.append((*futures[future], str(exc)[:160]))
            if i % 500 == 0:
                log(f"[{kind}] {i}/{len(pending)} in {time.time() - started:.0f}s, failed {len(failed)}")
    return failed


@lru_cache(maxsize=256)
def _cached_dem(lat: int, lon: int):
    path = earth_tile_path(lat, lon, ARCSEC, default_earth_source(lat))
    if not path.exists():
        return None
    with np.load(path) as data:
        return None if int(data["empty"]) else np.where(data["z"] == -32768, np.nan, data["z"].astype("float32")).astype("float32")


@lru_cache(maxsize=32)
def _cached_cover(lat: int, lon: int):
    path = landcover.tile_path(lat, lon)
    if not path.exists():
        return None
    with np.load(path) as data:
        return None if int(data["empty"]) else data["z"]


def _process_chunk(args) -> list[dict]:
    centers, window_m, scales = args
    half_m = window_m / 2 + margin_m(scales)
    rows = []
    for lat, lon in centers:
        patch = read_earth_patch(lat, lon, half_m, arcsec=ARCSEC, tile_loader=_cached_dem)
        if not np.isfinite(patch.z).any():
            continue
        sea = float(np.mean((patch.z == 0) | ~np.isfinite(patch.z)))
        if sea > MAX_SEA_FRAC:
            continue
        features = extract_features(patch, lat, lon, window_m, scales)
        if features.valid_frac < 0.95 or not np.isfinite(features.values["slope_mean_1km"]):
            continue
        cover = landcover.window_fractions(lat, lon, window_m, loader=_cached_cover)
        rows.append({"lat": lat, "lon": lon, "sea_frac": sea, "valid_frac": features.valid_frac, **features.values, **cover})
    return rows


def build_catalog(
    regions,
    output: Path,
    *,
    window_m: float = DEFAULT_WINDOW_M,
    stride_m: float = 20_000.0,
    scales_m=DEFAULT_SCALES_M,
    workers: int = 8,
    threads: int = 12,
    log=print,
) -> pd.DataFrame:
    started = time.time()
    centers = window_centers(regions, window_m, stride_m)
    log(f"{len(centers)} window centres")
    dem_tiles, cover_tiles = _tiles_for(centers, window_m)
    failed_dem = prefetch(dem_tiles, "dem", threads, log)
    failed_cover = prefetch(cover_tiles, "cover", threads, log)
    fetched = time.time()

    chunks = [centers[i : i + 400] for i in range(0, len(centers), 400)]
    rows: list[dict] = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(_process_chunk, (chunk, window_m, scales_m)) for chunk in chunks]
        for i, future in enumerate(as_completed(futures), 1):
            rows.extend(future.result())
            if i % 25 == 0:
                log(f"[features] {i}/{len(chunks)} chunks, {len(rows)} valid windows")

    frame = pd.DataFrame(rows).sort_values(["lat", "lon"]).reset_index(drop=True)
    measurements = frame.columns.difference(["lat", "lon"])
    frame[measurements] = frame[measurements].astype("float32")
    output.parent.mkdir(parents=True, exist_ok=True)
    frame.to_parquet(output, index=False)

    metadata = {
        "created_utc": utc_now(),
        "regions": regions,
        "window_m": window_m,
        "stride_m": stride_m,
        "scales_m": list(scales_m),
        "earth_grid_arcsec": ARCSEC,
        "n_centers": len(centers),
        "n_valid_windows": len(frame),
        "n_dem_tiles": len(dem_tiles),
        "failed_dem_tiles": failed_dem,
        "failed_cover_tiles": failed_cover,
        "feature_names": feature_names(scales_m),
        "method": METHOD,
        "data_sources": data_sources(),
        "software": software_versions(),
        "seconds_fetch": round(fetched - started),
        "seconds_features": round(time.time() - fetched),
    }
    output.with_suffix(".json").write_text(json.dumps(metadata, indent=2, ensure_ascii=False))
    log(f"Done: {len(frame)} windows -> {output}; failed tiles dem={len(failed_dem)} cover={len(failed_cover)}")
    return frame


def load_catalog(path: Path) -> tuple[pd.DataFrame, dict]:
    frame = pd.read_parquet(path)
    return frame, json.loads(path.with_suffix(".json").read_text())
