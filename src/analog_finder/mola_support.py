"""How much real laser data backs a MOLA grid window.

A 463 m cell may be interpolated from neighbouring orbit tracks. The MEGDR
"counts" product records how many laser shots fell in each 1/128° cell.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
import requests

from .bodies import MARS
from .dem import CACHE_DIR

BASE_URL = "https://pds-geosciences.wustl.edu/mgs/mgs-m-mola-5-megdr-l3-v1/mgsl_300x/meg128"
PPD = 128
TILE_LINES = 5632  # 44° of latitude
TILE_SAMPLES = 11520  # 90° of longitude


@dataclass
class ShotCoverage:
    measured_frac: float
    mean_shots: float
    product: str


def _tile_for(lat: float, lon_east: float) -> tuple[str, float, float]:
    """Tile name, its top latitude and its west longitude (0-360° East)."""
    for top in (88, 44, 0, -44):
        if top - 44 <= lat < top:
            break
    else:
        raise ValueError("Latitude outside the 128 ppd tiles.")
    west = int(lon_east // 90) * 90
    hemisphere = "s" if top < 0 else "n"
    return f"megc{abs(top):02d}{hemisphere}{west:03d}hb.img", float(top), float(west)


def shot_counts(lat: float, lon: float, window_m: float) -> tuple[np.ndarray, str]:
    """Laser shots per 1/128° cell over a window; row 0 is north.

    A window crossing a tile boundary is clipped to the tile holding its centre.
    """
    lon_east = lon % 360.0
    name, top, west = _tile_for(lat, lon_east)

    half_lat = math.degrees(window_m / 2 / MARS.radius_m)
    half_lon = half_lat / math.cos(math.radians(lat))
    line0 = max(int((top - lat - half_lat) * PPD), 0)
    line1 = min(int((top - lat + half_lat) * PPD) + 1, TILE_LINES)
    sample0 = max(int((lon_east - half_lon - west) * PPD), 0)
    sample1 = min(int((lon_east + half_lon - west) * PPD) + 1, TILE_SAMPLES)

    cache = CACHE_DIR / "mola_counts" / f"{name}_{line0}_{line1}.bin"
    if cache.exists():
        raw = cache.read_bytes()
    else:
        start = line0 * TILE_SAMPLES
        end = line1 * TILE_SAMPLES - 1
        response = requests.get(f"{BASE_URL}/{name}", headers={"Range": f"bytes={start}-{end}"}, timeout=180)
        if response.status_code != 206:
            raise RuntimeError(f"Range request refused for {name}: HTTP {response.status_code}")
        raw = response.content
        cache.parent.mkdir(parents=True, exist_ok=True)
        cache.write_bytes(raw)

    rows = np.frombuffer(raw, dtype="uint8").reshape(line1 - line0, TILE_SAMPLES)
    return rows[:, sample0:sample1], name


def shot_coverage(lat: float, lon: float, window_m: float) -> ShotCoverage:
    """Fraction of 1/128° cells in a window containing at least one laser shot."""
    counts, name = shot_counts(lat, lon, window_m)
    return ShotCoverage(float((counts > 0).mean()), float(counts.mean()), f"MEGDR counts {name}")
