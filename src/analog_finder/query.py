"""Compute a window's features for any body, using the catalog's method."""

from __future__ import annotations

from .bodies import Body
from .catalog import ARCSEC
from .dem import Patch, read_earth_patch, read_global_patch
from .terrain import DEFAULT_SCALES_M, DEFAULT_WINDOW_M, Features, extract_features, margin_m


def read_window(body: Body, lat: float, lon: float, window_m: float = DEFAULT_WINDOW_M, scales_m=DEFAULT_SCALES_M, arcsec: int = ARCSEC) -> Patch:
    half_m = window_m / 2 + margin_m(scales_m)
    if body.key == "earth":
        return read_earth_patch(lat, lon, half_m, arcsec=arcsec)
    return read_global_patch(body, lat, lon, half_m)


def window_features(body: Body, lat: float, lon: float, window_m: float = DEFAULT_WINDOW_M, scales_m=DEFAULT_SCALES_M, arcsec: int = ARCSEC) -> tuple[Features, Patch]:
    patch = read_window(body, lat, lon, window_m, scales_m, arcsec)
    return extract_features(patch, lat, lon, window_m, scales_m), patch
