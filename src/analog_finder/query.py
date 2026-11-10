"""Compute a window's features for any body, using the catalog's method."""

from __future__ import annotations

from .bodies import Body
from .catalog import ARCSEC
from .dem import Patch, read_earth_patch, read_global_patch
from .terrain import DEFAULT_SCALES_M, DEFAULT_WINDOW_M, Features, extract_features, margin_m, read_half_m


def read_window(body: Body, lat: float, lon: float, window_m: float = DEFAULT_WINDOW_M, scales_m=DEFAULT_SCALES_M, arcsec: int = ARCSEC) -> Patch:
    """Read the patch every feature needs.

    Next to the ±180° map edge the wider read is refused; the narrower one
    still supports every feature, but the radial profile then searches for
    its centre over a smaller area.
    """
    read = (lambda half: read_earth_patch(lat, lon, half, arcsec=arcsec)) if body.key == "earth" else (lambda half: read_global_patch(body, lat, lon, half))
    try:
        return read(read_half_m(window_m, scales_m))
    except ValueError:
        return read(window_m / 2 + margin_m(scales_m))


def window_features(body: Body, lat: float, lon: float, window_m: float = DEFAULT_WINDOW_M, scales_m=DEFAULT_SCALES_M, arcsec: int = ARCSEC) -> tuple[Features, Patch]:
    patch = read_window(body, lat, lon, window_m, scales_m, arcsec)
    return extract_features(patch, lat, lon, window_m, scales_m), patch
