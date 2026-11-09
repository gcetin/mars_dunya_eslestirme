"""Physically interpretable terrain features at matched physical scales.

The pipeline for every body is identical:

1. Low-pass filter the DEM so its spatial support equals the analysis scale,
   taking into account the support the source pixels already have.
2. Resample onto a local azimuthal-equidistant grid in metres, using the
   body's own radius. Geographic coordinates are never compared across bodies.
3. Compute features only from local, metric morphology.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy import ndimage

from .dem import Patch

DEFAULT_WINDOW_M = 60_000.0
DEFAULT_SCALES_M = (1_000.0, 2_000.0, 4_000.0)

# Gaussian sigma as a fraction of the analysis scale. At the Nyquist wavelength
# of the output grid this keeps about 29 % of the amplitude, limiting aliasing.
SIGMA_PER_SCALE = 0.5

# Largest filter reach needed around a window, in units of analysis scale.
MARGIN_PER_SCALE = 3 * SIGMA_PER_SCALE

# A scale is only meaningful if the source grid can represent it.
MIN_SCALE_PER_PIXEL = 1.9

PER_SCALE_FEATURES = ("slope_mean", "slope_std", "slope_p95", "tri", "rough")
WINDOW_FEATURES = ("relief", "elev_std", "hypso", "skew", "radial_r2", "basin_depth")


def margin_m(scales_m=DEFAULT_SCALES_M) -> float:
    return MARGIN_PER_SCALE * max(scales_m)


def feature_names(scales_m=DEFAULT_SCALES_M) -> list[str]:
    names = [f"{name}_{s / 1000:g}km" for s in scales_m for name in PER_SCALE_FEATURES]
    return names + list(WINDOW_FEATURES)


def aeqd_to_latlon(x, y, lat0: float, lon0: float, radius: float):
    """Inverse azimuthal-equidistant projection on a sphere (degrees)."""
    phi0 = math.radians(lat0)
    rho = np.hypot(x, y)
    c = rho / radius
    sin_c, cos_c = np.sin(c), np.cos(c)
    with np.errstate(invalid="ignore", divide="ignore"):
        y_over_rho = np.where(rho > 0, y / rho, 0.0)
    phi = np.arcsin(np.clip(cos_c * math.sin(phi0) + y_over_rho * sin_c * math.cos(phi0), -1, 1))
    lam = np.arctan2(x * sin_c, rho * math.cos(phi0) * cos_c - y * math.sin(phi0) * sin_c)
    return np.degrees(phi), lon0 + np.degrees(lam)


def _nan_gaussian(z: np.ndarray, sigma_px: tuple[float, float]) -> np.ndarray:
    """Gaussian filter that ignores NaN cells by normalised convolution."""
    valid = np.isfinite(z)
    if valid.all():
        return ndimage.gaussian_filter(z, sigma_px, mode="nearest", truncate=3.0)
    weights = ndimage.gaussian_filter(valid.astype("float64"), sigma_px, mode="nearest", truncate=3.0)
    values = ndimage.gaussian_filter(np.where(valid, z, 0.0), sigma_px, mode="nearest", truncate=3.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = values / weights
    out[weights < 0.5] = np.nan
    return out


def resample_to_scale(patch: Patch, lat0: float, lon0: float, window_m: float, scale_m: float) -> np.ndarray:
    """Return the DEM on a north-up metric grid with spacing ``scale_m``."""
    radius = patch.body.radius_m
    dy_m = math.radians(patch.dlat) * radius
    dx_m = math.radians(patch.dlon) * radius * math.cos(math.radians(lat0))
    if scale_m < MIN_SCALE_PER_PIXEL * max(dx_m, dy_m):
        raise ValueError(f"Scale {scale_m:.0f} m is too fine for {max(dx_m, dy_m):.0f} m pixels.")

    # A source pixel is a box of width d, with variance d²/12, so only the
    # remaining variance is applied. Coarse and fine sources converge on the
    # same total support.
    target_var = (SIGMA_PER_SCALE * scale_m) ** 2
    sigma_y = math.sqrt(max(target_var - dy_m**2 / 12, 0.0)) / dy_m
    sigma_x = math.sqrt(max(target_var - dx_m**2 / 12, 0.0)) / dx_m
    smooth = _nan_gaussian(patch.z, (sigma_y, sigma_x))

    n = int(round(window_m / scale_m))
    axis = (np.arange(n) - (n - 1) / 2) * scale_m
    x, y = np.meshgrid(axis, axis[::-1])
    lat, lon = aeqd_to_latlon(x, y, lat0, lon0, radius)
    rows = (patch.lat[0] - lat) / patch.dlat
    cols = (lon - patch.lon[0]) / patch.dlon
    return ndimage.map_coordinates(smooth, [rows, cols], order=1, mode="constant", cval=np.nan)


def slope_deg(z: np.ndarray, spacing_m: float) -> np.ndarray:
    dz_dy, dz_dx = np.gradient(z, spacing_m)
    return np.degrees(np.arctan(np.hypot(dz_dx, dz_dy)))


def _tri(z: np.ndarray) -> np.ndarray:
    """Terrain Ruggedness Index: mean absolute difference to the 8 neighbours."""
    diffs = []
    for dy in (-1, 0, 1):
        for dx in (-1, 0, 1):
            if dx or dy:
                diffs.append(np.abs(z[1:-1, 1:-1] - np.roll(np.roll(z, dy, 0), dx, 1)[1:-1, 1:-1]))
    return np.mean(diffs, axis=0)


def _detrend(z: np.ndarray, spacing_m: float) -> np.ndarray:
    """Remove the best-fitting plane so that regional tilt does not look like shape."""
    ny, nx = z.shape
    y, x = np.mgrid[:ny, :nx] * spacing_m
    design = np.column_stack([x.ravel(), y.ravel(), np.ones(z.size)])
    coef, *_ = np.linalg.lstsq(design, z.ravel(), rcond=None)
    return z - (design @ coef).reshape(z.shape)


def radial_profile(z: np.ndarray, spacing_m: float, bins: int = 12):
    """Azimuthal mean of plane-detrended elevation around the window centre.

    Returns radius (m), the mean profile, and the fraction of the elevation
    variance inside the inscribed circle that the profile explains.
    """
    detrended = _detrend(z, spacing_m)
    ny, nx = z.shape
    y, x = np.mgrid[:ny, :nx]
    r = np.hypot(x - (nx - 1) / 2, y - (ny - 1) / 2) * spacing_m
    r_max = (min(nx, ny) - 1) / 2 * spacing_m
    inside = r <= r_max
    index = np.minimum((r[inside] / r_max * bins).astype(int), bins - 1)
    values = detrended[inside]
    counts = np.bincount(index, minlength=bins)
    means = np.bincount(index, weights=values, minlength=bins) / np.maximum(counts, 1)
    centres = (np.arange(bins) + 0.5) / bins * r_max
    total = float(np.sum((values - values.mean()) ** 2))
    r2 = 0.0 if total == 0 else 1.0 - float(np.sum((values - means[index]) ** 2)) / total
    return centres, means, r2


def hypsometric_curve(z: np.ndarray, points: int = 41):
    """Share of the area above each normalised elevation."""
    low, high = np.percentile(z, [1, 99])
    levels = np.linspace(0, 1, points)
    span = max(high - low, 1e-6)
    area_above = np.array([np.mean(z >= low + level * span) for level in levels])
    return area_above, levels


def slope_cdf(z: np.ndarray, spacing_m: float, points: int = 61):
    """Slope values against the share of the area at or below them."""
    slope = np.sort(slope_deg(z, spacing_m)[1:-1, 1:-1].ravel())
    grid = np.linspace(0, 1, points)
    return np.quantile(slope, grid), grid


def scale_features(z: np.ndarray, spacing_m: float) -> dict[str, float]:
    """Slope and ruggedness statistics at one analysis scale."""
    slope = slope_deg(z, spacing_m)[1:-1, 1:-1]
    residual = z - ndimage.uniform_filter(z, 3, mode="nearest")
    return {
        "slope_mean": float(np.mean(slope)),
        "slope_std": float(np.std(slope)),
        "slope_p95": float(np.percentile(slope, 95)),
        "tri": float(np.mean(_tri(z))),
        "rough": float(np.std(residual[1:-1, 1:-1])),
    }


def window_features(z: np.ndarray, spacing_m: float) -> dict[str, float]:
    """Shape of the whole window, independent of absolute elevation."""
    low, high = np.percentile(z, [1, 99])
    relief = float(high - low)
    detrended = _detrend(z, spacing_m)
    std = float(np.std(detrended))
    skew = 0.0 if std == 0 else float(np.mean((detrended - detrended.mean()) ** 3) / std**3)
    hypso = 0.5 if relief == 0 else float(np.clip((np.mean(z) - low) / relief, 0, 1))

    radius, profile, r2 = radial_profile(z, spacing_m)
    r_max = radius[-1]
    centre = profile[radius <= 0.4 * r_max].mean()
    rim_zone = profile[radius >= 0.75 * r_max].mean()
    return {
        "relief": relief,
        "elev_std": std,
        "hypso": hypso,
        "skew": skew,
        "radial_r2": r2,
        "basin_depth": float(rim_zone - centre),
    }


@dataclass
class Features:
    values: dict[str, float]
    valid_frac: float
    grids: dict[float, np.ndarray]


def extract_features(
    patch: Patch,
    lat0: float,
    lon0: float,
    window_m: float = DEFAULT_WINDOW_M,
    scales_m=DEFAULT_SCALES_M,
) -> Features:
    """Compute the multi-scale feature vector of one window.

    ``valid_frac`` is the share of finest-scale cells backed by data. Values
    are NaN when it is below 0.95.
    """
    grids = {s: resample_to_scale(patch, lat0, lon0, window_m, s) for s in scales_m}
    finest = min(scales_m)
    valid_frac = float(np.isfinite(grids[finest]).mean())
    names = feature_names(scales_m)
    if valid_frac < 0.95 or not all(np.isfinite(g).mean() > 0.95 for g in grids.values()):
        return Features(dict.fromkeys(names, float("nan")), valid_frac, grids)

    filled = {s: _fill(g) for s, g in grids.items()}
    values: dict[str, float] = {}
    for s in scales_m:
        for name, value in scale_features(filled[s], s).items():
            values[f"{name}_{s / 1000:g}km"] = value
    values.update(window_features(filled[finest], finest))
    return Features(values, valid_frac, filled)


def _fill(z: np.ndarray) -> np.ndarray:
    """Replace the few remaining gaps with the nearest valid value."""
    mask = ~np.isfinite(z)
    if not mask.any():
        return z
    indices = ndimage.distance_transform_edt(mask, return_distances=False, return_indices=True)
    return z[tuple(indices)]


def hillshade(z: np.ndarray, spacing_m: float, azimuth: float = 315.0, altitude: float = 45.0) -> np.ndarray:
    """Shaded relief in the range 0 to 1, for display."""
    dz_dy, dz_dx = np.gradient(z, spacing_m)
    slope = np.arctan(np.hypot(dz_dx, dz_dy))
    # dz_dy is taken along the row index, which grows southwards.
    aspect = np.arctan2(dz_dy, -dz_dx)
    zenith = math.radians(90 - altitude)
    az = math.radians((360 - azimuth + 90) % 360)
    shade = math.cos(zenith) * np.cos(slope) + math.sin(zenith) * np.sin(slope) * np.cos(az - aspect)
    return np.clip(shade, 0, 1)
