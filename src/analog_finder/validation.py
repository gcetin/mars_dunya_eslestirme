"""Experiments measuring how far the method can be trusted."""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pandas as pd
from scipy import ndimage
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator

from .catalog import ARCSEC, _cached_dem
from .dem import Patch, read_earth_patch
from .similarity import DEFAULT_WEIGHTS, Normalizer, feature_weights, weighted_distance
from .terrain import DEFAULT_SCALES_M, DEFAULT_WINDOW_M, PROFILE_PREFIX, extract_features, read_half_m


def resample_mask(mask: np.ndarray, shape: tuple[int, int]) -> np.ndarray:
    """Nearest-neighbour resize of a boolean mask."""
    zoom = (shape[0] / mask.shape[0], shape[1] / mask.shape[1])
    return ndimage.zoom(mask.astype("uint8"), zoom, order=0, grid_mode=True, mode="nearest").astype(bool)


def with_track_gaps(patch: Patch, measured: np.ndarray) -> Patch:
    """Rebuild a DEM from only the cells a real orbit-track pattern would measure.

    Gaps are filled by linear interpolation, and outside the convex hull by the
    nearest measured value. MEGDR's own interpolator differs in detail, so this
    gives the size of the effect, not an exact replica of the product.
    """
    measured = resample_mask(measured, patch.z.shape) & np.isfinite(patch.z)
    rows, cols = np.nonzero(measured)
    points = np.column_stack([rows, cols])
    values = patch.z[measured]
    grid_rows, grid_cols = np.mgrid[: patch.z.shape[0], : patch.z.shape[1]]
    filled = LinearNDInterpolator(points, values)(grid_rows, grid_cols)
    holes = ~np.isfinite(filled)
    if holes.any():
        filled[holes] = NearestNDInterpolator(points, values)(grid_rows[holes], grid_cols[holes])
    return replace(patch, z=filled, source=f"{patch.source} + MOLA track gaps")


def self_retrieval(
    catalog: pd.DataFrame,
    names: list[str],
    measured: np.ndarray,
    *,
    n: int = 200,
    seed: int = 7,
    window_m: float = DEFAULT_WINDOW_M,
    scales_m=DEFAULT_SCALES_M,
    family_weights: dict[str, float] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Degrade known Earth windows to MOLA-like coverage and search for them.

    Each sampled window is rebuilt from a real MOLA track pattern, then used as
    a query. ``measured`` must cover ``2 * read_half_m(window_m, scales_m)``.
    Returns one row per window (rank of its true location among all windows)
    and the median absolute change of every feature, in percent. The radial
    profile is summarised in one row, ``rprof``: the RMS change of the whole
    profile relative to its RMS, since single rings can be close to zero.

    This is an upper bound on real performance: query and catalog come from
    the same planet, so it tests the method and sampling, not cross-planet validity.
    """
    rng = np.random.default_rng(seed)
    sample = catalog.iloc[rng.choice(len(catalog), size=min(n, len(catalog)), replace=False)]
    normalizer = Normalizer.fit(catalog, names)
    z_all = normalizer.apply_frame(catalog)
    weights = feature_weights(names, family_weights or DEFAULT_WEIGHTS)
    half_m = read_half_m(window_m, scales_m)
    profile_names = [n_ for n_ in names if n_.startswith(PROFILE_PREFIX)]
    other_names = [n_ for n_ in names if not n_.startswith(PROFILE_PREFIX)]

    rows, changes = [], []
    for index, row in sample.iterrows():
        try:
            patch = read_earth_patch(row["lat"], row["lon"], half_m, arcsec=ARCSEC, tile_loader=_cached_dem)
        except ValueError:  # next to ±180°
            continue
        if not np.isfinite(patch.z).all():
            continue
        gappy = with_track_gaps(patch, measured)
        degraded = extract_features(gappy, row["lat"], row["lon"], window_m, scales_m).values
        if not np.isfinite(list(degraded.values())).all():
            continue
        position = catalog.index.get_loc(index)
        distance = weighted_distance(normalizer.apply_values(degraded), z_all, weights)
        rank = int(np.sum(distance < distance[position])) + 1
        rows.append({"lat": row["lat"], "lon": row["lon"], "rank": rank, "percentile": 100.0 * (1 - rank / len(catalog)), "distance_to_self": float(distance[position])})
        change = {n_: 100.0 * abs(degraded[n_] - row[n_]) / max(abs(row[n_]), 1e-6) for n_ in other_names}
        if profile_names:
            before = np.array([row[n_] for n_ in profile_names], dtype="float64")
            after = np.array([degraded[n_] for n_ in profile_names])
            change["rprof"] = 100.0 * float(np.sqrt(np.mean((after - before) ** 2)) / max(np.sqrt(np.mean(before**2)), 1e-6))
        changes.append(change)

    changes_frame = pd.DataFrame(changes)
    summary = pd.DataFrame({"median_abs_change_pct": changes_frame.median(), "p90_abs_change_pct": changes_frame.quantile(0.9)})
    return pd.DataFrame(rows), summary


def support_bias(fine: dict[str, float], coarse: dict[str, float]) -> pd.DataFrame:
    """Relative change of each feature between two versions of the same window."""
    rows = [{"feature": name, "fine": fine[name], "coarse": coarse[name], "bias_pct": 100.0 * (coarse[name] - fine[name]) / max(abs(fine[name]), 1e-6)} for name in fine]
    return pd.DataFrame(rows)
