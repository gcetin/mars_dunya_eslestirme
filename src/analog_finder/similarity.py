"""Explainable ranking of Earth windows against a planetary reference."""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from .terrain import DEFAULT_WINDOW_M, PROFILE_BINS, PROFILE_PREFIX

FAMILIES = {"slope": "Eğim", "ruggedness": "Engebelilik", "relief": "Rölyef", "shape": "Biçim", "profile": "Radyal profil"}
DEFAULT_WEIGHTS = {"slope": 0.25, "ruggedness": 0.20, "relief": 0.15, "shape": 0.20, "profile": 0.20}

# Relative weight of a feature within its family, by analysis scale. Rebuilding
# Earth DEMs from real MOLA track patterns changes 1 km features by 8-17 %, but
# 2-4 km features by mostly 1-8 % (see notebooks/01), so 1 km counts half.
SCALE_WEIGHTS = {"1km": 0.5}

_LABELS = {
    "slope_mean": ("Ortalama eğim", "°"),
    "slope_std": ("Eğim std. sapması", "°"),
    "slope_p95": ("Eğim 95. yüzdelik", "°"),
    "tri": ("Arazi Engebelilik İndeksi", "m"),
    "rough": ("Yerel pürüzlülük", "m"),
    "relief": ("Yükselti farkı (P1–P99)", "m"),
    "elev_std": ("Eğimden arındırılmış yükseklik std.", "m"),
    "hypso": ("Hipsometrik integral", ""),
    "skew": ("Yükseklik çarpıklığı", ""),
    "radial_r2": ("Dairesel simetri R²", ""),
    "basin_depth": ("Havza derinliği (kenar − merkez)", "m"),
}


def _base(name: str) -> str:
    return name.rsplit("_", 1)[0] if name.endswith("km") else name


def is_profile(name: str) -> bool:
    return name.startswith(PROFILE_PREFIX)


def family_of(name: str) -> str:
    base = _base(name)
    if is_profile(base):
        return "profile"
    if base.startswith("slope"):
        return "slope"
    if base in ("tri", "rough"):
        return "ruggedness"
    if base in ("relief", "elev_std"):
        return "relief"
    return "shape"


def _km(value_m: float) -> str:
    return f"{value_m / 1000:.1f}".rstrip("0").rstrip(".").replace(".", ",")


def label_of(name: str) -> tuple[str, str]:
    """Human-readable name and unit of a feature.

    Profile rings are labelled for the catalog's 60 km window, whose disc has
    a radius of 29.5 km at 1 km spacing.
    """
    if name == "rprof":
        return "Radyal profil (tüm halkalar)", "σ"
    if is_profile(name):
        ring = int(name.removeprefix(PROFILE_PREFIX))
        width = (DEFAULT_WINDOW_M / 2 - 500.0) / PROFILE_BINS
        return f"Radyal profil, {_km(ring * width)}–{_km((ring + 1) * width)} km halkası", "σ"
    label, unit = _LABELS[_base(name)]
    if name.endswith("km"):
        label = f"{label} @ {name.rsplit('_', 1)[1].removesuffix('km')} km"
    return label, unit


def transform(name: str, values):
    """Map a feature to a space where differences mean comparable things."""
    base = _base(name)
    values = np.asarray(values, dtype="float64")
    if base == "basin_depth":
        return np.arcsinh(values / 100.0)
    if base in ("hypso", "skew", "radial_r2") or is_profile(base):
        return values
    return np.log1p(values)


@dataclass
class Normalizer:
    """Robust centring and scaling learned from the Earth window population.

    The rings of the radial profile share one scale, the median of their
    interquartile ranges, so the profile is compared as one curve in its own
    units rather than ring by ring.
    """

    names: list[str]
    center: np.ndarray
    scale: np.ndarray

    @classmethod
    def fit(cls, catalog: pd.DataFrame, names: list[str]) -> "Normalizer":
        matrix = np.column_stack([transform(n, catalog[n].to_numpy()) for n in names])
        q25, q50, q75 = np.nanpercentile(matrix, [25, 50, 75], axis=0)
        scale = q75 - q25
        rings = np.array([is_profile(n) for n in names])
        if rings.any():
            scale[rings] = np.median(scale[rings])
        scale[scale == 0] = 1.0
        return cls(names, q50, scale)

    def apply_frame(self, frame: pd.DataFrame) -> np.ndarray:
        matrix = np.column_stack([transform(n, frame[n].to_numpy()) for n in self.names])
        return (matrix - self.center) / self.scale

    def apply_values(self, values: dict[str, float]) -> np.ndarray:
        vector = np.array([float(transform(n, values[n])) for n in self.names])
        return (vector - self.center) / self.scale


def _scale_weight(name: str) -> float:
    return SCALE_WEIGHTS.get(name.rsplit("_", 1)[-1], 1.0)


def feature_weights(names: list[str], family_weights: dict[str, float]) -> np.ndarray:
    """Split each family's weight among its members; the total is 1."""
    families = [family_of(n) for n in names]
    relative = np.array([_scale_weight(n) for n in names])
    family_totals = {f: relative[[g == f for g in families]].sum() for f in set(families)}
    weights = np.array(
        [max(family_weights.get(f, 0.0), 0.0) * r / family_totals[f] for f, r in zip(families, relative)]
    )
    total = weights.sum()
    if total <= 0:
        raise ValueError("At least one family weight must be positive.")
    return weights / total


def weighted_distance(reference: np.ndarray, candidates: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Weighted RMS difference, in robust standard units of the Earth population."""
    return np.sqrt(np.sum(weights * (candidates - reference) ** 2, axis=-1))


def score_from_distance(distance):
    """0–100 display score. A convenience, not a probability of being an analog."""
    return 100.0 * np.exp(-np.asarray(distance))


@dataclass
class Filters:
    """Land-cover limits on the fraction of each window, in the range 0 to 1."""

    min_bare: float = 0.0
    max_tree: float = 1.0
    max_crop: float = 1.0
    max_built: float = 1.0
    max_water: float = 1.0
    max_snow: float = 1.0

    @classmethod
    def bare_ground(cls) -> "Filters":
        """Default for Mars and the Moon: exposed ground, little vegetation or water."""
        return cls(min_bare=0.5, max_tree=0.10, max_crop=0.15, max_built=0.03, max_water=0.15, max_snow=0.2)

    def mask(self, catalog: pd.DataFrame) -> np.ndarray:
        return (
            (catalog["lc_bare"] >= self.min_bare)
            & (catalog["lc_tree"] <= self.max_tree)
            & (catalog["lc_crop"] <= self.max_crop)
            & (catalog["lc_built"] <= self.max_built)
            & (catalog["lc_water"] <= self.max_water)
            & (catalog["lc_snow"] <= self.max_snow)
        ).to_numpy()


@dataclass
class Ranking:
    table: pd.DataFrame
    names: list[str]
    reference_z: np.ndarray
    normalizer: Normalizer
    weights: np.ndarray
    n_scanned: int
    n_after_filter: int
    sorted_distance: np.ndarray
    out_of_range: list[str] = field(default_factory=list)

    def locate(self, values: dict[str, float]) -> tuple[float, float]:
        """Distance and percentile of any feature vector against this reference."""
        z = self.normalizer.apply_values(values)
        distance = float(weighted_distance(self.reference_z, z, self.weights))
        farther = 1.0 - np.searchsorted(self.sorted_distance, distance, side="right") / len(self.sorted_distance)
        return distance, 100.0 * float(farther)


def great_circle_km(lat1, lon1, lat2, lon2, radius_km: float = 6371.0):
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dp, dl = p2 - p1, np.radians(lon2 - lon1)
    a = np.sin(dp / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    return 2 * radius_km * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


def rank_candidates(
    reference: dict[str, float],
    catalog: pd.DataFrame,
    names: list[str],
    *,
    family_weights: dict[str, float] | None = None,
    filters: Filters | None = None,
    top_n: int = 10,
    min_separation_km: float = 60.0,
    stability_draws: int = 300,
    seed: int = 0,
) -> Ranking:
    """Rank catalog windows by morphological distance to ``reference``.

    The result table carries, for each candidate:

    * ``distance``: weighted RMS difference, smaller is closer;
    * ``percentile``: share of all scanned windows that are farther away;
    * ``stability``: share of randomly perturbed weightings under which the
      candidate stays among the closest 0.1 % of eligible windows;
    * ``d_<family>``: the distance computed on that feature family alone.
    """
    family_weights = family_weights or DEFAULT_WEIGHTS
    normalizer = Normalizer.fit(catalog, names)
    z_all = normalizer.apply_frame(catalog)
    z_ref = normalizer.apply_values(reference)
    weights = feature_weights(names, family_weights)

    all_distance = weighted_distance(z_ref, z_all, weights)
    sorted_all = np.sort(all_distance)

    mask = filters.mask(catalog) if filters else np.ones(len(catalog), bool)
    eligible = np.nonzero(mask)[0]
    out_of_range = [n for n, v in zip(names, z_ref) if abs(v) > 3.0]
    empty = Ranking(catalog.iloc[0:0].copy(), names, z_ref, normalizer, weights, len(catalog), 0, sorted_all, out_of_range)
    if eligible.size == 0:
        return empty

    # Greedy non-maximum suppression so neighbouring windows do not fill the list.
    order = eligible[np.argsort(all_distance[eligible])]
    chosen: list[int] = []
    for index in order:
        row = catalog.iloc[index]
        if chosen:
            chosen_rows = catalog.iloc[chosen]
            separation = great_circle_km(row["lat"], row["lon"], chosen_rows["lat"].to_numpy(), chosen_rows["lon"].to_numpy())
            if separation.min() < min_separation_km:
                continue
        chosen.append(int(index))
        if len(chosen) == top_n:
            break

    table = catalog.iloc[chosen].copy()
    table["distance"] = all_distance[chosen]
    table["score"] = score_from_distance(all_distance[chosen])
    table["percentile"] = 100.0 * (1.0 - np.searchsorted(sorted_all, all_distance[chosen], side="right") / len(sorted_all))

    families = [family_of(n) for n in names]
    for family in FAMILIES:
        member = np.array([f == family for f in families])
        if member.any():
            table[f"d_{family}"] = np.sqrt(np.mean((z_all[chosen][:, member] - z_ref[member]) ** 2, axis=1))

    table["stability"] = _stability(z_all, eligible, chosen, z_ref, names, family_weights, stability_draws, seed)
    table.insert(0, "rank", np.arange(1, len(table) + 1))
    return Ranking(table.reset_index(drop=True), names, z_ref, normalizer, weights, len(catalog), int(eligible.size), sorted_all, out_of_range)


def _stability(z_all, eligible, chosen, z_ref, names, family_weights, draws: int, seed: int, max_windows: int = 20_000) -> np.ndarray:
    """How often each candidate remains near the top when the weights change.

    A random subset of eligible windows stands in for the full population, so
    the answer stays interactive on a global catalog.
    """
    rng = np.random.default_rng(seed)
    population = eligible
    if eligible.size > max_windows:
        population = np.union1d(rng.choice(eligible, max_windows, replace=False), chosen)
    z = z_all[population]
    chosen_positions = np.searchsorted(population, chosen)
    families = sorted(family_weights)
    alpha = np.array([max(family_weights[f], 1e-3) for f in families]) * 20.0
    top_k = max(int(round(z.shape[0] * 0.001)), 20)
    hits = np.zeros(len(chosen_positions))
    for draw in rng.dirichlet(alpha, draws):
        weights = feature_weights(names, dict(zip(families, draw)))
        distance = weighted_distance(z_ref, z, weights)
        threshold = np.partition(distance, min(top_k, distance.size - 1))[min(top_k, distance.size - 1)]
        hits += distance[chosen_positions] <= threshold
    return hits / draws


def explain(ranking: Ranking, reference: dict[str, float], candidate: pd.Series) -> pd.DataFrame:
    """Feature-by-feature comparison, largest disagreement first."""
    z_cand = ranking.normalizer.apply_values({n: candidate[n] for n in ranking.names})
    rows = []
    for i, name in enumerate(ranking.names):
        label, unit = label_of(name)
        difference = z_cand[i] - ranking.reference_z[i]
        rows.append(
            {
                "feature": name,
                "label": label,
                "family": FAMILIES[family_of(name)],
                "unit": unit,
                "reference": reference[name],
                "candidate": float(candidate[name]),
                "z_diff": float(difference),
                "contribution": float(ranking.weights[i] * difference**2),
            }
        )
    frame = pd.DataFrame(rows)
    frame["contribution"] /= max(frame["contribution"].sum(), 1e-12)
    return frame.sort_values("contribution", ascending=False).reset_index(drop=True)


def family_agreement(candidate: pd.Series) -> dict[str, float]:
    """Display score per family, derived from the family distances."""
    return {
        label: float(100.0 * math.exp(-float(candidate[f"d_{family}"])))
        for family, label in FAMILIES.items()
        if f"d_{family}" in candidate
    }
