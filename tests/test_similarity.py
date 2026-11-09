"""Ranking behaviour on a synthetic catalog where the right answer is known."""

import numpy as np
import pandas as pd
import pytest

from analog_finder.landcover import COLUMNS
from analog_finder.similarity import (
    DEFAULT_WEIGHTS,
    Filters,
    explain,
    feature_weights,
    rank_candidates,
)
from analog_finder.terrain import feature_names

NAMES = feature_names()


def make_catalog(n: int = 2000, seed: int = 1) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    data = {
        "lat": rng.uniform(-50, 50, n),
        "lon": rng.uniform(-170, 170, n),
    }
    for name in NAMES:
        if name in ("hypso", "radial_r2"):
            data[name] = rng.uniform(0, 1, n)
        elif name == "skew":
            data[name] = rng.normal(0, 1, n)
        elif name == "basin_depth":
            data[name] = rng.normal(0, 300, n)
        else:
            data[name] = rng.lognormal(1.5, 0.8, n)
    frame = pd.DataFrame(data)
    for column in COLUMNS:
        frame[column] = 0.0
    frame["lc_bare"] = 1.0
    return frame


def test_identical_window_ranks_first_with_zero_distance():
    catalog = make_catalog()
    reference = catalog.iloc[123][NAMES].to_dict()
    ranking = rank_candidates(reference, catalog, NAMES, top_n=5)
    assert ranking.table.iloc[0]["lat"] == pytest.approx(catalog.iloc[123]["lat"])
    assert ranking.table.iloc[0]["distance"] == pytest.approx(0.0)
    assert ranking.table.iloc[0]["score"] == pytest.approx(100.0)
    assert ranking.table.iloc[0]["percentile"] > 99.9
    assert list(ranking.table["distance"]) == sorted(ranking.table["distance"])


def test_absolute_units_do_not_dominate():
    """Equal family weights must hold even if families have many members."""
    weights = feature_weights(NAMES, {"slope": 1, "ruggedness": 1, "relief": 1, "shape": 1})
    families = {"slope": 0.0, "ruggedness": 0.0, "relief": 0.0, "shape": 0.0}
    from analog_finder.similarity import family_of

    for name, weight in zip(NAMES, weights):
        families[family_of(name)] += weight
    assert all(value == pytest.approx(0.25) for value in families.values())


def test_filter_removes_vegetated_windows():
    catalog = make_catalog()
    catalog.loc[123, ["lc_bare", "lc_tree"]] = [0.0, 1.0]
    reference = catalog.iloc[123][NAMES].to_dict()
    ranking = rank_candidates(reference, catalog, NAMES, filters=Filters.bare_ground(), top_n=5)
    assert not np.isclose(ranking.table["lat"], catalog.iloc[123]["lat"]).any()
    assert ranking.n_after_filter == len(catalog) - 1


def test_results_are_spatially_separated():
    catalog = make_catalog()
    twin = catalog.iloc[[7]].copy()
    twin["lat"] += 0.05
    catalog = pd.concat([catalog, twin], ignore_index=True)
    reference = catalog.iloc[7][NAMES].to_dict()
    ranking = rank_candidates(reference, catalog, NAMES, top_n=3, min_separation_km=60)
    assert not np.isclose(ranking.table.iloc[1]["lat"], catalog.iloc[-1]["lat"])


def test_empty_result_when_nothing_passes():
    catalog = make_catalog()
    catalog["lc_bare"] = 0.0
    ranking = rank_candidates(catalog.iloc[0][NAMES].to_dict(), catalog, NAMES, filters=Filters.bare_ground())
    assert ranking.table.empty
    assert ranking.n_after_filter == 0


def test_explanation_points_at_the_changed_feature():
    catalog = make_catalog()
    reference = catalog.iloc[50][NAMES].to_dict()
    candidate = catalog.iloc[50].copy()
    candidate["basin_depth"] += 3000.0
    ranking = rank_candidates(reference, catalog, NAMES, top_n=3)
    frame = explain(ranking, reference, candidate)
    assert frame.iloc[0]["feature"] == "basin_depth"
    assert frame["contribution"].sum() == pytest.approx(1.0)


def test_out_of_range_reference_is_flagged():
    catalog = make_catalog()
    reference = catalog.iloc[0][NAMES].to_dict()
    reference["relief"] = 1e7
    assert "relief" in rank_candidates(reference, catalog, NAMES, top_n=3).out_of_range


def test_stability_is_a_probability():
    catalog = make_catalog()
    ranking = rank_candidates(catalog.iloc[9][NAMES].to_dict(), catalog, NAMES, family_weights=DEFAULT_WEIGHTS, top_n=4)
    assert ranking.table["stability"].between(0, 1).all()
    assert ranking.table.iloc[0]["stability"] == pytest.approx(1.0)


def test_one_km_features_count_half_within_their_family():
    weights = dict(zip(NAMES, feature_weights(NAMES, DEFAULT_WEIGHTS)))
    assert weights["tri_1km"] == pytest.approx(weights["tri_2km"] / 2)
    assert sum(weights.values()) == pytest.approx(1.0)


def test_locate_reports_distance_for_arbitrary_window():
    catalog = make_catalog()
    reference = catalog.iloc[3][NAMES].to_dict()
    ranking = rank_candidates(reference, catalog, NAMES, top_n=3)
    assert ranking.locate(reference)[0] == pytest.approx(0.0)
    assert ranking.locate(reference)[1] > 99.9
    assert ranking.locate(catalog.iloc[10][NAMES].to_dict())[0] > 0
