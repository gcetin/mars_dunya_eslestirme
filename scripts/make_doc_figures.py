"""Draw the figures of docs/nasil-calisir.md from real data.

Each figure is written twice, for light and dark page themes:
docs/img/<name>.png and docs/img/<name>-dark.png. The worked-example numbers
quoted in the guide are printed at the end.

    python scripts/make_doc_figures.py
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from matplotlib.patches import Circle, Rectangle
from matplotlib.ticker import FuncFormatter, ScalarFormatter

from analog_finder.bodies import EARTH, MARS
from analog_finder.catalog import CATALOG_DIR, load_catalog
from analog_finder.dem import read_earth_patch, read_global_patch
from analog_finder.mola_support import shot_counts, shot_coverage
from analog_finder.query import window_features
from analog_finder.similarity import DEFAULT_WEIGHTS, Filters, explain, rank_candidates
from analog_finder.surface import surface_for_window
from analog_finder.terrain import (
    DEFAULT_SCALES_M,
    SIGMA_PER_SCALE,
    extract_features,
    hillshade,
    hypsometric_curve,
    radial_profile,
    read_half_m,
    slope_cdf,
    slope_deg,
)

OUT = ROOT / "docs" / "img"
WINDOW = 60_000.0
HALF = read_half_m(WINDOW, DEFAULT_SCALES_M)
JEZERO = (18.38, 77.58)
SALDA = (37.55, 29.68)
TOOTING = (23.17, -152.2)

# Reference palette, light and dark steps of the same hues.
RAMP = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7", "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]
THEMES = {
    "light": {
        "surface": "#fcfcfb", "ink": "#0b0b0b", "ink2": "#52514e", "muted": "#898781",
        "grid": "#e1e0d9", "axis": "#c3c2b7", "earth": "#2a78d6", "reference": "#eb6834",
        "neutral": "#c3c2b7", "off": "#ecebe6", "old": "#898781", "ramp": RAMP,
    },
    "dark": {
        "surface": "#1a1a19", "ink": "#f0efec", "ink2": "#c3c2b7", "muted": "#898781",
        "grid": "#2c2c2a", "axis": "#383835", "earth": "#3987e5", "reference": "#d95926",
        "neutral": "#4a4a46", "off": "#262624", "old": "#898781", "ramp": RAMP[::-1],
    },
}


def num(value: float, decimals: int = 1) -> str:
    """Number with a Turkish decimal comma."""
    return f"{value:.{decimals}f}".replace(".", ",")


def count(value: int) -> str:
    """Integer with a Turkish thousands separator."""
    return f"{value:,}".replace(",", ".")


TICKS = FuncFormatter(lambda v, _: (count(int(round(v))) if abs(v) >= 1000 else f"{v:g}".replace(".", ",")).replace("-", "−"))
# One contrast stretch for every shaded-relief panel, so panels stay comparable.
SHADE = {"cmap": "gray", "vmin": 0.4, "vmax": 0.95}


def style(t: dict) -> None:
    plt.rcParams.update({
        "font.family": "sans-serif", "font.size": 10, "figure.dpi": 120, "savefig.dpi": 120,
        "figure.facecolor": t["surface"], "axes.facecolor": t["surface"], "savefig.facecolor": t["surface"],
        "text.color": t["ink"], "axes.labelcolor": t["ink2"], "axes.edgecolor": t["axis"],
        "xtick.color": t["muted"], "ytick.color": t["muted"], "axes.titlecolor": t["ink"],
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.color": t["grid"], "grid.linewidth": 0.6, "axes.axisbelow": True,
        "axes.titlesize": 10.5, "axes.titleweight": "bold", "axes.titlelocation": "left",
        "lines.linewidth": 2, "legend.frameon": False, "legend.labelcolor": t["ink2"],
    })


def save(fig, name: str, theme: str) -> None:
    # Only numeric axes; a category axis keeps its labels.
    for ax in fig.axes:
        for axis in (ax.xaxis, ax.yaxis):
            if isinstance(axis.get_major_formatter(), ScalarFormatter):
                axis.set_major_formatter(TICKS)
    suffix = "" if theme == "light" else "-dark"
    fig.savefig(OUT / f"{name}{suffix}.png", bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)


def bare(ax) -> None:
    ax.grid(False)
    ax.set_xticks([]); ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def km_axes(ax, half_km: float = 30) -> None:
    ax.grid(False)
    ax.set_xticks([-half_km, 0, half_km]); ax.set_yticks([-half_km, 0, half_km])
    ax.set_xlabel("Doğu (km)"); ax.set_ylabel("Kuzey (km)")


def end_label(ax, x, y, text, t, dx=0.0, dy=0.0) -> None:
    ax.annotate(text, (x[-1], y[-1]), xytext=(4 + dx, dy), textcoords="offset points", color=t["ink2"], va="center", fontsize=9)


def raw_slope(patch, lat0: float) -> np.ndarray:
    """Slope on the source pixels over the central 60 km, as a naive comparison would do."""
    dy = math.radians(patch.dlat) * patch.body.radius_m
    dx = math.radians(patch.dlon) * patch.body.radius_m * math.cos(math.radians(lat0))
    ny, nx = int(WINDOW / 2 / dy), int(WINDOW / 2 / dx)
    cy, cx = patch.z.shape[0] // 2, patch.z.shape[1] // 2
    z = patch.z[cy - ny : cy + ny, cx - nx : cx + nx]
    gy, gx = np.gradient(z, dy, dx)
    return np.degrees(np.arctan(np.hypot(gx, gy)))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    catalog, meta = load_catalog(CATALOG_DIR / "earth_global.parquet")
    names = meta["feature_names"]

    jezero_patch = read_global_patch(MARS, *JEZERO, HALF)
    salda_patch = read_earth_patch(*SALDA, HALF, arcsec=1)
    jezero = extract_features(jezero_patch, *JEZERO)
    salda = extract_features(salda_patch, *SALDA)
    salda_15, _ = window_features(EARTH, *SALDA)

    ranked = rank_candidates(jezero.values, catalog, names, filters=Filters.bare_ground(), top_n=10)
    best = ranked.table.iloc[0]
    why = explain(ranked, jezero.values, best)
    salda_distance, salda_percentile = ranked.locate(salda_15.values)
    top_windows = [window_features(EARTH, float(r.lat), float(r.lon))[0] for r in ranked.table.head(5).itertuples()]
    counts, _ = shot_counts(*JEZERO, 2 * HALF)
    coverage = shot_coverage(*JEZERO, WINDOW)
    surface = surface_for_window(float(best["lat"]), float(best["lon"]), WINDOW)

    # Jezero seen from a window 10 km east of the crater centre.
    east_deg = 10_000 / (MARS.radius_m * math.radians(1) * math.cos(math.radians(JEZERO[0])))
    shifted_at = (JEZERO[0], JEZERO[1] + east_deg)
    shifted = extract_features(read_global_patch(MARS, *shifted_at, HALF), *shifted_at)
    fixed_r, fixed_p, fixed_r2 = radial_profile(shifted.grids[1000.0], 1000.0)
    # Fixed-centre profile in the same units as the searched one.
    fixed_p = fixed_p / np.std(shifted.grids[1000.0] - np.mean(shifted.grids[1000.0]))
    tooting, _ = window_features(MARS, *TOOTING)
    only_profile = dict.fromkeys(DEFAULT_WEIGHTS, 0.0) | {"profile": 1.0}
    crater_match = rank_candidates(tooting.values, catalog, names, family_weights=only_profile, top_n=3, stability_draws=50).table.iloc[0]
    crater_feats, _ = window_features(EARTH, float(crater_match["lat"]), float(crater_match["lon"]))

    raw_j, raw_s = raw_slope(jezero_patch, JEZERO[0]), raw_slope(salda_patch, SALDA[0])
    slope_j, slope_s = slope_deg(jezero.grids[1000.0], 1000.0)[1:-1, 1:-1], slope_deg(salda.grids[1000.0], 1000.0)[1:-1, 1:-1]
    extent = [-30, 30, -30, 30]

    for theme, t in THEMES.items():
        style(t)
        seq = LinearSegmentedColormap.from_list("seq", t["ramp"])

        # 1. Jezero and Salda at the same physical scale.
        fig, axes = plt.subplots(1, 2, figsize=(9, 4.6))
        for ax, feats, title in zip(axes, (jezero, salda), ("Jezero Krateri (Mars) · MGS MOLA", "Salda Gölü (Dünya) · NASADEM")):
            ax.imshow(hillshade(feats.grids[1000.0], 1000.0), **SHADE, extent=extent)
            ax.set_title(title); km_axes(ax)
        fig.tight_layout(); save(fig, "jezero-salda", theme)

        # 2. Support matching: raw pixels versus a common 1 km support.
        fig, axes = plt.subplots(1, 3, figsize=(12, 4.3), constrained_layout=True)
        panels = (
            (raw_s, f"Salda · ham 30 m piksel\nortalama eğim {num(raw_s.mean())}°"),
            (slope_s, f"Salda · 1 km destek\nortalama eğim {num(salda.values['slope_mean_1km'])}°"),
            (slope_j, f"Jezero · 1 km destek\nortalama eğim {num(jezero.values['slope_mean_1km'])}°"),
        )
        for ax, (grid, title) in zip(axes, panels):
            image = ax.imshow(grid, cmap=seq, extent=extent, vmin=0, vmax=20, interpolation="nearest" if grid.shape[0] < 100 else "antialiased")
            ax.set_title(title); km_axes(ax)
        bar = fig.colorbar(image, ax=axes, shrink=0.85, pad=0.02, extend="max")
        bar.set_label("Eğim (°)"); bar.outline.set_visible(False)
        save(fig, "destek-esleme", theme)

        # 3. The three analysis scales.
        fig, axes = plt.subplots(1, 3, figsize=(11, 4))
        for ax, scale in zip(axes, DEFAULT_SCALES_M):
            grid = jezero.grids[scale]
            ax.imshow(hillshade(grid, scale), **SHADE, extent=extent, interpolation="nearest")
            ax.set_title(f"{scale / 1000:.0f} km · {grid.shape[0]} × {grid.shape[1]} hücre"); km_axes(ax)
        fig.tight_layout(); save(fig, "olcekler", theme)

        # 4. Feature curves for Jezero and Salda.
        fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.8))
        g_j, g_s = jezero.grids[1000.0], salda.grids[1000.0]
        curves = (
            (slope_cdf(g_j, 1000.0), slope_cdf(g_s, 1000.0), "Eğim (°)", "Bu eğimde veya altında kalan alan", "Eğim dağılımı"),
            (hypsometric_curve(g_j), hypsometric_curve(g_s), "Daha yüksekte kalan alan payı", "Bağıl yükseklik", "Hipsometrik eğri"),
        )
        for ax, (a, b, xlabel, ylabel, title) in zip(axes[:2], curves):
            ax.plot(*a, color=t["reference"], label="Jezero (Mars)")
            ax.plot(*b, color=t["earth"], label="Salda (Dünya)")
            ax.set_xlabel(xlabel); ax.set_ylabel(ylabel); ax.set_title(title)
        pj, ps = jezero.profile, salda.profile
        ax = axes[2]
        ax.plot(pj.radius_m / 1000, pj.shape, color=t["reference"], label="Jezero (Mars)")
        ax.plot(ps.radius_m / 1000, ps.shape, color=t["earth"], label="Salda (Dünya)")
        end_label(ax, pj.radius_m / 1000, pj.shape, f"Jezero\nR² {num(pj.r2, 2)}", t)
        end_label(ax, ps.radius_m / 1000, ps.shape, f"Salda\nR² {num(ps.r2, 2)}", t)
        ax.set_xlim(0, 36)
        ax.set_xlabel("Profil merkezinden uzaklık (km)"); ax.set_ylabel("Halka ortalaması (σ)"); ax.set_title("Radyal profil")
        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 1.06))
        fig.tight_layout(); save(fig, "ozellikler", theme)

        # 4b. The searched centre and what it buys.
        fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.1), gridspec_kw={"width_ratios": [1, 1.25, 1.25]})
        ax = axes[0]
        ax.imshow(hillshade(shifted.grids[1000.0], 1000.0), **SHADE, extent=extent)
        found = shifted.profile
        ax.plot(0, 0, marker="+", markersize=14, markeredgewidth=2.5, color=t["ink"])
        ax.add_patch(Circle((found.dx_m / 1000, found.dy_m / 1000), 29.5, fill=False, edgecolor=t["reference"], linewidth=2))
        ax.plot(found.dx_m / 1000, found.dy_m / 1000, marker="x", markersize=9, markeredgewidth=2, color=t["reference"])
        ax.set_xlim(-30, 30); ax.set_ylim(-30, 30)
        ax.set_title("Pencere kraterden 10 km doğuda\n+ pencere merkezi · × profil merkezi"); km_axes(ax)
        ax = axes[1]
        ax.plot(jezero.profile.radius_m / 1000, jezero.profile.shape, color=t["reference"], label="Krater merkezli pencere")
        ax.plot(found.radius_m / 1000, found.shape, color=t["earth"], linestyle=(0, (4, 2)), label="Kaymış pencere, merkez aranarak")
        ax.plot(fixed_r / 1000, fixed_p, color=t["old"], label="Kaymış pencere, sabit merkez")
        ax.set_xlabel("Profil merkezinden uzaklık (km)"); ax.set_ylabel("Halka ortalaması (σ)")
        ax.set_title(f"Jezero'nun profili (R²: {num(jezero.profile.r2, 2)} · {num(found.r2, 2)} · {num(fixed_r2, 2)})")
        ax.legend(loc="upper left", fontsize=8.5)
        ax = axes[2]
        ax.plot(tooting.profile.radius_m / 1000, tooting.profile.shape, color=t["reference"], label="Tooting (Mars, 29 km)")
        ax.plot(crater_feats.profile.radius_m / 1000, crater_feats.profile.shape, color=t["earth"], label="En yakın Dünya profili: Mistastin")
        ax.set_xlabel("Profil merkezinden uzaklık (km)"); ax.set_ylabel("Halka ortalaması (σ)")
        ax.set_title("Kenarı yükseltilmiş krater")
        ax.legend(loc="lower right", fontsize=8.5)
        fig.tight_layout(); save(fig, "radyal-profil", theme)

        # 5. Where the catalog windows are.
        fig, ax = plt.subplots(figsize=(12, 5.6))
        mask = Filters.bare_ground().mask(catalog)
        ax.scatter(catalog["lon"], catalog["lat"], s=0.15, color=t["neutral"], linewidths=0, rasterized=True, label=f"Tüm pencereler ({count(len(catalog))})")
        ax.scatter(catalog["lon"][mask], catalog["lat"][mask], s=0.15, color=t["earth"], linewidths=0, rasterized=True, label=f"Çıplak zemin filtresinden geçen ({count(int(mask.sum()))})")
        top = ranked.table
        ax.scatter(top["lon"], top["lat"], s=70, facecolor="none", edgecolor=t["reference"], linewidths=2, label="Jezero'nun ilk 10 adayı", zorder=3)
        ax.set_xlim(-180, 180); ax.set_ylim(-60, 86); ax.set_aspect("equal")
        ax.set_xlabel("Boylam (°)"); ax.set_ylabel("Enlem (°)")
        ax.set_title("Dünya kataloğu: 60 km'lik pencerelerin merkezleri (20 km arayla)")
        legend = ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=3, handletextpad=0.4)
        for handle in legend.legend_handles[:2]:
            handle.set_sizes([30])
        fig.tight_layout(); save(fig, "katalog", theme)

        # 6. Distance of every Earth window to Jezero.
        fig, ax = plt.subplots(figsize=(10, 3.6))
        upper = np.percentile(ranked.sorted_distance, 99.9)
        ax.hist(ranked.sorted_distance[ranked.sorted_distance <= upper], bins=90, color=t["earth"], rwidth=0.85)
        ax.set_xlim(0, upper)
        ymax = ax.get_ylim()[1]
        for value, text in ((float(best["distance"]), f"En yakın aday\n{num(float(best['distance']), 2)}"), (salda_distance, f"Salda\n{num(salda_distance, 2)}")):
            ax.axvline(value, color=t["reference"], linewidth=2)
            ax.text(value + upper * 0.01, ymax * 0.93, text, color=t["ink"], va="top", fontsize=9)
        ax.set_xlabel("Jezero'ya morfolojik mesafe (küçük = daha benzer)"); ax.set_ylabel("Pencere sayısı")
        ax.set_title("Bütün Dünya pencerelerinin Jezero'ya uzaklığı (en uzak %0,1 gösterilmiyor)")
        fig.tight_layout(); save(fig, "mesafe-dagilimi", theme)

        # 7. Jezero and its five closest Earth windows.
        fig, axes = plt.subplots(2, 3, figsize=(10.5, 7.4))
        axes = axes.ravel()
        axes[0].imshow(hillshade(jezero.grids[1000.0], 1000.0), **SHADE, extent=extent)
        axes[0].set_title("Referans: Jezero (Mars)\nMGS MOLA")
        for ax, row, feats in zip(axes[1:], ranked.table.head(5).itertuples(), top_windows):
            ax.imshow(hillshade(feats.grids[1000.0], 1000.0), **SHADE, extent=extent)
            ns = "K" if row.lat >= 0 else "G"
            ew = "D" if row.lon >= 0 else "B"
            ax.set_title(f"#{row.rank} · {abs(row.lat):.2f}°{ns} {abs(row.lon):.2f}°{ew}\nbenzerlik {row.score:.0f} · kararlılık %{row.stability * 100:.0f}".replace(".", ",", 2))
        for ax in axes:
            bare(ax)
        fig.tight_layout(); save(fig, "adaylar", theme)

        # 8. Why the closest candidate was chosen.
        head = why.head(8).iloc[::-1]
        fig, ax = plt.subplots(figsize=(9, 3.9))
        ax.barh(head["label"], head["contribution"] * 100, color=t["earth"], height=0.62)
        for y, value in enumerate(head["contribution"] * 100):
            ax.text(value + 0.6, y, f"%{value:.0f}", va="center", color=t["ink2"], fontsize=9)
        ax.set_xlabel("Toplam farka katkı (%)"); ax.grid(axis="y", visible=False)
        ax.set_title("En yakın adayla Jezero arasındaki farkın kaynakları")
        ax.tick_params(axis="y", colors=t["ink2"], length=0)
        fig.tight_layout(); save(fig, "katki", theme)

        # 9. Where MOLA actually measured.
        fig, ax = plt.subplots(figsize=(5.4, 5.2))
        half_km = HALF / 1000
        ax.imshow(counts > 0, cmap=ListedColormap([t["off"], t["earth"]]), vmin=0, vmax=1, extent=[-half_km, half_km, -half_km, half_km], interpolation="nearest")
        ax.add_patch(Rectangle((-30, -30), 60, 60, fill=False, edgecolor=t["reference"], linewidth=2))
        ax.text(-29, 32, "60 km pencere", color=t["ink"], fontsize=9, bbox={"facecolor": t["surface"], "edgecolor": "none", "pad": 1.5})
        km_axes(ax, 30)
        ax.set_title(f"Jezero çevresinde lazer ölçümü olan\nMOLA hücreleri (pencerede %{coverage.measured_frac * 100:.1f})".replace(".", ","))
        fig.tight_layout(); save(fig, "mola-izleri", theme)

        # 10. Sentinel-2 view of the closest candidate.
        if surface is not None:
            summary, rgb, _ = surface
            fig, ax = plt.subplots(figsize=(5.4, 5.6))
            ax.imshow(rgb, extent=extent)
            km_axes(ax)
            ax.set_title(f"En yakın aday · Sentinel-2 L2A, {summary.date}\nmedyan NDVI {summary.ndvi_median:.2f} · bulut %{summary.cloud_pct:.1f}".replace(".", ",", 2))
            fig.tight_layout(); save(fig, "sentinel2", theme)

    # Numbers for the worked examples in the guide.
    print("== Destek eşleme örneği (1 km)")
    for label, d in (("MOLA 463 m", MARS.dem.native_m), ("NASADEM 30 m", 30.0)):
        target = SIGMA_PER_SCALE * 1000.0
        print(f"  {label}: hedef σ {target:.0f} m, piksel varyansı d²/12 = {d * d / 12:,.0f} m², uygulanan σ {math.sqrt(target**2 - d * d / 12):.1f} m")
    print(f"  Ham ortalama eğim: Jezero {raw_j.mean():.2f}°, Salda {raw_s.mean():.2f}° -> oran {raw_s.mean() / raw_j.mean():.2f}")
    print(f"  1 km eğim: Jezero {jezero.values['slope_mean_1km']:.2f}°, Salda {salda.values['slope_mean_1km']:.2f}° -> oran {salda.values['slope_mean_1km'] / jezero.values['slope_mean_1km']:.2f}")
    print("== Jezero ve Salda özellikleri")
    for name in ("slope_mean_1km", "relief", "hypso", "skew", "radial_r2", "basin_depth"):
        print(f"  {name}: Jezero {jezero.values[name]:.3f}, Salda {salda.values[name]:.3f}")
    print("== Sıralama örneği")
    print(f"  en yakın aday {best['lat']:.3f}, {best['lon']:.3f}: mesafe {best['distance']:.3f}, benzerlik {best['score']:.1f}, yüzdelik {best['percentile']:.4f}, kararlılık {best['stability']:.2f}")
    print(f"  aile uyumu: " + ", ".join(f"{f} {100 * math.exp(-best[f'd_{f}']):.0f}" for f in ('slope', 'ruggedness', 'relief', 'shape', 'profile')))
    print(f"  Salda: mesafe {salda_distance:.3f}, yüzdelik {salda_percentile:.1f}")
    print(f"  filtreden geçen {ranked.n_after_filter} / {ranked.n_scanned}")
    weights = dict(zip(ranked.names, ranked.weights))
    for _, row in why.head(4).iterrows():
        print(f"  {row['feature']}: Jezero {row['reference']:.3f}, aday {row['candidate']:.3f}, z farkı {row['z_diff']:+.2f}, ağırlık {weights[row['feature']]:.4f}, katkı %{row['contribution'] * 100:.1f}")
    print(f"  MOLA ölçümlü hücre payı: %{coverage.measured_frac * 100:.1f}")
    print("== Radyal profil")
    print(f"  Jezero: R² {jezero.profile.r2:.2f}, merkez kayması ({jezero.profile.dx_m / 1000:+.0f}, {jezero.profile.dy_m / 1000:+.0f}) km")
    print(f"  10 km kaymış pencere: bulunan kayma ({shifted.profile.dx_m / 1000:+.0f}, {shifted.profile.dy_m / 1000:+.0f}) km, R² {shifted.profile.r2:.2f}; sabit merkezde R² {fixed_r2:.2f}")
    print(f"  Tooting R² {tooting.profile.r2:.2f}; yalnızca profille en yakın: {crater_match['lat']:.2f}, {crater_match['lon']:.2f} (R² {crater_match['prof_r2']:.2f})")
    print(f"  Jezero profil halkaları: {np.round(jezero.profile.shape, 2)}")


if __name__ == "__main__":
    main()
