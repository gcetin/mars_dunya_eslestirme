"""Planetary Analog Finder dashboard.

Run with:  streamlit run app/streamlit_app.py
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import folium
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st
from streamlit_folium import st_folium

from analog_finder.bodies import BODIES, EARTH
from analog_finder.catalog import CATALOG_DIR, load_catalog
from analog_finder.mola_support import shot_coverage
from analog_finder.provenance import METHOD
from analog_finder.query import window_features
from analog_finder.similarity import (
    DEFAULT_WEIGHTS,
    FAMILIES,
    Filters,
    explain,
    family_agreement,
    label_of,
    rank_candidates,
)
from analog_finder.sites import KNOWN_ANALOGS, REFERENCE_SITES
from analog_finder.surface import PRODUCT as S2_PRODUCT, surface_for_window
from analog_finder.terrain import PROFILE_SHIFT_M, hillshade, hypsometric_curve, radial_profile, slope_cdf

st.set_page_config(page_title="Planetary Analog Finder", page_icon="🪐", layout="wide")

SCALE = 1000.0
PALETTE = {
    "light": {"earth": "#2a78d6", "reference": "#eb6834", "grid": "#e1e0d9", "muted": "#898781"},
    "dark": {"earth": "#3987e5", "reference": "#d95926", "grid": "#2c2c2a", "muted": "#898781"},
}
COLORS = PALETTE["dark" if st.context.theme.type == "dark" else "light"]

# Context imagery only; these layers never enter the similarity calculation.
# Each entry: label, service path, file extension, deepest zoom level offered.
TREK_TILES = "https://trek.nasa.gov/tiles/{path}/1.0.0/default/default028mm/{{z}}/{{y}}/{{x}}.{extension}"
PLANET_LAYERS = {
    "mars": (
        ("Viking renkli mozaik", "Mars/EQ/Mars_Viking_MDIM21_ClrMosaic_global_232m", "jpg", 7),
        ("MOLA renkli rölyef", "Mars/EQ/Mars_MGS_MOLA_ClrShade_merge_global_463m", "jpg", 7),
    ),
    "moon": (
        ("LROC WAC mozaik", "Moon/EQ/LRO_WAC_Mosaic_Global_303ppd_v02", "jpg", 8),
        ("LOLA renkli rölyef", "Moon/EQ/LRO_LOLA_ClrShade_Global_128ppd_v04", "png", 5),
    ),
}


@st.cache_resource(show_spinner="Dünya kataloğu yükleniyor…")
def get_catalog():
    # Prefer the global catalog; otherwise use any regional one that was built.
    for path in [CATALOG_DIR / "earth_global.parquet", *sorted(CATALOG_DIR.glob("*.parquet"))]:
        if path.exists():
            return *load_catalog(path), path.name
    return None, None, None


@st.cache_resource(show_spinner="Referans bölgenin DEM verisi okunuyor…")
def get_reference(body_key: str, lat: float, lon: float, window_m: float):
    features, _ = window_features(BODIES[body_key], lat, lon, window_m)
    return features


@st.cache_data(show_spinner=False)
def get_shot_coverage(lat: float, lon: float, window_m: float):
    try:
        return shot_coverage(lat, lon, window_m).measured_frac
    except Exception:
        return None


@st.cache_resource(show_spinner="Aday pencereler sıralanıyor…", max_entries=32)
def get_ranking(_catalog, _names, reference_key, _reference, weights, filter_values, top_n):
    return rank_candidates(
        _reference,
        _catalog,
        _names,
        family_weights=dict(weights),
        filters=Filters(*filter_values),
        top_n=top_n,
    )


@st.cache_resource(show_spinner=False, max_entries=128)
def get_earth_window(lat: float, lon: float, window_m: float):
    return window_features(EARTH, lat, lon, window_m)[0]


@st.cache_data(show_spinner="Bilinen analoglar değerlendiriliyor…", max_entries=32)
def known_analog_table(_ranking, reference_key, weights_key, window_m: float) -> pd.DataFrame:
    rows = []
    for site in KNOWN_ANALOGS:
        try:
            features = get_earth_window(site.lat, site.lon, window_m)
        except Exception:
            features = None
        if features is None or not np.isfinite(list(features.values.values())).all():
            rows.append({"Alan": site.name, "Bilinen analog türü": site.note, "Mesafe": None, "Yüzdelik": None})
            continue
        distance, percentile = _ranking.locate(features.values)
        rows.append({"Alan": site.name, "Bilinen analog türü": site.note, "Mesafe": distance, "Yüzdelik": percentile})
    return pd.DataFrame(rows).sort_values("Mesafe", na_position="last").reset_index(drop=True)


@st.cache_data(show_spinner="Sentinel-2 sahnesi okunuyor…")
def get_surface(lat: float, lon: float, window_m: float):
    return surface_for_window(lat, lon, window_m)


def count(value: int) -> str:
    """Integer with Turkish thousands separator."""
    return f"{value:,}".replace(",", ".")


def relief_image(grid: np.ndarray) -> np.ndarray:
    return np.kron(hillshade(grid, SCALE), np.ones((7, 7)))


def base_layout(height: int, xtitle: str, ytitle: str) -> dict:
    axis = {"gridcolor": COLORS["grid"], "zeroline": False, "linecolor": COLORS["grid"]}
    return {
        "height": height,
        "margin": {"l": 10, "r": 10, "t": 10, "b": 10},
        "paper_bgcolor": "rgba(0,0,0,0)",
        "plot_bgcolor": "rgba(0,0,0,0)",
        "xaxis": {**axis, "title": xtitle},
        "yaxis": {**axis, "title": ytitle},
        "legend": {"orientation": "h", "y": 1.12, "x": 0},
        "hovermode": "x unified",
    }


def two_line_chart(x_ref, y_ref, x_cand, y_cand, xtitle, ytitle, ref_name, hover_x):
    fig = go.Figure()
    fig.add_scatter(x=x_ref, y=y_ref, name=ref_name, mode="lines", line={"color": COLORS["reference"], "width": 2}, hovertemplate=hover_x)
    fig.add_scatter(x=x_cand, y=y_cand, name="Dünya adayı", mode="lines", line={"color": COLORS["earth"], "width": 2}, hovertemplate=hover_x)
    fig.update_layout(**base_layout(260, xtitle, ytitle))
    return fig


def window_bounds(radius_m: float, lat: float, lon: float, window_m: float) -> list[list[float]]:
    """South-west and north-east corners of a square window, in degrees."""
    half_lat = float(np.degrees(window_m / 2 / radius_m))
    half_lon = half_lat / float(np.cos(np.radians(lat)))
    return [[lat - half_lat, lon - half_lon], [lat + half_lat, lon + half_lon]]


def planet_map(body, lat: float, lon: float, window_m: float) -> folium.Map:
    """Zoomable basemap of Mars or the Moon with the analysis window outlined."""
    bounds = window_bounds(body.radius_m, lat, lon, window_m)
    window_deg = bounds[1][0] - bounds[0][0]
    # The tile grid has 512 * 2^z / 360 pixels per degree. Start with the window
    # about 90 px tall, so its surroundings fill the rest of the map.
    zoom = int(np.clip(round(np.log2(90 * 360 / (512 * window_deg))), 1, 7))
    # Leaflet's scale bar assumes Earth's size, so it stays off here.
    view = folium.Map(location=[lat, lon], zoom_start=zoom, min_zoom=1, max_zoom=10, tiles=None, crs="EPSG4326", control_scale=False)
    for index, (label, path, extension, deepest) in enumerate(PLANET_LAYERS[body.key]):
        folium.TileLayer(
            tiles=TREK_TILES.format(path=path, extension=extension),
            attr="NASA Solar System Treks",
            name=label,
            max_native_zoom=deepest,
            max_zoom=10,
            no_wrap=True,
            show=index == 0,
        ).add_to(view)
    folium.Rectangle(
        bounds,
        color=COLORS["reference"],
        weight=2,
        fill=False,
        tooltip=f"Analiz penceresi: {window_m / 1000:.0f} km × {window_m / 1000:.0f} km",
    ).add_to(view)
    folium.LayerControl(collapsed=False).add_to(view)
    return view


def world_map(table: pd.DataFrame, selected: int, window_m: float) -> folium.Map:
    earth_map = folium.Map(location=[20, 10], zoom_start=2, tiles="OpenStreetMap", control_scale=True)
    folium.TileLayer("Esri.WorldImagery", name="Uydu görüntüsü").add_to(earth_map)
    for row in table.itertuples():
        is_selected = int(row.rank) == selected
        color = COLORS["reference"] if is_selected else COLORS["earth"]
        folium.Rectangle(
            window_bounds(EARTH.radius_m, row.lat, row.lon, window_m),
            color=color,
            weight=2,
            fill=True,
            fill_opacity=0.25 if is_selected else 0.12,
            tooltip=f"#{int(row.rank)} · benzerlik {row.score:.0f}/100 · {row.lat:.2f}, {row.lon:.2f}",
        ).add_to(earth_map)
        folium.Marker(
            [row.lat, row.lon],
            icon=folium.DivIcon(
                html=(
                    f'<div style="background:{color};color:#fff;border-radius:50%;width:24px;height:24px;'
                    f'line-height:24px;text-align:center;font:600 12px system-ui">{int(row.rank)}</div>'
                ),
                icon_size=(24, 24),
                icon_anchor=(12, 12),
            ),
            tooltip=f"#{int(row.rank)} · {row.lat:.2f}, {row.lon:.2f}",
        ).add_to(earth_map)
    folium.LayerControl(collapsed=True).add_to(earth_map)
    if not table.empty:
        earth_map.fit_bounds([[table.lat.min() - 2, table.lon.min() - 2], [table.lat.max() + 2, table.lon.max() + 2]])
    return earth_map


catalog, meta, catalog_name = get_catalog()
if catalog is None:
    st.title("Planetary Analog Finder")
    st.error("Dünya kataloğu bulunamadı. Önce `python scripts/build_catalog.py --global --output data/catalog/earth_global.parquet` çalıştırın.")
    st.stop()

names = meta["feature_names"]
window_m = float(meta["window_m"])

with st.sidebar:
    st.header("Referans bölge")
    body_label = st.radio("Gezegen", ["Mars", "Ay"], horizontal=True)
    body = BODIES["mars" if body_label == "Mars" else "moon"]
    sites = REFERENCE_SITES[body.key]
    choice = st.selectbox("Alan", [s.name for s in sites] + ["Özel koordinat"])
    if choice == "Özel koordinat":
        lat = st.number_input("Enlem (°K)", -80.0, 80.0, 18.38, 0.1, format="%.2f")
        lon = st.number_input("Boylam (°D, −180…180)", -180.0, 180.0, 77.58, 0.1, format="%.2f")
        note = ""
    else:
        site = next(s for s in sites if s.name == choice)
        lat, lon, note = site.lat, site.lon, site.note

    st.header("Eşleşme ağırlıkları")
    st.caption("Görevin önemsediği özelliği öne çıkarın.")
    weights = {
        family: st.slider(label, 0, 100, int(DEFAULT_WEIGHTS[family] * 100), 5, key=f"w_{family}") / 100
        for family, label in FAMILIES.items()
    }
    if sum(weights.values()) == 0:
        st.warning("En az bir ağırlık sıfırdan büyük olmalı. Varsayılanlar kullanılıyor.")
        weights = DEFAULT_WEIGHTS

    st.header("Dünya yüzeyi filtresi")
    use_filter = st.toggle("Çıplak zemin ara (WorldCover)", value=True)
    defaults = Filters.bare_ground()
    if use_filter:
        min_bare = st.slider("En az çıplak zemin", 0.0, 1.0, defaults.min_bare, 0.05)
        max_tree = st.slider("En fazla ağaç örtüsü", 0.0, 1.0, defaults.max_tree, 0.05)
        max_crop = st.slider("En fazla tarım alanı", 0.0, 1.0, defaults.max_crop, 0.05)
        max_water = st.slider("En fazla su", 0.0, 1.0, defaults.max_water, 0.05)
        filter_values = (min_bare, max_tree, max_crop, defaults.max_built, max_water, defaults.max_snow)
    else:
        filter_values = (0.0, 1.0, 1.0, 1.0, 1.0, 1.0)
    top_n = st.slider("Gösterilecek aday sayısı", 3, 20, 10)

st.title("Planetary Analog Finder")
st.markdown(
    f"**{body.name} · {choice}** ({lat:.2f}°, {lon:.2f}°){' — ' + note if note else ''}  \n"
    f"{window_m / 1000:.0f} km × {window_m / 1000:.0f} km pencere, Dünya genelinde "
    f"{count(len(catalog))} pencere ile karşılaştırılıyor."
)

# Shown before the DEM read, which can take a minute on first use of a site.
st.subheader("Seçilen bölgenin geniş görünümü")
st.caption(
    f"Turuncu kare, karşılaştırılan {window_m / 1000:.0f} km × {window_m / 1000:.0f} km pencereyi gösterir. "
    "Uzaklaştırarak bölgenin çevresini ve gezegen üzerindeki konumunu görebilirsiniz."
)
st_folium(
    planet_map(body, float(lat), float(lon), window_m),
    height=460,
    use_container_width=True,
    returned_objects=[],
    key=f"planet_{body.key}_{float(lat):.3f}_{float(lon):.3f}",
)

try:
    reference = get_reference(body.key, float(lat), float(lon), window_m)
except Exception as exc:
    st.error(f"Referans bölge okunamadı: {exc}")
    st.stop()

if not np.isfinite(list(reference.values.values())).all():
    st.error(f"Bu pencerede yeterli yükseklik verisi yok (geçerli alan: %{reference.valid_frac * 100:.0f}).")
    st.stop()

reference_key = (body.key, round(float(lat), 3), round(float(lon), 3), window_m)
ranking = get_ranking(catalog, names, reference_key, reference.values, tuple(sorted(weights.items())), filter_values, top_n)

notices = []
if body.key == "mars":
    measured = get_shot_coverage(float(lat), float(lon), window_m)
    if measured is not None:
        notices.append(
            f"**Veri desteği:** Bu penceredeki 463 m MOLA hücrelerinin yalnızca %{measured * 100:.0f}'ünde gerçek lazer ölçümü var. "
            "Kalanı enterpolasyon olduğu için 1 km özellikleri yarım ağırlıkla sayılıyor."
        )
if ranking.out_of_range:
    labels = ", ".join(label_of(n)[0] for n in ranking.out_of_range)
    notices.append(f"**Dünya'da nadir:** Referansın şu özellikleri Dünya pencerelerinin olağan aralığının dışında: {labels}.")
for notice in notices:
    st.info(notice)

if ranking.table.empty:
    st.warning("Bu filtrelerle hiçbir Dünya penceresi kalmadı. Yüzey filtresini gevşetin.")
    st.stop()

table = ranking.table
map_column, list_column = st.columns(2, gap="medium")
with list_column:
    st.subheader("En yakın adaylar")
    st.caption(
        f"Filtreden geçen {count(ranking.n_after_filter)} pencere içinden. "
        "Benzerlik göreli bir göstergedir, olasılık değildir. Yüzdelik: taranan pencerelerin kaçından daha yakın. "
        "Kararlılık: ağırlıklar rastgele değiştiğinde en yakın %0,1 içinde kalma oranı."
    )
    shown = pd.DataFrame(
        {
            "#": table["rank"],
            "Enlem": table["lat"].round(2),
            "Boylam": table["lon"].round(2),
            "Benzerlik": table["score"].round(0).astype(int),
            "Yüzdelik": table["percentile"].round(3),
            "Kararlılık %": (table["stability"] * 100).round(0).astype(int),
            "Çıplak %": (table["lc_bare"] * 100).round(0).astype(int),
        }
    )
    st.dataframe(shown, hide_index=True, width="stretch", height=380)
    selected = st.selectbox("Ayrıntılı incelenecek aday", table["rank"].tolist(), format_func=lambda r: f"#{r} · {table.loc[r - 1, 'lat']:.2f}, {table.loc[r - 1, 'lon']:.2f}")

with map_column:
    st.subheader("Dünya üzerindeki konumlar")
    st_folium(world_map(table, selected, window_m), height=470, use_container_width=True, returned_objects=[], key=f"map_{selected}_{hash(reference_key)}")

candidate = table.loc[selected - 1]
candidate_features = get_earth_window(float(candidate["lat"]), float(candidate["lon"]), window_m)
reasons = explain(ranking, reference.values, candidate)
agreement = family_agreement(candidate)

st.divider()
st.subheader(f"Aday #{selected} · Karşılaştırma raporu")

best_family = max(agreement, key=agreement.get)
worst = reasons.iloc[0]
st.markdown(
    f"En iyi uyum **{best_family.lower()}** ailesinde. Toplam farkın en büyük kaynağı **{worst['label']}**: "
    f"referansta {worst['reference']:.2f}{worst['unit']}, adayda {worst['candidate']:.2f}{worst['unit']}. "
    "Benzerlik yalnızca arazi biçimini ölçer, jeoloji veya mineralojiyi değil."
)

image_a, image_b, bars = st.columns([1, 1, 1.1], gap="medium")
with image_a:
    st.image(relief_image(reference.grids[SCALE]), caption=f"{body.name} · {body.dem.id.split(' (')[0]} · kuzey yukarı", clamp=True, width="stretch")
with image_b:
    st.image(relief_image(candidate_features.grids[SCALE]), caption="Dünya adayı · NASADEM/Copernicus · aynı fiziksel ölçek", clamp=True, width="stretch")
with bars:
    fig = go.Figure(
        go.Bar(
            x=list(agreement.values()),
            y=list(agreement.keys()),
            orientation="h",
            marker={"color": COLORS["earth"], "cornerradius": 4},
            text=[f"{value:.0f}" for value in agreement.values()],
            textposition="outside",
            hovertemplate="%{y}: %{x:.0f}/100<extra></extra>",
        )
    )
    layout = base_layout(280, "Aile bazında uyum (0–100)", "")
    layout["xaxis"]["range"] = [0, 112]
    layout["yaxis"]["autorange"] = "reversed"
    layout["hovermode"] = "closest"
    fig.update_layout(**layout)
    st.plotly_chart(fig, width="stretch", config={"displayModeBar": False})

grid_ref = reference.grids[SCALE]
grid_cand = candidate_features.grids[SCALE]
chart_a, chart_b, chart_c = st.columns(3, gap="medium")
with chart_a:
    st.markdown("**Eğim dağılımı @ 1 km**")
    x_a, y_a = slope_cdf(grid_ref, SCALE)
    x_b, y_b = slope_cdf(grid_cand, SCALE)
    st.plotly_chart(two_line_chart(x_a, y_a, x_b, y_b, "Eğim (°)", "Alan payı", body.name, "%{x:.1f}°: %{y:.0%}<extra></extra>"), width="stretch", config={"displayModeBar": False})
with chart_b:
    st.markdown("**Hipsometrik eğri**")
    x_a, y_a = hypsometric_curve(grid_ref)
    x_b, y_b = hypsometric_curve(grid_cand)
    st.plotly_chart(two_line_chart(x_a, y_a, x_b, y_b, "Daha yüksekte kalan alan payı", "Bağıl yükseklik", body.name, "%{x:.0%}: %{y:.2f}<extra></extra>"), width="stretch", config={"displayModeBar": False})
with chart_c:
    st.markdown("**Radyal profil**")
    prof_a, prof_b = reference.profile, candidate_features.profile
    if prof_a is not None and prof_b is not None:
        st.plotly_chart(two_line_chart(prof_a.radius_m / 1000, prof_a.shape, prof_b.radius_m / 1000, prof_b.shape, "Profil merkezinden uzaklık (km)", "Halka ortalaması (σ)", body.name, "%{x:.1f} km: %{y:.2f} σ<extra></extra>"), width="stretch", config={"displayModeBar": False})

        def shift(profile) -> str:
            return f"{np.hypot(profile.dx_m, profile.dy_m) / 1000:.0f} km"

        st.caption(
            f"Merkez, pencere merkezinin {PROFILE_SHIFT_M / 1000:.0f} km çevresinde en simetrik nokta olarak arandı "
            f"(kayma: {body.name} {shift(prof_a)}, aday {shift(prof_b)}). σ: penceredeki yükseklik oynaklığı. "
            f"Simetri R²: {body.name} {prof_a.r2:.2f}, aday {prof_b.r2:.2f}."
        )
    else:
        r_a, p_a, _ = radial_profile(grid_ref, SCALE)
        r_b, p_b, _ = radial_profile(grid_cand, SCALE)
        st.plotly_chart(two_line_chart(r_a / 1000, p_a, r_b / 1000, p_b, "Merkezden uzaklık (km)", "Bağıl yükseklik (m)", body.name, "%{x:.0f} km: %{y:.0f} m<extra></extra>"), width="stretch", config={"displayModeBar": False})

with st.expander("Özellik tablosu: her özelliğin farka katkısı"):
    frame = reasons.copy()
    frame["Katkı %"] = (frame["contribution"] * 100).round(1)
    st.dataframe(
        frame.rename(columns={"label": "Özellik", "family": "Aile", "unit": "Birim", "reference": body.name, "candidate": "Dünya adayı"})[
            ["Özellik", "Aile", "Birim", body.name, "Dünya adayı", "Katkı %"]
        ].round(3),
        hide_index=True,
        width="stretch",
    )

st.subheader("Yüzey doğrulaması: Sentinel-2")
surface_key = f"surface_{candidate['lat']:.3f}_{candidate['lon']:.3f}"
if st.button("Bu aday için Sentinel-2 sahnesini getir", key=surface_key) or st.session_state.get(f"done_{surface_key}"):
    st.session_state[f"done_{surface_key}"] = True
    try:
        surface = get_surface(float(candidate["lat"]), float(candidate["lon"]), window_m)
    except Exception as exc:
        surface = None
        st.error(f"Sentinel-2 okunamadı: {exc}")
    if surface is None:
        st.warning("Bu aday için uygun, az bulutlu Sentinel-2 sahnesi bulunamadı.")
    else:
        summary, rgb, ndvi = surface
        rgb_col, stats_col = st.columns([1, 1.4], gap="medium")
        with rgb_col:
            st.image(rgb, caption=f"{summary.scene_id} · {summary.date}", width="stretch")
        with stats_col:
            m1, m2, m3 = st.columns(3)
            m1.metric("Medyan NDVI", f"{summary.ndvi_median:.2f}")
            m2.metric("Yeşil alan (NDVI > 0,4)", f"%{summary.green_frac * 100:.0f}")
            m3.metric("Kapsama", f"%{summary.coverage * 100:.0f}")
            m4, m5, m6 = st.columns(3)
            m4.metric("Medyan NDMI", f"{summary.ndmi_median:.2f}")
            m5.metric("SWIR B11/B12", f"{summary.swir_ratio_median:.2f}")
            m6.metric("Sahne bulutu", f"%{summary.cloud_pct:.1f}")
            st.caption("Bu oranlar yüzey malzemesi hakkında ipucu verir. Belirli bir minerali tek başına tanımlamaz.")

st.subheader("Literatürdeki analoglar bu referansa ne kadar yakın?")
st.caption("Her alan kendi tam koordinatında ölçüldü. Yüzdelik, Dünya penceresinin kaçından daha yakın olduğunu gösterir.")
known = known_analog_table(ranking, reference_key, tuple(sorted(weights.items())), window_m)
st.dataframe(known.assign(Mesafe=known["Mesafe"].round(2), Yüzdelik=known["Yüzdelik"].round(1)), hide_index=True, width="stretch", height=35 * (len(known) + 1) + 3)

with st.expander("Bilimsel izlenebilirlik"):
    st.markdown(f"**Yöntem.** {METHOD}")
    st.markdown(
        f"**Referans DEM.** {body.dem.id}, yaklaşık {body.dem.native_m:.0f} m ızgara. {body.dem.reference}.  \n"
        f"**Katalog.** `{catalog_name}`, {meta['created_utc']} tarihinde üretildi. "
        f"Pencere {meta['window_m'] / 1000:.0f} km, adım {meta['stride_m'] / 1000:.0f} km, "
        f"ölçekler {', '.join(f'{s / 1000:.0f}' for s in meta['scales_m'])} km.  \n"
        f"**Başarısız karolar.** DEM {len(meta['failed_dem_tiles'])}, WorldCover {len(meta['failed_cover_tiles'])}.  \n"
        f"**Yüzey.** {S2_PRODUCT}.  \n"
        "**Geniş görünüm.** NASA Solar System Treks haritaları yalnızca görsel bağlam içindir, benzerlik hesabına girmez."
    )
    st.json({"veri_kaynaklari": meta["data_sources"], "yazilim": meta["software"]})
    st.markdown(
        "**Sınırlamalar.** Benzerlik puanı olasılık değildir ve yalnızca arazi biçimini karşılaştırır. "
        "Eşik altı ölçekler (yaklaşık 1 km'den küçük kaya, rover rotası) temsil edilmez. "
        "Yerçekimi ve aşınma farkı nedeniyle aynı biçimler farklı süreçlerden doğmuş olabilir."
    )
