"""Records of what data, versions and methods produced a result."""

from __future__ import annotations

import platform
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version

import rasterio

from . import __version__
from .bodies import MARS, MOON
from .landcover import PRODUCT as WORLDCOVER_PRODUCT

EARTH_DEM_NOTES = {
    "nasadem": "NASADEM_HGT v001 (NASA JPL; Microsoft Planetary Computer COG), 56°G–60°K",
    "copdem": "Copernicus DEM GLO-90 (AWS Open Data), 60°K üstü ve 56°G altı",
}

METHOD = (
    "Her gövde için aynı işlem: kaynak pikselin desteğini hesaba katan Gauss alçak geçiren filtre "
    "(sigma = ölçek × 0,5), kendi yarıçapıyla yerel azimutal eşit uzaklık ızgarasına yeniden örnekleme, "
    "ardından yalnızca metre cinsinden yerel morfoloji özellikleri. Mutlak yükseklik ve koordinat "
    "gövdeler arasında karşılaştırılmaz."
)


def software_versions() -> dict[str, str]:
    packages = {}
    for name in ("numpy", "scipy", "pandas", "rasterio"):
        try:
            packages[name] = version(name)
        except PackageNotFoundError:
            packages[name] = "unknown"
    return {
        "analog_finder": __version__,
        "python": platform.python_version(),
        "gdal": rasterio.__gdal_version__,
        **packages,
    }


def data_sources() -> dict[str, str]:
    return {
        "mars_dem": f"{MARS.dem.id} — {MARS.dem.url}",
        "moon_dem": f"{MOON.dem.id} — {MOON.dem.url}",
        "earth_dem": "; ".join(EARTH_DEM_NOTES.values()),
        "land_cover": WORLDCOVER_PRODUCT,
    }


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")
