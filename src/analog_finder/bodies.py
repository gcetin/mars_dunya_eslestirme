"""Planetary bodies and the DEM products used to represent them."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class DemProduct:
    id: str
    url: str
    native_m: float
    reference: str


@dataclass(frozen=True)
class Body:
    key: str
    name: str
    radius_m: float
    dem: DemProduct | None = None


MARS = Body(
    key="mars",
    name="Mars",
    radius_m=3_396_190.0,
    dem=DemProduct(
        id="MGS MOLA MEGDR 128 ppd (USGS global mosaic)",
        url="https://planetarymaps.usgs.gov/mosaic/Mars_MGS_MOLA_DEM_mosaic_global_463m.tif",
        native_m=463.09,
        reference="Smith et al. (2001), JGR 106(E10); USGS Astrogeology mosaic",
    ),
)

MOON = Body(
    key="moon",
    name="Ay",
    radius_m=1_737_400.0,
    dem=DemProduct(
        id="LRO LOLA LDEM 256 ppd (USGS global mosaic, Mar 2014)",
        url="https://planetarymaps.usgs.gov/mosaic/Lunar_LRO_LOLA_Global_LDEM_118m_Mar2014.tif",
        native_m=118.45,
        reference="Smith et al. (2010), GRL 37; USGS Astrogeology mosaic",
    ),
)

# Earth is treated as a sphere of authalic radius. Across a 60 km window the
# error against WGS84 radii of curvature stays below 0.4 %, far under DEM error.
EARTH = Body(key="earth", name="Dünya", radius_m=6_371_007.0)

BODIES = {body.key: body for body in (MARS, MOON, EARTH)}
