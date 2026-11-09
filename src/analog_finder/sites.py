"""Reference sites and literature-known Earth analogs."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Site:
    name: str
    lat: float
    lon: float
    note: str = ""


MARS_SITES = (
    Site("Jezero Krateri", 18.38, 77.58, "Perseverance; eski göl ve delta"),
    Site("Gale Krateri", -5.37, 137.81, "Curiosity; tabakalı merkezi höyük"),
    Site("Gusev Krateri", -14.57, 175.48, "Spirit"),
    Site("Meridiani Planum", -1.95, -5.53, "Opportunity; düz ova"),
    Site("Oxia Planum", 18.28, -24.63, "Rosalind Franklin hedefi; killi düzlük"),
    Site("Holden Krateri", -26.40, -34.00, "Fan ve göl çökelleri"),
)

MOON_SITES = (
    Site("Apollo 11 – Tranquillitatis", 0.67, 23.47, "Bazaltik mare"),
    Site("Apollo 15 – Hadley", 26.13, 3.63, "Mare, Apennin dağları ve rille"),
    Site("Apollo 17 – Taurus-Littrow", 20.19, 30.77, "Yüksek arazi vadisi"),
    Site("Copernicus Krateri", 9.62, -20.08, "Genç kompleks krater"),
    Site("Tycho Krateri", -43.31, -11.36, "Genç kompleks krater"),
)

REFERENCE_SITES = {"mars": MARS_SITES, "moon": MOON_SITES}

# Sites discussed in the analog literature. The note says what they are an
# analog *of*: many are chosen for chemistry or climate, not landform.
KNOWN_ANALOGS = (
    Site("Salda Gölü", 37.55, 29.68, "Mineraloji: Jezero karbonatları"),
    Site("Atacama – Yungay", -24.08, -69.99, "İklim: aşırı kuraklık"),
    Site("Rio Tinto", 37.72, -6.56, "Kimya: asidik ve demirli su"),
    Site("Death Valley", 36.50, -117.00, "Kurak havza ve evaporitler"),
    Site("Mauna Kea", 19.82, -155.47, "Volkanik regolit"),
    Site("Askja", 65.05, -16.75, "Volkanik arazi; Apollo eğitim alanı"),
    Site("Haughton Krateri", 75.38, -89.68, "Morfoloji: kutupsal çarpma krateri, 23 km"),
    Site("McMurdo Kuru Vadileri", -77.50, 161.80, "İklim: soğuk ve kuru"),
    Site("Ries Krateri", 48.88, 10.57, "Morfoloji: çarpma krateri, 24 km"),
    Site("Bosumtwi", 6.50, -1.41, "Morfoloji: çarpma krateri gölü, 10,5 km"),
    Site("Gosses Bluff", -23.82, 132.31, "Morfoloji: aşınmış çarpma yapısı, 22 km"),
    Site("Mistastin", 55.88, -63.31, "Morfoloji: çarpma krateri gölü, 28 km"),
    Site("Manicouagan", 51.38, -68.70, "Morfoloji: çarpma yapısı, yaklaşık 85 km"),
    Site("Vredefort", -27.00, 27.50, "Morfoloji: çok aşınmış çarpma yapısı"),
    Site("Richat Yapısı", 21.12, -11.40, "Morfoloji: dairesel erozyon yapısı, çarpma değil"),
)
