# Planetary Analog Finder

NASA Space Apps Challenge projesi. Mars veya Ay üzerinde seçilen bir bölgeye **arazi biçimi açısından** en çok benzeyen Dünya bölgelerini bulur ve her eşleşmenin nedenini, veri kaynağını ve belirsizliğini gösterir.

Benzer görünen iki uydu görüntüsü bilimsel olarak benzer iki yüzey anlamına gelmez. Bu yüzden sistem görüntüleri değil, metre cinsinden hesaplanan ve fiziksel olarak yorumlanabilen morfoloji özelliklerini karşılaştırır.

## Hızlı başlangıç

```bash
python3 -m venv .venv
.venv/bin/pip install -e ".[dev]"

# Etkileşimli panel
.venv/bin/streamlit run app/streamlit_app.py

# Testler (ağ gerektirmez)
.venv/bin/python -m pytest
# Gerçek uzak verilerle testler
.venv/bin/python -m pytest -m network
```

Hazır küresel katalog `data/catalog/earth_global.parquet` dosyasındadır, panel doğrudan onu kullanır. Bilimsel analiz [notebooks/01_jezero_salda.ipynb](notebooks/01_jezero_salda.ipynb) içindedir; `.venv` çekirdeğiyle açın.

Kataloğu yeniden üretmek için (hesap gerekmez; yaklaşık 45 dakika, 1,1 GB önbellek):

```bash
.venv/bin/python scripts/build_catalog.py --global --output data/catalog/earth_global.parquet
# Yalnızca bir bölge (güney kuzey batı doğu)
.venv/bin/python scripts/build_catalog.py --region 36 42 26 45 --output data/catalog/turkiye.parquet
```

## İlk sonuçlar: Jezero Krateri

| Soru | Sonuç |
|---|---|
| Salda, Jezero'ya topografik olarak benziyor mu? | Hayır. Dünya pencerelerinin yalnızca %11,6'sından daha yakın. Salda'nın analogluğu mineralojiktir. |
| Jezero'nun Dünya'daki en yakın eşi | Çin'in Sincan bölgesinde kapalı bir havza (43,47°K, 91,62°D). İlk 10 adayın 4'ü Atacama/Puna'da. |
| Katalog | 324.859 pencere (60 km, 20 km adım); 48.763'ü çıplak zemin filtresinden geçiyor. |
| MOLA veri desteği | Jezero penceresindeki hücrelerin %32,9'unda gerçek lazer ölçümü var. |
| Kör arama (399 pencere, gerçek MOLA iziyle bozulmuş) | %85,7 ilk sırada, %95,7 ilk 10'da. |

Adaylar dairesel krater değil, kapalı havzalardır. Metrikler dönmeden bağımsız özet istatistikler olduğu için yöntem "kenarları yüksek, tabanı düz havza" biçimini bulur. Kör arama yöntemin üst sınırını gösterir, gezegenler arası geçerliliği kanıtlamaz.

## Yöntem

Her gövde için aynı adımlar uygulanır:

1. **Destek eşleme.** DEM, kaynak pikselin sahip olduğu desteği hesaba katan bir Gauss filtresiyle analiz ölçeğine indirilir. Böylece 30 m ve 463 m kaynaklar aynı fiziksel desteğe gelir.
2. **Yerel metrik ızgara.** Pencere, gövdenin kendi yarıçapıyla azimutal eşit uzaklık ızgarasına örneklenir. Gövdeler arasında mutlak yükseklik ve koordinat karşılaştırılmaz.
3. **Özellikler.** 1, 2 ve 4 km ölçeklerinde eğim istatistikleri, TRI ve pürüzlülük; pencere düzeyinde rölyef, hipsometrik integral, çarpıklık, dairesel simetri ve havza derinliği.
4. **Sıralama.** Özellikler Dünya pencerelerinin dağılımına göre ölçeklenir ve dört aile (eğim, engebelilik, rölyef, biçim) ağırlıklandırılır. MOLA iz deseni deneyi 1 km özelliklerini en çok etkilediği için bu ölçek yarım ağırlıkla sayılır.
5. **Eleme ve yüzey.** WorldCover ile bitki örtüsü, tarım, su ve yerleşim elenir. Seçilen adaylar Sentinel-2 ile incelenir.
6. **Açıklama.** Her aday için özellik bazında fark katkısı, ağırlık değişimine kararlılık ve Dünya pencereleri içindeki yüzdelik verilir.

## Veri kaynakları

| Katman | Ürün | Erişim | Hesap |
|---|---|---|---|
| Mars topografyası | MGS MOLA MEGDR 128 ppd, yaklaşık 463 m | USGS global mozaik | Gerekmez |
| Mars veri desteği | MOLA MEGDR ölçüm sayısı (counts) | NASA PDS Geosciences | Gerekmez |
| Ay topografyası | LRO LOLA LDEM 256 ppd, yaklaşık 118 m | USGS global mozaik | Gerekmez |
| Dünya topografyası (56°G–60°K) | NASADEM HGT v001, 1 yay saniyesi | Microsoft Planetary Computer COG | Gerekmez |
| Dünya topografyası (diğer enlemler) | Copernicus DEM GLO-90 | AWS Open Data | Gerekmez |
| Arazi örtüsü | ESA WorldCover 2021 v200 | AWS Open Data | Gerekmez |
| Yüzey | Sentinel-2 L2A | Earth Search / AWS Open Data | Gerekmez |
| Geniş görünüm (yalnızca görsel) | Viking MDIM 2.1 ve MOLA renkli rölyef (Mars), LROC WAC ve LOLA renkli rölyef (Ay) | NASA Solar System Treks | Gerekmez |

### Araştırma raporundan farklar

- NASADEM, Earthdata hesabı gerektirmeyen Planetary Computer kopyasından okunur. Sentinel-2 için de Copernicus girişi yerine Earth Search kullanılır.
- Rapor NASADEM'i arazi yüksekliği, Copernicus DEM'i yüzey modeli olarak ayırır. SRTM kökenli NASADEM de radar tabanlı bir yüzey modelidir ve bitki örtüsünden etkilenir. Çıplak bölgelerde ve 1 km ölçeğinde bu fark küçüktür.
- Katalog, MOLA ızgarasıyla yaklaşık aynı aralıkta olan 15 yay saniyelik Dünya ızgarası üzerinde kurulur. 30 m kaynakla farkı notebook'ta ölçülmüştür: 1 km'de en çok %4,6.

## Proje yapısı

```text
src/analog_finder/
  bodies.py        Gövdeler ve DEM ürünleri
  dem.py           Uzak DEM okuma, karo birleştirme, önbellek
  terrain.py       Destek eşleme ve morfoloji özellikleri
  similarity.py    Sıralama, açıklama, kararlılık
  catalog.py       Küresel Dünya kataloğu
  landcover.py     WorldCover payları
  surface.py       Sentinel-2 göstergeleri
  mola_support.py  MOLA gerçek ölçüm yoğunluğu
  validation.py    Kör arama ve yanlılık deneyleri
  provenance.py    Yöntem, sürüm ve kaynak kayıtları
app/               Streamlit paneli
notebooks/         Jezero–Salda analizi
scripts/           Katalog üretimi
tests/             Analitik ve veri katmanı testleri
```

## Sınırlamalar

- Benzerlik puanı olasılık değildir ve jeoloji veya mineraloji içermez.
- MOLA, 1 km'den küçük yapıları güvenilir biçimde temsil etmez. Jezero deltası gibi birkaç km'lik yapılar için HRSC gerekir.
- Pencereler 20 km adımla kayar; bir yapının merkezi 14 km'ye kadar kaçabilir.
- Mars ve Ay'da yerçekimi, atmosfer ve aşınma farklıdır. Benzer biçimler farklı süreçlerden doğmuş olabilir.
- 80°'den yüksek enlemlerde Mars ve Ay için kutupsal ürünler henüz entegre değildir.

## Sonraki adımlar

1. Radyal profili doğrudan karşılaştıran, krater biçimine duyarlı bir özellik.
2. Birden çok pencere boyutu (ör. 15, 30, 120 km).
3. USGS Ay jeoloji haritası ile jeolojik bağlam.
4. Seçilmiş Mars bölgeleri için HRSC, ardından CRISM/OMEGA ile mineraloji.
