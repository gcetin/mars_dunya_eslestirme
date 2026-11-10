# Planetary Analog Finder

NASA Space Apps Challenge projesi. Mars veya Ay üzerinde seçilen bir bölgeye **arazi biçimi açısından** en çok benzeyen Dünya bölgelerini bulur ve her eşleşmenin nedenini, veri kaynağını ve belirsizliğini gösterir.

Benzer görünen iki uydu görüntüsü bilimsel olarak benzer iki yüzey anlamına gelmez. Bu yüzden sistem görüntüleri değil, metre cinsinden hesaplanan ve fiziksel olarak yorumlanabilen morfoloji özelliklerini karşılaştırır.

**Canlı uygulama:** https://mars-dunya-eslestirme.streamlit.app

Yazılımın baştan sona nasıl çalıştığı, hangi veriyi ne için kullandığı, hesapların ayrıntısı ve uygulamanın kullanımı [docs/nasil-calisir.md](docs/nasil-calisir.md) dosyasında anlatılır. Katkı vermek için [CONTRIBUTING.md](CONTRIBUTING.md) dosyasına bakın.

## Hızlı başlangıç

Python 3.11 veya üstü gerekir; panel 3.13 ile test edilmiştir.

### macOS ve Linux

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

### Windows

PowerShell'de:

```powershell
# Python ve Git yoksa
winget install Python.Python.3.13 Git.Git

git clone https://github.com/gcetin/mars_dunya_eslestirme.git
cd mars_dunya_eslestirme
py -3.13 -m venv .venv
.venv\Scripts\python -m pip install -e ".[dev]"

# Etkileşimli panel
.venv\Scripts\streamlit run app\streamlit_app.py

# Testler (ağ gerektirmez)
.venv\Scripts\python -m pytest
```

Bu bölümdeki diğer komutlar da Windows'ta aynıdır; yalnızca `.venv/bin/` yerine `.venv\Scripts\` yazılır. Satır sonlarının gereksiz fark üretmemesi için bir kez `git config --global core.autocrlf true` çalıştırın. WSL2 kullanıyorsanız macOS ve Linux komutları olduğu gibi çalışır. Windows kurulumu henüz gerçek bir makinede denenmedi; sorun çıkarsa bildirin.

### Veri ve katalog

Hazır küresel katalog `data/catalog/earth_global.parquet` dosyasındadır, panel doğrudan onu kullanır. Bilimsel analiz [notebooks/01_jezero_salda.ipynb](notebooks/01_jezero_salda.ipynb) içindedir; `.venv` çekirdeğiyle açın.

Kataloğu yeniden üretmek için (hesap gerekmez; sıfırdan yaklaşık bir saat, 1,1 GB önbellek; önbellek doluysa 15 dakika):

```bash
.venv/bin/python scripts/build_catalog.py --global --output data/catalog/earth_global.parquet
# Yalnızca bir bölge (güney kuzey batı doğu)
.venv/bin/python scripts/build_catalog.py --region 36 42 26 45 --output data/catalog/turkiye.parquet
```

## İlk sonuçlar: Jezero Krateri

| Soru | Sonuç |
|---|---|
| Salda, Jezero'ya topografik olarak benziyor mu? | Hayır. Dünya pencerelerinin yalnızca %9,8'inden daha yakın. Salda'nın analogluğu mineralojiktir. |
| Jezero'nun Dünya'daki en yakın eşi | Çin'in Sincan bölgesinde kapalı bir havza (43,47°K, 91,62°D). İlk 10 adayın 4'ü Atacama/Puna'da. |
| Katalog | 324.859 pencere (60 km, 20 km adım), pencere başına 33 özellik; 48.763'ü çıplak zemin filtresinden geçiyor. |
| MOLA veri desteği | Jezero penceresindeki hücrelerin %32,9'unda gerçek lazer ölçümü var. |
| Kör arama, sorgu pencere merkezinde (396 pencere, gerçek MOLA iziyle bozulmuş) | %88,6 ilk sırada, %91,4 ilk 10'da. Radyal profil olmadan %84,8 ve %95,7. |
| Kör arama, sorgu 10 km'ye kadar kaymış (398 nokta) | %20,1 ilk sırada, %33,9 ilk 10'da. Radyal profil olmadan %9,3 ve %23,6. |
| Krater biçimi (Tooting, Mars, 29 km; yalnızca radyal profil) | İlk aday Kanada'daki Mistastin çarpma krateri (28 km). |

Jezero'nun adayları dairesel krater değil, kapalı havzalardır. Jezero'nun kenarı aşınmış ve 60 km'lik pencerede dışarı doğru alçalmadığı için radyal profili de böyle bir havzanınkine benziyor; profil eklendikten sonra liste değişmedi. Kenarı belirgin genç kraterlerde ise profil gerçek kraterleri bulabiliyor (Tooting → Mistastin).

Kaymış sorgu deneyi daha gerçekçidir, çünkü bir Mars bölgesi hiçbir Dünya pencere merkezine denk gelmez. Pencere merkezindeki kör arama yöntemin başarısını olduğundan yüksek gösterir. Radyal profil kaymış sorguda isabeti iki katına çıkarıyor, pencere merkezindeki sorguda ise ilk 10 oranını biraz düşürüyor. Her iki deney de sorgu ve katalog aynı gezegenden geldiği için bir üst sınırdır, gezegenler arası geçerliliği kanıtlamaz.

## Yöntem

Her gövde için aynı adımlar uygulanır:

1. **Destek eşleme.** DEM, kaynak pikselin sahip olduğu desteği hesaba katan bir Gauss filtresiyle analiz ölçeğine indirilir. Böylece 30 m ve 463 m kaynaklar aynı fiziksel desteğe gelir.
2. **Yerel metrik ızgara.** Pencere, gövdenin kendi yarıçapıyla azimutal eşit uzaklık ızgarasına örneklenir. Gövdeler arasında mutlak yükseklik ve koordinat karşılaştırılmaz.
3. **Özellikler.** 1, 2 ve 4 km ölçeklerinde eğim istatistikleri, TRI ve pürüzlülük; pencere düzeyinde rölyef, hipsometrik integral, çarpıklık, dairesel simetri ve havza derinliği. Ayrıca 12 halkalık bir radyal profil: merkezi pencere merkezinin 15 km çevresinde en simetrik nokta olarak aranır ve yükseklik oynaklığına bölünerek birimsiz yazılır, böylece yerçekimine bağlı krater derinliği karşılaştırmaya girmez.
4. **Sıralama.** Özellikler Dünya pencerelerinin dağılımına göre ölçeklenir ve beş aile ağırlıklandırılır: eğim 0,25, engebelilik 0,20, rölyef 0,15, biçim 0,20, radyal profil 0,20. MOLA iz deseni deneyi 1 km özelliklerini en çok etkilediği için bu ölçek yarım ağırlıkla sayılır.
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
  terrain.py       Destek eşleme, morfoloji özellikleri, radyal profil
  similarity.py    Sıralama, açıklama, kararlılık
  catalog.py       Küresel Dünya kataloğu
  landcover.py     WorldCover payları
  surface.py       Sentinel-2 göstergeleri
  mola_support.py  MOLA gerçek ölçüm yoğunluğu
  validation.py    Kör arama (pencere merkezinde ve kaymış), yanlılık deneyleri
  provenance.py    Yöntem, sürüm ve kaynak kayıtları
app/               Streamlit paneli
notebooks/         Jezero–Salda analizi
scripts/           Katalog üretimi, belge şekilleri
tests/             Analitik ve veri katmanı testleri
```

## Sınırlamalar

- Benzerlik puanı olasılık değildir ve jeoloji veya mineraloji içermez.
- MOLA, 1 km'den küçük yapıları güvenilir biçimde temsil etmez. Jezero deltası gibi birkaç km'lik yapılar için HRSC gerekir.
- Pencereler 20 km adımla kayar; bir yapının merkezi 14 km'ye kadar kaçabilir. Radyal profil merkezini aradığı için bundan etkilenmez, diğer özellikler etkilenir.
- Pencere boyu tektir (60 km); çok küçük veya büyük yapılar bu ölçeğe sığmaz.
- Mars ve Ay'da yerçekimi, atmosfer ve aşınma farklıdır. Benzer biçimler farklı süreçlerden doğmuş olabilir.
- 80°'den yüksek enlemlerde Mars ve Ay için kutupsal ürünler henüz entegre değildir.

## Sonraki adımlar

1. Birden çok pencere boyutu (ör. 15, 30, 120 km) ve konum kaymasına daha dayanıklı arama.
2. USGS Ay jeoloji haritası ile jeolojik bağlam.
3. Seçilmiş Mars bölgeleri için HRSC, ardından CRISM/OMEGA ile mineraloji.

## Lisans

Kod [MIT lisansıyla](LICENSE) yayımlanır. Lisans yalnızca kodu kapsar: `data/catalog/` altındaki katalog ve `docs/img/` altındaki şekiller aşağıdaki verilerden türetilmiştir ve o verilerin kullanım koşullarına tabidir.

## Veri atıfları

- **NASADEM:** NASA JPL (2020). NASADEM Merged DEM Global 1 arc second V001. NASA EOSDIS Land Processes DAAC. doi:10.5067/MEaSUREs/NASADEM/NASADEM_HGT.001. Microsoft Planetary Computer üzerinden okunur.
- **Copernicus DEM GLO-90:** "© DLR e.V. 2010–2014 and © Airbus Defence and Space GmbH 2014–2018 provided under COPERNICUS by the European Union and ESA; all rights reserved." AWS Open Data üzerinden okunur.
- **ESA WorldCover 2021 v200:** "© ESA WorldCover project 2021 / Contains modified Copernicus Sentinel data (2021) processed by ESA WorldCover consortium." CC BY 4.0. Zanaga, D. ve diğ. (2022), doi:10.5281/zenodo.7254221.
- **Sentinel-2 L2A:** "Contains modified Copernicus Sentinel data (2025)." Element 84 Earth Search ve AWS Open Data üzerinden okunur.
- **MGS MOLA MEGDR:** Smith, D. E. ve diğ. (2001), JGR 106(E10). USGS Astrogeology global mozaiği; ölçüm sayıları NASA PDS Geosciences Node'dan.
- **LRO LOLA LDEM:** Smith, D. E. ve diğ. (2010), GRL 37. USGS Astrogeology global mozaiği.
- **NASA Solar System Treks:** Mars ve Ay geniş görünüm katmanları, yalnızca görsel bağlam için.
- **Uygulamadaki Dünya haritası:** © OpenStreetMap katkıcıları (ODbL); Esri World Imagery (Esri, i-cubed, USDA, USGS, AEX, GeoEye, Getmapping, Aerogrid, IGN, IGP, UPR-EGP ve GIS User Community).
