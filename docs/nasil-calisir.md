# Planetary Analog Finder: Nasıl Çalışır?

Bu belge yazılımı baştan sona anlatır: hangi veriyi ne için kullandığımızı, hangi hesapları yaptığımızı, kodun nasıl düzenlendiğini ve uygulamanın nasıl kullanılacağını. Kurulum komutları ve ilk bilimsel sonuçlar [README](../README.md) dosyasındadır.

## İçindekiler

1. [Ne yapıyoruz?](#1-ne-yapıyoruz)
2. [Genel akış](#2-genel-akış)
3. [Veri kaynakları ve görevleri](#3-veri-kaynakları-ve-görevleri)
4. [Bir pencerenin işlenmesi](#4-bir-pencerenin-işlenmesi)
5. [Dünya kataloğu](#5-dünya-kataloğu)
6. [Benzerlik ve sıralama](#6-benzerlik-ve-sıralama)
7. [Yüzey ve güvenilirlik katmanları](#7-yüzey-ve-güvenilirlik-katmanları)
8. [Uygulamayı kullanma](#8-uygulamayı-kullanma)
9. [Notebook, komut satırı ve Python](#9-notebook-komut-satırı-ve-python)
10. [Kod yapısı](#10-kod-yapısı)
11. [Geliştirme](#11-geliştirme)
12. [Sınırlamalar](#12-sınırlamalar)
13. [Yol haritası](#13-yol-haritası)
14. [Terimler](#14-terimler)

## 1. Ne yapıyoruz?

Mars veya Ay üzerinde seçilen 60 km × 60 km'lik bir bölgeye **arazi biçimi bakımından** en çok benzeyen Dünya bölgelerini buluyoruz. Her eşleşme için neden seçildiğini, hangi veriye dayandığını ve sonucun ne kadar sağlam olduğunu da gösteriyoruz.

Temel fikir şu: iki uydu görüntüsünün birbirine benzemesi, iki yüzeyin bilimsel olarak benzer olduğu anlamına gelmez. Güneş açısı, renk, çözünürlük ve atmosfer görüntüyü değiştirir. Bu yüzden görüntüleri değil, yükseklik modellerinden (DEM) metre cinsinden hesaplanan ve fiziksel anlamı olan biçim ölçülerini karşılaştırıyoruz: eğim, engebelilik, rölyef ve havza biçimi gibi.

Her adımda üç kural geçerlidir:

- **Mutlak yükseklik ve koordinat gövdeler arasında karşılaştırılmaz.** Mars'ın sıfır yüksekliği ile Dünya'nın deniz seviyesi aynı şey değildir.
- **Her veri önce aynı fiziksel çözünürlüğe getirilir.** 30 m'lik bir Dünya DEM'i ile 463 m'lik bir Mars DEM'i doğrudan karşılaştırılamaz.
- **Her sonuç kaynağı, sürümü ve belirsizliğiyle birlikte verilir.**

## 2. Genel akış

```mermaid
flowchart TD
    A["Mars / Ay DEM<br/>MOLA, LOLA"] --> C["Destek eşleme<br/>1, 2 ve 4 km"]
    B["Dünya DEM<br/>NASADEM, Copernicus"] --> C
    C --> D["Yerel metrik ızgara<br/>60 km pencere"]
    D --> E["21 morfoloji özelliği"]
    E -->|"Dünya: bir kez, önceden"| F[("Küresel katalog<br/>324.859 pencere")]
    E -->|"Referans: her sorguda"| G["Benzerlik ve sıralama"]
    F --> G
    H["ESA WorldCover"] -->|"yüzey filtresi"| G
    G --> I["Aday listesi, açıklama, kararlılık"]
    I --> J["Sentinel-2 yüzey kontrolü"]
    I --> K["Bilinen analoglarla karşılaştırma"]
```

İş iki aşamada yapılır:

- **Önceden (çevrimdışı):** Dünya'nın karasal yüzeyi 60 km'lik pencerelere bölünür. Her pencerenin özellikleri bir kez hesaplanıp bir katalog dosyasına yazılır (`scripts/build_catalog.py`).
- **Sorguda (çevrimiçi):** Kullanıcı bir Mars veya Ay bölgesi seçtiğinde yalnızca o pencerenin DEM'i okunur ve özellikleri hesaplanır. Ardından katalogdaki bütün pencerelerle bellekte karşılaştırılır. Bu birkaç saniye sürer.

## 3. Veri kaynakları ve görevleri

| Veri | Ne için kullanıyoruz | Kod | Erişim |
|---|---|---|---|
| MGS MOLA MEGDR 128 ppd (yaklaşık 463 m) | Mars referans penceresinin yüksekliği | `bodies.py`, `dem.py` | USGS global mozaik |
| MOLA MEGDR ölçüm sayısı (counts) | Mars penceresinde gerçek lazer ölçümü olan hücrelerin payı; iz deseni deneyi | `mola_support.py` | NASA PDS Geosciences |
| LRO LOLA LDEM 256 ppd (yaklaşık 118 m) | Ay referans penceresinin yüksekliği | `bodies.py`, `dem.py` | USGS global mozaik |
| NASADEM HGT v001, 1 yay saniyesi | Dünya yüksekliği, 56°G ile 60°K arası | `dem.py` | Microsoft Planetary Computer |
| Copernicus DEM GLO-90 | Dünya yüksekliği, 60°K üstü ve 56°G altı | `dem.py` | AWS Open Data |
| ESA WorldCover 2021 v200, 10 m | Dünya pencerelerinin arazi örtüsü payları; bitki örtüsü, tarım, su ve yerleşimin elenmesi | `landcover.py` | AWS Open Data |
| Sentinel-2 L2A | Seçilen adayın yüzeyine bakmak ve basit spektral göstergeler | `surface.py` | Earth Search STAC ve AWS |
| NASA Solar System Treks katmanları | Mars ve Ay haritasında geniş görünüm (yalnızca görsel) | `app/` | trek.nasa.gov |
| OpenStreetMap, Esri World Imagery | Dünya haritasında adayların konumu (yalnızca görsel) | `app/` | Folium karo katmanları |

Bu kaynakların hiçbiri hesap gerektirmez. Dosyalar Cloud Optimized GeoTIFF biçimindedir ve HTTP aralık istekleriyle okunur, yani yalnızca istenen bölgenin baytları indirilir. İndirilen her şey `data/cache/` altında saklanır. Bu klasör git'e girmez; küresel katalog için yaklaşık 1,1 GB tutar ve istenildiğinde silinebilir. Konumu `ANALOG_CACHE` ortam değişkeniyle değiştirilebilir.

Yükseklik verilerinin Mars ve Ay tarafı NASA görevlerinden gelir. Dünya tarafının büyük kısmı ise Avrupa'nın Copernicus programına aittir (Copernicus DEM, WorldCover, Sentinel-2).

## 4. Bir pencerenin işlenmesi

Aşağıdaki adımlar Mars, Ay ve Dünya pencereleri için birebir aynıdır. Kodun büyük kısmı `terrain.py` dosyasındadır.

### 4.1 DEM okuma (`dem.py`)

- **Mars ve Ay:** Global mozaikten pencere ile kenar payını kapsayan alan okunur. LOLA'nın ölçek çarpanı (0,5 m/DN) uygulanır. 80°'den yüksek enlemlerde kutup ürünü gerektiği için okuma reddedilir. ±180° boylam çizgisini kesen pencereler de henüz desteklenmez.
- **Dünya:** Veri 1° × 1°'lik karolar hâlinde okunur: 56°G ile 60°K arasında NASADEM, diğer enlemlerde Copernicus DEM. Her karo 15 yay saniyesine (ekvatorda yaklaşık 463 m, MOLA ızgarasına yakın) ortalama alınarak indirgenir ve tam metre olarak önbelleğe yazılır. Sunucuda bulunmayan bir karo (HTTP 404) açık deniz sayılır. Başka hatalar, örneğin hız sınırlaması, hata olarak bildirilir. Böylece geçici bir ağ sorunu karayı sessizce denize çevirmez.

### 4.2 Destek eşleme

Eğim, üzerinde ölçüldüğü mesafeye bağlıdır. 30 m'lik piksellerden ölçülen eğim, 463 m'lik piksellerden ölçülenden çok daha yüksek çıkar. Notebook'taki örnekte Salda'nın ortalama eğimi Jezero'nunkinin ham piksellerde 3,3 katı, aynı 1 km destekte ise 1,7 katıdır.

Bu yüzden her DEM, analiz ölçeğine (s = 1, 2 veya 4 km) bir Gauss filtresiyle yumuşatılır:

- Hedef yumuşatma genişliği σ = 0,5 × s.
- d genişliğindeki bir kaynak pikselin zaten d²/12 varyansı vardır. Bu yüzden yalnızca kalan kısım uygulanır: σ_uygulanan = √((0,5 s)² − d²/12). Böylece 30 m ve 463 m kaynaklar aynı toplam desteğe ulaşır.
- Filtre boş hücreleri yok sayar. Verinin yarıdan azını gören hücreler boş bırakılır.
- Ölçek, kaynak pikselin en az 1,9 katı olmalıdır; değilse hesap reddedilir. MOLA için en küçük ölçek yaklaşık 880 m olduğundan 1 km uygundur.

### 4.3 Yerel metrik ızgara

Yumuşatılmış DEM, pencere merkezli ve kuzeyi yukarıda olan bir metre ızgarasına örneklenir. Bunun için gövdenin kendi yarıçapıyla azimutal eşit uzaklık izdüşümü kullanılır: Mars 3.396,19 km, Ay 1.737,4 km, Dünya 6.371,007 km. Dünya küre kabul edilir; 60 km içinde bunun hatası %0,4'ün altındadır ve DEM hatasından çok küçüktür.

Hücre aralığı ölçeğe eşittir: 1 km'de 60 × 60, 2 km'de 30 × 30, 4 km'de 15 × 15 hücre. 1 km ızgaranın en az %95'i veriyle dolu değilse pencere geçersiz sayılır. Kalan küçük boşluklar en yakın değerle doldurulur.

### 4.4 Özellikler

Her pencere 21 sayıyla özetlenir. Her ölçekte 5 özellik hesaplanır (3 ölçekte toplam 15). Buna pencerenin bütününü 1 km ızgarada anlatan 6 özellik eklenir.

| Özellik | Nasıl hesaplanır | Ne anlatır | Aile |
|---|---|---|---|
| `slope_mean`, `slope_std`, `slope_p95` | Eğim = arctan(\|∇z\|), merkezi farklarla; ortalama, standart sapma ve 95. yüzdelik | Yüzey ne kadar dik, diklik ne kadar değişken | Eğim |
| `tri` | Her hücrenin 8 komşusuna ortalama mutlak yükseklik farkı (m) | Kısa mesafede engebe | Engebelilik |
| `rough` | Yüksekliğin 3 × 3 ortalamasından sapmalarının standart sapması (m) | Yerel pürüzlülük | Engebelilik |
| `relief` | Yüksekliğin 99. ve 1. yüzdelikleri arasındaki fark (m) | Pencerenin toplam yükselti aralığı | Rölyef |
| `elev_std` | En uygun düzlem çıkarıldıktan sonra yüksekliğin standart sapması (m) | Bölgesel eğimden arındırılmış oynaklık | Rölyef |
| `hypso` | Hipsometrik integral: (ortalama − P1) / (P99 − P1), 0 ile 1 arası | Alanın çoğu alçakta mı yüksekte mi | Biçim |
| `skew` | Düzlemden arındırılmış yüksekliğin çarpıklığı | Negatif: dar ve derin çukurlar; pozitif: tepeler | Biçim |
| `radial_r2` | Pencere merkezi etrafındaki 12 halkanın ortalama profilinin açıkladığı varyans payı | Dairesel simetri; krater, havza ve kubbede yüksek | Biçim |
| `basin_depth` | Dış halkaların (r ≥ 0,75 R) ortalaması eksi merkezin (r ≤ 0,4 R) ortalaması (m) | Pozitifse merkez kenarlardan alçak | Biçim |

Pencereye dik bir düzlem uydurulup çıkarılır, çünkü bölgesel eğim biçim gibi görünmemelidir. Testler, eğimli bir düzlem eklemenin ve yüksekliği topluca kaydırmanın biçim özelliklerini değiştirmediğini doğrular.

## 5. Dünya kataloğu

`scripts/build_catalog.py`, `catalog.py` modülünü çalıştırır:

1. **Pencere merkezleri:** 56°G ile 83°K arası, −179° ile 179° boylamları ve ayrıca McMurdo Kuru Vadileri (Antarktika) taranır. Merkezler yerde 20 km arayla dizilir; boylam adımı enleme göre genişler. Toplam 1.155.795 merkez çıkar.
2. **İndirme:** Gereken 50.801 DEM karosu ve WorldCover karoları 12 paralel bağlantıyla indirilir.
3. **Eleme:** Yüzde 30'dan fazlası deniz olan (yükseklik 0 veya veri yok) ya da %95'ten azı geçerli olan pencereler atılır. Geriye 324.859 pencere kalır.
4. **Özellikler:** Her pencerenin 21 özelliği 8 paralel süreçte hesaplanır. Pencerenin WorldCover sınıf payları da eklenir: ağaç, çalı, çayır, tarım, yerleşim, çıplak zemin, kar, su, sulak alan, mangrov ve yosun. Bu paylar WorldCover'ın 1/64'lük önizlemesinden (yaklaşık 640 m örnekler, pencere başına yaklaşık 10.000 örnek) okunur.
5. **Çıktı:** `data/catalog/earth_global.parquet` dosyasında her satır bir penceredir. Yanındaki `.json` dosyası üretim zamanını, bölgeleri, pencere ve adım boyunu, ölçekleri, özellik listesini, yöntemi, veri kaynaklarını, yazılım sürümlerini, başarısız karoları ve süreleri kaydeder.

Sıfırdan üretim yaklaşık 45 dakika sürer. Önbellek doluysa birkaç dakika yeter. Uygulama özellik listesini kataloğun `.json` dosyasından okur, bu yüzden yeni özellikler eklenmiş bir katalog uygulamada kod değişikliği gerektirmeden çalışır.

15 yay saniyelik ızgaranın etkisi notebook'ta ölçülmüştür: 30 m kaynağa göre fark 1 km ölçekte en çok %4,6, 2–4 km ölçekte çoğunlukla %1,5'in altındadır.

## 6. Benzerlik ve sıralama

Kod `similarity.py` dosyasındadır.

1. **Dönüşüm.** Farkların karşılaştırılabilir anlamlar taşıması için özellikler dönüştürülür. Eğim, TRI, pürüzlülük, rölyef ve `elev_std` için log(1 + x) alınır, çünkü bu büyüklükler çarpımsal yayılır. Havza derinliği için arcsinh(derinlik / 100 m) kullanılır; bu dönüşüm işareti korur, sıfır çevresinde doğrusal, uçlarda logaritmiktir. `hypso`, `skew` ve `radial_r2` olduğu gibi kalır.
2. **Ölçekleme.** Her özellik, bütün Dünya pencerelerinin medyanına göre ortalanır ve çeyrekler arası aralığına (IQR) bölünür. Böylece birim "Dünya'daki tipik yayılım" olur ve aykırı değerler sonucu ele geçiremez. Referansın bir özelliği bu birimde 3'ten büyük çıkarsa uygulama "Dünya'da nadir" uyarısı gösterir.
3. **Ağırlıklar.** Varsayılan aile ağırlıkları şöyledir: eğim 0,30, engebelilik 0,25, rölyef 0,15, biçim 0,30. Her ailenin ağırlığı üyelerine bölünür, bu yüzden çok üyeli bir aile tek başına baskın olamaz. Bir aile içinde 1 km özellikleri yarım ağırlık alır. Bunun nedeni MOLA iz deseni deneyidir: 1 km özellikleri %8–17, 2–4 km özellikleri %1–8 değişiyor. Toplam ağırlık 1'dir. Kullanıcı kaydırıcılarla aile ağırlıklarını değiştirebilir.
4. **Mesafe.** d = √(Σ wᵢ (z_aday,ᵢ − z_referans,ᵢ)²). Bu, Dünya yayılımı biriminde ağırlıklı ortalama karesel farktır. 0 birebir aynı demektir.
5. **Yüzey filtresi.** Mars'ta ve Ay'da bitki örtüsü ve açık su olmadığı için varsayılan olarak çıplak Dünya yüzeyi aranır: çıplak zemin en az %50; ağaç en çok %10, tarım %15, yerleşim %3, su %15, kar %20.
6. **Seçim.** Komşu pencereler birbiriyle örtüştüğü için aynı yer listeyi doldurabilir. Bunu önlemek için en yakından başlanarak seçilmiş bir adaya 60 km'den yakın olan pencereler atlanır.
7. **Her adayın göstergeleri:**
   - **Benzerlik (0–100)** = 100 × e^(−d). Yalnızca okumayı kolaylaştırır, olasılık değildir.
   - **Yüzdelik:** Taranan bütün pencerelerin (filtreden önce) kaçının referansa bu adaydan daha uzak olduğu.
   - **Kararlılık:** Ağırlıklar seçilen değerlerin çevresinde 300 kez rastgele değiştirilir (Dirichlet dağılımı). Adayın bu denemelerin kaçında uygun pencerelerin en yakın %0,1'i içinde kaldığı ölçülür; bu dilim en az 20 penceredir. Hız için en çok 20.000 pencerelik rastgele bir örnek kullanılır. Yüksek kararlılık, sonucun ağırlık seçimine bağlı olmadığını gösterir.
   - **Aile uyumu:** Her ailenin kendi fark değeri d_f, 100 × e^(−d_f) olarak gösterilir.
   - **Özellik katkısı:** wᵢ Δzᵢ² / Σ. Farkın en çok hangi özellikten geldiğini söyler.

## 7. Yüzey ve güvenilirlik katmanları

### 7.1 Sentinel-2 yüzey kontrolü (`surface.py`)

Yalnızca kullanıcı seçili aday için düğmeye bastığında çalışır.

- Earth Search STAC kataloğunda pencerenin merkezini içeren sahneler aranır. Koşullar: Haziran 2024 ile Eylül 2026 arası, bulut %15'ten az, eksik piksel %25'ten az. En az bulutlu 5 sahne denenir. Kullanılabilir alanı %85'i geçen ilk sahne alınır; hiçbiri geçmezse en iyisi kullanılır.
- Pencere için kırmızı, yakın kızılötesi (NIR), SWIR1 (B11), SWIR2 (B12), sahne sınıflandırması (SCL) ve gerçek renkli görüntü 512 × 512 olarak okunur. SCL ile veri yok, doygun, gölge, bulut ve sirrus pikselleri maskelenir.
- Göstergeler:
  - NDVI = (NIR − kırmızı) / (NIR + kırmızı), medyan ve 90. yüzdelik. Yeşil alan payı NDVI > 0,4 olan piksellerdir.
  - NDMI = (NIR − SWIR1) / (NIR + SWIR1). Nem göstergesidir.
  - SWIR B11/B12 oranı. Kil ve karbonat gibi bazı mineraller SWIR2'de soğurma yaptığı için bir ipucu verir, ama tek başına mineral tanımlamaz.
  - Kapsama: pencerenin kullanılabilir piksel payı.

### 7.2 MOLA veri desteği (`mola_support.py`)

MOLA'nın 463 m'lik hücrelerinin çoğu, yörünge izleri arasından enterpolasyonla doldurulmuştur. MEGDR "counts" ürünü her hücreye kaç lazer atımı düştüğünü söyler. Uygulama, seçilen Mars penceresinde gerçek ölçümü olan hücrelerin payını gösterir; Jezero'da bu oran %32,9'dur. 1 km özelliklerinin yarım ağırlık almasının nedeni budur.

### 7.3 Güvenilirlik deneyleri (`validation.py`, notebook)

- **İz deseni deneyi:** Bir Dünya penceresi, gerçek MOLA izlerinin ölçeceği hücrelerden yeniden kurulur; boşluklar doğrusal enterpolasyonla doldurulur. Ardından her özelliğin ne kadar değiştiği ölçülür.
- **Kör arama:** Bu şekilde bozulmuş 399 rastgele Dünya penceresi katalogda aranır. Kendi yerini %85,7 oranında ilk sırada, %95,7 oranında ilk 10'da buluyor. Sorgu ve katalog aynı gezegenden geldiği için bu bir üst sınırdır, gezegenler arası geçerliliği kanıtlamaz.
- **Çözünürlük yanlılığı:** Salda için 30 m ve 15 yay saniyelik kaynaklardan hesaplanan özellikler karşılaştırılır.

### 7.4 Bilinen analoglar (`sites.py`)

Literatürde geçen 15 Dünya analogu, neyin analoğu oldukları notuyla birlikte listelenir: mineraloji, iklim, kimya, volkanizma veya morfoloji. Her biri kendi tam koordinatında aynı yöntemle ölçülür ve seçilen referansa göre mesafesi ve yüzdeliği gösterilir. Bu, morfolojileri nedeniyle seçilmiş analogların kimya veya iklim nedeniyle seçilenlerden daha yakın çıkıp çıkmadığını sınamayı sağlar.

## 8. Uygulamayı kullanma

Canlı sürüm Streamlit Community Cloud'da çalışır. Yerelde `streamlit run app/streamlit_app.py` ile açılır (kurulum için README'ye bakın). Bir alan ilk kez seçildiğinde DEM okunurken bir dakikaya kadar beklemek gerekebilir; sonraki açılışlar önbellekten gelir.

### Sol panel

1. **Gezegen:** Mars veya Ay.
2. **Alan:** Hazır alanlardan biri ya da "Özel koordinat". Mars için Jezero, Gale, Gusev, Meridiani Planum, Oxia Planum ve Holden; Ay için Apollo 11, 15 ve 17 bölgeleri, Copernicus ve Tycho var. Özel koordinatta enlem −80° ile 80°, boylam −180° ile 180° arasında olmalı; doğu pozitiftir.
3. **Eşleşme ağırlıkları:** Dört aile için 0–100 arası kaydırıcılar. Örneğin bir rover rotası için eğim ailesi öne çıkarılabilir.
4. **Dünya yüzeyi filtresi:** Açılıp kapatılabilir; çıplak zemin, ağaç, tarım ve su eşikleri ayarlanabilir.
5. **Aday sayısı:** 3 ile 20 arası.

### Ana sayfa, yukarıdan aşağıya

1. **Geniş görünüm:** Seçilen bölgenin yakınlaştırılabilir haritası. Turuncu kare, karşılaştırılan 60 km'lik penceredir.
2. **Uyarılar:** Mars için gerçek lazer ölçümü payı; referansın Dünya'da nadir görülen özellikleri.
3. **En yakın adaylar:** Sıra, konum, benzerlik, yüzdelik, kararlılık ve çıplak zemin payı tablosu. Yanında Dünya haritası var; uydu görüntüsüne katman menüsünden geçilir.
4. **Karşılaştırma raporu:** Seçilen aday için şunları içerir:
   - özet cümle,
   - iki pencerenin aynı fiziksel ölçekte, kuzeyi yukarıda gölgelendirilmiş rölyefi,
   - aile uyumu çubukları,
   - üç grafik: 1 km'de eğim dağılımı, hipsometrik eğri, radyal profil,
   - her özelliğin farka katkısını gösteren tablo.
5. **Sentinel-2:** Düğmeye basınca sahne, NDVI, NDMI ve SWIR göstergeleri gelir.
6. **Literatürdeki analoglar:** 15 alanın bu referansa uzaklığı.
7. **Bilimsel izlenebilirlik:** Yöntem, referans DEM, katalog tarihi ve parametreleri, başarısız karolar, veri kaynakları, yazılım sürümleri ve sınırlamalar.

### Sonuçları okurken

- Yüksek yüzdelik ve yüksek kararlılık birlikte sağlam bir eşleşme gösterir. Kararlılığı düşük bir aday ağırlık seçimine duyarlıdır.
- Özellik tablosu "neden" sorusunu yanıtlar. Toplam farkın büyük kısmı tek bir özellikten geliyorsa o özelliğe bakın.
- Benzerlik yalnızca arazi biçimini ölçer. Jeoloji, mineraloji ve oluşum süreci aynı olmayabilir.

## 9. Notebook, komut satırı ve Python

**Notebook.** [notebooks/01_jezero_salda.ipynb](../notebooks/01_jezero_salda.ipynb) bilimsel analizin tamamını içerir: veri okuma, aynı ölçekte görünüm, ham ve eşlenmiş eğim karşılaştırması, özellik tablosu, MOLA veri desteği, yanlılık ve kör arama deneyleri, küresel arama, bilinen analoglar, açıklama, yüzey katmanı ve bulgular. `.venv` çekirdeğiyle açılır.

**Katalog üretimi:**

```bash
.venv/bin/python scripts/build_catalog.py --global --output data/catalog/earth_global.parquet
.venv/bin/python scripts/build_catalog.py --region 36 42 26 45 --output data/catalog/turkiye.parquet
```

Seçenekler: `--window-km` (varsayılan 60), `--stride-km` (20), `--workers` (8 süreç), `--threads` (12 bağlantı).

**Testler:** `pytest` ağ gerektirmeyen testleri, `pytest -m network` gerçek uzak verilerle çalışan testleri koşturur. `tests/test_terrain.py` analitik yanıtı bilinen yapay arazileri, `tests/test_similarity.py` doğru cevabı bilinen yapay bir kataloğu, `tests/test_data_layers.py` karo birleştirme, deniz karoları, WorldCover ve MOLA okumayı sınar.

**Python'dan kullanım.** Proje kökünden `PYTHONPATH=src` ile:

```python
from analog_finder.bodies import MARS
from analog_finder.catalog import CATALOG_DIR, load_catalog
from analog_finder.query import window_features
from analog_finder.similarity import Filters, explain, rank_candidates

catalog, meta = load_catalog(CATALOG_DIR / "earth_global.parquet")
reference, _ = window_features(MARS, 18.38, 77.58)  # Jezero
ranking = rank_candidates(reference.values, catalog, meta["feature_names"], filters=Filters.bare_ground())
print(ranking.table[["rank", "lat", "lon", "score", "percentile", "stability"]])
print(explain(ranking, reference.values, ranking.table.iloc[0]).head())
```

## 10. Kod yapısı

| Yer | Görevi |
|---|---|
| `src/analog_finder/bodies.py` | Gövdeler (yarıçap) ve Mars/Ay DEM ürünleri |
| `src/analog_finder/dem.py` | Uzak DEM okuma, Dünya karoları, birleştirme ve önbellek |
| `src/analog_finder/terrain.py` | Destek eşleme, metrik ızgara ve özellikler |
| `src/analog_finder/catalog.py` | Dünya kataloğunun üretimi ve okunması |
| `src/analog_finder/query.py` | Herhangi bir gövdede tek pencerenin özellikleri |
| `src/analog_finder/similarity.py` | Dönüşüm, ölçekleme, ağırlık, sıralama, kararlılık ve açıklama |
| `src/analog_finder/landcover.py` | WorldCover payları |
| `src/analog_finder/surface.py` | Sentinel-2 sahne arama ve göstergeler |
| `src/analog_finder/mola_support.py` | MOLA gerçek ölçüm yoğunluğu |
| `src/analog_finder/validation.py` | İz deseni, kör arama ve yanlılık deneyleri |
| `src/analog_finder/sites.py` | Referans alanları ve bilinen analoglar |
| `src/analog_finder/provenance.py` | Yöntem metni, veri kaynakları ve sürüm kayıtları |
| `app/streamlit_app.py` | Web arayüzü |
| `scripts/build_catalog.py` | Katalog üretim komutu |
| `notebooks/` | Bilimsel analiz |
| `tests/` | Testler |
| `data/catalog/` | Hazır katalog (git'te) |
| `data/cache/` | İndirilen veriler (git'te değil) |

Başlıca kütüphaneler: numpy ve scipy (filtreler, enterpolasyon), pandas ve pyarrow (katalog), rasterio (GDAL üzerinden GeoTIFF ve HTTP aralık okuma), requests (STAC ve PDS istekleri), Streamlit (arayüz), folium (haritalar), plotly (grafikler), matplotlib (notebook).

## 11. Geliştirme

- **Dal düzeni:** `main` dalı canlı uygulamadır. Streamlit Community Cloud her push'ta uygulamayı yeniden kurar. Yeni işler ayrı bir dalda yapılır ve hazır olunca `main`'e birleştirilir.
- **Bağımlılıklar:** Yeni bir paket hem `pyproject.toml` dosyasına hem de sürümü sabitlenmiş olarak `requirements.txt` dosyasına eklenir. Streamlit Cloud ikincisini kullanır.
- **Yeni bir özellik eklemek:**
  1. Hesabı `terrain.py` dosyasına yazın ve adını `PER_SCALE_FEATURES` ya da `WINDOW_FEATURES` listesine ekleyin.
  2. `similarity.py` dosyasında etiketini (`_LABELS`), ailesini (`family_of`) ve dönüşümünü (`transform`) tanımlayın.
  3. Yanıtı bilinen yapay bir araziyle test yazın.
  4. Kataloğu yeniden üretin.
  5. Katalog dosyasını commit edin. GitHub 50 MB üstündeki dosyalar için uyarı verir; dosya başına kesin sınır 100 MB'tır.
- **Yeni bir veri kaynağı eklemek:**
  - Mars veya Ay için `bodies.py` dosyasına yeni bir `DemProduct` eklenir. Ürün, HTTP aralık isteğiyle okunabilen global, basit silindirik bir GeoTIFF olmalıdır. Çözünürlüğün farklı olması sorun değildir, destek eşleme bunu karşılar. Yalnızca en küçük ölçeğin piksel boyunun 1,9 katından büyük olması gerekir.
  - Dünya için `dem.fetch_earth_tile` ve `default_earth_source` güncellenir.
  - Her iki durumda kaynak `provenance.py` dosyasına kaydedilir.
- **Platformlar:** macOS, Linux ve Windows. Komutlar README'de.

## 12. Sınırlamalar

- Benzerlik puanı olasılık değildir ve jeoloji veya mineraloji içermez.
- MOLA, 1 km'den küçük yapıları güvenilir biçimde temsil etmez. Jezero deltası gibi birkaç km'lik yapılar için HRSC gerekir.
- Pencereler 20 km adımla dizilir; bir yapının merkezi en yakın pencere merkezinden 14 km'ye kadar uzakta olabilir.
- Pencere boyu tektir (60 km). Çok daha küçük veya büyük yapılar bu ölçeğe sığmaz; örneğin Manicouagan pencereden büyüktür.
- Mevcut biçim ölçüleri dönmeden bağımsız özet istatistiklerdir. Bu yüzden yöntem dairesel kraterleri değil, "kenarları yüksek, tabanı düz" kapalı havzaları bulur.
- Mars ve Ay'da yerçekimi, atmosfer ve aşınma farklıdır. Benzer biçimler farklı süreçlerden doğmuş olabilir.
- 80°'den yüksek enlemlerde Mars ve Ay için kutup ürünleri henüz eklenmedi.
- Sentinel-2 özeti tek sahneye dayanır ve MGRS karo kenarlarında pencerenin bir kısmını kapsamayabilir.

## 13. Yol haritası

1. **Krater biçimine duyarlı radyal profil.** `feature/radial-profile` dalında devam ediyor. Profil merkezi, pencere merkezinin 15 km çevresinde en simetrik nokta olarak aranır. Böylece 20 km'lik pencere adımının yarattığı kayma giderilir ve kraterler kapalı havzalardan ayrılır.
2. **Birden çok pencere boyutu**, örneğin 15, 30 ve 120 km.
3. **USGS Ay jeoloji haritası** ile jeolojik bağlam.
4. **Seçilmiş Mars bölgeleri için HRSC.** ESA'nın Mars Express HRSC ve MOLA birleşik modeli yaklaşık 200 m çözünürlüktedir. Ardından CRISM ve OMEGA ile mineraloji.

Değerlendirilen diğer kaynaklar: Ay için SLDEM2015 (LOLA ve Kaguya birleşimi, ±60° içinde yaklaşık 59 m), kutup ürünleri, Chang'e-2'nin 20 m'lik Ay DEM'i. Chang'e-2 verisi kayıt ve kullanım koşulları gerektirdiği için şimdilik dışarıda.

## 14. Terimler

| Terim | Anlamı |
|---|---|
| DEM | Sayısal yükseklik modeli: her hücrede bir yükseklik değeri bulunan ızgara |
| COG | Cloud Optimized GeoTIFF; parça parça okunabilen GeoTIFF |
| HTTP aralık isteği | Dosyanın yalnızca istenen bayt aralığını indirmek |
| STAC | Uydu sahnelerini aramak için ortak katalog standardı |
| Destek | Bir değerin temsil ettiği fiziksel alan; yumuşatmanın genişliği |
| Azimutal eşit uzaklık izdüşümü | Merkezden her yöne uzaklıkları doğru koruyan harita izdüşümü |
| TRI | Terrain Ruggedness Index, arazi engebelilik indeksi |
| IQR | Çeyrekler arası aralık: 75. ve 25. yüzdelikler arasındaki fark |
| Hipsometrik integral | Yüksekliğin, en düşük ile en yüksek arasında nerede yoğunlaştığı |
| NDVI, NDMI | Bitki örtüsü ve nem için normalize fark indeksleri |
| SCL | Sentinel-2 sahne sınıflandırması (bulut, gölge, su vb.) |
| MOLA, LOLA | Mars ve Ay'ın lazer altimetre ölçümleri |
| MEGDR | MOLA'nın ızgaralanmış yükseklik ürünü |
| ppd | Derece başına piksel |
| Yay saniyesi | Derecenin 1/3600'ü; Dünya'da ekvatorda yaklaşık 31 m |
