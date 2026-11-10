# Katkı rehberi

Planetary Analog Finder'a katkı vermek istediğiniz için teşekkürler. Bu belge, bir değişikliğin fikirden canlı uygulamaya nasıl gittiğini anlatır.

Başlamadan önce okuyun:

- [README](README.md): kurulum (macOS, Linux, Windows) ve ilk sonuçlar.
- [docs/nasil-calisir.md](docs/nasil-calisir.md): yazılımın nasıl çalıştığı, hangi veriyi ne için kullandığı ve hesapların ayrıntısı.
- Canlı uygulama: https://mars-dunya-eslestirme.streamlit.app

## En önemli kural: `main` canlı uygulamadır

`main` dalına giren her değişiklik birkaç dakika içinde canlı uygulamaya yansır. Bu yüzden `main`'e doğrudan push edilmez; her değişiklik bir pull request (PR) üzerinden gelir. GitHub bunu zorunlu tutar.

## Çalışma akışı

### 1. Repoya erişim

- **Ekip üyeyseniz** proje sahibine GitHub kullanıcı adınızı verin. Size collaborator daveti gelir; e-postadaki bağlantıdan ya da github.com/notifications sayfasından kabul edin.
- **Ekip dışındaysanız** repoyu fork edin ve kendi kopyanızda çalışın. PR'ı fork'unuzdan açarsınız.

### 2. Dal açın

Güncel `main`'den, ne yaptığını anlatan kısa bir adla yeni bir dal açın:

```bash
git switch main
git pull
git switch -c feature/pencere-boyutlari
```

Dal adı önekleri: `feature/` (yeni özellik), `fix/` (hata düzeltme), `docs/` (yalnızca belge), `data/` (veri kaynağı veya katalog).

### 3. Değişikliği yapın

Bir PR tek bir işi yapsın. Bir özellik ile ilgisiz bir düzeltme aynı PR'a girmesin; küçük PR'lar daha hızlı incelenir ve sorun çıkarsa kolay geri alınır.

### 4. PR açmadan önce kontrol listesi

- [ ] **Testler geçiyor:** `pytest`. Veri okuma kodunu (`dem.py`, `landcover.py`, `mola_support.py`, `surface.py`) değiştirdiyseniz `pytest -m network` de.
- [ ] **Uygulama açılıyor:** `streamlit run app/streamlit_app.py` ile yerelde bir Mars ve bir Ay bölgesi deneyin.
- [ ] **Yeni hesap test edildi:** Yeni bir özellik ya da hesap için yanıtı önceden bilinen yapay bir araziyle test yazın (örnekler `tests/test_terrain.py` içinde).
- [ ] **Özellik değiştiyse katalog yenilendi:** `scripts/build_catalog.py --global` ile yeniden üretin ve `data/catalog/` dosyalarını commit edin. GitHub 50 MB üstündeki dosyalar için uyarı verir, dosya başına kesin sınır 100 MB'tır; katalog bu sınıra yaklaşıyorsa önce konuşalım.
- [ ] **Sonuçlar değiştiyse belgeler güncellendi:**
  - notebook'u yeniden çalıştırın,
  - README'deki ve rehberdeki sayıları düzeltin,
  - rehberin şekillerini `scripts/make_doc_figures.py` ile yeniden üretin.
- [ ] **Bağımlılık eklendiyse iki yerde:** hem `pyproject.toml` hem de sürümü sabitlenmiş olarak `requirements.txt`. Streamlit Cloud ikincisini kullanır.

### 5. PR açın

PR açıklamasında şunları yazın:

- ne değişti ve neden,
- nasıl denediniz,
- sonuçlar değiştiyse eski ve yeni sayılar.

Onay zorunlu değildir, ama birleştirmeden önce bir ekip arkadaşının göz atması önerilir. Bilimsel sonucu değiştiren her PR'ı en az bir kişi incelemeli.

### 6. Birleştirme ve kontrol

PR'ı "Squash and merge" ile birleştirin; böylece `main`'in geçmişi her PR için tek ve okunaklı bir commit olarak kalır. Birleştirdikten birkaç dakika sonra canlı uygulamayı açıp çalıştığını kontrol edin. Uygulama hata verirse proje sahibine haber verin; Streamlit'teki "Reboot app" yalnızca onun hesabından yapılabilir.

## Bilimsel kurallar

Projenin değeri sonuçlarının güvenilir olmasında. Bu yüzden:

- **Gövdeler arasında mutlak yükseklik, koordinat ya da ham piksel karşılaştırmayın.** Her veri önce aynı fiziksel desteğe getirilir (`terrain.resample_to_scale`).
- **Yeni bir veri kaynağı:**
  - mümkünse hesap gerektirmeden erişilebilmeli,
  - lisansı açık olmalı,
  - atfı README'nin "Veri atıfları" bölümüne eklenmeli,
  - `provenance.py` dosyasına kaydedilmeli.
- **Olumsuz sonuçları da yazın.** "Salda, Jezero'ya topografik olarak benzemiyor" ya da "kaymış sorguda kör arama isabeti %20" gibi bulgular projenin en değerli çıktılarındandır; gizlenmez, belgelenir.
- **Benzerlik puanı olasılık değildir.** Arayüzde ve belgelerde öyleymiş gibi anlatmayın.

## Kod ve yazım düzeni

- Python 3.11 veya üstü.
- Kod tanımlayıcıları, docstring'ler ve kod yorumları İngilizce; arayüz metinleri, belgeler ve commit mesajları Türkçe.
- Commit mesajı: ilk satırda kısa bir özet, gerekirse boş bir satırdan sonra neden ve nasıl.
- `data/cache/` git'e girmez; indirilen veriler oradan her zaman yeniden üretilebilir.

## Hata ve fikirler

Bir hata bulduysanız ya da bir fikriniz varsa GitHub'da bir issue açın. Hata için şunları yazın:

- hangi bölge ve ayarlarla oluştuğu,
- ne beklediğiniz ve ne gördüğünüz,
- varsa hata mesajı.

## E-posta gizliliği

Commit'lerinizde kişisel e-posta adresiniz herkese açık görünür. Gizlemek isterseniz GitHub → Settings → Emails sayfasında "Keep my email addresses private" ayarını açın ve orada gösterilen `noreply` adresini bu repoda kullanın:

```bash
git config user.email "KULLANICI_ID+kullanici-adi@users.noreply.github.com"
```

## Lisans

Katkılarınız projenin [MIT lisansı](LICENSE) altında yayımlanır.
