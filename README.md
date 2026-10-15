# Identifikace potenciálních chyb VMB pomocí spektrální odlehlosti
Spolehlivost dat z mapování biotopů, publikovaných ve formě vrstvy mapování biotopů (VMB), je komplikována jejich velkým datovým a geografickým rozsahem, který prakticky znemožňuje vizuální přehled a zpětnou kontrolu v terénu. Navíc data v publikované vrstvě jsou až 15 let stará, takže od jejich pořízení často došlo ke změnám v pokryvu, jejichž odhalení dosud vyžadovalo nové terénní mapování.
Jednou z možných metod pro odhalení potenciálních atributových či polohových anomálií, ať už způsobených neaktuálností nebo technickou chybou, je využití spektrální odlehlosti. Princip spočívá v identifikaci míst, která mají odlišný spektrální projev oproti segmentům téže biotopové kategorie, a to na základě analýzy družicových dat. Tento postup představuje optimalizovanou a inovovanou analýzu spektrální odlehlosti vyvinutou v rámci projektu Sentimap (Šímová et al. 2023; Moravec et al. 2023).

Analýza je provedena na území celého Česka. Z původního datasetu VMB jsou odstraněny segmenty s kategoriemi *X, –1 a moz.* Na zbývající segmenty je aplikován negativní buffer, tj. geometrie je zmenšena směrem dovnitř o 15 metrů. Důvodem tohoto kroku je minimalizace vlivu krajních pixelů a zároveň odstranění velmi malých segmentů. Pro každý segment jsou spočítány průměr, medián a směrodatná odchylka. Analýza je prováděna samostatně pro každou biotopovou kategorii — například segmenty typu L3.1 jsou hodnoceny výhradně vzhledem k ostatním segmentům kategorie L3.1 z celého území Česka.

Aplikace zahrnuje dvě analytické úrovně:

- **Polygonovou analýzu**, kde jsou využity dvě komplementární metody — robustní Mahalanobisova vzdálenost, která měří odlehlost pozorování od centra typických hodnot při zohlednění jejich rozptylu, a lokální faktor odlehlosti (LOF), který detekuje body v řídkých oblastech vícerozměrného prostoru. Vstupem jsou agregované statistiky jednotlivých polygonů (průměr, medián a směrodatná odchylka).

- **Pixelovou mapu odlehlosti**, která je založena na výpočtu robustní Mahalanobisovy vzdálenosti jednotlivých pixelů vůči vzorku pixelů dané biotopové kategorie. Vstupem jsou přímo hodnoty jednotlivých pixelů.

## Struktura repozitáře
```plaintext
spektralni-odlehlost/
├── data/
│   ├── s1_mosaic/            # Podpůrná vrstva dlaždic ČR pro stahování S1 (tiles_czechia.*)
│   └── s2_mosaic/            # Hranice ČR pro oříznutí S2 (czechia.geojson)
├── scripts/
│   ├── analysis/             # R skripty pro vlastní analýzu
│   └── s1_mosaic/            # Jupyter Notebook pro běh v Google Colab
└── README.md
```
## Data a jejich získání

**VMB:** Stáhněte aktuální vrstvu *Aktualni_Segment* z agenturního disku K.
   
**Sentinel-1:** Pro získání mozaik využijte skript `scripts/s1_mosaic/S1_tiles.ipynb`, který je určený pro běh v prostředí Google Colab. K běhu vyžaduje vrstvu dlaždic uloženou ve složce `data/s1_mosaic/tiles_czechia.*`. Po dokončení výpočtu budou výsledné měsíční mediánové mozaiky uloženy na vašem Google Disku, odkud je stáhněte.

**Sentinel-2:** Data získáte prostřednictvím platformy Sentinel-2 Global Mosaic na adrese https://s2gm.land.copernicus.eu/mosaic-hub. Pro oříznutí zájmového území nahrajte soubor s hranicemi Česka ze složky `data/s2_mosaic/czechia.geojson` a stáhněte bezešvou medoidovou mozaiku pokrývající vegetační období od 1. května do 31. srpna pro všechna dostupná pásma s prostorovým rozlišením 10 m/px a projekcí UTM.

## Nastavení lokálního prostředí
Pro samotný výpočet si na lokálním disku připravte pracovní adresář. V něm vytvořte složku `data/`, do které nakopírujte stažená družicová data a vrstvu VMB. Dále vytvořte složku `spektralni_odlehlost_R/`, kam nakopírujete všechny *.R* skripty. Vše podle následující struktury: 

```plaintext
spektralni_odlehlost/
├── data/
│   ├── S1_may_asc/            
│   ├── S1_may_desc/           
│   ├── S1_june_asc/          
│   ├── S1_june_desc/         
│   ├── S1_july_asc/           
│   ├── S1_july_desc/         
│   ├── S1_aug_asc/           
│   ├── S1_aug_desc/          
│   ├── s2gm_copernicus/      
│   └── VMB/                  
└── spektralni_odlehlost_R/    # Zdrojové R skripty
    ├── main.R                 # Hlavní spouštěcí skript
    ├── vmb_filt.R
    ├── sentinel_to_stack.R
    ├── compute_stats.R
    ├── polygon_error.R
    └── pixel_error.R
```

## Spuštění analýzy

1. Ujistěte se, že jsou všechna vstupní data stažena a umístěna v příslušných podsložkách uvnitř `data/` podle schématu výše.
2. Otevřete RStudio a nastavte pracovní adresář do kořenové složky projektu.
3. Otevřete skripty ze složky `spektralni_odlehlost_R/`.
4. Spusťte skript `main.R`, který zajistí celý proces analýzy: od filtrace VMB vrstvy přes tvorbu rastrového stacku až po finální výpočet polygonové a pixelové odlehlosti.
5. V kořenové složce projektu vzniknou automaticky dvě nové složky `interim/` pro ukládání meziproduktů výpočtů a `outputs/` pro výsledné vrstvy.
6. Výsledné vrstvy budou po dokončení analýzy připraveny ve složce `outputs/` pro další práci v GIS.

## Interpretace výsledků

Výstupy polygonové a pixelové analýzy je vhodné vnímat jako komplementární. Vzhledem k odlišnostem ve statistickém přístupu se nemusí v konkrétních místech stoprocentně překrývat. Například segment vyhodnocený jako odlehlý nemusí být tvořen odlehlými pixely a naopak shluky odlehlých pixelů se mohou nacházet i uvnitř polygonu, který jako celek kritéria odlehlosti nesplnil. Metoda obecně dobře funguje při detekci výrazných zásahů do území (např. těžba dřeva, rozšiřování lomů či nová výstavba) a spolehlivě odhaluje také technické chyby ve vstupních datech, například nechtěné sloučení dvou odlišných biotopů do jednoho segmentu.

## Zdroje


- *Šímová, P., Barták, V., Moravec, D., et al. (2023) Možnosti aktualizace Vrstvy mapování biotopů s využitím dálkového průzkumu Země. Metodika, Katedra prostorových věd, Fakulta životního prostředí, Česká zemědělská univerzita v Praze.*
 
- *Moravec, D., Barták, V., Šímová P., (2023) Finding Misclassified Natura 2000 Habitats by Applying Outlier Detection to Sentinel-1 and Sentinel-2 Data, Remote Sensing 15, no. 18: 4409. https://doi.org/10.3390/rs15184409*

--- 

V případě dotazů nebo problémů se spuštěním mě kontaktujte na e-mailu: vit.jezek@aopk.gov.cz
