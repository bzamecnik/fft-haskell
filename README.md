# fft-haskell

FFT implementation and a simple spectral analyzer in Haskell

This is the original project which inspirated [HarmonEye](http://harmoneye.com), a real-time
music pitch analyzer and visualizer several years later.

This was a semester project in [Non-procedural Programming](http://is.cuni.cz/studium/predmety/index.php?do=predmet&kod=NPRG005) course
at the [Faculty of Mathematics and Physics of Charles University in Prague](http://www.mff.cuni.cz/).

Created: April 2008

[Bohumír Zámečník](http://zamecnik.me)

Original README in Czech:

## FFT a jeho aplikace (spektralní analyzér)

Použitý jazyk: Haskell

### Implementované algoritmy:

* Rekurzivní FFT (Fast Fourier Transform)
	* (dopředná a inverzní) na seznamech komplexních čísel
* Iterativní FFT
	* (dopředná a inverzní) na polích komplexních čísel
	* (dopředná a inverzní) na seznamech komplexních čísel 
* Naivní DFT (Discrete Fourier Transform) podle vzorce
	* (dopředná a inverzní) na seznamech

### Praktická aplikace algoritmů:

#### Spektrální analyzér signálu

Ze vstupního signálu pomocí FFT spočítá a pak zobrazí frekvenční spektrum.

Na vstupu dostane signál jako seznam reálných čísel. Čísla bere po oknech zvolené velikosti. Jednotlivá okna prožene přes FFT a ze získaných komplexních čísel spočítá výkon pro každé pásmo, a tím získá frekvenční spektrum. To poté přeškáluje do logaritmického měřítka a podle zadané velikosti okna. Spektrogram je pak pomocí sloupců hvězdiček zobrazen na výstup.

Původně bylo v plánu použít pro vstup signálu zvukové soubory nebo zařízení. Buď pomocí SDL-Audio, OpenAL nebo ALSA. Ukázalo se ale, že bindingy těchto knihoven jsou v Haskellu zatím na takové úrovni, že vstup surového zvukového signálu nebyl možný. SDL a OpelAL sice umožňují pracovat se zvukovými soubory či zařízeními, ale není možné sahat přímo na holý signál. Bining pro ALSU by to měl umožnovat, ale bohužel je zatím v alpha stavu, balíček neexistuje a ani přes usilovnou snahu toto knihovnu nainstalovat se mi to nepodařilo. Proto jsem se rozhodl udělat pouze textový vstup a neinteraktivní zpracování.  

### Použití

#### Algoritmy

##### Rekurzivní FFT

* `fft` - FFT seznamu komplexních čísel rekurzivně
* `ifft` - inverzní FFT seznamu komplexních čísel rekurzivně
* `fftR` - FFT seznamu reálných čísel rekurzivně
* `ifftR` - inverzní FFT seznamu reálných čísel rekurzivně

##### Iterativní FFT
 
* `fftiterArray` - FFT pole komplexních čísel iterativně
* `ifftiterArray` - inverzní FFT pole komplexních čísel iterativně
* `fftiterList` - FFT seznamu komplexních čísel iterativně
* `ifftiterList` - inverzní FFT seznamu komplexních čísel iterativně
* `fftiter` - zkratka za fftiterList
* `ifftiter` - zkratka za ifftiterList

##### DFT

* `dft` - DFT jednoduše podle vzorce, ale nepříliš efektivně
* `idft` - inverzní DFT podle vzorce

#### Spektrální analyzér

`processSignal` - spektrální analyzér

* parametry
	* velikost okna
	* vstupní signál jako seznam reálných čísel
* výstup - spektrogram na stdout

#### Testovací sada

Viz poslední sekce zdrojového kódu. Je tam mnoho testovacích funkcí pro jednotlivé části programu. Jejich název začíná prefixem 'test'.

Více informací je možno nalézt v komentářích k jednotlivým funkcím.

#### TODO aneb, co už se nevešlo:

* Iterativní FFT na polích moc žere paměť. Např. ve WinHugsu s 7 MB heapu spolkne seznam délky nejvýš 256, pak už se nevejde.
	* ŘEŠENÍ: iterativní FFT na listech 
* Vstup přes stdin nebo nejlépe přes zvukovou knihovnu (závisí na problému s bindingy).
* Opravdu grafický výstup.
* Interaktivní ovládání.

## Použitá literatura

Cormen, Leiserson, Rivest, Stein: Introduction to Algorithms, 2nd Ed, MIT Press 2001
