Projekt vypracovala skupina Boháč, Bubeník, Koudelka

Tento projekt se zabývá numerickým řešením stacionárního trojrozměrného vedení tepla v homogenním tělese, ve kterém je umístěn vnitřní tepelný zdroj. Cílem je simulovat ustálené rozložení teploty v doméně, tedy takový stav, kdy se teplota v čase již nemění.

**Matematický model**
Fyzikálním základem úlohy je Poissonova rovnice pro stacionární vedení tepla:
<img width="331" height="95" alt="image" src="https://github.com/user-attachments/assets/d6c62021-0e05-45f2-ab23-0ff31953c8d4" />

V našem modelu je vnitřní zdroj umístěn do definované centrální oblasti tělesa.
Na povrchu tělesa jsou uvažovány konvekční (Robinovy) okrajové podmínky, které popisují přestup tepla mezi povrchem a okolním prostředím. Tím je simulována situace, kdy těleso není tepelně izolované, ale odevzdává teplo do okolí.

**Numerický postup**
Prostorová doména je rozdělena na pravidelnou 3D mřížku. V každém jejím uzlu se numericky aproximuje Laplaceův operátor pomocí metody konečných rozdílů.
Výsledný diskrétní systém rovnic je řešen pomocí Gauss–Seidelovy iterační metody, která postupně aktualizuje odhad teploty ve všech vnitřních uzlech mřížky. Výpočet končí v okamžiku, kdy změna teploty mezi dvěma iteracemi klesne pod zadaný práh tolerance.
Součástí řešení je i výpočet rezidua (residuals.csv), které umožňuje sledovat průběh konvergence
<k vykreslení reziduí slouží gnuplot: **gnuplot ./plot_residuals.gnuplot**>

**Geometrie a fyzikální parametry**
Těleso má definované rozměry a počty uzlů v jednotlivých směrech. Z těchto parametrů se automaticky odvodí krok sítě (dx, dy, dz). Každý uzel obsahuje:
                                                                                                                                                - Teplotu
                                                                                                                                                - Případnou hodnotu zdroje
                                                                                                                                                - Fyzikální vlastnosti materiálu

Konvekční koeficienty mohou být různé na různých stranách tělesa, což umožňuje simulovat například kontakt se vzduchem, podložkou nebo jiným prostředím.

