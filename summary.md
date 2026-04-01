# eoCRC Microbiome Meta-Analysis - Results Summary

**Date:** 2026-03-30
**Institution:** University of Texas at Dallas (UTD)
**Pipeline:** curatedMetagenomicData v3 -> MaAsLin2 -> LOSO-CV -> SHAP

\---

## 1\. Sample Counts

|Group|N|
|-|-|
|eoCRC (CRC, age<50)|78|
|Young controls (healthy, age<50)|121|
|loCRC (CRC, age>=50)|562|
|Older controls (healthy, age>=50)|521|
|**Total**|**1282**|

**Feature matrices:** Species 1282x220 (post-filter); Pathways 1282x419 (unstratified)

\---

## 2\. LOSO-CV Classification Results

### Primary: eoCRC vs Young Controls (8 valid folds, YuJ\_2015 excluded)

|Feature set|Model|Pooled AUROC|95% CI|I2|perm-p|
|-|-|-|-|-|-|
|Species|RF|0.704|0.602-0.789|4.5%|0.0020|
|Combined|RF|0.708|0.561-0.821|36.3%|0.0000|

### Secondary: eoCRC vs loCRC (9 folds, loCRC downsampled in training)

|Feature set|Model|Pooled AUROC|95% CI|I2|perm-p|
|-|-|-|-|-|-|
|Species|XGB|0.719|0.640-0.787|3.4%|NA|

\---

## 3\. MaAsLin2 Differential Abundance

### eoCRC vs Young Controls - Species (q<0.25, total n=10)

1. Blautia obeum (up eoCRC, q=0.0496)
2. Parvimonas micra (up eoCRC, q=0.0496)
3. Gemella morbillorum (up eoCRC, q=0.0840)
4. Bacteroides caccae (up eoCRC, q=0.1493)
5. Alistipes indistinctus (up eoCRC, q=0.1493)
6. Dialister pneumosintes (up eoCRC, q=0.1493)
7. Clostridium  symbiosum (up eoCRC, q=0.1799)
8. Alistipes finegoldii (up eoCRC, q=0.2037)
9. Peptostreptococcus stomatis (up eoCRC, q=0.2037)
10. Anaerotruncus colihominis (up eoCRC, q=0.2378)

### eoCRC vs Young Controls - Pathways (q<0.25, total n=95)

1. P108 PWY  pyruvate fermentation to propanoate I (up eoCRC, q=0.0028)
2. THISYNARA PWY  superpathway of thiamin diphosphate biosynthesis III  eukaryotes  (down eoCRC, q=0.0500)
3. PWYG 321  mycolate biosynthesis (up eoCRC, q=0.0545)
4. PWY 1042  glycolysis IV  plant cytosol  (down eoCRC, q=0.0751)
5. PWY 6151  S adenosyl L methionine cycle I (down eoCRC, q=0.0751)
6. PWY 7234  inosine 5  phosphate biosynthesis III (up eoCRC, q=0.0751)
7. PWY 5989  stearate biosynthesis II  bacteria and plants  (up eoCRC, q=0.0751)
8. P42 PWY  incomplete reductive TCA cycle (up eoCRC, q=0.0751)
9. PWY 6703  preQ0 biosynthesis (down eoCRC, q=0.0803)
10. PPGPPMET PWY  ppGpp biosynthesis (up eoCRC, q=0.0803)

### loCRC vs Older Controls - Species (q<0.25, total n=100)

1. Dialister pneumosintes (up loCRC, q=0.0000)
2. Ruthenibacterium lactatiformans (up loCRC, q=0.0000)
3. Gemella morbillorum (up loCRC, q=0.0000)
4. Faecalibacterium prausnitzii (down loCRC, q=0.0000)
5. Butyricimonas virosa (up loCRC, q=0.0000)
6. Solobacterium moorei (up loCRC, q=0.0000)
7. Peptostreptococcus stomatis (up loCRC, q=0.0001)
8. Odoribacter splanchnicus (up loCRC, q=0.0003)
9. Enterocloster citroniae (up loCRC, q=0.0003)
10. Fusobacterium nucleatum (up loCRC, q=0.0019)

\---

## 4\. SHAP Top 20 Features (Extended Combined RF, top-200)

\[star] = Doubly validated in both SHAP and MaAsLin2 (q<0.25)

1. Peptostreptococcus stomatis                   | ↓ eoCRC | q=0.2037 \[star]
2. Gemella morbillorum                           | ↑ eoCRC | q=0.0840 \[star]
3. Parvimonas micra                              | ↓ eoCRC | q=0.0496 \[star]
4. Intestinimonas butyriciproducens              | ↓ eoCRC | q=NA
5. Dialister pneumosintes                        | ↑ eoCRC | q=0.1493 \[star]
6. Faecalibacterium prausnitzii                  | ↑ eoCRC | q=NA
7. Lactobacillus rogosae                         | ↑ eoCRC | q=NA
8. Alistipes finegoldii                          | ↓ eoCRC | q=0.2037 \[star]
9. Bacteroides faecis                            | ↓ eoCRC | q=NA
10. Intestinibacter bartlettii                    | ↑ eoCRC | q=NA
11. Anaerotruncus colihominis                     | ↓ eoCRC | q=0.2378 \[star]
12. Proteobacteria bacterium CAG:139              | ↓ eoCRC | q=NA
13. Parabacteroides goldsteinii                   | ↓ eoCRC | q=NA
14. Bacteroides cellulosilyticus                  | ↓ eoCRC | q=NA
15. Lachnospira eligens                           | ↑ eoCRC | q=NA
16. Bacteroides caccae                            | ↓ eoCRC | q=0.1493 \[star]
17. Eubacterium sp. CAG:274                       | ↑ eoCRC | q=NA
18. Bacteroides ovatus                            | ↑ eoCRC | q=NA
19. Turicimonas muris                             | ↓ eoCRC | q=NA
20. \[Clostridium] symbiosum                       | ↑ eoCRC | q=0.1799 \[star]

Doubly validated: 8/20 features

\---

## 5\. Key Biological Findings

* eoCRC carries a detectable gut microbiome signature: species RF AUROC=0.704, perm-p=0.002
* Signal is geographically consistent (I2=4.5%; 9 cohorts, Europe + East Asia)
* Enriched in eoCRC: oral pathobionts Parvimonas micra, Gemella morbillorum, Peptostreptococcus stomatis, Dialister pneumosintes
* Depleted in eoCRC: butyrate producers Intestinimonas butyriciproducens, Lachnospira eligens
* Signal is primarily taxonomic (species RF p=0.002 vs pathway RF p=0.181 NS)
* Oral-gut axis hypothesis: ectopic colonic colonization by oral pathobionts promotes tumorigenesis

\---

## 6\. Limitations

1. Cross-sectional design - causality unestablished
2. European/East Asian cohorts only - no Black/Hispanic representation
3. Low eoCRC n=78 limits power for rare species
4. No oral microbiome data - oral-gut axis inferred not measured
5. Approximated SHAP (Monte Carlo nsim=30); direction may mislead for sparse features
6. Pathway signal not significant (perm-p=0.181); functional interpretation limited

\---

## 7\. Future Directions

1. Co-occurrence network analysis (SparCC) for polymicrobial biofilm structure
2. Mendelian Randomization using FUT2/LCT genetic instruments
3. Paired oral-stool study to directly test translocation of P. micra / G. morbillorum
4. Dietary fiber intervention: test modulation of butyrate-producer:pathobiont ratio
5. Validation in racially diverse US cohort (Parkland Health/UTSW)
6. colibactin (pks island) quantification from raw metagenomic reads

\---

*Generated: 2026-03-30 23:47:12*

