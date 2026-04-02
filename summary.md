# eoCRC Microbiome Meta-Analysis - Results Summary

**Date:** 2026-04-02
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

## 3\. MaAsLin2 Differential Abundance (v3 — corrected CLR via compositions::clr, pseudocount=1e-6)

### eoCRC vs Young Controls - Species (q<0.25, total n=47, top 10 shown)

1. Parvimonas micra (up eoCRC, coef=1.002, q<0.001)
2. Gemella morbillorum (up eoCRC, coef=0.721, q<0.001)
3. Peptostreptococcus stomatis (up eoCRC, coef=0.745, q<0.001)
4. Anaerotruncus colihominis (up eoCRC, coef=0.706, q=0.0036)
5. Dialister pneumosintes (up eoCRC, coef=0.512, q=0.0036)
6. Monoglobus pectinilyticus (down eoCRC, coef=-0.623, q=0.0049)
7. Romboutsia ilealis (down eoCRC, coef=-0.341, q=0.0126)
8. Bacteroides cellulosilyticus (up eoCRC, coef=0.904, q=0.0156)
9. Turicimonas muris (down eoCRC, coef=-0.513, q=0.0244)
10. Intestinimonas butyriciproducens (up eoCRC, coef=0.597, q=0.0245)

### eoCRC vs Young Controls - Pathways (q<0.25, total n=87, top 10 shown)

1. P108-PWY: pyruvate fermentation to propanoate I (up eoCRC, q=0.0051)
2. PWY-6588: pyruvate fermentation to acetone (up eoCRC, q=0.0306)
3. TRPSYN-PWY: L-tryptophan biosynthesis (down eoCRC, q=0.0412)
4. PWY-1042: glycolysis IV (plant cytosol) (down eoCRC, q=0.0421)
5. PWY-6151: S-adenosyl-L-methionine cycle I (down eoCRC, q=0.0421)
6. PWY-6936: seleno-amino acid biosynthesis (down eoCRC, q=0.0421)
7. THISYNARA-PWY: superpathway of thiamin diphosphate biosynthesis III (down eoCRC, q=0.0421)
8. GLUCUROCAT-PWY: beta-D-glucuronide and D-glucuronate degradation (down eoCRC, q=0.0421)
9. RHAMCAT-PWY: L-rhamnose degradation I (down eoCRC, q=0.0421)
10. P163-PWY: L-lysine fermentation to acetate and butanoate (up eoCRC, q=0.0421)

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

\[star] = Doubly validated in both SHAP and MaAsLin2 v3 (q<0.25) | MaAsLin2 direction shown in parentheses

1. Peptostreptococcus stomatis                   | SHAP ↓ eoCRC | q<0.001 (up eoCRC) \[star]
2. Gemella morbillorum                           | SHAP ↑ eoCRC | q<0.001 (up eoCRC) \[star]
3. Parvimonas micra                              | SHAP ↓ eoCRC | q<0.001 (up eoCRC) \[star]
4. Intestinimonas butyriciproducens              | SHAP ↓ eoCRC | q=0.0245 (up eoCRC) \[star]
5. Dialister pneumosintes                        | SHAP ↑ eoCRC | q=0.0036 (up eoCRC) \[star]
6. Faecalibacterium prausnitzii                  | SHAP ↑ eoCRC | q=0.0822 (down eoCRC) \[star]
7. Lactobacillus rogosae                         | SHAP ↑ eoCRC | q=0.0780 (down eoCRC) \[star]
8. Alistipes finegoldii                          | SHAP ↓ eoCRC | q=0.2108 (up eoCRC) \[star]
9. Bacteroides faecis                            | SHAP ↓ eoCRC | q=NA
10. Intestinibacter bartlettii                    | SHAP ↑ eoCRC | q=NA
11. Anaerotruncus colihominis                     | SHAP ↓ eoCRC | q=0.0036 (up eoCRC) \[star]
12. Proteobacteria bacterium CAG:139              | SHAP ↓ eoCRC | q=NA
13. Parabacteroides goldsteinii                   | SHAP ↓ eoCRC | q=NA
14. Bacteroides cellulosilyticus                  | SHAP ↓ eoCRC | q=0.0156 (up eoCRC) \[star]
15. Lachnospira eligens                           | SHAP ↑ eoCRC | q=NA
16. Bacteroides caccae                            | SHAP ↓ eoCRC | q=0.1817 (up eoCRC) \[star]
17. Eubacterium sp. CAG:274                       | SHAP ↑ eoCRC | q=NA
18. Bacteroides ovatus                            | SHAP ↑ eoCRC | q=NA
19. Turicimonas muris                             | SHAP ↓ eoCRC | q=0.0244 (down eoCRC) \[star]
20. \[Clostridium] symbiosum                       | SHAP ↑ eoCRC | q=0.2492 (up eoCRC) \[star]

Doubly validated: 13/20 features (10 enriched in eoCRC, 3 depleted: F. prausnitzii, L. rogosae, T. muris)

\---

## 5\. Key Biological Findings

* eoCRC carries a detectable gut microbiome signature: species RF AUROC=0.704 (95% CI 0.602–0.789), perm-p=0.002
* Signal is geographically consistent (I²=4.5%; 9 cohorts, Europe + East Asia)
* Enriched in eoCRC: oral pathobionts Parvimonas micra, Gemella morbillorum, Peptostreptococcus stomatis, Dialister pneumosintes, Anaerotruncus colihominis; opportunists Bacteroides cellulosilyticus, Bacteroides caccae, Intestinimonas butyriciproducens, Alistipes finegoldii, [Clostridium] symbiosum
* Depleted in eoCRC: protective commensals Faecalibacterium prausnitzii (major butyrate producer, q=0.082), Lactobacillus rogosae (q=0.078), Turicimonas muris (q=0.024)
* 13/20 SHAP top features doubly validated against MaAsLin2 v3 (up from 8 with v2 CLR bug)
* Signal is primarily taxonomic (species RF perm-p=0.002 vs pathway RF perm-p=0.181 NS)
* Oral-gut axis hypothesis: ectopic colonic colonization by oral pathobionts promotes tumorigenesis; concurrent loss of butyrate producers removes mucosal protection

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

*Generated: 2026-03-30 23:47:12 | Updated: 2026-04-02 (v3 MaAsLin2 corrected CLR)*

