# Supplementary Table S2: Zebrafish Heart Layer Estimation Rationale

**System:** Zebrafish (*Danio rerio*) ventricular apex regeneration after resection  
**Outcome:** Complete regeneration via cardiomyocyte dedifferentiation and proliferation  
**Key Reference:** Jopling et al., Nature 2010 (DOI: 10.1038/nature09391)

## Layer Estimates

| Layer | Value | Justification |
|-------|-------|---------------|
| **P** | 0.80 | Gata4-dependent enhancer activation opens chromatin at regenerative genes; cardiomyocytes retain epigenetic plasticity absent in mammals |
| **C** | 0.70 | Positional identity maintained but less critical than in appendage regeneration; heart regenerates apex without precise patterning |
| **S** | 0.85 | Coordinated FGF, PDGF, Notch, and Wnt signaling; epicardial-derived factors synchronize cardiomyocyte proliferation |
| **B** | 0.75 | Transient fibrin clot rapidly resolved; TGF-β3 present; no permanent scar formation unlike mammalian myocardial infarction |
| **T** | 0.80 | Cardiomyocytes re-enter cell cycle safely; no evidence of malignant transformation; controlled proliferation terminates when mass restored |
| **M** | 0.70 | Metabolic shift from fatty acid oxidation to glycolysis supports proliferation; less extreme than liver regeneration |

## Calculation

Using ARC v4.0 weights:
Φ = σ(1.45×0.80 + 0.63×0.70 + 1.42×0.85 + 1.78×0.75 + 1.20×0.80 + 0.93×0.70 + 0.85×(0.85×0.75) − 4.52)
Φ = σ(1.16 + 0.441 + 1.207 + 1.335 + 0.96 + 0.651 + 0.542 − 4.52)
Φ = σ(6.296 − 4.52)
Φ = σ(1.776)
Φ ≈ 0.871


**Prediction:** Φ > 0.5 → Regeneration  
**Observed:** Complete regeneration ✅

## Why This System Was Chosen for Independent Validation

1. **Phylogenetic distance:** Zebrafish (teleost fish) diverged from mammalian lineage ~450 million years ago
2. **Tissue novelty:** Cardiac tissue not represented in training set (which included liver, skin, limb, fin, lens, tentacle, antler, head)
3. **Clinical relevance:** Heart regeneration is a major unsolved problem in human medicine
4. **Well-characterized:** Molecular basis extensively studied, enabling confident layer estimation
5. **Binary outcome:** Clear regenerative success vs. failure (unlike partial regeneration systems)

Successful prediction on this out-of-sample system provides strong evidence that ARC captures fundamental principles of regeneration rather than overfitting to specific taxa or tissues.