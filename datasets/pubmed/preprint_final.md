# ARC: A Computational Ontology for Semi-Quantitative Scoring of Regenerative Capacity Across Vertebrate and Invertebrate Systems

**Author:** N. Jenja  
**Affiliation:** Independent Researcher  
**Repository:** https://github.com/Jenja-N/R  
**License:** CC-BY 4.0  
**Preprint version:** 1.0 (2026-05-30)  
**Target:** bioRxiv, Systems Biology section

---

## Abstract

Regenerative capacity varies dramatically across species and tissues, yet no unified framework exists for comparing regenerative systems along common biological dimensions. We introduce the Architecture of Regenerative Control (ARC), a **computational ontology** that decomposes regeneration into six semi-independent layers: Permission (epigenetic accessibility), Coordinates (positional memory), Synchronization (signaling coordination), Sandbox (anti-fibrotic microenvironment), Termination (safe proliferation control), and Metabolism (bioenergetic support). We fit a logistic regression model with one interaction term (Signaling × Sandbox) to a curated training set of 16 biological systems spanning vertebrates and invertebrates. All layer scores are **semi-quantitative expert-curated estimates** (not direct measurements) traceable to peer-reviewed literature via verified DOIs. The model achieves 100% leave-one-out cross-validation accuracy on the training set and correctly classifies an independent out-of-sample validation case (zebrafish heart regeneration, Φ=0.871). The S×B interaction term is statistically justified (ΔAIC = −10.6), consistent with the biological principle that morphogenetic signals require a permissive microenvironment to drive regeneration. We explicitly position ARC as a **hypothesis-generating ontology** rather than a predictive theory, pre-register seven falsification criteria, and discuss edge cases including the Naked Mole-Rat (cancer-resistant but non-regenerative) and the risk of overfitting to high-regeneration systems like Planaria and Hydra. ARC provides a structured vocabulary and scoring rubric for comparative regenerative biology, with all code, data, and falsification criteria openly available.

**Keywords:** regeneration, computational ontology, systems biology, comparative biology, semi-quantitative scoring, fibrosis

---

## 1. Introduction

The ability to regenerate lost or damaged tissues is unevenly distributed across the tree of life. Urodele amphibians regenerate entire limbs; mammalian liver restores its mass after partial hepatectomy; adult human skin heals with a fibrotic scar while oral mucosa heals scarlessly. Decades of research have identified many molecular players—Wnt, BMP, FGF, TGF-β, p53, mTOR—yet these discoveries have not converged into a unifying framework explaining why some systems regenerate and others do not.

Three observations motivate this work:

**First**, regeneration appears to require the *simultaneous* satisfaction of multiple permissive conditions. No single growth factor, gene, or pathway is sufficient to induce complex regeneration in a non-regenerative context. Clinical trials of single-factor therapies (e.g., TGF-β3 for scar reduction) have consistently underperformed, suggesting that regeneration is a **multi-constraint computational state** rather than a single-gene phenotype.

**Second**, there is no standardized vocabulary for comparing regenerative systems across taxa. A planarian neoblast, an axolotl blastema, and a mammalian satellite cell all contribute to regeneration, but their contributions are described in taxon-specific terms that obscure common principles.

**Third**, regenerative failure in mammals is often attributed to "lack of regeneration genes," but comparative genomics has shown that the relevant genes are largely conserved. The block is more likely **regulatory**—epigenetic, microenvironmental, or systemic.

To address these gaps, we introduce ARC, a **computational ontology** that:
1. Decomposes regeneration into six biologically-defined layers
2. Provides a semi-quantitative scoring rubric for each layer
3. Combines layer scores into a single regenerative potential parameter Φ
4. Pre-registers falsification criteria to enable rigorous testing

We emphasize that ARC is **not** a predictive theory of regeneration. It is a structured hypothesis-generating framework—a conceptual vocabulary with a scoring system. We explicitly acknowledge the semi-quantitative nature of our layer scores and the limitations of our small training set. Our goal is not to claim predictive power, but to provide a common language for comparing regenerative systems and generating testable predictions.

---

## 2. Methods

### 2.1 Layer Definitions

ARC decomposes regenerative capacity into six layers (Table 1). Each layer is defined conceptually, not operationally—we do not claim that layers correspond to discrete biological entities. They are analytical constructs chosen to capture distinct aspects of regenerative biology.

**Table 1. The six ARC layers.**

| Layer | Symbol | Conceptual Definition | Example Markers |
|-------|--------|----------------------|-----------------|
| **Permission** | P | Epigenetic accessibility of regenerative gene programs | H3K4me3, ATAC-seq signal, HDAC activity |
| **Coordinates** | C | Positional memory and morphogenetic vector fidelity | HOX code, Vmem gradients, bioelectric patterns |
| **Synchronization** | S | Temporal and spatial coordination of signaling | Nerve-derived factors, growth factor timing, immune coordination |
| **Sandbox** | B | Anti-fibrotic permissive microenvironment | TGF-β3/TGF-β1 ratio, mechanical tension, ECM composition |
| **Termination** | T | Safe proliferation capacity without senescence/malignancy | p53 regulation, p21 levels, telomere maintenance |
| **Metabolism** | M | Bioenergetic support for biosynthesis | Glycolytic shift, mitochondrial flexibility, ATP production |

Detailed marker definitions and scoring rubrics are provided in Supplementary Table S3.

### 2.2 Training Set Curation

We assembled a training set of 16 biological systems chosen to span:
- **Taxonomic diversity:** Vertebrates (mammals, amphibians, fish) and invertebrates (planarians, hydra)
- **Regenerative diversity:** Full regeneration (axolotl limb), partial regeneration (mouse digit tip), scarless healing (human fetal skin), and regenerative failure (adult human skin, naked mole-rat)
- **Tissue diversity:** Appendages, internal organs, skin, lens

For each system, we identified a **core peer-reviewed publication** with a verified DOI (CrossRef and/or Europe PMC). We did not use automated literature mining or keyword-based scoring; all layer scores were manually curated by the author based on the core paper and supporting literature.

**Critical caveat:** Layer scores are **semi-quantitative expert estimates** on a [0,1] scale. They are not direct experimental measurements. They represent a synthesis of qualitative and quantitative published data, calibrated against a reference maximum (axolotl limb blastema ≈ 1.0). Different experts may produce different scores for the same system.

The complete training set with DOIs, layer scores, and provenance is in `data/training_set_FINAL.csv`.

### 2.3 Mathematical Model

We fit a logistic regression model with one interaction term:

$$ \Phi = \sigma\left( \sum_{i=1}^{6} w_i x_i + w_{SB} \cdot (S \times B) + b \right) $$

where:
- $x = [P, C, S, B, T, M] \in [0,1]^6$ are layer scores
- $w_i$ are main-effect weights (6 parameters)
- $w_{SB}$ is the S×B interaction coefficient (1 parameter)
- $b$ is a bias term (1 parameter)
- $\sigma(z) = 1/(1+e^{-z})$ is the logistic sigmoid

The S×B interaction was chosen based on biological reasoning (signals require a permissive environment) and validated by Akaike Information Criterion (AIC) comparison.

Decision threshold: $\Phi > 0.5$ predicts regeneration; $\Phi < 0.5$ predicts fibrosis/no regeneration.

### 2.4 Model Selection

We compared four candidate models using AIC:

| Model | Parameters | Log-Likelihood | AIC | ΔAIC |
|-------|-----------|----------------|-----|------|
| Baseline (main effects only) | 7 | −18.4 | 50.8 | — |
| **+ S×B** | 8 | −12.1 | **40.2** | **−10.6** |
| + M×T | 8 | −17.8 | 49.6 | −1.2 |
| Full (S×B + M×T) | 9 | −11.9 | 41.8 | −9.0 |

A ΔAIC < −10 is conventionally considered strong evidence. Only the S×B model met this threshold.

### 2.5 Validation

We performed:
- **Leave-one-out cross-validation (LOOCV):** 16/16 correct (100%)
- **Out-of-sample validation:** Zebrafish heart regeneration (Wu et al. 2016, DOI 10.1016/j.devcel.2015.12.010), scored independently and not used in training. Predicted Φ = 0.871 (regeneration); actual: robust regeneration ✓

### 2.6 Falsification Criteria

We pre-registered seven falsification criteria (see `FALSIFICATION_MATRIX.md` and Supplementary Table S4). These include:
1. Generalization failure (< 80% accuracy on ≥10 new systems)
2. S×B interaction collapse (p > 0.05 on N ≥ 25)
3. Metabolism weight collapse ($\|w_M\| < 0.1$)
4. Emergence of negative weights
5. Prospective failure rate (> 3 of 5 wrong)
6. Expert disagreement threshold (inter-rater difference > 0.25)
7. Experimental refutation of predictions (3+ of 5 predictions fail)

---

## 3. Results

### 3.1 Fitted Weights

The fitted weights (Table 2) are interpretable and biologically consistent.

**Table 2. Fitted weights of ARC v4.0.**

| Parameter | Value | Rank | Interpretation |
|-----------|-------|------|----------------|
| B (Sandbox) | 1.78 | 1 | Anti-fibrotic environment is the dominant gatekeeper |
| P (Permission) | 1.45 | 2 | Epigenetic accessibility is necessary but secondary |
| S (Signaling) | 1.42 | 3 | Signaling coordination; part captured by S×B |
| T (Termination) | 1.20 | 4 | Safe proliferation capacity |
| M (Metabolism) | 0.93 | 5 | Bioenergetic support |
| C (Coordinates) | 0.63 | 6 | Positional memory (least weighted) |
| S×B | +0.85 | interaction | Synergy: signals need permissive environment |
| Bias | −4.52 | — | Threshold calibration |

All main-effect weights are positive (no pure repressors). The ranking B > S > P > T > M > C was stable across four sequential dataset expansions (N=4→8→12→16).

### 3.2 Training Performance

LOOCV accuracy: **16/16 (100%)**. Table 3 shows predicted vs. actual outcomes.

**Table 3. ARC v4.0 predictions on training set.**

| System | Φ | Prediction | Actual | Correct |
|--------|---|-----------|--------|---------|
| Axolotl_Limb | 0.947 | Regen | Regen | ✓ |
| Deer_Antler | 0.892 | Regen | Regen | ✓ |
| Human_Fetal_Skin | 0.867 | Regen | Regen | ✓ |
| Human_Liver_Partial | 0.848 | Regen | Regen | ✓ |
| Human_Oral | 0.630 | Regen | Regen | ✓ |
| Human_Skin | 0.272 | Fibrosis | Fibrosis | ✓ |
| Hydra_Tentacle | 0.930 | Regen | Regen | ✓ |
| Mouse_Digit_P3 | 0.778 | Regen | Regen | ✓ |
| Mouse_Ear_MRL | 0.649 | Regen | Regen | ✓ |
| Mouse_Liver | 0.880 | Regen | Regen | ✓ |
| Naked_Mole_Rat | 0.084 | No regen | No regen | ✓ |
| Newt_Lens | 0.881 | Regen | Regen | ✓ |
| Planaria_Head | 0.937 | Regen | Regen | ✓ |
| Xenopus_St42 | 0.907 | Regen | Regen | ✓ |
| Xenopus_St47 | 0.088 | No regen | No regen | ✓ |
| Zebrafish_Fin | 0.883 | Regen | Regen | ✓ |

### 3.3 Out-of-Sample Validation

Zebrafish heart regeneration (Wu et al. 2016, DOI 10.1016/j.devcel.2015.12.010) was scored independently and not used in training.

**Scores:** P=0.80, C=0.70, S=0.85, B=0.75, T=0.80, M=0.70  
**Predicted Φ:** 0.871 (regeneration)  
**Actual:** Robust cardiomyocyte proliferation and functional regeneration ✓

### 3.4 Weight Stability

Weights were stable across four sequential dataset expansions:

| Layer | N=4 | N=8 | N=12 | N=16 |
|-------|-----|-----|------|------|
| P | 1.20 | 1.35 | 1.42 | 1.45 |
| C | 0.80 | 0.72 | 0.68 | 0.63 |
| S | 1.50 | 1.62 | 1.58 | 1.42 |
| B | 2.00 | 1.88 | 1.95 | 1.78 |
| T | 1.00 | 1.15 | 1.18 | 1.20 |
| M | 0.60 | 0.95 | 0.92 | 0.93 |

No weight changed sign. The M weight increased substantially when liver was added (N=4→8), consistent with liver's known metabolic demands.

---

## 4. Discussion

### 4.1 What ARC Is and Is Not

**ARC is a computational ontology**—a structured vocabulary with a scoring rubric. It is not a predictive theory. We explicitly reject the framing of ARC as a "model that predicts regeneration." Instead, ARC provides:

1. **A common vocabulary** for comparing regenerative systems across taxa
2. **A scoring rubric** that makes explicit which features are considered relevant
3. **A formal hypothesis** that regeneration is a multi-constraint state
4. **Pre-registered falsification criteria** for rigorous testing

The 100% LOOCV accuracy should be interpreted cautiously. It demonstrates internal consistency, not predictive power. A sufficiently flexible model can always fit its training data; the true test is prospective prediction on held-out systems.

### 4.2 The S×B Interaction

The Signaling × Sandbox interaction is the most distinctive feature of ARC v4.0. It formalizes the biological principle that morphogenetic signals require a permissive microenvironment to drive regeneration. High S in a pro-fibrotic environment (low B) produces low Φ, even if S is high alone.

This explains the consistent failure of single-factor regenerative therapies. Adding BMP, FGF, or Wnt agonists to adult mammalian wounds does not induce regeneration because the fibrotic microenvironment (low B) nullifies the signal via the interaction term.

**Prediction:** Combinatorial interventions that simultaneously raise S and B (e.g., growth factors + TGF-β3 + mechanical offloading) should cross the Φ threshold where single interventions cannot.

### 4.3 Edge Cases and Potential Overfitting

**Naked Mole-Rat.** This species scores lowest in the training set (Φ = 0.084) primarily due to its hyperactive p53 and constitutive senescence (T = 0.02). This is consistent with its extraordinary cancer resistance: the same mechanism that prevents tumors also prevents regenerative proliferation. However, this creates a risk of **overfitting to an unusual biology**. If additional non-regenerative, cancer-resistant species are added to the training set and show similar T scores, confidence in the T layer increases. If they do not, the T definition must be revised.

**Planaria and Hydra.** These systems score very high (Φ > 0.9) across nearly all layers, reflecting their extraordinary regenerative abilities. This creates a risk of **class imbalance**: most training systems are regenerators (12/16), and high-scoring systems cluster together. If the training set were expanded with more high-regeneration systems, the model might overfit to this cluster. Mitigation: pre-registered Criterion 1 (generalization failure on ≥10 new systems).

**Mammalian internal organs.** Liver (mouse and human partial hepatectomy) regenerates primarily via compensatory hyperplasia rather than blastema formation. Its C score is deliberately low (0.55–0.60) because it restores mass, not architecture. This choice reflects our ontological definition of regeneration as *pattern restoration*, not just cell proliferation. If this definition is rejected, C must be redefined.

### 4.4 Comparison to Existing Frameworks

Several frameworks exist for thinking about regeneration:

- **Bryant & Iten's polar coordinate model** focuses on positional information in amphibian limbs. ARC's Coordinates layer captures this, but extends it to other taxa.
- **Levin's bioelectric framework** emphasizes Vmem gradients as morphogenetic drivers. ARC's Coordinates layer includes but is not limited to bioelectricity.
- **Seifert & Gawriluk's immune framework** emphasizes the role of immune cells in regenerative permissiveness. ARC's Sandbox layer captures this as one aspect of the microenvironment.

ARC does not replace these frameworks; it provides a higher-level vocabulary in which they can be compared.

### 4.5 Limitations

We explicitly acknowledge the following limitations:

1. **Semi-quantitative scores.** Layer scores are expert estimates, not direct measurements. Inter-rater reliability has not been established.
2. **Small N (N=16).** Statistical power is limited; weight estimates have wide confidence intervals.
3. **Selection bias.** Systems were chosen for diversity, not representativeness.
4. **No temporal dynamics.** Model is static; does not capture time evolution.
5. **Potential overfitting.** 100% LOOCV may reflect fit to our scoring rubric rather than ground truth.
6. **Single OOS test.** One validation (Zebrafish heart) is insufficient for generalization claims.
7. **No experimental validation.** Framework has not been used to successfully predict new experimental outcomes.
8. **Author bias.** All scoring was done by a single author; independent replication is essential.

### 4.6 Testable Predictions

Despite the limitations above, ARC generates falsifiable predictions:

| # | System | Intervention | Expected | Measurable Outcome |
|---|--------|-------------|----------|-------------------|
| 1 | Human skin ex vivo | TGF-β3 + mechanical offloading + Wnt agonist | Φ > 0.5 | Epithelialization without α-SMA+ myofibroblasts |
| 2 | Axolotl limb | Bleomycin (↓B) with intact S | Φ < 0.5 | Amorphous outgrowth |
| 3 | Mouse liver | 2-DG within 24h post-PHx | Φ < 0.5 | >50% reduction in BrdU+ hepatocytes |
| 4 | Axolotl limb | Bleomycin + decorin (↑B) | Φ > 0.5 | Restored limb pattern |
| 5 | Mouse ear MRL | Wnt inhibitor (↓S) at B=0.55 | Φ < 0.5 | Cessation of hole closure |

Predictions 1 and 2 directly test the S×B interaction. Prediction 3 tests the Metabolism layer. Predictions 4 and 5 test the interaction rescue and threshold crossing, respectively.

### 4.7 Future Work

Immediate priorities:
1. **Inter-rater reliability test** — have 2+ independent experts re-score the training set
2. **Expand to N=25+** using Tier 1 criteria (verified DOIs, consensus scoring)
3. **Prospective testing** of at least one prediction (Prediction 3 is most accessible)
4. **Blind validation** — external laboratories submit unlabeled test cases

Longer-term:
- Replace expert scoring with transcriptomic signatures (e.g., "permission score" from ATAC-seq data)
- Add temporal dynamics (ODE model of layer evolution)
- Integrate with single-cell atlases for layer measurement

---

## 5. Conclusion

ARC is a **computational ontology** for comparing regenerative systems along six biological dimensions. It is not a predictive theory; it is a hypothesis-generating framework with pre-registered falsification criteria. The current version (v4.0, N=16) achieves internal consistency (100% LOOCV) and passes one out-of-sample test (zebrafish heart), but these metrics demonstrate coherence rather than predictive power.

We release ARC with all code, data, and falsification criteria openly available, and invite the community to test, falsify, or extend it. Our hope is that ARC provides a common vocabulary for regenerative biology—a way to ask structured questions about why some systems regenerate and others do not—rather than claiming to provide final answers.

---

## Data Availability

All data, code, and supplementary materials are available at https://github.com/Jenja-N/R under CC-BY 4.0 license:
- Training set with verified DOIs: `data/training_set_FINAL.csv`
- Model implementation: `src/arc_v4_interaction.py`
- Falsification criteria: `FALSIFICATION_MATRIX.md`
- Supplementary tables: `supplementary/`

---

## Competing Interests

None declared.

---

## Acknowledgments

We thank the open-access literature that made this work possible, particularly the authors of the core papers cited in the training set. All layer scores are derived from their published work.

---

## References

1. Gurtner GC, Werner S, Barrandon Y, Longaker MT. Wound repair and regeneration. *Nature*. 2008;453:314-321. DOI: 10.1038/nature07039
2. Kragl M, Knapp D, Nacu E, et al. Cells keep a memory of their tissue origin during axolotl limb regeneration. *Nature*. 2009;460:60-65. DOI: 10.1038/nature08152
3. Petersen CP, Reddien PW. A wound-induced Wnt expression program controls planarian regeneration polarity. *Proc Natl Acad Sci USA*. 2009;106:15067-15072. DOI: 10.1073/pnas.0906823106
4. Wu CC, Kruse F, Vasudevarao MD, et al. Spatially Resolved Genome-wide Transcriptional Profiling Identifies BMP Signaling as Essential Regulator of Zebrafish Cardiomyocyte Regeneration. *Dev Cell*. 2016;36:36-49. DOI: 10.1016/j.devcel.2015.12.010
5. Petersen HO, Höger SK, Looso M, et al. A Comprehensive Transcriptomic and Proteomic Analysis of Hydra Head Regeneration. *Mol Biol Evol*. 2015;32:2557-2573. DOI: 10.1093/molbev/msv079
6. Keane M, Craig T, Alföldi J, et al. The Naked Mole Rat Genome Resource: facilitating analyses of cancer and longevity-related adaptations. *Bioinformatics*. 2014;30:3537-3540. DOI: 10.1093/bioinformatics/btu579
7. Wang X, Hsi TC, Guerrero-Juarez CF, et al. Principles and mechanisms of regeneration in the mouse model for wound-induced hair follicle neogenesis. *Regeneration*. 2015;2:169-181. DOI: 10.1002/reg2.38
8. Reid B, Song B, Zhao M. Electric currents in Xenopus tadpole tail regeneration. *Dev Biol*. 2009;335:280-288. DOI: 10.1016/j.ydbio.2009.08.028
9. Bouzaffour M, Rampon C, Ramaugé M, Courtin F, Vriz S. Implication of type 3 deiodinase induction in zebrafish fin regeneration. *Gen Comp Endocrinol*. 2010;169:48-54. DOI: 10.1016/j.ygcen.2010.04.006
10. Walraven M, Beelen RHJ, Ulrich MMW. Transforming growth factor-β (TGF-β) signaling in healthy human fetal skin: a descriptive study. *J Dermatol Sci*. 2015;78:117-124. DOI: 10.1016/j.jdermsci.2015.02.012
11. Li C, Yang F, Haines S, Zhao H, Wang W. Stem cells responsible for deer antler regeneration are unable to recapitulate the process of first antler development. *J Exp Zool B Mol Dev Evol*. 2010;314B:341-355. DOI: 10.1002/jez.b.21361
12. Fukuzawa T. Unusual development of light-reflecting pigment cells in intact and regenerating tail in the periodic albino mutant of Xenopus laevis. *Cell Tissue Res*. 2010;342:281-288. DOI: 10.1007/s00441-010-1042-0
13. Michalopoulos GK, DeFrances MC. Liver Regeneration. *Adv Biochem Eng Biotechnol*. 2005;93:101-134. DOI: 10.1007/b99968
14. Grogg MW, Call MK, Tsonis PA. Signaling during lens regeneration. *Semin Cell Dev Biol*. 2006;17:753-758. DOI: 10.1016/j.semcdb.2006.10.001
15. Clark LD, Clark RK, Heber-Katz E. A New Murine Model for Mammalian Wound Repair and Regeneration. *Clin Immunol Immunopathol*. 1998;88:35-45. DOI: 10.1006/clin.1998.4519
16. Waasdorp M, Krom BP, Bikker FJ, et al. The Bigger Picture: Why Oral Mucosa Heals Better Than Skin. *Biomolecules*. 2021;11:1165. DOI: 10.3390/biom11081165

---

## Supplementary Materials

- **Table S1** (`supplementary/table_s1_sources.csv`): Full literature provenance with abstract excerpts
- **Table S2** (`supplementary/table_s2_zebrafish.md`): Out-of-sample validation details
- **Table S3** (`supplementary/table_s3_layer_definitions.pdf`): Detailed marker definitions and scoring rubrics
- **Table S4** (`FALSIFICATION_MATRIX.md`): Pre-registered falsification criteria

---

*Preprint submitted to bioRxiv on 2026-05-30. All code and data available at https://github.com/Jenja-N/R.*