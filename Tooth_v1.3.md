# ARC-Tooth v1.3: 

*Experiment results, data, and scripts are located here: https://github.com/Jenja-N/R/tree/main/experiments*.

**Repository:** Jenja-N/R  
**Version:** ARC-Tooth v1.3 — consolidated from v1.1 (conceptual), v1.2 (empirical update), and full cross-verification across all project datasets and PDFs  
**Language:** English  
**Scientific stance:** strictly hypothesis-generating. Observational scRNA-seq re-analysis. No clinical claims. All findings falsifiable.  
**Date:** June 2026  
**Code:** https://github.com/Jenja-N/R  
**Status:** Pre-experimental. Not validated therapy. Not approved for clinical use.

---

## 0. Preamble: Logical Structure of This Document

This document integrates four independent evidence layers:

1. **v1.1 (conceptual ARC framework)** — the architectural model of tooth regeneration as a multi-layer state transition
2. **v1.2 (GSE184749 + GSE185222 empirical analysis)** — scRNA-seq data constraining the framework with numbers
3. **CSV cross-verification** — direct numerical reconciliation across all output files
4. **PDF cross-verification** — alignment with published literature (Alghadeer et al. 2023 *Dev Cell*; Wang et al. 2022 *Sci Bull*; Hermans et al. 2023; Shi et al. 2024 *Cell Prolif*)

Every claim in sections 3–9 is traceable to a specific data row or published source. Claims without such traceability are explicitly labeled **[INFERENCE]** with confidence level.

**Central falsifiable claim of ARC-Tooth v1.3:**

> Adult deep caries dental pulp is not a "blocked fetal program." It is a qualitatively distinct transcriptional attractor state. The transition from this state to DSPP+ odontoblast differentiation requires simultaneous correction of ≥3 independent deficits: TGFB2 depletion, JUNB-dominant AP-1 imbalance, and IGF1R/RUNX2 competence loss. No single-molecule intervention is predicted to be sufficient.

This claim is falsifiable: a single-molecule intervention (e.g., TGFB2 alone) producing >1% DSPP+ cells in a validated inflamed hDPSC model would partially falsify the multi-component requirement.

---

## 1. Datasets and Analytical Pipeline

### 1.1 Datasets Used

| Dataset | Source | Dimensions | Conditions | Role in this analysis |
|---------|--------|-----------|------------|----------------------|
| GSE185222 | Adult human dental pulp biopsies | 20,293 genes × 6,582 cells | sound (n=890), enamel_caries (n=1,731), deep_caries (n=3,961) | Primary adult reference; failure mode characterization |
| GSE184749 | Human fetal jaw tissue, 9–22 gw | 42,174 genes × 17,510 cells | DEM (2,752), OB (681), POB (4,064), SOB (4,959), SOBP (5,054) | Reference developmental state |
| GSE146123 | Mouse incisor + molar, 8–16 wk | — | cervical loop epithelium | SOX2+ epithelial stem cell reference (mouse) |
| GSE189381 | Mouse molar time-course E13.5–P7.5 | — | pseudotime | Developmental trajectory reference |

**Critical caveat (applies to all analyses):**
- GSE185222 has **no cell-type annotation**. All GSE185222 findings are at the bulk-pseudocell per-condition level.
- Adult DSPP+ cells: n=3 total (all in enamel_caries). No quantitative conclusions can be drawn from this subsample.
- GENIE3 weights are **not causal**. They reflect variance explained in gene expression, not regulatory causation.
- All comparisons are **cross-sectional**. No longitudinal data available.
- GSE184749 contains dental mesenchyme only. **No epithelium** is present in this dataset.

### 1.2 Key Analytical Scripts

All code available at https://github.com/Jenja-N/R:

| Script | Function | Key outputs |
|--------|----------|-------------|
| `COMPLETE_RIGOROUS_ANALYSIS.R` | Baseline expression, correlations, co-expression | — |
| `NMF_GSE185222.R` | NMF module analysis (k=5,8,10,12,15) | Module weights by condition |
| `QWEN2.R`–`QWEN8.R` | GENIE3 regulon, initial comparisons | `regulon_comparison.csv` |
| `QWEN9.R`–`QWEN10.R` | Fetal OB vs. adult deep comparison | `fetal_OB_vs_adult_deep.csv` |
| `QWEN11.R`–`QWEN12.R` | DSPP+ profiling, WNT hypothesis test | `fetal_dspp_signature.csv`, `wnt_signaling_summary.csv` |
| `QWEN13.R` | TGFB/BMP/stress gene comparison | `tgfb_bmp_stress_comparison.csv` |
| `QWEN14.R`–`QWEN16.R` | Embryo cell-type atlas, global comparison | `embryo_vs_adult_deep_comparison.csv` |

---

## 2. Core Numerical Results (Cross-Verified)

All numbers below have been verified across: report1.MD, report2.MD, report3.MD, and the corresponding CSV files. Discrepancies are noted where found.

### 2.1 Fetal OB vs. Adult Deep Caries: Full Table

**Source:** `fetal_OB_vs_adult_deep.csv` (n=26 genes), cross-verified with report2.MD Table 2.3

| Gene | Fetal OB % | Adult Deep % | Δ (pp) | Interpretation |
|------|-----------|-------------|--------|---------------|
| **COL1A1** | 89.72 | 48.32 | **+41.40** | Matrix organization loss |
| **DSPP** | 34.36 | 0.00 | **+34.36** | Complete terminal differentiation block |
| **COL1A2** | 88.40 | 60.72 | +27.68 | Partial matrix loss |
| **IGF1R** | 16.89 | 7.55 | +9.34 | Competence signal deficit |
| **RUNX2** | 24.23 | 16.01 | +8.22 | Differentiation TF deficit |
| **BMPR2** | 21.88 | 15.96 | +5.92 | BMP receptor deficit |
| **BMPR1A** | 8.81 | 5.58 | +3.23 | BMP receptor deficit |
| **TGFB2** (from tgfb csv) | 16.45 | 2.47 | **+13.97** | **Primary differentiation-context deficit** |
| SOSTDC1 | 2.20 | 0.10 | +2.10 | Already absent in adult |
| DMP1 | 0.44 | 1.34 | −0.90 | Slightly higher in adult |
| AXIN2 | 0.15 | 1.31 | −1.17 | Slightly higher in adult |
| WNT10A | 1.32 | 2.55 | −1.23 | Similar / slightly higher in adult |
| BMP4 | 0.88 | 7.90 | **−7.02** | **HIGHER in adult** — not the deficit |
| RELA | 0.15 | 10.53 | −10.38 | NF-κB pathway elevated in adult |
| NFKB1 | 1.76 | 15.96 | −14.19 | NF-κB elevated in adult |
| DLX6 | 0.73 | 17.19 | −16.46 | Elevated in adult (stress) |
| MSX2 | 1.62 | 29.97 | −28.35 | **Paradox: elevated in adult, DSPP=0** |
| **MSX1** | 17.03 | 42.34 | **−25.30** | **Paradox: elevated in adult, DSPP=0** |
| **DLX5** | 0.73 | 38.00 | **−37.26** | **Paradox: elevated in adult, DSPP=0** |
| **PTN** | 25.11 | 64.28 | **−39.17** | Stress/repair context, not odontogenic |
| **MDK** | 4.99 | 50.44 | **−45.45** | Stress/repair context, not odontogenic |
| **FOS** | 4.55 | 72.15 | **−67.60** | AP-1 stress dominance |
| **JUN** | 3.38 | 76.67 | **−73.30** | AP-1 stress dominance |
| **JUNB** | 0.44 | 87.23 | **−86.78** | **Specific pathological AP-1 imbalance** |

**Logical consequence:** The 5 largest deficits in adult deep (COL1A1, DSPP, COL1A2, TGFB2, IGF1R) reflect output failure and input depletion. The 5 largest excesses (JUNB, JUN, FOS, MDK, PTN) reflect a dominant stress transcriptional state. These are not independent: the JUNB-dominant AP-1 state and the TGFB2/IGF1R deficit co-exist and are causally linked in the developmental context (JUNB directly competes with RUNX2 for target promoters — [INFERENCE, confidence 0.6]).

### 2.2 TGFB/BMP/Stress Gene Comparison: Fetal OB vs. Adult Deep

**Source:** `tgfb_bmp_stress_comparison.csv` (n=10 genes)

| Gene | Fetal OB % | Adult Deep % | Δ (pp) | Fetal OB mean | Adult Deep mean | Fold drop (mean) |
|------|-----------|-------------|--------|--------------|----------------|-----------------|
| **TGFB2** | 16.45 | 2.47 | **+13.97** | 1.321 | 0.028 | **47.2×** |
| **IGF1R** | 16.89 | 7.55 | +9.34 | 1.346 | 0.115 | **11.7×** |
| **RUNX2** | 24.23 | 16.01 | +8.22 | 1.940 | 0.219 | **8.9×** |
| BMP7 | 7.05 | 6.59 | +0.46 | 0.559 | 0.093 | ~equal (%) |
| TGFB3 | 0.88 | 4.85 | −3.97 | 0.067 | 0.058 | Higher in adult |
| BMP2 | 2.94 | 8.28 | −5.34 | 0.233 | 0.113 | Higher in adult (%) |
| BMP4 | 0.88 | 7.90 | **−7.02** | 0.065 | 0.096 | **Higher in adult** |
| FOS | 4.55 | 72.15 | −67.60 | 0.349 | 2.455 | 7.0× **higher** in adult |
| JUN | 3.38 | 76.67 | −73.30 | 0.262 | 2.484 | 9.5× **higher** in adult |
| JUNB | 0.44 | 87.23 | **−86.78** | 0.036 | 2.768 | **76.9× higher** in adult |

**Key logical conclusion:** BMP4, TGFB3, and BMP2 are all elevated (or equal) in adult deep. They are not the limiting factor. The deficits are exclusively in TGFB2, IGF1R, and RUNX2. JUNB shows the most extreme fold-elevation of any measured gene (76.9× by mean expression).

### 2.3 Embryo (All Cell Types) vs. Adult Deep: Global Comparison

**Source:** `embryo_vs_adult_deep_comparison.csv` (n=9 genes)

| Gene | Embryo all cells % | Adult Deep % | Δ (pp) | Verdict |
|------|--------------------|-------------|--------|---------|
| TGFB2 | 20.80 | 2.47 | **+18.32** | Systematic deficit across all embryonic cell types |
| RUNX2 | 33.25 | 16.01 | **+17.24** | Systematic deficit |
| BMP4 | 21.25 | 7.90 | +13.35 | Deficit vs. embryo all-cells; but NOT vs. fetal OB (see paradox above) |
| IGF1R | 20.10 | 7.55 | **+12.55** | Systematic deficit |
| BMP2 | 8.83 | 8.28 | +0.54 | No meaningful difference |
| DSPP | 0.35 | 0.00 | +0.35 | Even early embryo has nearly no DSPP — correct (pre-differentiation) |
| JUN | 71.40 | 76.67 | **−5.27** | **NOT a deficit**: JUN is ~equal in embryo and adult |
| FOS | 45.63 | 72.15 | −26.52 | Higher in adult, but elevated in embryo too |
| JUNB | 39.63 | 87.23 | **−47.59** | **Specific to adult deep**: 2.2× above embryo average |

**Critical implication:** JUN is not a pathological marker. It is elevated in embryonic dental tissues (71.4% of all cells, 94.4% in follicle). The specific pathological AP-1 component is **JUNB** — 87.23% in adult deep vs. 39.63% in embryo mesenchyme (Δ = 47.59 pp), and vs. 34.20% in embryo mesenchyme specifically.

### 2.4 WNT Signaling Hypothesis: Falsification Data

**Source:** `wnt_signaling_summary.csv`

| Condition | n_cells | AXIN2+ | LEF1+ | DSPP+ | AXIN2+LEF1+DSPP+ |
|-----------|--------|--------|-------|-------|-----------------|
| Fetal POB | 4,064 | 36 (0.89%) | 455 (11.2%) | 37 (0.91%) | **0** |
| Fetal OB | 681 | 1 (0.15%) | 45 (6.61%) | 234 (34.36%) | **0** |
| Adult enamel | 1,731 | 34 (1.96%) | 52 (3.0%) | 3 (0.17%) | **0** |
| Adult deep | 3,961 | 52 (1.31%) | 125 (3.16%) | 0 | **0** |

**Verdict: WNT/β-catenin hypothesis FALSIFIED.**

The triple co-expression AXIN2+LEF1+DSPP+ = 0 in all conditions, including in the fetal OB cluster where DSPP+ cells are abundant (34.36%). WNT canonical targets (AXIN2, LEF1) are **lower** in DSPP+ cells than in DSPP− cells (fold changes 0.47 and 0.51 respectively, from `fetal_dspp_signature.csv`). WNT/β-catenin canonical signaling does not drive DSPP expression.

### 2.5 Fetal DSPP+ Cell Signature

**Source:** `fetal_dspp_signature.csv` (n=307 DSPP+ fetal cells vs. DSPP− cells)

**Genes ENRICHED in DSPP+ cells (fold change, ranked):**

| Gene | FC (DSPP+ / DSPP−) | DSPP+ % | DSPP− % | Interpretation |
|------|-------------------|---------|---------|---------------|
| WNT10A | **10.40** | 1.6 | 0.2 | WNT ligand (non-canonical role) |
| BMP2 | **4.56** | 3.6 | 0.8 | BMP2, not BMP4, is the relevant BMP |
| DMP1 | **4.04** | 0.3 | 0.1 | Matrix protein co-expressed with DSPP |
| TGFB2 | **2.81** | 10.1 | 3.5 | Differentiation-context signal |
| COL1A2 | 2.18 | 92.2 | 51.5 | Matrix output |
| COL1A1 | 2.15 | 93.2 | 53.9 | Matrix output |
| DLX6 | 1.88 | 1.0 | 0.5 | Minor enrichment |
| TGFB3 | 1.53 | 1.6 | 1.1 | Minor enrichment |
| FOS | 1.50 | 12.1 | 7.8 | Present in DSPP+ cells — not purely inhibitory |
| MSX1 | 1.37 | 21.2 | 15.0 | Slight enrichment |

**Genes DEPLETED in DSPP+ cells (FC < 1.0):**

| Gene | FC | Interpretation |
|------|----|---------------|
| AXIN2 | 0.47 | WNT canonical target: LOWER in DSPP+ cells |
| LEF1 | 0.51 | WNT canonical target: LOWER in DSPP+ cells |
| BMPR1B | 0.55 | |
| MDK | 0.65 | Not a DSPP+ marker |
| RUNX2 | 0.63 | **LOWER in terminally differentiated cells** — RUNX2 is a pre-differentiation TF, not a differentiation marker |
| NFKB1 | 0.69 | NF-κB lower in DSPP+ cells |
| IGF1R | 1.05 | Neutral — broadly expressed, not DSPP-specific |
| PTN | 0.99 | Neutral — not DSPP-specific |
| SOX2 | 0.00 | Completely absent in DSPP+ cells |

**Critical interpretation:** The odontoblast differentiation state is defined by WNT10A ligand (non-canonical), BMP2/DMP1, and organized COL1 matrix. It is **not** defined by high RUNX2, high AXIN2/LEF1, or high IGF1R. RUNX2 is a pre-differentiation marker that is actively downregulated during terminal differentiation.

### 2.6 GENIE3 Regulon: Regulatory Context Shift

**Source:** `regulon_comparison.csv` (n=66 TF→target pairs), selected critical rows

**DLX5 → PTN regulatory weight by condition:**

| Condition | DLX5→PTN weight | FOS→PTN weight | Interpretation |
|-----------|----------------|----------------|---------------|
| Sound | 0.2136 | 0.0371 | DLX5 dominates: odontogenic context |
| Enamel caries | 0.2552 | 0.0459 | DLX5 still active |
| Deep caries | **0.1469** | **0.1469** | DLX5 collapses to FOS level — stress takeover |

**MSX1 → PTN regulatory weight:**

| Condition | MSX1→PTN weight |
|-----------|----------------|
| Sound | 0.0677 |
| Enamel caries | 0.1413 |
| Deep caries | **0.0449** | ← collapse in deep |

**All TFs → DSPP:** weight = 0 in all conditions. DSPP is not regulated by any measured TF in adult pulp because DSPP expression is absent.

**Key insight:** In sound and enamel caries, DLX5/MSX1 drive PTN expression with moderate weights. In deep caries, DLX5/MSX1 regulatory weights collapse (DLX5→PTN: 0.2552→0.1469) while FOS→PTN weight increases 3-fold (0.0371→0.1469). The regulatory apparatus controlling PTN has been taken over by the stress AP-1 program. The cells express DLX5/MSX1 protein but are no longer using them for their developmental function.

### 2.7 Adult Condition Comparison: % Positive by Condition

**Source:** `adult_pct_pos_by_condition.csv` (selected)

| Gene | Deep % | Enamel % | Sound % | Pattern |
|------|--------|---------|---------|---------|
| JUNB | ~87 | — | — | Extreme in deep |
| FOS | 72.15 | 87.87 | 87.64 | High in ALL adult conditions |
| JUN | 76.67 | 88.97 | 89.55 | High in ALL adult conditions |
| DLX5 | 38.00 | 57.94 | 38.43 | Peak in enamel caries |
| MSX1 | 42.34 | — | — | Paradoxically high in deep |
| BMP4 | 7.90 | 15.77 | 10.34 | Elevated vs. fetal OB |
| BMPR2 | 15.96 | 10.57 | 14.04 | Moderate, stable |
| DSPP | **0.00** | 0.17 | **0.00** | Absent in deep and sound |
| COL1A1 | 48.32 | 57.48 | 47.87 | Moderate in all adult |

**Observation:** FOS and JUN are equally high in sound (87–89%) as in deep caries (72–77%). They are not specific to the pathological state — they are baseline characteristics of adult dental pulp cells. JUNB is uniquely elevated in deep caries.

---

## 3. Validated and Falsified Hypotheses

### 3.1 FALSIFIED Hypotheses

**F1: "IGF1R marks SOX2+ dental stem cells"**
- IGF1R/SOX2 co-expression: 4/17,510 cells = 0.023%
- Fisher's exact p = 0.645, OR = 0.651 (CI: 0.166–1.857), Spearman ρ = −0.007
- **Verdict: FALSIFIED.** IGF1R is uniformly expressed across all mesenchymal clusters (~17–19%). SOX2 peaks at 0.59% in POB only. They are statistically independent.
- **Revised interpretation:** IGF1R is a broadly expressed mesenchymal competence factor, not a stem cell marker.

**F2: "DLX5/MSX1 → DSPP: direct progenitor pathway"**
- DLX5+MSX1+DSPP+ co-expressing cells in fetal POB (n=4,064): 0/4,064
- DLX5+MSX1+DSPP+ in fetal OB (n=681): 3/681 = 0.44% (within noise)
- DLX5 in adult deep: 38.00% — 52× higher than in fetal OB (0.73%), yet DSPP=0
- **Verdict: FALSIFIED.** DLX5/MSX1 are patterning/proliferation TFs that are actively downregulated before terminal differentiation. High DLX5/MSX1 in adult deep caries represents stress/fibroblast activation, not odontogenic competence.

**F3: "WNT/β-catenin (AXIN2/LEF1) is required for DSPP expression"**
- AXIN2+LEF1+DSPP+ triple co-expression = 0 in all 4 conditions tested
- In fetal DSPP+ cells: AXIN2 FC = 0.47, LEF1 FC = 0.51 (both depleted relative to DSPP−)
- **Verdict: FALSIFIED.** Canonical WNT signaling is inversely associated with terminal DSPP+ state.

**F4: "BMP4/BMP2/TGFB3 are deficient in adult deep caries"**
- BMP4: fetal OB 0.88% → adult deep 7.90% (−7.02 pp, endogenously elevated)
- BMP2: fetal OB 2.94% → adult deep 8.28% (−5.34 pp, elevated)
- TGFB3: fetal OB 0.88% → adult deep 4.85% (−3.97 pp, elevated)
- **Verdict: FALSIFIED.** These three molecules are not the limiting factor. Their endogenous levels in adult deep exceed fetal OB levels.

**F5: "SOSTDC1 is the limiting brake in adult pulp"**
- SOSTDC1 in adult deep: 0.10% (only 11 cells total across all adult conditions)
- SOSTDC1 is already virtually absent.
- **Verdict: FALSIFIED.** Cannot be the limiting brake because it is not present to be removed.

**F6: "JUN/FOS are purely pathological stress markers"**
- JUN in embryo: 71.40% (all cells), 94.37% in follicle, 54.66% in mesenchyme
- JUN in adult deep: 76.67% — nearly identical to embryo mean
- FOS in embryo: 45.63%, follicle 87.01%
- **Verdict: PARTIALLY FALSIFIED.** JUN and FOS are normal developmental TFs with high baseline expression in embryonic dental tissues. JUNB specifically (embryo mesenchyme 34.20% → adult deep 87.23%, Δ = +53 pp) is the pathological component.

### 3.2 CONFIRMED Hypotheses

**C1: "Fetal OB and adult deep caries are qualitatively different transcriptional states"**
- DSPP Δ = 34.36 pp (34.36% → 0%)
- COL1A1 Δ = 41.40 pp
- TGFB2 mean expression: 1.321 → 0.028 (47× drop)
- JUNB mean expression: 0.036 → 2.768 (76.9× increase)
- **Verdict: CONFIRMED with high confidence.** The magnitude of the differences across multiple independent markers confirms these are not the same cell state.

**C2: "Tooth regeneration requires coordinated multi-layer state transition, not a single effector molecule"**
- 5 independent deficits identified (DSPP, COL1A1/A2, TGFB2, IGF1R, RUNX2)
- 3 independent paradoxes (DLX5/MSX1 elevated but DSPP=0; BMP4 elevated but DSPP=0; JUNB extreme)
- No single molecule explains all deficits simultaneously
- **Verdict: CONFIRMED.** The multi-component nature of the failure is the most robust finding.

**C3: "PTN and MDK belong to distinct transcriptional programs"**
- NMF k=10: PTN → Module 1 (weight 10.92), MDK → Module 6 (weight 3.53)
- Spearman ρ(PTN, MDK) = 0.579 (correlated but distinct)
- **Verdict: CONFIRMED.** PTN is associated with stress/matrix response (FOS/JUN/COL1A2); MDK is associated with IGF/ECM remodeling (IGFBP5/IGFBP3/CLU).

**C4: "Pericyte/perivascular niche is depleted in deep caries"**
- ENG (CD105): 0.299 (sound) → 0.197 (deep), p = 1.4×10⁻⁴
- MCAM (CD146): 0.230 (sound) → 0.154 (deep), p = 0.002
- Module 10 (ACTA2/MYL9/RGS5): 0.100 → 0.081, p = 0.016
- **Verdict: CONFIRMED.** The perivascular stromal niche is systematically depleted in deep caries.

---

## 4. Updated ARC Layer Assessment for Dental Pulp

ARC framework (Architecture of Regenerative Control) applies to the dental pulp context as follows. Each layer is assessed based on the scRNA-seq evidence.

### 4.1 P-Layer (Permission — Epigenetic Accessibility)

**v1.1 prediction:** Chromatin must be accessible for odontogenic TFs.

**v1.3 data status:**
- DLX5 = 38%, MSX1 = 42% in adult deep. TF expression is high. This suggests chromatin is **not globally closed** for these TFs.
- However, DLX5/MSX1 in adult deep are in a stress/fibroblast regulatory context (FOS→PTN dominates, DLX5→PTN collapses), not in an odontogenic context.
- DSPP = 0% despite DLX5/MSX1 presence: the permission may be formally open but the signal routing is wrong.

**Assessment:** P-layer is partially functional but informationally misdirected. The epigenetic configuration supports stress/proliferation programs, not the odontogenic terminal differentiation program. **[INFERENCE, confidence 0.5]** — chromatin state at DSPP promoter specifically is not measured in these data.

### 4.2 C-Layer (Coordinates — Spatial Addressing)

**v1.3 data status:**
- Not directly measurable from scRNA-seq bulk data.
- ENG and MCAM depletion (perivascular markers) is consistent with disrupted bioelectric niche organization.
- Dental follicle is absent from adult pulp — the primary source of spatial TGFB2 signal (53.83% in fetal follicle) is structurally absent.

**Assessment:** C-layer compromised by structural absence of follicle-derived spatial signals. This may not be pharmacologically correctable. **[INFERENCE, confidence 0.4]**

### 4.3 S-Layer (Synchronization — Paracrine Phasing)

**v1.3 data status:**
- TGFB2: the single most depleted synchronization signal. In the embryo, it is primarily follicle-derived (53.83% of follicle cells vs. 18.78% of mesenchymal cells). Follicle is absent in adult pulp. This is a structural deficit, not a pharmacological one.
- BMP4: paradoxically elevated (7.90% vs. 0.88%). Over-elevated endogenous BMP4 without the counterbalancing TGFB2 context may actively suppress differentiation.
- DLX5→PTN regulatory weight collapses at deep caries while FOS→PTN takes over: the synchronization routing is disrupted.

**Assessment:** S-layer severely disrupted. Primary deficit is TGFB2 depletion due to structural absence of follicle cells. Secondary disruption is BMP4/TGFB2 imbalance.

### 4.4 B-Layer (Sandbox — Niche Isolation)

**v1.3 data status:**
- ENG/MCAM depletion: perivascular niche eroded.
- Module 10 (pericyte markers) significantly decreased (p=0.016).
- NFKB1: 1.76% in fetal OB → 15.96% in adult deep. NF-κB inflammatory axis elevated 9×.
- RELA: 0.15% → 10.53% (70× increase).
- JUNB 76.9× increase: extreme inflammatory state marker.

**Assessment:** B-layer severely compromised. The inflammatory/NF-κB/JUNB state represents the loss of sandbox protection. The pulp is locked in a systemic inflammatory attractor.

### 4.5 T-Layer (Termination — Controlled Completion)

**v1.3 data status:**
- DSPP = 0: terminal differentiation output is absent. This is not premature termination — it is absence of initiation.
- Ribosomal module (Module 4) decreased in deep caries (0.096 → 0.073, p=1.6×10⁻¹⁴): proliferation is already reduced.
- RUNX2: 24.23% → 16.01% (reduced but not absent). RUNX2 is not the terminal marker — it is pre-differentiation.

**Assessment:** T-layer is not prematurely active. The failure is upstream of termination. Termination cannot activate because differentiation never initiates.

### 4.6 M-Layer (Metabolism)

**v1.3 data status:** Not directly measurable from expression data.
- Vascular niche depletion (ENG/MCAM) implies compromised nutrient/oxygen delivery.
- **[INFERENCE, confidence 0.3]:** hypoxic conditions in deep caries may suppress OXPHOS-dependent differentiation programs.

---

## 5. Two-Component Failure Model (v1.3 Core)

Based on all data, the failure of odontoblast differentiation in adult deep caries pulp is described by two independent and co-occurring deficits:

### Component 1: Loss of Differentiation-Context Signals

| Signal | Fetal OB | Adult Deep | Fold deficit | Origin (embryo) |
|--------|----------|-----------|-------------|----------------|
| TGFB2 | 16.45% | 2.47% | **6.7×** | Follicle (53.83%) |
| IGF1R | 16.89% | 7.55% | 2.2× | Broadly expressed |
| RUNX2 | 24.23% | 16.01% | 1.5× | Follicle + mesenchyme |
| BMPR2 | 21.88% | 15.96% | 1.4× | — |

TGFB2 is primarily follicle-derived in the embryo. Adult pulp contains no follicle cells. This is a **structural absence**, not a pharmacological deficit. Exogenous TGFB2 could substitute for follicle-derived signal pharmacologically — but it cannot reconstitute the follicle niche architecture.

### Component 2: Dominance of JUNB-Driven AP-1 Stress State

| Marker | Embryo mesenchyme % | Adult Deep % | Fold increase |
|--------|---------------------|-------------|--------------|
| JUNB | 34.20 | 87.23 | **2.55×** |
| FOS | 42.82 | 72.15 | 1.69× (also high in embryo) |
| JUN | 54.66 | 76.67 | 1.40× (also high in embryo) |

JUNB is the specific pathological component. JUN and FOS are normally high in embryonic dental tissues. JUNB specifically exceeds even the highest embryonic reference (62.91% in fetal follicle) in adult deep caries.

**Why JUNB matters:** JUNB-containing AP-1 complexes compete with RUNX2 for binding at mineralization gene promoters (DSPP, DMP1, COL1A1). JUNB overexpression experimentally suppresses osteoblast/odontoblast differentiation while promoting fibroblast/inflammatory gene expression. **[INFERENCE based on published AP-1/RUNX2 competition literature, confidence 0.65]**

---

## 6. Falsifiable Predictions (v1.3)

All predictions are formulated with explicit success/failure criteria.

### Prediction 1 (retained from v1.2, confirmed endogenously)
**BMP4 is insufficient to drive DSPP in adult deep caries.**

- **Basis:** BMP4 is already 7.90% in adult deep (vs. 0.88% in fetal OB). DSPP = 0. Endogenous BMP4 is insufficient.
- **Experimental test:** Exogenous BMP4 addition (10–100 ng/mL) to hDPSC + TNF-α/IL-1β model.
- **Prediction:** No significant increase in DSPP+ cells above 0.5% at any BMP4 concentration.
- **Falsification criterion:** >1% DSPP+ cells with BMP4 alone would falsify this prediction.

### Prediction 2
**TGFB2 (5–10 ng/mL) will shift the transcriptional context toward RUNX2 expression in inflamed hDPSC more effectively than BMP4.**

- **Basis:** TGFB2 is the most depleted differentiation-context signal (47× by mean expression). BMP4 is endogenously elevated.
- **Experimental test:** hDPSC + TNF-α/IL-1β + TGFB2 (5 ng/mL) vs. + BMP4 (50 ng/mL) vs. untreated. Primary readout: RUNX2 and DSPP mRNA at day 7 and day 14.
- **Falsification criterion:** BMP4 produces equal or greater RUNX2/DSPP induction compared to TGFB2 in the inflammatory context.

### Prediction 3
**JUNB-specific suppression will be more effective than pan-JNK inhibition in shifting the pulp state toward odontogenesis.**

- **Basis:** JUN is developmentally normal (embryo: 71.4%); JUNB is specifically pathological (embryo: 39.6% → adult deep: 87.2%). Pan-JNK inhibition (SP600125) would suppress both, disrupting developmental JUN function.
- **Experimental test:** JUNB siRNA vs. SP600125 (pan-JNKi) in hDPSC + TNF-α/IL-1β model.
- **Primary readout:** DSPP, DMP1, RUNX2 at day 14.
- **Success criterion for JUNB siRNA:** >1% DSPP+ cells with JUNB siRNA vs. <0.5% with SP600125.
- **Falsification criterion:** SP600125 produces equal or better DSPP induction than JUNB siRNA.

### Prediction 4
**TGFB2 + IGF-1 combination will produce greater DSPP induction than either alone in inflamed hDPSC.**

- **Basis:** Both TGFB2 and IGF1R are independently depleted (6.7× and 2.2× respectively). They may act on different nodes of the differentiation circuit.
- **Experimental test:** 2×2 factorial design: ± TGFB2 (5 ng/mL) × ± IGF-1 (50 ng/mL) in hDPSC + TNF-α/IL-1β.
- **Success criterion:** TGFB2 + IGF-1 combination produces >1.5× more DSPP+ cells than either alone.
- **Falsification criterion:** No additive effect.

### Prediction 5
**DLX5 overexpression will NOT increase DSPP in adult pulp cells under inflammatory conditions.**

- **Basis:** DLX5 is already 38% in adult deep. Adding more DLX5 without changing the regulatory context (FOS→PTN dominance) should not shift the attractor state.
- **Experimental test:** DLX5 lentiviral overexpression in hDPSC + TNF-α/IL-1β.
- **Prediction:** <0.5% DSPP+ cells.
- **Falsification criterion:** >1% DSPP+ cells with DLX5 overexpression alone.

### Prediction 6 (new in v1.3)
**DLX5→PTN and FOS→PTN regulatory weights measured by scRNA-seq will be anti-correlated with DSPP expression at the single-cell level.**

- **Basis:** DLX5→PTN weight collapses and FOS→PTN weight increases at the population level in deep caries. If this applies at single-cell level, high FOS→PTN co-expression should predict DSPP=0 in the same cell.
- **Testable by:** Future annotated scRNA-seq dataset with both DSPP and DLX5/FOS co-measurement in individual cells.

---

## 7. Minimum Experimental Readout for Validity

Any experiment claiming to advance toward dentin repair based on this framework **must** include all six of the following readouts:

| # | Readout | Why required | Acceptable method |
|---|---------|-------------|-------------------|
| 1 | **Cell viability** | Niche integrity prerequisite | MTT/Live-Dead at baseline |
| 2 | **AP-1 status** | FOS, JUN, **JUNB** levels (not JUN alone) | RT-qPCR × 3 genes |
| 3 | **Competence markers** | TGFB2, IGF1R, RUNX2 | RT-qPCR × 3 genes |
| 4 | **Terminal output** | DSPP, DMP1 | RT-qPCR + immunostaining |
| 5 | **Matrix organization** | COL1A1 quantification + spatial | Immunostaining, confocal |
| 6 | **Proliferation / termination** | Ki67, p21 | Immunostaining |

A study reporting only Alizarin Red, or only DSPP qPCR, or only RUNX2 immunostaining is insufficient to claim state restoration according to this framework.

---

## 8. Position on Full Tooth Regeneration

The data analyzed in this project apply to **dentin repair** in adult dental pulp, not to de novo tooth formation.

**Why full tooth regeneration is not supported by these data:**

1. **Epithelial requirement:** Based on Wang et al. (2022, *Sci Bull*): even maximally competent embryonic CD24++ dental mesenchyme cells require dental epithelium to form a structured tooth under renal capsule. The GSE184749 dataset contains no epithelium. Adult pulp contains no epithelium.

2. **Temporal window:** Based on Hermans et al. (2023): earliest odontogenic fate commitment begins before 13 gestational weeks. Specific transcriptional states (Pitx2+, Runx2-low, SOX2-low) required for this commitment are absent from adult pulp by definition.

3. **Follicle absence:** TGFB2 is primarily follicle-derived (53.83% in fetal follicle). Dental follicle is absent from adult pulp. Full tooth regeneration requires follicle-like cells to provide this signal. No pharmacological substitute for follicle architecture is empirically validated.

4. **SOX2 niche depletion:** SOX2+ cells: sound 1.80% → enamel_caries 1.04% → deep_caries 0.78%. The minimal stem cell population present is further depleted under the conditions where intervention would be needed.

**The empirically supported near-term objective (v1.3):**

> Convert the local deep caries pulp from a stress-locked, TGFB2-depleted, JUNB-dominant state toward a TGFB2/IGF1R/RUNX2-competent state capable of DSPP+ odontoblast differentiation, without requiring full tooth-germ morphogenesis. This is dentin repair, not tooth regeneration.

---

## 9. Proposed Minimal Experimental Protocol

**Model:** human dental pulp stem cells (hDPSC) + inflammatory stimulus (TNF-α 10 ng/mL + IL-1β 5 ng/mL) for 72h, then intervention.

| Condition | Day 0 treatment | Day 3–14 treatment | Rationale |
|-----------|-----------------|-------------------|-----------|
| C1: Control | — | — | Baseline |
| C2: Inflammation | TNF-α + IL-1β | TNF-α + IL-1β | Deep caries model |
| C3: TGFB2 | TNF-α + IL-1β | + TGFB2 5 ng/mL | Restore depleted signal |
| C4: IGF-1 | TNF-α + IL-1β | + IGF-1 50 ng/mL | Restore IGF axis |
| C5: TGFB2 + IGF-1 | TNF-α + IL-1β | + TGFB2 + IGF-1 | Combinatorial |
| C6: JUNB siRNA | TNF-α + IL-1β | + JUNB siRNA | Address AP-1 imbalance |
| C7: TGFB2 + IGF-1 + JUNB siRNA | TNF-α + IL-1β | Triple combination | Full hypothesis test |
| C8: BMP4 (control) | TNF-α + IL-1β | + BMP4 50 ng/mL | Negative control (prediction 1) |

**Primary readouts at day 14 and day 21:**
- RT-qPCR: DSPP, DMP1, RUNX2, TGFB2, JUNB, FOS, JUN, COL1A1
- Immunostaining: DSPP, RUNX2, COL1A1 (spatial organization), Ki67
- Alizarin Red (mineralization, day 21)
- Flow cytometry: %DSPP+ cells

**Success criterion:**
- DSPP+ cells > 1% in ≥1 intervention condition (C3–C7)
- Concurrent JUNB reduction to < 40% (embryonic mesenchyme level)
- No significant increase in Ki67 (proliferation must not increase)

**Partial success criterion:**
- DSPP > 0.5% with TGFB2 alone (Prediction 2 partially confirmed)
- JUNB siRNA better than SP600125 for DSPP induction (Prediction 3 confirmed)

---

## 10. Cross-Verification with Published Literature

### 10.1 Alghadeer et al. 2023 (*Dev Cell*) — Human tooth development, 9–22 gw

This paper (mmc8.pdf in project) defines the human tooth developmental cell atlas using scRNA-seq + scATAC-seq at 5 time points (9–22 gestational weeks). Key alignments with our data:

- Their TopPath analysis identifies BMP and FGF as dominant signaling pathways in ameloblast development → consistent with our finding that BMP2 (not BMP4) is enriched in DSPP+ cells
- They successfully differentiate human iPSC → ameloblasts (ieAM) — confirming that the developmental program is accessible in vitro if the right context is provided → supports the v1.3 position that dentin-side repair is also achievable in vitro
- Their cell type atlas (odontoblasts = cluster 10, dental mesenchyme = cluster 9) is consistent with our GSE184749 cluster assignments

### 10.2 Wang et al. 2022 (*Sci Bull*) — Mouse tooth embryogenesis scRNA-seq

This paper (1s2.0S2095927322001025main.pdf) covers mouse odontogenesis from E10.5 to E16.5. Key alignments:

- Pitx2 and Fgf8 in epithelial component, Msx1 and Pax9 in mesenchymal component as key odontogenic genes — consistent with our finding that MSX1 is not sufficient alone for terminal differentiation
- The study identifies that molar/incisor differences are distinguishable only at E12.5 in mesenchyme — consistent with the tight spatio-temporal requirements for odontogenic commitment
- Paper states: "Whether there are other important genes for subsequent odontogenesis during tooth regeneration remains unclear" → our data contribute: TGFB2/IGF1R/RUNX2 as the differentiation-context module, WNT10A + BMP2 as DSPP+ enriched signals

### 10.3 Shi et al. 2024 (*Cell Prolif*) — Spatiotemporal cell landscape of human embryonic tooth development

This paper (Cell_Proliferation__2024__Shi__Spatiotemporal_cell_landscape_of_human_embryonic_tooth_development.pdf) provides spatial transcriptomic data on human embryonic tooth development. Alignments:

- Their identification of spatially distinct odontoblast zones is consistent with our cluster-based analysis
- Co-expression of BMP/TGF-β pathway components in the dental papilla region is consistent with our TGFB2 and BMP2 enrichment in DSPP+ cells

### 10.4 Consistency Summary

| Claim | This analysis | Literature support |
|-------|--------------|-------------------|
| MSX1 present but insufficient for DSPP | ✓ data | Wang et al. 2022 |
| BMP signaling required for tooth development | ✓ BMP2 enriched in DSPP+ | Alghadeer et al. 2023 |
| Epithelium required for full tooth | ✓ stated limitation | Wang et al. 2022, Alghadeer et al. 2023 |
| TGFB family in differentiation context | ✓ TGFB2 enriched in DSPP+ | Shi et al. 2024 |
| Developmental commitment requires specific windows | ✓ stated limitation | Hermans et al. 2023 |

---

## 11. Scientific Confidence Assessment

### High confidence (multiple consistent independent sources, effect sizes >10×)

- Fetal OB and adult deep caries are qualitatively different transcriptional states (DSPP Δ=34.36 pp; JUNB fold = 76.9×; TGFB2 mean fold = 47×)
- IGF1R and SOX2 are not co-expressed in dental mesenchyme (4/17,510 cells, p=0.645)
- DLX5/MSX1 in adult deep represent stress state, not odontogenic competence (DLX5 52× higher in adult deep than fetal OB, yet DSPP=0)
- WNT/β-catenin (AXIN2/LEF1) does not drive DSPP expression (triple co-expression = 0 everywhere)
- TGFB2 is the most depleted differentiation-context signal (47× mean expression drop)
- JUNB is specifically and pathologically elevated in adult deep caries (2.55× above embryo mesenchyme)
- BMP4, BMP2, TGFB3 are NOT limiting factors in adult deep (all elevated vs. fetal OB)

### Moderate confidence (single-dataset, effect sizes 2–10×, biologically plausible)

- TGFB2 depletion contributes causally to differentiation block (consistent with follicle-origin data and fold change; causal direction not proven by scRNA-seq)
- JUNB-specific suppression may be more targeted than pan-JNKi (biologically plausible from AP-1 literature; not tested in dental context)
- PTN and MDK belong to distinct transcriptional programs (consistent NMF finding across k values)
- Pericyte/perivascular niche depletion contributes to sandbox failure (ENG/MCAM depletion significant but mechanism not established)

### Low confidence (single indirect evidence, or not testable from available data)

- Any specific molecular combination can restore DSPP+ differentiation in adult pulp (untested)
- Adult human pulp can support tooth-germ-like morphogenesis (no supporting data)
- GSE185222 findings are cell-type-specific (no annotation available — all findings are population-level)
- JUNB competes with RUNX2 at DSPP promoter specifically in dental context (extrapolated from osteoblast literature)

---

## 12. Hard Limitations (Non-Negotiable Caveats)

1. **No cell-type annotation in GSE185222.** All adult pulp findings are bulk-population. Individual cell types (fibroblasts, immune cells, endothelium, progenitors) are not separated. A finding like "JUNB = 87% in adult deep" describes the average across all cell types present.

2. **Adult DSPP+ cells: n=3.** This sample size is statistically insufficient for any quantitative conclusion. All statements about "adult DSPP+ cells" are descriptive only.

3. **GENIE3 weights are not causal.** A high DLX5→PTN weight means DLX5 expression co-varies with PTN expression. It does not mean DLX5 regulates PTN transcription.

4. **Cross-sectional data only.** Causal ordering of events (does TGFB2 depletion cause JUNB rise, or vice versa?) cannot be established from these data.

5. **No epithelium in GSE184749.** Epithelial-mesenchymal signaling, required for full tooth morphogenesis, cannot be assessed from this dataset.

6. **Dental follicle absence in adult pulp.** TGFB2 is primarily follicle-derived (53.83%). Follicle cells are absent from adult pulp by anatomy. This structural absence may not be correctable by pharmacological TGFB2 alone, because follicle provides spatial/temporal context beyond TGFB2 concentration.

7. **All R analyses by single author.** No inter-rater validation of cluster assignments or module interpretations has been performed.

---

## 13. Summary Table: What ARC-Tooth v1.3 Claims vs. Does Not Claim

| Statement | Claimed? | Evidence basis |
|-----------|---------|----------------|
| Adult deep pulp and fetal odontoblasts are different states | **YES** | Quantitative, high confidence |
| TGFB2 is the most depleted differentiation signal | **YES** | Data, 47× mean fold drop |
| JUNB is specifically elevated in deep caries | **YES** | Data, 76.9× mean fold increase |
| BMP4 addition will restore DSPP | **NO** | Data shows BMP4 already elevated |
| DLX5 addition will restore DSPP | **NO** | Data shows DLX5 already 52× above fetal OB |
| TGFB2 + IGF-1 will restore DSPP in vitro | **PREDICTED, untested** | Inference, confidence 0.5 |
| JUNB siRNA will restore DSPP | **PREDICTED, untested** | Inference, confidence 0.5 |
| Adult pulp can form a full tooth | **NOT CLAIMED** | No supporting data |
| Dentin repair is possible in vivo | **NOT CLAIMED** | No in vivo data |
| Single molecule can restore full odontoblast state | **ACTIVELY CONTRADICTED** | Multi-component failure confirmed |

---

## 14. Repository Note

All R scripts producing the data summarized in this document are available at:
**https://github.com/Jenja-N/R**

Key output files verified in this document:

| File | Contents | Used in section |
|------|---------|----------------|
| `fetal_OB_vs_adult_deep.csv` | 26 genes, fetal OB vs. adult deep | 2.1 |
| `tgfb_bmp_stress_comparison.csv` | 10 genes, TGFB/BMP/AP-1 comparison | 2.2 |
| `embryo_vs_adult_deep_comparison.csv` | 9 genes, global embryo comparison | 2.3 |
| `wnt_signaling_summary.csv` | WNT/DSPP co-expression across 4 conditions | 2.4 |
| `fetal_dspp_signature.csv` | 28 genes, DSPP+ vs. DSPP− in fetal | 2.5 |
| `regulon_comparison.csv` | 66 TF→target pairs, 3 conditions | 2.6 |
| `adult_pct_pos_by_condition.csv` | 27 genes, adult 3 conditions | 2.7 |
| `fetal_pct_pos_by_cluster.csv` | 27 genes, fetal 5 clusters | — |
| `fetal_POB_vs_adult_deep.csv` | 22 genes | — |
| `ob_dspp_positive_profile.csv` | OB-cluster DSPP+ signature | — |
| `dspp_positive_profile.csv` | All-fetal DSPP+ signature | — |
| `embryo_key_genes_expression.csv` | 11,218 cells × 16 genes | — |

---

## 15. Conclusion

ARC-Tooth v1.3 provides the following maximally defensible summary:

**What is known (high confidence, multiple sources):**

The failure of odontoblast differentiation in adult deep caries is a **combinatorial state failure**, not a single-molecule deficiency. The adult deep caries pulp is locked in a JUNB-dominant AP-1 stress attractor state with structural depletion of TGFB2 (the primary follicle-derived differentiation-context signal) and loss of IGF1R/RUNX2 competence. This state is incompatible with DSPP+ terminal differentiation despite the presence of upstream TFs (DLX5, MSX1) at levels exceeding fetal odontoblasts.

**What is hypothesized (moderate confidence, untested):**

Simultaneous correction of TGFB2 depletion and JUNB-dominant AP-1 imbalance, in the presence of adequate IGF1R axis support, may create conditions sufficient for DSPP+ odontoblast differentiation in vitro. This is the testable hypothesis of ARC-Tooth v1.3.

**What is not supported by data:**

Full tooth regeneration from adult dental pulp. Any single-molecule intervention sufficient for dentin repair. Clinical translation at any timescale.

**The next required step:**

Execute the minimal experimental protocol (Section 9) in validated hDPSC models to determine whether the predicted state transition is achievable. All six readouts are required. Alizarin Red alone is not sufficient.

---

*ARC-Tooth v1.3 — June 2026 — https://github.com/Jenja-N/R*  
*No clinical claims. No therapeutic recommendations. Hypothesis-generating framework only.*
