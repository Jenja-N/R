# ARC 7.0 — Architecture of Regenerative Control

> **Status:** Working research document — not a clinical protocol, not a validated theory  
> **Scope:** Comparative regenerative biology, developmental biology, systems biology  
> **License:** CC-BY 4.0  
> **Zenodo:** [10.5281/zenodo.20526831](https://doi.org/10.5281/zenodo.20526831)

---

## Core principle

No claim in this document should be stronger than the evidence it rests on.

ARC proposes that adult mammalian regeneration is best understood not as the firing of a single "regeneration gene," and not as raw proliferation, but as a **controlled transition of damaged tissue into a coordinated morphogenetic state**. Reaching that state requires the joint satisfaction of six functional tasks: **P, C, S, B, M, T**.

These are **functional categories**, not anatomical structures or genes. A given molecule or tissue can serve more than one gate; different species can realize the same gate through different mechanisms.

---

## What ARC Is (and Is Not)

**ARC IS:**
- A falsifiable functional meta-model of regenerative control
- An operational vocabulary precise enough to generate testable predictions
- A structured way to compare regenerative systems across taxa
- A research program capable of testing necessity, interactions, and tissue-specific realizations of its gates

**ARC IS NOT:**
- A predictive theory (it is hypothesis-generating, not yet prospectively validated)
- A mechanistic simulation (it does not model dynamics)
- A clinical tool (not validated for medical decision-making)
- A claim that any single intervention is sufficient or safe

---

## The Six Gates — Operational Definitions

| Gate | Symbol | Function | Includes | Is NOT |
|:---|:---|:---|:---|:---|
| **Permission / Competence** | P | Ability of relevant cells to enter a plastic, proliferative, differentiation-competent state | Suitable progenitors; accessible regulatory chromatin; reversible cell-state switching; intact organ-specific networks; absence of irreversible senescent block | Global epigenetic de-repression (risks identity loss or transformation) |
| **Coordinates / Positional Information** | C | Specifying what, where, in what orientation, topology, and to what size something should be built | Stromal positional memory; regional transcriptional codes; epithelial polarity; geometry of residual niche; signaling gradients; mechanical fields; bioelectric states; ECM spatial organization | A HOX code alone, or membrane voltage alone, or tissue mass alone |
| **Synchronization** | S | Coordinating cell populations, signals, and phases in time and space | Neural, endocrine, vascular, paracrine, immune, electrical channels — used in different combinations by different systems | A single universal channel (e.g., "nerve is required for all regeneration" is false) |
| **Boundary / Permissive Microenvironment** | B | Keeping the local environment compatible with morphogenesis rather than rapid scarring, infection, or mechanical collapse | Mechanical offloading; geometric constraint; moisture; matrix stiffness/composition; controlled proteolysis; immune–matrix interaction; protection from premature epithelialization | Total suppression of inflammation or TGF-β (human fingertip regenerates with substantial inflammatory activity) |
| **Metabolic & Material Support** | M | Supplying oxygen, substrates, minerals, vascular access, redox control, biosynthetic and energetic capacity | Vascular contribution; stage-specific metabolic regimes; mineral supply | A single Warburg-type switch |
| **Termination & Stabilization** | T | Stopping growth, stabilizing identity, restoring barrier function, transitioning to homeostasis | Cell-cycle control; mechanosensing/contact inhibition; Hippo/YAP–TAZ; p53/Rb; differentiation; apoptosis; immune resolution; vascular maturation; ECM remodeling | The hypothesis that "the organism will stop growth on its own" — this is unproven |

**Key distinction:** One mechanism can serve several gates simultaneously (vasculature contributes to M, C, and S). This functional overlap is a real property of biology, not evidence that the categories are statistically independent.

---

## Retraction Ledger (v4.0 → v7.0)

Earlier versions of this project (v1.0–v4.0, archived in commit history) contained claims that did not survive independent verification. They are listed here so they are not silently reintroduced.

| # | Retracted claim | Correction |
|:---|:---|:---|
| 1 | Fixed logistic-regression weights: B=1.78, P=1.45, S=1.42, T=1.20, M=0.93, C=0.63, S×B=+0.85, bias=−4.52 | Weights were single-author expert judgment with no independent inter-rater validation; N=16 is far too small for this claim; not a biological constant. **Retracted in full.** |
| 2 | "100% LOOCV accuracy on 16 systems" | Explaining known phenomena post hoc is not proof. The same phenomena are independently explained by named mechanisms published before and independently of ARC. **Retracted.** |
| 3 | "Out-of-sample PASS on zebrafish heart regeneration (Φ=0.871)" | No corresponding computation exists in the project files. **Retracted.** |
| 4 | S×B interaction term with ΔAIC=−10.6 on a "curated 16-system, 6-species matrix" | Neither the underlying analysis nor the AIC comparison could be located in the project's data, scripts, or prior documents. The cited weights are the same unvalidated numbers retracted in item 1. **Retracted in full.** |
| 5 | Specific illustrative per-tissue gate scores (e.g., adult skin vs. oral mucosa, Φ=0.128 vs. 0.772) presented as worked examples from real data | No such per-tissue gate scoring exists in the project's actual scripts or outputs; this reads as an invented illustration, not a measurement. **Retracted.** |

**Current honest formal statement:**

```
R(t) = F_tissue [ P(t), C(t), S(t), B(t), M(t), T(t), interactions(t) ]
```

where `R(t)` is a multidimensional outcome vector (volume, geometry, tissue composition, vascularity, innervation, mechanical function, physiological function, stability, safety) and `F_tissue` is an **unknown, tissue- and species-specific function**. A threshold or logistic approximation is appropriate only after each gate is measured by independent operational markers and the model is tested out of sample.

---

## Three Independent Convergences (2025–2026)

The following papers were published without knowledge of ARC. Their findings are compatible with the six-gate architecture:

### 1. P — Permission/Competence
**"Reactivation of mammalian regeneration by turning on an evolutionarily disabled genetic switch"**  
*Science*, 26 June 2025. DOI: [10.1126/science.adp0176](https://doi.org/10.1126/science.adp0176)  
Lin W. et al., NIBS China

- Rabbits regenerate ear pinna; mice and rats do not.
- Cause: insufficient retinoic acid (RA) production due to evolutionary inactivation of *Aldh1a2* enhancers in murine rodents.
- Switching on *Aldh1a2* or RA supplementation **reactivated regeneration in mice**.
- **ARC reading:** Regeneration in mammals is not genetically lost; it is blocked at the regulatory level (P-gate).

### 2. S × B — Synchronization × Boundary
**"Digit regeneration in mice is stimulated by sequential treatment with FGF2 and BMP2"**  
*Nature Communications*, 17 April 2026. DOI: [10.1038/s41467-026-72066-8](https://doi.org/10.1038/s41467-026-72066-8)  
Yu L., Muneoka K. et al., Texas A&M University

- Non-regenerative amputation level in mice was converted to epimorphic regeneration by **sequential** FGF2 → BMP2 treatment.
- FGF2 induces blastema formation; BMP2 drives morphogenesis and differentiation.
- **ARC reading:** Regenerative failure is not limited by competent cells (P), but by the absence of correctly sequenced signals in the wound environment (S × B).

### 3. P — Permission/Competence (intercellular block)
**"Restoration of retinal regenerative potential of Müller glia by disrupting intercellular Prox1 transfer"**  
*Nature Communications*, 25 March 2025. DOI: [10.1038/s41467-025-58290-8](https://doi.org/10.1038/s41467-025-58290-8)  
Lee E.J., Kim J.W. et al., Korea

- Mammalian Müller glia fail to regenerate retina because **Prox1 protein is transferred from neighboring neurons** into glia, suppressing reprogramming.
- Blocking this transfer restores retinal progenitor potential.
- **ARC reading:** A specific intercellular suppression mechanism operates in mammals that is absent in regenerating species (P-gate block).

**None of these groups cite ARC.** This is consilience — independent lines of evidence converging on a common structure.

---

## Epistemic Standard

Every claim in ARC is tagged by evidentiary strength:

| Tier | Symbol | Criterion |
|:---|:---|:---|
| **Established fact** | F | Directly observed, ideally replicated or obtained interventionally |
| **Strict inference** | I | Logical consequence of facts, requiring no new biological assumption |
| **Working hypothesis** | H | Consistent with the facts, generates a testable prediction |
| **Engineering guess** | E | A specific proposed intervention, not yet shown sufficient or safe |

"Consistent with the data" means just that — compatibility, not causal proof.

---

## How to Falsify ARC

ARC is **not** falsified by the failure of any single intervention cocktail. It is falsified by results such as:

1. Full regeneration at B = 0 — organized regeneration in a stable contractile, collagen-I-rich, α-SMA-high fibrotic environment.
2. Full anatomically correct regeneration at C = 0 — without any positional/coordinate hub.
3. S × B fails to improve prediction in an independent, preregistered dataset compared to a main-effects-only model.
4. S alone reliably produces full, safe, organized form without a permissive boundary or coordinate input.
5. B + C + S does not outperform simpler interventions across several independent injury models.

See `FALSIFICATION_MATRIX.md` for pre-registered test designs.

---

## Safety

ARC is **not** grounds for self-administered or unsupervised clinical use of HDAC inhibitors, senolytics, mTOR activators, Wnt/BMP agonists, TGF-β modulators, growth factors, electrical/bioelectric stimulation devices, or anti-inflammatory cocktails.

Principal risks that any ARC-derived protocol must treat as primary, multidimensional safety endpoints: neoplasia, dysplasia, odontoma/ectopic mineralization, chronic infection, alveolar bone resorption, uncontrolled proteolysis, ischemia, mis-innervation and chronic pain, ankylosis, non-functional tissue mass, incorrect form, and unstable termination.

---

## Criteria for ARC to Graduate from Heuristic to Validated Theory

ARC moves from "functional meta-model" to "supported theory" only when, jointly:

1. Stable operational definitions exist for P/C/S/B/M/T, fixed before analysis of new data.
2. Multiple independent raters can score them with high inter-rater agreement.
3. ARC generates a prediction for a new system **in advance** of seeing the outcome.
4. That prediction is confirmed out of sample.
5. ARC outperforms simpler alternative models on the same data.
6. At least one intervention experiment confirms a structural prediction unique to ARC.
7. The result is independently replicated.
8. Negative results and the limits of applicability are published, not omitted.
9. A concrete experiment exists that could falsify the architecture itself.
10. Terminology and gate definitions are not changed after the fact to rescue the model.

**Until these are met, the correct status of ARC is:** A falsifiable functional meta-model of regenerative control, with strong comparative-biology motivation, but without completed prospective validation.

---

## Citation

If you use ARC in your work, please cite:

> ARC 7.0 — Architecture of Regenerative Control. A falsifiable framework for adult mammalian regeneration. Zenodo. https://doi.org/10.5281/zenodo.20526831

---

## Provenance

This document is the product of an iterative, adversarial verification process: successive drafts were checked against primary sources, fabricated statistics and misattributed citations were identified and removed, and genuine improvements were retained. Every new numeric claim should specify whether it is a fact, a strict inference, a hypothesis, or an engineering guess, and every new citation should be checked against the actual source before being added.

---

*This is a research notebook, not medical guidance. Nothing in this document authorizes self-experimentation or unsupervised clinical application of any compound, device, or protocol.*
