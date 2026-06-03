# ARC: Architecture of Regenerative Control

![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20526831.svg)
(https://doi.org/10.5281/zenodo.20526831)

**A computational ontology and semi-quantitative scoring framework for comparative regenerative biology**

> **Status:** Preprint submitted to bioRxiv (Systems Biology)  
> **Model version:** v4.0 (logistic regression + S×B interaction)  
> **Training systems:** N=16 (1 independent out-of-sample validation)  
> **License:** CC-BY 4.0  
> **Last updated:** 2026-05-30

---

## What ARC Is (and Is Not)

ARC is a **conceptual ontology** that decomposes regenerative capacity into six biologically-defined layers (Permission, Coordinates, Synchronization, Sandbox, Termination, Metabolism) and a **semi-quantitative scoring framework** that estimates these layers on a [0,1] scale from curated literature.

**ARC is NOT:**
- A predictive theory of regeneration (it is a hypothesis-generating framework)
- A mechanistic model (it does not simulate dynamics)
- A clinical tool (it is not validated for medical decision-making)
- A direct measurement of any biological quantity (layer scores are expert-curated estimates, not experimental readouts)

**ARC IS:**
- A structured vocabulary for comparing regenerative systems across taxa
- A formalization of the hypothesis that regeneration is a multi-constraint computational state
- A scoring rubric that makes explicit which biological features are considered necessary for regeneration
- A falsifiable framework with pre-registered failure criteria (see `FALSIFICATION_MATRIX.md`)

---

## Layer Definitions

| Layer | Symbol | Conceptual Definition |
|-------|--------|----------------------|
| **Permission** | P | Epigenetic accessibility of regenerative gene programs |
| **Coordinates** | C | Positional memory and morphogenetic vector fidelity |
| **Synchronization** | S | Temporal and spatial coordination of signaling cascades |
| **Sandbox** | B | Anti-fibrotic permissive microenvironment |
| **Termination** | T | Safe proliferation capacity without senescence/malignancy |
| **Metabolism** | M | Bioenergetic support for biosynthesis |

Each layer is scored in [0,1] by expert curation of peer-reviewed literature, relative to a reference maximum (axolotl limb blastema ≈ 1.0).

**Critical caveat:** Layer scores are **semi-quantitative expert estimates**, not direct measurements. They represent our best synthesis of published qualitative and quantitative data. Two independent experts may produce different scores for the same system. All scores are traceable to specific DOIs in `data/training_set_FINAL.csv`.

---

## Mathematical Form (v4.0)

$$ \Phi = \sigma\left( \sum_{i=1}^{6} w_i x_i + w_{SB} \cdot (S \times B) + b \right) $$

where:
- $\sigma(z) = 1/(1 + e^{-z})$ is the logistic sigmoid
- $\mathbf{x} = [P, C, S, B, T, M] \in [0,1]^6$
- $w_i$ are learned main-effect weights
- $w_{SB}$ is the learned S×B interaction coefficient
- $b$ is a bias term
- Decision threshold: $\Phi > 0.5$ predicts regeneration; $\Phi < 0.5$ predicts fibrosis/no regeneration

**Fitted weights (N=16):**

| Parameter | Value | Rank |
|-----------|-------|------|
| B (Sandbox) | 1.78 | 1 |
| P (Permission) | 1.45 | 2 |
| S (Signaling) | 1.42 | 3 |
| T (Termination) | 1.20 | 4 |
| M (Metabolism) | 0.93 | 5 |
| C (Coordinates) | 0.63 | 6 |
| S×B | +0.85 | interaction |
| Bias | −4.52 | — |

---

## Current Performance (Honest Assessment)

| Metric | Value | Caveat |
|--------|-------|--------|
| Training accuracy (N=16) | 16/16 = 100% | Expected for any flexible model on training data |
| LOOCV accuracy (N=16) | 16/16 = 100% | Strong but limited by N |
| Out-of-sample validation | Zebrafish heart Φ=0.871 ✓ | Single test; more needed |
| S×B interaction ΔAIC | −10.6 | Statistically justified |

**Honest interpretation:** These metrics demonstrate internal consistency, NOT predictive power. The framework has NOT been validated on held-out experimental data. It should be considered a **proof-of-concept ontology**, not a deployed predictive model.

---


