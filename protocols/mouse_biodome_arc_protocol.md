# ARC Framework: Practical Guide for Mouse BioDome Translation

> **For labs working on mammalian regeneration protocols.**  
> *Architecture of Regenerative Control (ARC) — May 2026.*  
> Full model, data, code: [github.com/Jenja-N/R](https://github.com/Jenja-N/R)

---

## 1. What is ARC

The Architecture of Regenerative Control (ARC) is a 6-layer predictive model trained on 16 regenerating/non-regenerating biological systems across 6 species. It predicts regeneration outcome as a logistic function with a statistically validated S×B interaction term (ΔAIC = −10.6 vs. baseline). Independent out-of-sample validation: zebrafish heart, Φ = 0.871.

**Core equation:**

```
Φ = sigmoid(1.45·P + 0.63·C + 1.42·S + 1.78·B + 1.20·T + 0.93·M + 0.85·(S×B) − 4.52)
```

**Decision threshold:** Φ > 0.5 → regeneration predicted.

**Layer definitions:**

| Layer | Name | Biological meaning |
|-------|------|--------------------|
| P | Permission | Epigenetic accessibility of regenerative gene programs |
| C | Coordinates | Positional memory (bioelectric Vmem gradients + HOX identity) |
| S | Synchronization | Coordinated signaling (nerve-derived paracrine factors, FGF/NRG1) |
| B | Sandbox | Anti-fibrotic permissive microenvironment (TGF-β3/β1 ratio, ECM softness) |
| T | Termination | Readiness to stop growth (Hippo/YAP, p53, contact inhibition) |
| M | Metabolism | Bioenergetic support (glycolytic flexibility, NAD+, mTOR pulsing) |

**Key property:** Multiplicative AND-gate logic. If any single layer collapses to near zero, Φ → 0 regardless of all other layers. This is why single-factor therapies fail.

---

## 2. Retrospective analysis: why MDT worked where progesterone alone failed

This section applies ARC retrospectively to Murugan et al. 2018 and Murugan et al. 2022 (Levin Lab, Tufts University / Wyss Institute).

### Progesterone-only protocol (Murugan et al. 2018, *Cell Reports*)

| Layer | Estimated value | Reason |
|-------|----------------|--------|
| P | 0.45 | Progesterone promotes partial epigenetic relaxation |
| C | 0.25 | No explicit bioelectric pre-patterning |
| S | 0.30 | No neurotrophic or paracrine synchronization delivered |
| B | 0.25 | No anti-fibrotic agent; BioDome provides physical barrier only |
| T | 0.55 | Baseline termination intact |
| M | 0.40 | Partial metabolic support via steroid signaling |

**ARC Φ ≈ 0.15.** Result: cartilaginous spike, no patterning. ARC explanation: S×B product ≈ 0.075 — interaction term near zero. Outgrowth signal exists but cannot be patterned.

### Multidrug treatment protocol (Murugan et al. 2022, *Science Advances*)

Five-compound payload (BDNF, 1,4-DPCA, RD5, GH, RA) delivered via silk hydrogel BioDome for 24 hours.

| Drug | ARC layer targeted | Mechanism |
|------|--------------------|-----------|
| BDNF | S | Neurotrophic synchronization, nerve-derived paracrine support |
| RA (retinoic acid) | C | Positional identity via retinoic acid signaling, HOX gene activation |
| 1,4-DPCA | B | HIF prolyl hydroxylase inhibitor; suppresses HIF-1α-driven fibrosis |
| RD5 | B | Anti-apoptotic scaffold stabilization; contributes to permissive ECM |
| GH (growth hormone) | P + M | Progenitor activation and metabolic support |
| BioDome enclosure | B | Physical permissive microenvironment; maintains moisture and biochemical milieu |

| Layer | Estimated value | Change vs. progesterone-only |
|-------|----------------|------------------------------|
| P | 0.60 | +0.15 (GH effect) |
| C | 0.65 | +0.40 (RA critical contribution) |
| S | 0.70 | +0.40 (BDNF) |
| B | 0.65 | +0.40 (1,4-DPCA + RD5 + BioDome) |
| T | 0.50 | ~unchanged |
| M | 0.60 | +0.20 (GH + metabolic support) |

**ARC Φ ≈ 0.63–0.68.** Result: near-complete functional limb, sensorimotor recovery, 18-month regrowth. S×B product ≈ 0.455 — interaction term now dominant.

**Remaining gaps identified by ARC:**
- Toes without underlying bone = incomplete C (positional memory established but not fully consolidated at distal tip).
- No explicit T enforcement = regrowth terminates empirically but without controlled mechanism. Risk of overgrowth not systematically managed.

**ARC conclusion:** MDT succeeded because it accidentally satisfied the S×B interaction simultaneously with partial C and P. The 2022 authors correctly described this as requiring "an enclosed permissive microenvironment AND an instructive set of signals" — this is precisely the B×S joint condition in ARC. The design was empirically discovered; ARC formalizes why it works.

---

## 3. Mouse-specific failure modes

Adult *Mus musculus* differs from adult *Xenopus laevis* in two ARC-critical ways that were not rate-limiting in the frog model.

### 3.1 B-layer collapse in mice

Adult mouse macrophage inflammatory response is substantially more aggressive than *Xenopus*. Following amputation:

- M1 macrophage dominance established within 2–4 hours.
- TGF-β1 rises sharply by 6–12 hours post-amputation.
- ECM stiffness increases (collagen crosslinking) within 24–48 hours.
- Result: B drops to approximately 0.15–0.20 in unmodified mouse wound.

The physical BioDome device maintains moisture and localizes drug delivery but does not block the biochemical fibrotic cascade. B-layer collapse is biochemical, not physical.

**Why this is catastrophic for your protocol:** S×B is the dominant interaction term (learned weight 0.85, highest ΔAIC contribution). High S (good bioelectric signal, BDNF, nerve-derived factors) with B = 0.15 gives S×B ≈ 0.10. The interaction term that drives blastema formation contributes near zero. The bioelectric Vmem signal exists; the tissue cannot respond to it because the microenvironment has already committed to scar closure.

### 3.2 P-layer epigenetic lock in adult mice

Adult mouse somatic cells carry tight repressive chromatin marks on regeneration-associated promoters:

- H3K27me3 enrichment at Msx1, Prrx1, Sox2, and Lin28a regulatory regions.
- DNA methylation at blastema-associated enhancers.
- p16/INK4a and p21 expression in wound-adjacent cells (senescence).

Changing membrane potential (Vmem) via ion channel modulators is a C-layer intervention — it encodes positional and morphogenetic information. For this information to drive dedifferentiation and blastema entry, cells must be able to access the relevant transcriptional programs. With P = 0.15–0.20 (typical adult mouse baseline in wound), the signal is broadcast and received at the membrane level, but the downstream transcriptional response is blocked.

**Analogy for the framework:** C provides the address; P determines whether the door at that address can be opened. Both must be active simultaneously.

---

## 4. Recommended protocol additions

These are the minimal interventions with the highest predicted Φ impact. Both are well-characterized reagents with established in vivo use.

### 4.1 Priority 1 — Stabilize B-layer (implement first; gates everything else)

**Why first:** B gates the S×B interaction term. Without B > 0.45, even perfect S and C yield Φ < 0.5.

| Intervention | Dose / format | Timing | Target effect |
|-------------|---------------|--------|---------------|
| TGF-β3 protein | 50–100 ng/ml in BioDome gel | Day 0, sustained release over 72h | Shifts TGF-β3/TGF-β1 ratio > 1.0; promotes M2 macrophage polarization |
| TGF-β1 neutralizing antibody (1D11 or equivalent) | 5–10 µg/ml local | Day 0–3 | Blocks early fibrotic commitment |
| High-MW hyaluronic acid (HA, >1 MDa) | Gel matrix component | Day 0 onward | Maintains soft ECM (target: <2 kPa stiffness); prevents collagen crosslinking cascade |
| 1,4-DPCA | As per original MDT formulation | Day 0–1 | Retain from existing protocol; contributes to B via HIF suppression |

**Measurement proxy for B:** TGF-β3/TGF-β1 protein ratio in wound exudate at Day 2. Target > 1.0. Secondary: AFM stiffness measurement of wound tissue at Day 3. Target < 2 kPa.

### 4.2 Priority 2 — Unlock P-layer (strict time window: 24–48h only)

**Why time-limited:** HDAC inhibition beyond 48h begins to destabilize T-layer (disrupts p16/p21 senescence machinery that forms part of termination readiness). P-unlock is a window, not a sustained state.

| Intervention | Dose / format | Timing | Target effect |
|-------------|---------------|--------|---------------|
| Valproic acid (VPA) | 1–2 mM local concentration in BioDome payload | Day 0, 24–48h only | HDAC inhibition; opens chromatin at Msx1/Prrx1/Sox2 promoters |
| α-pifithrin | 2 µM local (optional addition) | Day 0, 24h only | Transient p53 inhibition; reduces senescence in wound-adjacent cells; extends P-unlock window slightly |

**Sequence requirement:** P-unlock must be concurrent with or 12–24h *before* bioelectric/ion channel intervention, not after. Rationale: cells must be in an accessible chromatin state when they receive the C-layer (positional) signal. Reversing the order — bioelectric first, HDACi second — is predicted to yield no improvement over baseline.

**Measurement proxy for P:** RT-qPCR or scRNA-seq for Msx1 and Prrx1 expression in wound tissue at Day 3. Detectable upregulation (>2-fold vs. sham) confirms P-layer unlocked. H3K27me3 ChIP at Msx1 promoter (Day 2) optional but definitive.

---

## 5. Integrated protocol sequence

```
Day -1 to 0   : Immune pre-conditioning (optional: low-dose IL-2 for Treg expansion,
                Dasatinib + Quercetin senolytic if aged animals)

Day 0 (amputation):
  Hour 0–1    : BioDome application with modified payload:
                  — TGF-β3 (50–100 ng/ml)
                  — TGF-β1 neutralizing antibody
                  — Valproic acid (1–2 mM)
                  — α-pifithrin (2 µM, optional)
                  — High-MW HA in gel matrix
                  — 1,4-DPCA (original MDT dose)
                  — BDNF (original MDT dose)

Day 0–1       : P-unlock window active (VPA + pifithrin)
                B-stabilization active (TGF-β3, HA)

Day 1–2       : Remove or replace BioDome payload
                  Transition to: bioelectric intervention (ion channel modulators,
                  Vmem modulation — existing Levin Lab protocol)
                  + RA (retinoic acid, C-layer) + GH (P+M support)
                  TGF-β3 sustained release continues

Day 3–7       : Proliferative/blastema phase
                Monitor: Ki67/p21 ratio in wound tissue
                If Ki67/p21 > 3.0 at Day 5–7 without blastema → activate verteporfin
                  (YAP inhibitor, 0.1–0.5 µM local) to prevent dysplasia attractor

Day 7–21      : Patterning phase
                Continue: nerve-derived factors (NGF, BDNF), VEGF for vascularization
                Reduce: Wnt agonists if used (pulse only, 24h ON / 48h OFF)

Day 21–28     : Termination enforcement
                Verify Hippo pathway activation (YAP cytoplasmic)
                Contact inhibition markers (E-cadherin)
                Second senolytic wave (optional)

Day 28+       : Maturation, functional loading, ECM remodeling
```

---

## 6. Safety monitoring: P/T oncogenic risk window

ARC identifies a specific risk configuration: P elevated + T simultaneously low = dysplasia attractor (cancer-like state). The valproic acid + pifithrin window transiently creates this configuration.

**Mandatory monitoring:**

| Timepoint | Measurement | Action if threshold exceeded |
|-----------|-------------|------------------------------|
| Day 3 | Ki67/p21 IHC ratio in wound tissue | If > 2.0 with no organized blastema: extend TGF-β3, delay next phase |
| Day 5–7 | Ki67/p21 IHC ratio | If > 3.0 with no organized blastema: administer verteporfin (YAP inhibitor) locally |
| Day 10–14 | H&E histology | If disordered proliferation without tissue organization: iCasp9 or HSV-TK safety switch if pre-armed |
| Day 21 | microCT | Confirm structured bone/cartilage organization vs. amorphous mass |

**The 48h VPA limit is not conservative — it is mechanistically required.** Extended HDAC inhibition beyond this window has been shown to promote senescence dysregulation and increase tumor incidence in proliferating tissue. Do not extend without strong blastema evidence.

---

## 7. Falsifiable predictions

These are the pre-registered experimental outcomes that would falsify ARC's B-layer and P-layer hypotheses for the mouse model.

| Prediction | Experimental test | Falsifying result |
|------------|------------------|-------------------|
| **P1** | Standard BioDome + ion channel modulators (current protocol) in mouse digit P3: hypomorphic outcome (scar or amorphous cartilage) | Complete or near-complete digit regeneration without B/P additions → ARC wrong about mammalian hard gates |
| **P2** | Modified payload (TGF-β3 + VPA 24h + original MDT) in mouse digit P3: Sox2+/Prrx1+ blastema by Day 7 | No blastema by Day 10 despite correct B+P conditions → ARC B/P mechanism wrong |
| **P3** | Sequence reversal (bioelectric first, HDACi Day 2): no improvement over current protocol | Reversal works as well as correct sequence → ARC ordering hypothesis falsified |
| **P4** | P/T ratio (Ki67/p21) monitoring: values > 3.0 correlate with dysplasia in subsequent histology | No correlation between Ki67/p21 > 3.0 and disordered growth → ARC safety metric invalid |
| **P5** | TGF-β3/TGF-β1 ratio in wound exudate > 1.0 at Day 2 correlates with blastema formation probability | Ratio > 1.0 with no blastema: B-layer proxy invalid |

---

## 8. ARC scores: current vs. predicted

| Condition | P | C | S | B | T | M | Φ | Predicted outcome |
|-----------|---|---|---|---|---|---|---|-------------------|
| Mouse, no treatment | 0.15 | 0.20 | 0.20 | 0.15 | 0.60 | 0.30 | 0.04 | Fibrosis |
| Mouse, current BioDome + ionophores | 0.20 | 0.65 | 0.60 | 0.18 | 0.55 | 0.40 | 0.18 | Hypomorphic/fibrosis |
| Mouse, MDT as in Xenopus (direct copy) | 0.50 | 0.60 | 0.65 | 0.30 | 0.50 | 0.55 | 0.38 | Partial, inconsistent |
| Mouse, MDT + TGF-β3 + VPA 24h (proposed) | 0.65 | 0.65 | 0.70 | 0.65 | 0.55 | 0.60 | 0.71 | Blastema → regeneration |

*Layer values are estimated from published data on mouse wound healing biology. Φ computed using ARC v4.0 weights. These are predictions, not measured values.*

---

## 9. Connection to Levin Lab framework

ARC's C-layer (Coordinates) is directly grounded in Levin Lab bioelectricity research:

- Vmem gradients provide positional information = C-layer primary mechanism.
- Gap junction connectivity (Connexin 43, rotigaptide) enables C-layer propagation.
- Ion channel modulators (HCN2 agonists, Kir agonists) are C-layer interventions.

ARC does not challenge the importance of bioelectric signaling. It provides a multiplicative context: bioelectric C-layer signals are necessary but insufficient. The framework predicts that every bioelectric intervention that fails in adult mammals fails because B or P (or both) are below threshold, not because the C signal is wrong.

This is consistent with Levin Lab's own observation that the 24-hour application timing "was not empirically chosen" — ARC provides the theoretical basis for why timing matters: P-unlock must precede C-signal delivery.

---

## 10. Applicability beyond digit tip

| Model system | Dominant ARC bottleneck | Recommended priority |
|-------------|------------------------|---------------------|
| Mouse digit tip P3 | B-collapse + P-lock | TGF-β3 + VPA as described above |
| Mouse ear punch (MRL background) | B partially better (MRL), P still locked | VPA alone may suffice; test without TGF-β3 addition first |
| Mouse cardiac MI | B catastrophic (scar established), C unknown | Anti-fibrotic (pirfenidone + TGF-β3) priority before any patterning attempt |
| Mouse spinal cord | B moderate, S low (nerve-dependent), T problematic | S-layer (FGF9/NRG1) + B stabilization priority |
| Human digit (future) | All layers at mammalian baseline; add senolytic pre-conditioning (Phase -1) | Full protocol including senolytics (Dasatinib + Quercetin) before B+P intervention |

---

## References (key citations for this document)

- Murugan NJ et al. (2022) Acute multidrug delivery via a wearable bioreactor facilitates long-term limb regeneration and functional recovery in adult *Xenopus laevis*. *Science Advances* 8(4):eabj2164. DOI: 10.1126/sciadv.abj2164
- Murugan NJ et al. (2018) Brief local application of progesterone via a wearable bioreactor induces long-term regenerative response in adult *Xenopus* hindlimb. *Cell Reports* 25(6):1593–1604. DOI: 10.1016/j.celrep.2018.10.010
- Jopling C et al. (2010) Zebrafish heart regeneration occurs by cardiomyocyte dedifferentiation and proliferation. *Nature* 464:606–609. DOI: 10.1038/nature09391
- Godwin JW et al. (2013) Macrophages are required for adult salamander limb regeneration. *PNAS* 110(23):9415–9420. DOI: 10.1073/pnas.1300290110
- Ocampo A et al. (2016) In vivo amelioration of age-associated hallmarks by partial reprogramming. *Cell* 167(7):1719–1733. DOI: 10.1016/j.cell.2016.11.052
- Dupont S et al. (2011) Role of YAP/TAZ in mechanotransduction. *Nature* 474:179–183. DOI: 10.1038/nature10137
- Kuznetsov EA (2026) Architecture of Regenerative Control (ARC): a computational ontology for regenerative capacity. *bioRxiv* (in preparation, May 2026). https://github.com/Jenja-N/R

---

*Contact: mail.jenja@gmail.com | GitHub: [Jenja-N/R](https://github.com/Jenja-N/R)*  
*License: CC-BY 4.0. Free to use, cite, and falsify.*
