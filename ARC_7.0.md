# ARC 7.0 — Architecture of Regenerative Control

### A falsifiable framework for adult mammalian regeneration

**Status:** working research document — not a clinical protocol, not a validated theory
**Scope:** comparative regenerative biology, developmental biology, systems biology
**Core principle:** no claim in this document should be stronger than the evidence it rests on.

---

## 0. Abstract

ARC proposes that adult mammalian regeneration is best understood not as the firing of a single "regeneration gene," and not as raw proliferation, but as a controlled transition of damaged tissue into a coordinated morphogenetic state. Reaching that state requires the joint satisfaction of six functional tasks: **P** (cell competence), **C** (positional/coordinate information), **S** (synchronization of signals across cell populations), **B** (a permissive, anti-fibrotic boundary/microenvironment), **M** (metabolic and material support), and **T** (correct termination and stabilization).

These are functional categories, not anatomical structures or genes. A given molecule or tissue can serve more than one gate; different species can realize the same gate through different mechanisms (the clearest example: limb regeneration in salamanders is nerve-dependent, antler regeneration in deer is not — both are valid solutions to the *synchronization* problem).

ARC 7.0 separates four epistemic tiers for every claim — **established fact**, **strict inference**, **working hypothesis**, **engineering guess** — and explicitly retracts numerical artifacts from earlier drafts (fixed regression weights, an interaction-term significance claim, a 16-system validation matrix) that were never independently substantiated. The strongest forward-looking claim this document defends is not that regeneration has been engineered, but that ARC supplies an operational vocabulary precise enough to generate falsifiable predictions — and a registry of exactly which prior predictions failed.

---

## 1. Epistemic standard

| Tier | Symbol | Criterion |
|---|---|---|
| Established fact | **F** | Directly observed, ideally replicated or obtained interventionally |
| Strict inference | **I** | Logical consequence of facts, requiring no new biological assumption |
| Working hypothesis | **H** | Consistent with the facts, generates a testable prediction |
| Engineering guess | **E** | A specific proposed intervention, not yet shown sufficient or safe |

"Consistent with the data" means just that — compatibility, not causal proof. The presence of a protein across a process does not establish that adding it exogenously reproduces the process. An observed natural sequence of phases does not establish that the same sequence is the optimal order for therapeutic intervention.

Evidentiary strength is ranked, not by how appealing an explanation is, but by study design:

```
descriptive profile
< spatiotemporal association
< loss-of-function
< rescue
< factorial causal experiment
< independent cross-lab replication
< prospective out-of-sample prediction
```

---

## 2. Logical foundation

1. **F** — A multicellular organism builds a stable three-dimensional form from a single cell via morphogenesis.
2. **F** — Morphogenesis depends on the interaction of genetic networks, signaling, mechanics, metabolism, cell–cell communication, and tissue spatial history.
3. **F** — Regeneration restores part of a lost structure and its function.
4. **I** — Regeneration reuses at least part of the machinery of morphogenesis, patterning, differentiation, and tissue homeostasis.
5. **F** — Across vertebrates, many developmental genes and regulatory modules are deeply conserved; at the same time, lineage-specific components, regulatory rewiring, and genuine evolutionary losses exist.
6. **I** — The absence of full regeneration in an adult mammal cannot automatically be explained by the absence of all necessary genetic components.
7. **F** — Extensive regeneration of complex structures has independently persisted, arisen, been constrained, or been lost across different phylogenetic lineages.
8. **I** — Complex regeneration is biologically realizable as a class of process; it does not follow that it is automatically realizable for every organ in every species, including humans.
9. **F** — Adult mammalian tissues frequently respond to major injury with repair, fibrosis, or incomplete restoration rather than regeneration.
10. **F** — The outcome along the repair–fibrosis–regeneration axis depends on context: immune dynamics, extracellular matrix, mechanical tension, cell competence, and the scale/geometry of the injury.
11. **F** — Regeneration and tumor growth share overlapping executional machinery: plasticity, proliferation, migration, matrix remodeling, angiogenesis.
12. **I/H** — The core engineering problem is not only triggering growth, but recovering a controlled trajectory:

```
competence → coordinates → synchronization → construction
→ differentiation → termination → homeostasis
```

### 2.1 What does not follow

From the above, it does **not** follow that:

- the adult human genome contains a literal, fully accessible blueprint for any lost organ;
- the mechanisms required for human regeneration are identical to those in axolotl;
- removing fibrosis automatically triggers organogenesis;
- activating Wnt, BMP, FGF, YAP, or MYC alone will produce correct form;
- regeneration and cancer differ only by the presence of a single stop signal;
- evolution layered six ARC-relevant constraints in one linear, universal order;
- reversing a presumed evolutionary constraint is a validated therapeutic algorithm;
- regeneration of one human organ implies readiness of all adult organs for full regeneration.

---

## 3. Definition

**ARC — Architecture of Regenerative Control**: a functional meta-model describing the minimal classes of tasks whose joint failure makes complete, safe anatomical restoration unlikely.

```
ARC = { P, C, S, B, M, T }
```

ARC does not equate a gate with a single gene or protein. A gate is a *function*; proteins, cells, tissues, electrical states, and physical fields are possible *realizations* of that function. One mechanism can serve several gates at once (vasculature, for instance, contributes to M, C, and S simultaneously) — this functional overlap is a real property of biology, not evidence that the categories are statistically independent variables.

---

## 4. The six gates — operational definitions

### 4.1 P — Permission / Competence
**Function:** the ability of relevant cells to enter the plastic, proliferative, differentiation-competent state required.
Includes: availability of suitable cells or local progenitors; accessible regulatory chromatin; reversible cell-state switching; intact organ-specific regulatory networks; absence of an irreversible senescent or genomically unstable block.
P is *not* "open chromatin in general" — global epigenetic de-repression is not organ-specific and risks identity loss, dysplasia, or transformation.
**Falsifier of universal necessity:** full regeneration documented with experimentally confirmed absence of competent cells and no possible state change. In practice this is hard to establish, so P remains a strong functional hypothesis, not a proven universal law.

### 4.2 C — Coordinates / Positional information
**Function:** specifying *what*, *where*, in what orientation, topology, and to what size something should be built.
Includes: stromal positional memory; regional transcriptional codes; epithelial polarity; the geometry of the residual niche; signaling gradients; mechanical fields; bioelectric states; cell lineage of origin; ECM spatial organization.
C is not reducible to a HOX code or to membrane voltage alone. Bioelectric mechanisms participate in patterning in several model systems, but no universal electrical "blueprint" of a human organ has been decoded.
**Key distinction:** increased tissue mass ≠ restored C. C is supported by recovery of correct axis, boundaries, tissue interfaces, topology, regional identity, and functional geometry — not by volume alone.

### 4.3 S — Synchronization
**Function:** coordinating cell populations, signals, and phases in time and space.
S is best modeled as a vector, not a single channel:

```
S = { S_neural, S_endocrine, S_vascular, S_paracrine, S_immune, S_electrical, ... }
```

Different systems use different combinations. "A nerve is required for any regeneration" is false as a universal law — antler regrowth proceeds without dependence on its sensory innervation, while salamander limb regeneration is classically nerve-dependent (see §8.3, §8.5). This supports modularity of S, not free interchangeability of every channel with every other.

### 4.4 B — Boundary / Permissive microenvironment
**Function:** keeping the local environment within a range compatible with morphogenesis rather than rapid scarring, infection, or mechanical collapse.
Includes: mechanical offloading; geometric constraint; moisture; matrix stiffness/composition; controlled proteolysis; immune–matrix interaction; protection from premature epithelialization; protection from contamination.
B is *not* total suppression of inflammation, TGF-β, or proteases — human fingertip regeneration proceeds with substantial inflammatory and proteolytic activity (§7). What matters is regime, duration, and balance, not minimization.

### 4.5 M — Metabolic and material support
**Function:** supplying oxygen, substrates, minerals, vascular access, redox control, biosynthetic and energetic capacity.
M is not reducible to a single Warburg-type switch; different cells and stages use different metabolic regimes. Vascular contribution is a clear case of cross-gate overlap with C and S.

### 4.6 T — Termination and stabilization
**Function:** stopping growth, stabilizing identity, restoring barrier function, transitioning to homeostasis.
Includes: cell-cycle control, mechanosensing/contact inhibition, Hippo/YAP–TAZ, p53/Rb, differentiation, apoptosis, immune resolution, vascular maturation, ECM stabilization.
**Important correction relative to earlier drafts of this project:** "T doesn't need to be touched — the organism will stop growth on its own" is a hypothesis, not an established fact. Intact tumor suppressors do not guarantee organ-specific termination will fire automatically when growth is artificially extended.

---

## 5. Gate interaction logic

### 5.1 Why an AND-intuition persists

None of the following is sufficient on its own: competent cells without coordinates; coordinates without material; signals delivered into a fibrotic or destructive environment; growth without termination; a permissive environment without an organ-specific program. ARC therefore retains, as a **qualitative hypothesis**, that complete regeneration requires the joint satisfaction of P, C, S, B, M, and T. This is *not yet* a proven Boolean conjunction — the real variables are continuous, phase-dependent, partially overlapping, and may sometimes compensate for one another.

### 5.2 The S×B hypothesis — status, and what this document does not claim

A recurring idea across drafts of this project is that the *interaction* between signaling (S) and a permissive boundary (B) matters more than either alone: the same growth signal can drive organized morphogenesis in a permissive matrix and drive myofibroblast expansion / scarring in a fibrotic one. This is biologically motivated by several independent observations assembled in this document:

- adult skin (fibrotic-biased) vs. oral mucosa (permissive-biased) heal the same kind of wound very differently;
- the repeated clinical failure of single-growth-factor therapies (e.g., topical VEGF or BMP) in chronic, fibrotic wounds;
- wound-contraction biology, where blocking contraction with a dermal regeneration template can shift outcome toward regeneration-like repair (§8.6);
- Acomys vs. Mus, where the same early ERK activation event diverges in duration and outcome depending on tissue context (§8.1).

**This document explicitly does not claim** that this interaction has been statistically validated. An earlier internal draft asserted a specific logistic-regression result — fixed weights (B = 1.78, P = 1.45, S = 1.42, T = 1.20, M = 0.93, C = 0.63, an S×B interaction term of +0.85, bias = −4.52), fit on "16 systems across 6 species," with a reported ΔAIC of −10.6 favoring the interaction term, and an out-of-sample "pass" on zebrafish heart regeneration. **None of this could be located in the project's underlying data or scripts when checked**, and the actual exploratory scripts in the project show contradictory, unresolved attempts on a *single* dataset (a mouse digit-crush time series), including formula sign errors, memory failures, and zero-gene-matched runs. The numeric claim is retracted in full (see §13, items 1 and 17–19). S×B remains a **working hypothesis (H)**, motivated by real comparative biology, not a measured result.

### 5.3 Honest formal statement

Instead of fixed weights:

```
R(t) = F_tissue [ P(t), C(t), S(t), B(t), M(t), T(t), interactions(t) ]
```

where `R(t)` is a multidimensional outcome vector (volume, geometry, tissue composition, vascularity, innervation, mechanical function, physiological function, stability, safety) and `F_tissue` is an unknown, tissue- and species-specific function. A threshold or logistic approximation is appropriate only *after* each gate is measured by independent operational markers and the model is tested out of sample — not before.

---

## 6. Human anchor: fingertip regeneration

**Source (F):** Schultz et al., *npj Regenerative Medicine*, 2025 (NCT03089060, PMID 41193504). Longitudinal wound-fluid mass-spectrometry proteomics, 22 patients aged 2–72, 974 proteins detected.

### 6.1 What is established (F)

- Four clinically and proteomically distinct phases: **coagulation → hypergranulation → proliferation → epithelialization**, supported by unbiased clustering on phase, not on age/sex/injury severity.
- Duration and outcome are **not** statistically predicted by age (2–72 years), sex, bone involvement, or injury severity in this cohort (multiple regression).
- Epithelialization is the **last** phase, not the first — the wound stays open and moist for weeks, the reverse order from axolotl/mouse models.
- A classical blastema or apical epidermal cap (AEC) is **not identified or confirmed** in this human dataset; the authors call this an open question, not a settled negative.
- Predicted upstream regulators (by pathway analysis, not direct measurement): two waves of immune activation (TNF, IL6, IL1B, IFNG), TGF-β1 across phase transitions, the calcium channel RYR1, MYC/TP53/OSM, a metabolic shift.
- A protease/anti-protease paradox: MMP8, MMP9, and neutrophil elastase rise alongside SERPINA1, SERPINB, and A2M — a profile resembling a chronic non-healing wound that nonetheless resolves without scarring. The authors explicitly flag this mechanism as unresolved.
- The authors state plainly that their method (spectral counting) is not designed for absolute quantitative thresholds, and that the role of physical/electrical/Wnt/FGF signals is *future work*, not a result extracted from this dataset.

### 6.2 What is not established — explicit boundary

This study does **not** show: a fixed day-by-day calendar for any phase; a decoded "molecular code of form"; that age has no effect on human regeneration in general; that the four fingertip phases transfer directly onto a tooth socket; that the natural temporal order is the optimal *therapeutic* order.

### 6.3 LTF as a cross-species candidate (H)

Lactotransferrin (LTF) is abundant across all four human phases. Separately, Simkin et al. (*Developmental Cell*, 2024, PMID 38228141) show that LTF and VEGFC are specifically secreted by tissue-resident macrophages in regenerating (not scarring) Acomys ear tissue. The co-occurrence in two independent systems is a reasonable basis for prioritizing LTF as a candidate marker of a permissive niche. It is **not** established as a sufficient or causal driver of complex regeneration in either species, and it is explicitly *not* treated here as a seventh gate or a proven morphogen.

---

## 7. Comparative models

### 7.1 Acomys (African spiny mouse)

- (**F**) Acomys regenerates hair follicles, dermis, and cartilage in ear-punch injuries without scarring (Seifert et al., *Nature*, 2012, DOI 10.1038/nature11499).
- (**F**) Macrophage depletion at the time of injury abolishes the regenerative response (Simkin, Gawriluk, Gensel & Seifert, *eLife*, 2017) — macrophages are functionally necessary, not merely correlated.
- (**F**) Resident macrophages in regenerating Acomys ear tissue specifically secrete LTF and VEGFC; this is absent in non-regenerating Mus wounds (Simkin et al., *Developmental Cell*, 2024, PMID 38228141).
- (**F**) Fibroblasts from highly regenerative mammals (including Acomys) are more resistant to ROS-induced senescence than fibroblasts from non-regenerative species (Saxena, Vekaria, Sullivan & Seifert, *Nature Communications*, 2019, DOI 10.1038/s41467-019-12398-w).
- (**F**) Sustained ERK activation is associated with regeneration in Acomys, transient activation with scarring in Mus; this has been shown bidirectionally — inhibiting ERK in Acomys shifts the outcome toward scar, activating it in Mus partially unlocks regenerative features (*Science Advances*, 2023).
- **ARC interpretation (H):** these findings support the functional importance of B (permissive microenvironment), P/T (resistance to a senescent block), S_immune, and M (redox/mitochondrial resilience) — but they do not by themselves demonstrate the completeness of the six-gate architecture.

### 7.2 Deer antler

Antlers are the only mammalian organ that fully regrows annually after being shed. Denervation studies show that nerve supply affects size, quality, and rate, but is **not** a universally required initiator of antler growth — vascular, endocrine, and local paracrine mechanisms carry strong independent support (synthesized in Dolan et al., reviewed in *Cell Regeneration*, 2023, PMC10368610; original stem-cell characterization in Qin et al., *Science*, 2023, DOI 10.1126/science.add0488, and independently in Wang et al., *Science*, 2023, DOI 10.1126/science.adg9968).

- **Strict inference (I):** S_neural is not a necessary condition for antler regeneration.
- **Working hypothesis (H):** the S-function in antler is realized through a redundant, non-neural architecture.
- **Explicitly retracted claim from earlier drafts:** "all six ARC gates are fully causally proven in deer." Several gates have spatial, transcriptomic, or associative support, but their ARC-category assignment was made *post hoc*; this is compatibility with the model, not independent validation of it.

### 7.3 Axolotl / newt — nerve dependence and its molecular basis

Nerve dependence of salamander limb regeneration is classically established and has a specific molecular mechanism: nAG, a secreted ligand of the Prod1 receptor, is sequentially expressed in the regenerating nerve and wound epidermis after amputation; its expression at both sites is abolished by denervation, and local nAG delivery by electroporation is sufficient to rescue a denervated blastema (Kumar, Godwin, Gates, Garza-Garcia & Brockes, *Science*, 2007, DOI 10.1126/science.1147710).

- **ARC interpretation (H):** S can be neural-dominant in this system.
- **Not established:** that a cocktail of VEGF/IGF-1/FGF2/ANGPT1 can substitute for the nerve in this or any mammalian system. A negative result in such an experiment would not by itself falsify ARC as a whole — it could equally reflect a wrong factor combination, dose, timing, carrier, or a confound of S with C or B at the denervation site (see §11 for how to design a discriminating test).

### 7.4 Mouse digit tip

Two papers, published in close succession, established that the mouse digit-tip blastema is composed of lineage-restricted progenitor cells, not a multipotent pluripotent pool:

- Lehoczky, Robert & Tabin, *PNAS*, 2011, "Mouse digit tip regeneration is mediated by fate-restricted progenitor cells," DOI 10.1073/pnas.1118017108.
- Rinkevich et al., *Nature*, 2011, "Germ-layer and lineage-restricted stem/progenitors regenerate the mouse digit tip," *Nature* 476:409–413.

Regeneration in this model strictly requires the overlying nail organ (Takeo et al., *Nature*, 2013 — see §7.5) and is bone-growth-dependent, not purely soft-tissue-dependent.

### 7.5 Nail organ and Wnt — a direct C-gate support

Takeo et al. (*Nature*, 2013, *Nature* 499:228–232) showed that Wnt activation in the nail epithelium couples nail growth to digit-tip regeneration, and is required both for nail regeneration itself and for attracting the nerve fibers that support mesenchymal blastema growth. This is one of the more direct mechanistic links available between a positional/coordinate hub (the nail organ) and a synchronization channel (nerve attraction) in a mammalian system.

### 7.6 Wound contraction

Yannas and colleagues (synthesized in Yannas et al., *npj Regenerative Medicine*, 2021) document that in several mammalian injury models, scar-free outcomes correlate with near-absence of myofibroblast (α-SMA+) activity, and that physically blocking wound contraction — for example with a dermal regeneration template/scaffold — can shift the outcome toward a regeneration-like trajectory rather than dense scar. This is a direct empirical anchor for the B gate: signal presence does not guarantee regeneration if the mechanical/matrix environment is actively contracting and scarring.

### 7.7 A general mechanism for why B matters: matrix elasticity

Engler et al. (*Cell*, 2006, DOI 10.1016/j.cell.2006.06.044) showed that substrate stiffness alone — independent of soluble signaling — directs mesenchymal stem cell lineage commitment. This is not a regeneration study, but it is a clean, frequently-replicated mechanistic basis for why ARC treats B as a gate with causal weight of its own, not merely a passive backdrop for signaling.

---

## 8. Evolutionary hypothesis: a network, not a ladder

It cannot be defended as established science that evolution laid down a single linear sequence of constraints (e.g. "Wnt → BMP → FGF → SHH → Notch → TGF-β/immunity") that must be removed in the exact reverse order to restore regeneration. Evolution is branching, tissue-specific, repeatedly reuses the same ancient circuits for different purposes, and does not design top-down. Most ARC-relevant mechanisms predate vertebrates and have been redeployed many times.

A defensible, weaker evolutionary statement is that selection has tuned the **weights, connections, time windows, and redundancy** of these gates differently across lineages and tissues — trading regenerative capacity against other fitness costs:

| Gain | Plausible cost |
|---|---|
| Fast wound closure | Less time available for morphogenesis |
| Strong anti-infective response | Damage to the permissive niche |
| Tumor-suppressive safety margin | Reduced plasticity |
| Stiff ECM | Reduced migration / re-patterning capacity |
| Energy economy | Limits on sustained organogenesis |
| Early reproductive investment | Weaker late-life somatic regeneration |
| Rapid mechanical stabilization | Scar instead of restored geometry |

These trade-offs arose independently in different lineages, not in one order — so any evolutionary map for ARC should be drawn as a network, not a ladder.

**Working safety principle, not a reconstruction of evolutionary history:** stabilize the external environment (B) and resource supply (M) before amplifying signaling (S); confirm coordinates (C) before expanding competence (P); termination (T) must be actively managed before, during, and after growth, not "switched off first." Globally removing T before C and B are restored is the path to dysplasia or neoplasia, not regeneration.

---

## 9. ARC applied to teeth

### 9.1 What is known

The adult jaw retains periodontal ligament cells, bone progenitor populations, perivascular cells, fibroblasts, immune cells, residual epithelial populations (epithelial rests of Malassez), and intact vascular/neural structures. Dental pulp contains cells with reparative potential, and several developmental tooth-signaling pathways are partially retained in adult tissue.

USAG-1 (SOSTDC1) tonically suppresses both BMP and Wnt signaling around supernumerary tooth primordia.
- (**F**, animal) A single anti-USAG-1 antibody dose triggers growth of additional functional teeth in mice, ferrets, and dogs, replicated across labs.
- (**F**, human, safety only) TRG035 (Toregem BioPharma) is in a Phase I human safety trial (~30 adults, age 30–64, Kyoto University Hospital, since autumn 2024); no published human efficacy data exist as of mid-2026.

### 9.2 The critical boundary

Growing a supernumerary tooth from a *retained developmental primordium* and *de novo* regenerating a tooth after ordinary adult extraction are different problems. The first may require lifting a single dominant brake. The second requires creating or restoring an epithelial–mesenchymal organogenetic center, correct crown–root geometry, pulp, periodontium, an apical foramen, vasculature, innervation, cementum, dentin, and either an enamel-forming epithelium or an engineered substitute for it. **USAG-1 data cannot be used as evidence that an empty adult socket is ready to regrow a tooth from a single injection.**

### 9.3 Tooth-specific gate realizations

| Gate | Possible realization | Status |
|---|---|---|
| P | DPSCs, PDLSCs, epithelial rests of Malassez, dental-network competence | partially known; organ-level sufficiency unproven |
| C | socket geometry, PDL, epithelial–mesenchymal interfaces, regional memory | critical unknown |
| S | Wnt/BMP/FGF/SHH/Notch + immune/vascular/neural/paracrine channels | components known; order and sufficiency unestablished |
| B | infection control, epithelialization control, ECM, mechanical offloading | engineerable; optimal regime unknown |
| M | early diffusion, then vascularization, minerals, biosynthesis | necessary; parameters unknown |
| T | size/shape control, proliferation arrest, maturation | mandatory; cannot be assumed automatic |

### 9.4 Internal pilot data (GSE185222) — hedged

A within-project analysis of a public single-cell dataset of human dental pulp (sound vs. carious third molars) found, in the adult deep-caries subgroup compared to a fetal odontoblast-like reference:

- DSPP-positive cells fall from ~34% (fetal odontoblast-like) to **0%** (adult deep caries);
- TGFB2-positive fraction falls from ~16% to ~2.5%;
- JUNB (not JUN/FOS) rises sharply and tracks pathology — an explicit, documented correction of an earlier in-project hypothesis that had implicated FOS/JUN;
- a BMP4 paradox: BMP4-positive fraction is roughly nine-fold *higher* in the adult group, while DSPP output is zero — elevated BMP4 alone does not rescue odontoblast output;
- WNT10A is enriched in DSPP-positive cells but functions as a ligand, not as evidence of canonical pathway activation: the AXIN2+/LEF1+/DSPP+ triple-positive population is **zero cells** across all four groups (n = 10,437 cells total).

**Honest limitation statement:** this is one exploratory, non-replicated, pilot single-cell dataset with no independent cell-type annotation; specific subgroups are small (some n in the single digits). These results generate testable hypotheses about lost odontoblast competence and network-level (not single-ligand) failure; they are not, and should not be reported as, "mathematical certainty" about human biology.

### 9.5 Minimal realistic first product

Rather than targeting a complete tooth with enamel, a more defensible near-term engineering target is: (1) a viable, vascularized dentin–pulp complex; (2) a root-like structure with a periodontal interface; (3) functional fixation in the alveolus; (4) controlled mineralization; (5) an artificial crown covering, given that enamel-forming ameloblasts are terminally differentiated and lost after eruption even in species that do regenerate teeth (sharks, crocodilians).

---

## 10. Experimental falsification program

ARC is not falsified by the failure of any single intervention cocktail. It *is* falsified, as an architecture, by results such as:

1. **Full regeneration at B = 0** — organized regeneration documented in a stable, contractile, collagen-I-rich, α-SMA-high fibrotic environment.
2. **Full anatomically correct regeneration at C = 0** — without any positional/coordinate hub.
3. **S×B fails to improve prediction** in an independent, preregistered dataset compared to a main-effects-only model.
4. **`S only` reliably produces full, safe, organized form** without a permissive boundary or coordinate input — this would mean ARC over-weights the necessity of functional decomposition.
5. **`B + C + S` does not outperform simpler interventions** across several independent injury models — this would sharply reduce ARC's practical value, even if it remains descriptively true.

### 10.1 A registered test design

Before applying ARC to a new system:
1. Fix operational definitions for P/C/S/B/M/T, specific markers, measurement methods, an aggregation rule, the expected outcome, and the error criterion — *in advance*.
2. Apply ARC to systems **not used to develop it**.
3. Compare against simpler alternative models (fibrotic-environment-only; competence-only; signaling-only; an unstructured multivariate baseline; ARC).
4. ARC earns independent value only if it improves out-of-sample prediction, calibration, interpretability, or transfer between systems — not by fitting better post hoc.

### 10.2 Specific factorial tests

- **Necessity of each gate:** in a reproducible model, selectively null one gate while holding the others maximally intact; primary endpoint must include shape, tissue composition, function, stability, and safety — not length or mass alone.
- **Modularity of S:** factorial design crossing S_neural × S_vascular × S_paracrine × S_endocrine, measuring main effects and interactions; possible outcomes range from a single bottleneck channel to full redundancy to phase-dependent dynamics — all are informative, none assumed in advance.
- **B×S interaction:** apply the same signaling stimulus in a permissive vs. fibrotic vs. mechanically offloaded vs. mechanically loaded environment; if the effect of S depends strongly on environment, this is direct support for the interaction (this is the correct way to test the S×B hypothesis — not by reusing an unvalidated regression).
- **Tooth C-gate:** compare an intact periodontal/epithelial-rest niche, a selectively ablated niche, a reconstructed epithelial–mesenchymal niche, a geometric scaffold without positional memory, a cellular niche without scaffold, and the full combination; primary endpoint is organized dentin–pulp–PDL architecture with a reproducible axis and functional interface, not mineralization volume.

---

## 11. Safety

ARC is not grounds for self-administered or unsupervised clinical use of HDAC inhibitors, senolytics, mTOR activators, Wnt/BMP agonists, TGF-β modulators, growth factors, electrical/bioelectric stimulation devices, or anti-inflammatory cocktails, in any context — dental or otherwise.

Principal risks that any ARC-derived protocol must treat as primary, multidimensional safety endpoints (not absence-of-visible-tumor in a short observation window): neoplasia, dysplasia, odontoma/ectopic mineralization, chronic infection, alveolar bone resorption, uncontrolled proteolysis, ischemia, mis-innervation and chronic pain, ankylosis, non-functional tissue mass, incorrect form, and an unstable regenerate.

---

## 12. Retraction ledger

This is a record of specific claims that appeared in earlier drafts of this project and did not survive independent verification. They are listed here, with corrections, precisely so they are not silently reintroduced in a future version.

| # | Earlier claim | Correction |
|---|---|---|
| 1 | Fixed equation Φ = σ(1.78B + 1.45P + 1.42S + 1.20T + 0.93M + 0.63C + 0.85·S×B − 4.52), with "100% LOOCV on 16 systems" | Weights are single-author expert judgment with no independent inter-rater validation; N=16 is far too small for this claim; not a biological constant. Retracted. |
| 2 | "ARC has already been proven through its explanatory power over paradoxes" | Explaining known phenomena post hoc is not proof. The same phenomena are independently explained by named mechanisms (ERK switching, M1/M2 macrophage balance, positional memory) published before and independently of ARC. |
| 3 | "Age does not affect human regeneration; the old view that regenerative competence declines is a lie; the adult genome is fully competent" | Non-significance of age as a predictor in a single n=22 regression is weak evidence of *no large effect in that cohort*, not proof of equivalence between a 2-year-old and a 72-year-old, and not a general claim about all human tissues. |
| 4 | A fixed day-by-day intervention calendar for human fingertip regeneration (Day 1–3 / 2–7 / 7–21 / 21–42+) | The source study explicitly states duration is not predicted by the variables tested and that its method is not designed for absolute quantitative thresholds. No such calendar exists in the data. |
| 5 | "The molecular code of form (Wnt/FGF/bioelectric) has already been decoded from wound-fluid proteomics" | No Wnt/FGF/SHH proteins were detected in the wound fluid; the source authors explicitly list physical/electrical signaling as future work. |
| 6 | "Blastema-like mesenchymal condensation is confirmed" in human fingertip regeneration; invented phase names ("Proliferative/Blastema-like," "Patterning/Differentiation") | The source study explicitly states a blastema/AEC was not identified in humans and treats this as an open question. The real phases are coagulation, hypergranulation, proliferation, epithelialization. |
| 7 | Macrophage-depletion citation given as "Seifert et al., 2021" | Correct citation is Simkin, Gawriluk, Gensel & Seifert, *eLife*, 2017. |
| 8 | "TGF-β1/β3 balance (Simkin et al., 2024)" | Simkin 2024 is about LTF/VEGFC secretion by macrophages, not about TGF-β isoform balance. The TGF-β1/β3 balance literature is older (fetal scarless-healing work, early 1990s) and should be cited separately. |
| 9 | A proposed universal evolutionary patch order ("Wnt → BMP → FGF → SHH → Notch → TGF-β/immunity") and a matching reverse "shutdown order" for teeth | Unsupported extrapolation; the human fingertip dataset says nothing about teeth or about the evolutionary order in which signaling pathways arose. |
| 10 | A bioreactor intervention protocol (senolytics, exogenous LTF, mTOR activators, anti-α-SMA agents) presented as derived "from the data" | The source study is purely observational; none of these interventions were tested in it. A list of proteins present during a natural process is not an intervention protocol. |
| 11 | References to non-existent local files/paths (e.g. `sandbox:///mnt/agents/output/...`) as already-created downloadable deliverables | No such paths exist in the working environment; no such files were created. |
| 12 | "DSPP = 0... Confidence: ABSOLUTE / mathematical certainty" for the internal GSE185222 analysis | A correct measurement in one small, non-replicated pilot dataset whose own limitations (small subgroup n, no independent cell-type annotation) were stated elsewhere in the same document. Incompatible with an "absolute certainty" label. |
| 13 | "All six ARC gates are fully causally proven in deer antler" | Several gates have spatial/transcriptomic/associative support, but category assignment to ARC was made post hoc. Compatible with the model; not independent validation of it. |
| 14 | A wearable-bioreactor citation given as "Blackiston et al., 2021, *Nature Communications*, DOI 10.1038/s41467-021-25337-8" | This paper/DOI does not exist. The real wearable-bioreactor work on adult *Xenopus* hindlimb regeneration is Herrera-Rincon et al., *Cell Reports*, 2018 (progesterone) and Murugan et al., *Science Advances*, 2022 (multidrug cocktail). |
| 15 | A deer-antler citation given as "Li et al., 2023, *Science*, 'A stem cell-based mechanism for deer antler regeneration,' DOI 10.1126/science.adg6332" | This title/author/DOI combination does not exist. Two real, independent 2023 *Science* papers cover this topic: Qin et al. (DOI add0488) and Wang et al. (DOI adg9968). |
| 16 | A mouse-digit-tip citation given as "Lehoczky et al., 2011, *Nature*, 'Mouse digit tip regeneration is mediated by a blastema'" | Wrong title and wrong journal. The real paper is Lehoczky, Robert & Tabin, *PNAS*, 2011, "Mouse digit tip regeneration is mediated by **fate-restricted progenitor cells**" (DOI 10.1073/pnas.1118017108). The DOI given for the false citation (10.1038/nature10579) belongs to an unrelated *Drosophila* cell-cycle paper. |
| 17 | A DOI given for Seifert et al., 2012 (Acomys) as 10.1038/nature11214 | Wrong DOI. Correct DOI is 10.1038/nature11499. |
| 18 | "S×B = +0.85 with ΔAIC = −10.6 on a curated 16-system, 6-species matrix" attributed as already established "in the uploaded ARC document" | Neither the underlying analysis nor the AIC comparison could be located in the project's data, scripts, or prior documents. The cited weights are the same unvalidated numbers retracted in item 1; this re-presentation added a fabricated statistical justification on top of an already-retracted number. |
| 19 | An out-of-sample "PASS" on zebrafish heart regeneration (Φ = 0.871) used to defend the S×B weight | No corresponding computation exists in the project. Retracted along with item 18. |
| 20 | Specific illustrative numbers for "adult skin vs. oral mucosa" (P/C/S/B/T/M values to two decimal places, Φ = 0.128 vs. 0.772) presented as a worked example from real data | No such per-tissue gate scoring exists in the project's actual scripts or outputs; this reads as an invented illustration, not a measurement. |

---

## 13. Minimal defensible core

> Complete restoration of a complex structure requires more than cell proliferation: it requires the coordinated satisfaction of several functional tasks — competent cells, spatial information, synchronization across cell populations, a permissive microenvironment, metabolic/material support, and correct termination. These tasks are realized by different molecular and cellular mechanisms in different regenerating systems. ARC offers an operational decomposition of these tasks and a research program capable of testing their necessity, their interactions, and their tissue-specific realization.

This statement is compatible with developmental biology and comparative regeneration research, does not depend on any single molecular mechanism, allows for genuine cross-species differences, makes no therapeutic promise, generates falsifiable predictions, and can be disproven by the experiments in §10.

---

## 14. Criteria for ARC to graduate from heuristic to validated theory

ARC moves from "functional meta-model" to "supported theory" only when, jointly:

1. Stable operational definitions exist for P/C/S/B/M/T, fixed *before* analysis of new data.
2. Multiple independent raters can score them with high inter-rater agreement.
3. ARC generates a prediction for a new system *in advance* of seeing the outcome.
4. That prediction is confirmed out of sample.
5. ARC outperforms simpler alternative models on the same data.
6. At least one intervention experiment confirms a structural prediction unique to ARC (not just "more growth," but a specific change in geometry/composition/function predicted by the gate framework).
7. The result is independently replicated.
8. Negative results and the limits of applicability are published, not omitted.
9. A concrete experiment exists that could falsify the architecture itself, not just one intervention cocktail.
10. Terminology and gate definitions are not changed after the fact to rescue the model from a disconfirming result.

Until these are met, the correct status of ARC is:

> **A falsifiable functional meta-model of regenerative control, with strong comparative-biology motivation, but without completed prospective validation.**

---

## 15. References

Only sources independently verified during the development of this document are listed. Where a claim above is hedged as hypothesis, no citation implies it is established fact.

1. Zoonomia Consortium. Comparative genomics across 240 placental mammals. *Science*, 2023.
2. Schultz J., Patel P.A., Aires R. et al. Human fingertip regeneration follows clinical phases with distinct proteomic signatures. *npj Regenerative Medicine*, 2025. PMID: 41193504. NCT03089060.
3. Simkin J., Gawriluk T.R., Gensel J.C., Seifert A.W. Macrophages are necessary for epimorphic regeneration in African spiny mice. *eLife*, 2017.
4. Simkin J., Aloysius A., Adam M. et al. Tissue-resident macrophages specifically express lactoferrin and VEGFC during ear regeneration in *Acomys*. *Developmental Cell*, 2024. PMID: 38228141.
5. Saxena S., Vekaria H., Sullivan P.G., Seifert A.W. Connective tissue fibroblasts from highly regenerative mammals are refractory to ROS-induced cellular senescence. *Nature Communications*, 2019. DOI: 10.1038/s41467-019-12398-w.
6. Sustained ERK activity as a switch between regeneration and scarring in Acomys vs. Mus. *Science Advances*, 2023.
7. Seifert A.W., Kiama S.G., Seifert M.G., Goheen J.R., Palmer T.M., Maden M. Skin shedding and tissue regeneration in African spiny mice (*Acomys*). *Nature*, 2012. DOI: 10.1038/nature11499.
8. Dolan C.P. et al. (synthesis). Deer antler renewal gives insights into mammalian epimorphic regeneration. *Cell Regeneration*, 2023. PMC10368610.
9. Qin T. et al. A population of stem cells with strong regenerative potential discovered in deer antlers. *Science*, 2023. DOI: 10.1126/science.add0488.
10. Wang D. et al. Stem cells drive antler regeneration. *Science*, 2023. DOI: 10.1126/science.adg9968.
11. Kumar A., Godwin J.W., Gates P.B., Garza-Garcia A.A., Brockes J.P. Molecular basis for the nerve dependence of limb regeneration in an adult vertebrate. *Science*, 2007. DOI: 10.1126/science.1147710.
12. Lehoczky J.A., Robert B., Tabin C.J. Mouse digit tip regeneration is mediated by fate-restricted progenitor cells. *PNAS*, 2011. DOI: 10.1073/pnas.1118017108.
13. Rinkevich Y., Lindau P., Ueno H., Longaker M.T., Weissman I.L. Germ-layer and lineage-restricted stem/progenitors regenerate the mouse digit tip. *Nature*, 2011. 476:409–413.
14. Takeo M. et al. Wnt activation in nail epithelium couples nail growth to digit regeneration. *Nature*, 2013. 499:228–232.
15. Yannas I.V., Tzeranis D.S. et al. Wound contraction and the mechanism by which mammals fail to spontaneously regenerate organs. *npj Regenerative Medicine*, 2021.
16. Engler A.J., Sen S., Sweeney H.L., Discher D.E. Matrix elasticity directs stem cell lineage specification. *Cell*, 2006. DOI: 10.1016/j.cell.2006.06.044.
17. Herrera-Rincon C. et al. Brief local application of progesterone via a wearable bioreactor induces long-term regenerative response in adult *Xenopus* hindlimb. *Cell Reports*, 2018.
18. Murugan N.J. et al. Acute multidrug delivery via a wearable bioreactor facilitates long-term limb regeneration and functional recovery in adult *Xenopus laevis*. *Science Advances*, 2022.
19. USAG-1 (SOSTDC1) antibody and supernumerary tooth formation — animal studies and ongoing Phase I human safety trial (TRG035, Toregem BioPharma; Kyoto University Hospital, from autumn 2024). Human efficacy data not yet published as of mid-2026.
20. Internal pilot single-cell analysis of GSE185222 (human dental pulp, sound vs. carious third molars). Project-internal; not independently replicated; hypothesis-generating only.

---

## 16. Provenance note

This document is the product of an iterative, adversarial verification process: successive drafts were checked against primary sources, fabricated statistics and misattributed citations were identified and removed, and genuine improvements were retained. §12 exists specifically so that future revisions of this project do not silently reintroduce claims that have already been checked and found unsupported. Readers extending this document are asked to preserve that discipline: every new numeric claim should specify whether it is a fact, a strict inference, a hypothesis, or an engineering guess, and every new citation should be checked against the actual source before being added to §15.

**License / use:** this is a research notebook, not medical guidance. Nothing in this document authorizes self-experimentation or unsupervised clinical application of any compound, device, or protocol mentioned.
