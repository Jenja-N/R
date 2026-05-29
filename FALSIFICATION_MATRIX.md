# ARC Falsification Matrix – Version 2.0 (Complete, Overdetermined)

Each row provides a **pre‑registered experimental test** and the **quantitative condition** that would falsify the corresponding ARC layer or the entire framework.

## A. Falsification of individual layers

| Layer | Hypothesis | Experimental test | Falsification condition | Empirical anchor |
|-------|------------|-------------------|--------------------------|-------------------|
| **S** | Paracrine phasing (FGF9/NRG1), not electrical impulses | Full protocol + TTX (tetrodotoxin, impulse blocker) vs. Full protocol + anti‑FGF9/NRG1 neutralising antibodies | TTX abolishes regeneration **while** anti‑FGF9 does **not** (or vice versa) → then electrical oscillation would be required, contradicting literature | Xu 2026; Singer 1952; Kumar & Brockes 2012 |
| **C** | Vmem gradient is required for patterned outgrowth | Full protocol minus electrode array (C=0) but with retinoic acid gradient alone | Regenerated digit is structurally correct (bone volume >70% by microCT) → then Vmem is dispensable | Levin 2021; Pai 2015 |
| **C** | HOX positional memory is required | Full protocol minus retinoic acid (HOX de‑repression prevented) | Digit still forms correct proximo‑distal pattern → then HOX memory not necessary | Rinn 2008; Nogueira 2022 |
| **B** | Physical sandbox is required | Full protocol but PLGA scaffold omitted; only TGF‑β3 perimeter | Regeneration proceeds without any scaffold (no deformity) → then mechanical confinement not needed | Denis 2016; Badylak 2005‑2022 |
| **P** | Transient p53 inhibition window is critical | Full protocol with p53i extended from 72h to 7 days (Group 5) | No dysplasia / tumours observed after 90 days → then p53 is not a termination brake | Hartl 2024; Villiard 2007 |
| **T** (endogenous) | p53 terminates regenerative state | Conditional Trp53 knockout in regenerating tissue (e.g., intestinal crypts after injury) | Cells do **not** lock into hyperproliferative, Wnt‑independent glycolysis → then p53 not required | Hartl 2024 Science Advances |
| **T** (engineered) | iCasp9 suicide switch eliminates >99% of vector‑bearing cells | Deliver AP1903 (10 nM, days 28 & 35) to iCasp9‑edited cells | After two pulses, >5% of vector‑bearing cells remain detectable (flow cytometry) → safety layer insufficient | Straathof 2005 Blood; Zhou 2015 |
| **M** (hypothesised) | Anaerobic glycolysis required | Full ARC + 2‑deoxyglucose (500 mg/kg i.p., days 0‑7) | Bone volume at day 42 >70% (no difference from Group 4) → M not independent | Pfefferli 2011; Shyh‑Chang 2013 |

## B. Falsification of the multiplicative gate (whole framework)

| Test | Falsification condition | Rationale |
|------|-------------------------|-----------|
| **Single‑layer omission** (Groups 2,3, and additional single‑layer knockouts) | Any omission yields **complete, patterned epimorphic regeneration** (bone volume >70%, correct digit architecture) | Multiplicative necessity disproven: that layer is not required |
| **Out‑of‑sample prediction** | Fit exponents α‑ε‑ζ on Groups 1‑3 (incomplete layers) then predict Group 4 (full protocol) outcome. If predicted Φ(t) < Φ* but actual regeneration occurs **OR** predicted Φ(t) > Φ* but no regeneration, the model is descriptive, not predictive | Shows that the multiplicative gate is not a true order parameter |
| **Additive model superiority** | A linear model (Φ_add = w₁P+…+w₆M) gives lower AIC and higher cross‑validated R² than the multiplicative model on the same data | Multiplicative interaction is unnecessary |

## C. Stopping rules for preclinical trials (pre‑registered)

| Rule | Trigger | Action |
|------|---------|--------|
| **Early success** | Group 4 meets all three primary endpoints (blastema, bone volume >70%, no dysplasia) after first replicate (n=10) | Stop, confirm in two additional replicates; publish as proof‑of‑concept |
| **Futility** | After two replicates (n=20), Group 4 bone volume <30% **or** dysplasia in ≥2 animals | Stop; redesign protocol (dose, timing, scaffold) |
| **Safety** | Any animal in Group 4 or Group 5 develops macroscopic tumour before day 42 | Stop that group; pathological evaluation; revise T‑layer |
| **Technical failure** | >20% animals in any group lose scaffold or electrode before day 14 | Exclude, redesign stabilisation |

## D. What would falsify the entire ARC framework (the “killer experiment”)

> **If a single animal in Group 4 (full protocol) regenerates a complete, anatomically correct digit after **any** single‑layer ablation** (e.g., denervation, electrode removal, extended p53i without dysplasia, or 2‑DG without effect), the multiplicative necessity claim is falsified and ARC is invalid.

No such result has ever been reported.

## E. Pre‑registration recommendation

All tests and stopping rules should be deposited on OSF (Open Science Framework) **before** data collection. Any deviation must be reported in the final manuscript.
