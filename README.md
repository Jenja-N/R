# Architecture of Regenerative Control (ARC) / Computational Morphogenic Runtime (CMR)

**Version 2.0** – CC‑BY‑4.0 – Last update: 2026‑05‑29

## Core Thesis

Regeneration in adult mammals is a **latent computational state transition**, not an absent biological property. It requires the simultaneous orchestration of a 6‑dimensional phase space. Failure of any single layer collapses the regenerative order parameter **Φ(t)** below the critical threshold **Φ***, driving the system into **Fibrosis** (default) or **Dysplasia** (cancer).

## The 6 Layers of ARC (5 verified + 1 hypothesised)

The framework operates as a **strict multiplicative gate**:

**Φ(t) = P^α · C^β · S^γ · B^δ · T^ε · M^ζ** ;  Regeneration ⇔ Φ(t) > Φ* (preliminary Φ* ≈ 0.35)

### Verified Layers (Diamond‑hard empirical basis)

| Layer | Biological role | Operational definition (threshold) | Key references (necessity & sufficiency) |
|-------|----------------|--------------------------------------|---------------------------------------------|
| **P – Permission** | Transient epigenetic de‑repression | ATAC‑seq OCR ratio >2.0× contralateral (Day 1) | Lu et al. Nature 2020 (OSK); Yang et al. Cell 2023; Villiard et al. 2007 (p53i); Wang et al. 2019 (HDACi) |
| **C – Coordinates** | Vmem gradient + HOX positional memory | Spatial correlation r >0.7 with target Vmem map (Day 0) | Levin 2021; Rinn 2008; Pai et al. 2015; Nogueira et al. 2022 (ATAC‑seq HoxA13) |
| **S – Synchronisation** | Phasic paracrine nerve signals (FGF9, NRG1, BMP, SHH) | Order parameter r(t) >0.6 (Day 0) from MEA or bulk RNAseq phasing | Singer 1952; Farkas 2016; Makanae 2014 (BMP/FGF replaces nerve); Xu et al. Science 2026 (DRG phasic FGF9) |
| **B – Sandboxed Growth** | Physical + chemical confinement (PLGA scaffold + TGF‑β3 perimeter) | Boundary co‑localisation <500 µm >80% (Day 0) | Badylak ECM; Ferguson 2009 (Juvista); Denis 2016 (TGF‑β inhibition blocks regeneration) |
| **T – Timed Termination** | p53 recovery + Hippo/YAP‑TAZ + iCasp9 suicide switch | log₂(p53⁺/Ki‑67⁺) >1.0 (Day 28,42) | Hartl 2024 Sci Adv; Gong 2023; Straathof 2005 Blood; TEAD inhibitor ODM‑212 (Phase 2) |

### Hypothesised Sixth Layer – Metabolic Priming (M)

**Statement:** A transient shift to anaerobic glycolysis is required for blastema survival and proliferation.

**Operational definition:** Extracellular lactate/pyruvate ratio >1.5 in wound fluid (Day 1‑7) OR 2‑DG sensitivity test.

**Falsification test (Group 6):** Full ARC protocol + 2‑deoxyglucose (2‑DG, 500 mg/kg i.p. daily, Day 0‑7).  
 – If regeneration is blocked → M is independent → sixth law confirmed.  
 – If regeneration proceeds → M is derivative of P/S → keep 5‑dimensional model.

**Supportive but not yet conclusive:** Pfefferli & Wicky 2011 (hypoxia accelerates planarian regeneration); Shyh‑Chang & Daley 2013 (Lin28/let‑7 controls glycolytic switch).

## Mathematical Formalism (6D Dynamical System)

See **`arc_dynamical_system.py`** for full ODE implementation.

**State vector:** **x(t) = (P, C, S, B, T, M) ∈ [0,1]⁶**

**Order parameter:** Φ(t) = P^α·C^β·S^γ·B^δ·T^ε·M^ζ

**Attractors:**  
– Fibrosis (x*_F) : P→0, C→0 (stability lock intact)  
– Dysplasia (x*_D) : P→1, T→0 (permission without termination)  
– Regenerative corridor (x*_R) : transient trajectory with Φ>Φ*, exiting to normal state (P→0, T→1, M→0)

## How to Use This Repository

1. **Validate the theory** with the provided falsification matrix (see `FALSIFICATION_MATRIX.md`).  
2. **Reproduce the mathematical model** by running `arc_dynamical_system.py` (Python 3.9+ with numpy, scipy).  
3. **Design your own preclinical study** using the operational definitions and thresholds.  
4. **Contribute data** – we invite independent labs to test the axolotl layer‑ablation experiments and the mouse P3 full protocol.

## License

CC‑BY‑4.0 – You are free to share and adapt, provided you give appropriate credit.

## Contact & Collaboration

Evgeny A. Kuznetsov – github.com/Jenja-N/R – Open for replication and joint validation.
