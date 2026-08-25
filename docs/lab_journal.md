# Quantum Bio Recovery Engine — Lab Journal

**Researcher:** Rishabh Kumar

**Project:** Quantum Bio Recovery Engine

**Focus:** Hybrid biological–probabilistic–combinatorial–quantum optimization for therapeutic combination discovery

**Start Date:** 2026-01-16 (Foundation Phase, Days 1-9)

**Last updated:** 2026-08-25

---

## Project Vision

To design a reproducible, extensible research platform that models disease biology as a networked system, simulates pharmacological perturbations, formulates therapeutic combination discovery as a Hamiltonian (QUBO) optimization problem, and benchmarks exact, heuristic, and quantum-ready solvers across combinatorial scaling regimes.

Long-term goal: enable hybrid quantum–classical approaches for exploring high-order therapeutic combination spaces where classical exact solvers become computationally impractical.

---

## Day 1 — Foundation & Architecture

### Objectives

* Define a full-stack conceptual pipeline for quantum-assisted drug discovery
* Establish research framing and reproducibility mindset
* Create a lab-style workflow and logging structure

### Conceptual Pipeline

```
Disease / Protein Set
        ↓
Protein Interaction Network (Healthy vs Diseased)
        ↓
Drug Input (SMILES / Known Compounds)
        ↓
Classical Docking & Pre-Filtering
        ↓
Quantum Optimization Layer
        ↓
Systems-Level Network Simulation
        ↓
Bayesian Probability Engine
        ↓
Therapeutic Success & Risk Scores
```

### Notes

* The system is designed to treat biological recovery as a **network perturbation problem** rather than a single-target optimization.
* Optimization layer framed as a **Hamiltonian minimization / QUBO maximization** task, allowing future deployment on quantum hardware.

### Reflection

> Established a hybrid architecture integrating biological modeling, probabilistic inference, and Hamiltonian-based optimization. The system is framed as a reproducible research platform rather than a black-box predictor, enabling benchmarking across classical, heuristic, and quantum solvers.

---

## Day 2 — Disease Network Construction

### Objectives

* Build a protein–protein interaction (PPI) network for a reference disease system (Multiple Sclerosis)
* Validate biological plausibility and pipeline integrity

### Implementation

* Used curated protein set related to immune signaling and myelin biology
* Queried interaction confidence scores from external biological knowledge base (STRING-style model)
* Constructed a weighted graph representation

### Results

**Network Summary:**

* Nodes: 10
* Edges: 35

**Sample Interactions:**

```
IFNG ↔ MBP | confidence=0.667
IFNG ↔ VCAM1 | confidence=0.748
IFNG ↔ MOG | confidence=0.823
IFNG ↔ HLA-DRB1 | confidence=0.830
IFNG ↔ STAT3 | confidence=0.907
IFNG ↔ CD40 | confidence=0.927
IFNG ↔ CXCL10 | confidence=0.959
IFNG ↔ IL2RA | confidence=0.967
IFNG ↔ TNF | confidence=0.991
MBP ↔ HLA-DRB1 | confidence=0.851
```

### Observations

* IFNG emerges as a dominant hub, consistent with immune-centric MS pathology
* Network density validates system coherence and biological relevance

### Reflection

> Successfully constructed a high-confidence MS PPI network and validated hub dominance of IFNG, aligning with known immune-mediated disease mechanisms. This reference network serves as the baseline system state for downstream perturbation and optimization experiments.

---

## Day 3 — Disease Distance Metric

### Objectives

* Define a quantitative measure of system deviation between healthy and perturbed network states

### Metric

**System Disease Distance Score**

* Measures aggregate edge-weight deviation under simulated perturbations
* Interpreted as network-level biological stress

### Sample Outputs

```
System Disease Distance Score: 0.1221
System Disease Distance Score: 1.2708
```

### Interpretation

* Low score → near-baseline state
* High score → strongly perturbed / diseased system

### Reflection

> Introduced a network-level scalar metric capturing global biological deviation, enabling consistent downstream mapping between drug perturbations and system recovery potential.

---

## Day 4 — Drug Recovery Scoring

### Objectives

* Simulate pharmacological perturbations
* Map drug effects onto network recovery space

### Example

**Drug:** Fingolimod
**Recovery Score:** 0.4451

### Interpretation

* Recovery score represents estimated normalization of network topology relative to disease state

### Reflection

> Established a functional mapping from drug perturbations to network recovery metrics, enabling formal comparison between pharmacological candidates.

---

## Day 5 — Probabilistic Modeling

### Objectives

* Convert recovery scores into interpretable success probabilities

### Bayesian Output

```
Probability of Therapeutic Success:
Mean: 0.412
95% CI: (0.255, 0.579)
```

### Notes

* Enables uncertainty-aware ranking of drugs
* Provides statistical framing rather than deterministic scoring

### Reflection

> Added probabilistic interpretation layer, enabling uncertainty quantification and confidence intervals for therapeutic success estimates.

---

## Day 6 — Drug Panel Screening

### Objectives

* Screen multiple drugs under unified probabilistic and recovery framework

### Results

```
=== MS Drug Screening Results ===
1. Fingolimod           | Recovery: 0.629 | P(Success): 0.667 [0.537, 0.785]
2. Dimethyl Fumarate  | Recovery: 0.380 | P(Success): 0.463 [0.333, 0.596]
3. Natalizumab        | Recovery: 0.237 | P(Success): 0.185 [0.094, 0.298]
4. Interferon-beta   | Recovery: 0.166 | P(Success): 0.130 [0.055, 0.230]
```

### Reflection

> Demonstrated rank consistency between network recovery and probabilistic success estimates, validating internal coherence of the modeling stack.

---

## Day 7 — Synergy Discovery

### Objectives

* Identify combinatorial effects between drug pairs

### Results

```
=== MS Combination Synergy Results ===
1. Fingolimod + Dimethyl Fumarate | Recovery: 0.881 | Synergy: 0.252 | P(Success): 0.87
2. Fingolimod + Natalizumab      | Recovery: 0.785 | Synergy: 0.156 | P(Success): 0.704
3. Dimethyl Fumarate + Natalizumab | Recovery: 0.617 | Synergy: 0.237 | P(Success): 0.537
```

### Interpretation

* Synergy metric captures non-linear improvement beyond additive recovery

### Reflection

> Identified dominant therapeutic synergies and established a formal combinatorial scoring framework for drug pair selection.

---

## Day 8 — Hamiltonian Formulation & Optimization

### Objectives

* Encode drug selection as a constrained optimization problem
* Validate Hamiltonian formulation

### Formulation

**Objective:**

Maximize:

```
Σ recovery[i] * x_i + Σ synergy[i][j] * x_i * x_j
```

Subject to:

```
Σ x_i = k
```

### Implementation

* Encoded as QUBO / QuadraticProgram
* Solved using exact minimum eigensolver (NumPy baseline)

### Result

```
Selected drugs:
- Fingolimod
- Dimethyl Fumarate
Objective value: 1.261
```

### Reflection

> Validated Hamiltonian formulation by independently rediscovering the same optimal combination identified by network and synergy layers, confirming cross-layer consistency.

---

## Day 9 — Scaling Benchmarks (k = 3)

### Objectives

* Benchmark optimization hardness
* Compare exact, greedy, and random solvers

### Results

| N  | Exact | Greedy | Random |
| -- | ----- | ------ | ------ |
| 4  | 3.158 | 3.158  | 3.158  |
| 6  | 3.180 | 3.180  | 3.180  |
| 8  | 3.118 | 3.118  | 3.118  |
| 10 | 3.634 | 3.634  | 3.634  |
| 12 | 3.842 | 3.575  | 3.842  |

### Interpretation

* Onset of heuristic failure observed at N ≥ 12
* Greedy solver trapped in local optimum

### Reflection

> Identified transition from polynomially tractable to combinatorially hard regime for k=3 selection, validating need for advanced optimization strategies.

---

## Day 9 — Scaling Benchmarks (k = 4)

### Results

| N  | Exact Time (s) | Exact Value | Greedy Value | Random Value |
| -- | -------------- | ----------- | ------------ | ------------ |
| 8  | 0.0064         | 5.188       | 5.188        | 5.188        |
| 12 | 0.0095         | 5.118       | 5.118        | 4.998        |
| 16 | 0.0260         | 5.065       | 5.065        | 5.034        |
| 20 | 0.1879         | 5.226       | 5.143        | 5.033        |
| 24 | 3.3665         | 5.520       | 5.494        | 5.124        |

### Observations

* Exact solver runtime increased by >500× from N=8 to N=24
* Greedy and random solvers remained constant-time but showed growing optimality gaps

### Reflection

> Demonstrated combinatorial scaling behavior and identified a computational regime where classical exact solvers become impractical and heuristics sacrifice solution quality, motivating hybrid quantum–classical optimization strategies.

---

## Public Release

### Repository

* Name: `quantum-bio-recovery-engine`
* Status: Public

### Research Layer Added

* README.md (scientific front door)
* MIT License
* CITATION.cff
* Curated requirements.txt
* One-command `setup.sh`
* Figures and reproducibility instructions

### Reflection

> Released the full research platform as a public, reproducible repository with experimental logs, figures, and setup automation, establishing academic-grade transparency and collaboration readiness.

---

---

# Phase 2 — Publication-Grade MS Combination Screen

*Entries below are date-stamped. The foundation phase above (Days 1–9) was
committed on **2026-01-16**; the work in this phase was carried out on
**2026-08-25**.*

---

## 2026-08-25 — Session A: Conservative re-scoping of the MS screen (v2)

### Motivation

The Day 6–8 screen ranked 4 drugs and reported "P(Success)" values to three
decimals. Those numbers were not calibrated against anything, and the synergy
metric treated a network-distance improvement as evidence of pharmacological
synergy. Both would fail review.

### Changes

* Rebuilt the MS workflow around an interpretable, multi-criteria
  **prioritisation score** rather than a success probability.
* Expanded the panel to **22 candidates** with explicit evidence tiers,
  mechanism classes, directional target hypotheses, and broad safety classes.
* Expanded the signature to **25 signed genes**.
* Excluded shared-mechanism and shared-safety-class pairs from the primary
  ranking, retaining them for audit.
* Added a seeded bootstrap interval over curated target-effect uncertainty.

### Reflection

> Replaced an uncalibrated probability with a transparent, auditable
> prioritisation score, and made the screen's own limitations explicit rather
> than implicit. Correct in kind, but still too small to support a publication:
> 22 candidates, 6 parameters, no significance testing, no controls.

---

## 2026-08-25 — Session B: Extension to a publishable screen (v3)

### Objectives

Scale the candidate space, broaden the evidence base beyond signature reversal,
and add the statistical apparatus a reviewer will require.

### 1. Inputs

| Input | v2 | v3 |
| --- | --- | --- |
| Candidates | 22 | **74** (55 mechanism classes) |
| Signature genes | 25 | **112** across 10 pathways |
| Scoring parameters | 6 | **20** |
| Interactome | none | STRING v12 — 261 nodes, 14,784 edges |

Evidence tiers: 24 approved · 12 phase 3 ·
28 phase 2 · 10 preclinical.

The panel is regenerated by `tools/curate_ms_panel_v3.py`, which **fails the
build** if any target gene is absent from the signature. That check caught two
curation errors on first run (an undeclared risk domain, a stray gene symbol).

### 2. Methodological decisions

**Directional signature.** Raw logFC is the wrong optimisation target for
compensatory transcripts. `HMOX1`, `SOD2`, `TREM2`, and `PDGFRA` rise in MS as a
*protective* response, so a naive "reverse the disease signature" objective
scores dimethyl fumarate — an approved therapy — as harmful. Genes are now
annotated `pathogenic` (73), `protective_deficit` (35), or `compensatory` (4),
and scoring targets the curated therapeutic direction. This annotation is the
assumption most open to challenge and is flagged as such in the protocol.

**Bliss independence replaces clipping.** Same-direction effects combine as
`a + b − ab`; opposing effects add, so antagonism is represented rather than
hidden by a hard clip.

**Double normalisation of reversal.** Absolute reversal over a 112-gene
signature spans only [0, 0.1], because any one agent engages a handful of genes.
Left uncorrected this let low-efficacy, low-risk agents win the composite on the
bonus terms alone — the first full run ranked **High-dose biotin**, a failed
phase 3 negative control, in first place. Adding `reversal_efficiency`
(therapeutic movement ÷ engaged signal) restored a comparable [0,1] scale and
fixed the face-validity failure.

**Target-family redundancy.** The efficacy/safety figure exposed
`Firategrast + Natalizumab-biosimilar` as a top Pareto point. Both block
α4-integrin, but their `mechanism_class` labels differ by modality, so the
string-equality exclusion missed them. Redundancy is now declared by explicit
`target_family` plus a target-overlap rule.

### 3. Statistical apparatus

Permutation null (400 draws/pair, size- and magnitude-preserving) → empirical
*p* → BH-FDR; a weight-free Pareto front; a ±50% weight-jitter sweep; bootstrap
resampling for score **and rank** intervals; and mechanism-stratum enrichment.

### 4. Results

**Headline — mechanism strata (the stable unit of inference):**

| axis pair | n | median score | bootstrap 95% CI | q |
| --- | --- | --- | --- | --- |
| cns innate + remyelination | 79 | 1.226 | [1.149, 1.307] | <1e-05 |
| immunomodulation + remyelination | 299 | 1.190 | [1.129, 1.274] | <1e-05 |
| cns innate + metabolic repair | 147 | 1.153 | [1.105, 1.236] | 1e-05 |
| metabolic repair + remyelination | 92 | 1.151 | [1.079, 1.284] | 9e-05 |
| neuroprotection + remyelination | 105 | 1.150 | [1.080, 1.264] | 9e-05 |
| cns innate + neuroprotection | 160 | 1.142 | [1.087, 1.224] | 0.0001 |
| cns innate + cns innate | 56 | 1.137 | [1.052, 1.237] | 0.11 |
| immunomodulation + neuroprotection | 619 | 1.135 | [1.099, 1.192] | <1e-05 |
| cns innate + immunomodulation | 421 | 1.126 | [1.088, 1.194] | 0.00015 |
| remyelination + remyelination | 17 | 1.115 | [0.947, 1.243] | 0.64 |
| immunomodulation + metabolic repair | 527 | 1.099 | [1.057, 1.162] | 1 |
| metabolic repair + neuroprotection | 185 | 1.077 | [1.037, 1.177] | 1 |
| neuroprotection + neuroprotection | 103 | 1.050 | [0.958, 1.156] | 1 |
| metabolic repair + metabolic repair | 81 | 1.038 | [0.942, 1.145] | 1 |
| immunomodulation + immunomodulation | 754 | 1.021 | [0.975, 1.054] | 1 |

**Top pairs (diversity-capped at 3 appearances per agent):**

| # | combination | score | q | separation |
| --- | --- | --- | --- | --- |
| 1 | Clemastine + Minocycline | 1.511 | 0.0083 | 0.28 |
| 2 | Minocycline + Opicinumab | 1.502 | 0.0083 | 0.28 |
| 3 | Minocycline + Ocrelizumab | 1.488 | 0.0083 | 0.22 |
| 4 | Clemastine + Ocrelizumab | 1.467 | 0.0083 | 1.00 |
| 5 | Clemastine + Ibudilast | 1.446 | 0.0083 | 0.22 |
| 6 | Opicinumab + Tolebrutinib | 1.427 | 0.0083 | 0.47 |
| 7 | Ibudilast + Ocrelizumab | 1.426 | 0.0083 | 0.22 |
| 8 | Methylprednisolone + Siponimod | 1.423 | 0.0083 | 0.05 |
| 9 | Methylprednisolone + Opicinumab | 1.422 | 0.0083 | 0.54 |
| 10 | Ibudilast + Ofatumumab | 1.413 | 0.0241 | 0.35 |

**Controls, declared before ranking:**

Redundant positive control (Natalizumab + Natalizumab-biosimilar):
excluded as expected — **True**.

| safety control | best primary rank | in top 30 |
| --- | --- | --- |
| Daclizumab | 299 | no |
| Cyclophosphamide | 317 | no |
| Mitoxantrone | 385 | no |

| negative-efficacy control | best primary rank | in top 30 |
| --- | --- | --- |
| Ustekinumab | 16 | yes |
| High-dose biotin | 13 | yes |
| Opicinumab | 2 | yes |
| Evobrutinib | 242 | no |

**Robustness:**

* Weight sensitivity: Spearman **0.9742** mean, 0.9174 min under ±50% jitter.
* Pairs at q < 0.05: **1,315 / 1,793**.
* Pareto-optimal pairs: 361 of 2,016.
* Top-30 bootstrap Jaccard: **0.2458**.

### 5. The negative result, stated plainly

**Individual pair ranks are not reliable.** Under the curated target-effect
uncertainty the top-30 bootstrap Jaccard is 0.2458, and rank
95% CIs span hundreds of positions. A leaderboard of named pairs is therefore
*not* a defensible output of this screen, and the runner prints an explicit
warning whenever that statistic falls below 0.5.

Stratum medians behave differently: their bootstrap intervals are narrow
relative to the spread between strata, and the ordering survived the
target-family redundancy fix unchanged. The defensible finding is therefore the
**ordering of mechanism combinations**, not any specific pair.

Three negative-efficacy controls (Opicinumab, High-dose biotin, Ustekinumab)
still rank inside the top 30.
This is reported rather than filtered: the screen scores mechanism and
transcriptional direction, not trial outcome, and this bounds what it can claim.

### 6. Interpretation

Pairing a **CNS-innate agent with a remyelination agent** is the top-ranked
stratum (1.226), while pairing **two peripheral
immunomodulators** ranks last (1.021) despite being
the largest stratum (n=754) and the one most represented in
current practice. Every stratum pairing remyelination or CNS-innate action with
a second axis outranks every immunomodulation-only combination.

This is consistent with compartmentalised inflammation behind an intact
blood-brain barrier as the reason peripheral immunosuppression plateaus in
progressive MS. It is a hypothesis the screen was built to generate, not
evidence that it is true.

### 7. Artefacts

```
core/biology/ms_scoring.py           20-parameter pair scoring
core/biology/network_proximity.py    STRING proximity / separation
core/biology/screen_statistics.py    nulls, FDR, Pareto, stability, strata
experiments/ms/run_extended_screen.py
experiments/ms/plot_screen.py
tools/curate_ms_panel_v3.py          panel curation + validation
tools/fetch_string_network.py        cached interactome fetch
tests/                               19 tests, all passing
figures/fig_mechanism_strata.png
figures/fig_rank_stability.png
figures/fig_efficacy_safety_tradeoff.png
```

Runs are deterministic under a stated seed and carry a provenance block
recording git revision, Python/NumPy versions, platform, and the exact command.

### Reflection

> Scaled the MS screen from 22 to 74 candidates and 6 to 20
> parameters, and added the significance, robustness, and control apparatus a
> manuscript requires. The most important output was not the ranking but the
> demonstration that the ranking is unstable at the pair level — which relocated
> the claim to the mechanism stratum, where it survives bootstrap resampling, a
> ±50% weight sweep, and a redundancy fix. Two curation errors and one
> face-validity failure were caught by the validation and control machinery
> rather than by inspection, which is the argument for building it.

### Next steps

1. Replace the illustrative signature with discovery + validation cohorts,
   stratified by RRMS / active SPMS / PPMS.
2. Re-derive `target_effects` from measured pharmacology — the rank instability
   is driven by `target_uncertainty`, so this is the highest-value fix.
3. Wire the Hamiltonian/QAOA layer to the v3 scores to select higher-order
   (k > 2) combinations under the same safety and redundancy constraints.
4. Validate the top stratum in glia/oligodendrocyte co-culture with full
   dose matrices (Bliss, HSA, Loewe, ZIP).
5. Preprint.
