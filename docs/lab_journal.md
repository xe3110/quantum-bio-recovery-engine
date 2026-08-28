# Quantum Bio Recovery Engine — Lab Journal

**Researcher:** Rishabh Kumar

**Project:** Quantum Bio Recovery Engine

**Focus:** Hybrid biological–probabilistic–combinatorial–quantum optimization for therapeutic combination discovery

**Start Date:** 2026-01-16 (Foundation Phase, Days 1-9)

**Last updated:** 2026-08-29 (Phase 6)

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

---

---

# Phase 3 — De Novo Molecular Design

*Carried out on **2026-08-26**.*

---

## 2026-08-26 — Session A: From selecting molecules to designing one

### Motivation

Phase 2 ended with a screen that ranks pairs drawn from a fixed panel of 74
agents. That question has a ceiling built into it: **the answer can only ever
be a molecule someone has already made.** In MS every approved agent is a
peripheral immunomodulator, and the screen's own headline finding was that
peripheral-only combinations rank last. The screen could identify the gap and
could not, in principle, fill it.

So the question changed from *"which two existing agents pair best?"* to *"what
should a molecule do, and what would such a molecule look like?"*

That required three things the repository did not have: a specification of what
a single molecule must do, a chemistry layer that can represent and measure a
structure, and a generator that searches chemical space against the
specification.

### Architectural decision: build it disease-agnostic

Phase 2's vocabularies — `MS_DISEASE_PATHWAYS`, `THERAPEUTIC_AXES`,
`RISK_DOMAINS` — were module constants inside the MS scoring code. Fine for one
disease and wrong for two: adding Parkinson's would have meant editing scoring
logic. Phase 3 introduces a **disease model** (`core/models/disease.py`) loaded
from a registry entry under `data/diseases/`. Everything downstream reads the
context and never names a disease.

### 1. The chemistry model

RDKit has no wheel for the interpreter this project pins, so `core/chemistry/`
was written from scratch: SMILES **parser and writer** (a generator must emit
structures, not only read them), valence model, ring perception, Ertl TPSA,
Wildman–Crippen logP, Lipinski/Veber, Wager CNS MPO, a declared
synthetic-tractability proxy, structural alerts as readable graph predicates,
and ECFP-style circular fingerprints.

Validation state, all test-enforced:

| Quantity | Status |
| --- | --- |
| SMILES round-trip | 15-molecule corpus, composition preserved exactly |
| Formula and mass | all 42 curated structures match published values |
| TPSA | matches published Ertl values exactly |
| cLogP | reduced Crippen typing; ~0.4 MAE vs experiment |
| Stereochemistry | absent — parsed and discarded |
| Aromaticity | trusted as written, not re-perceived |

Sanity check on real MS drugs behaved correctly on first run: dimethyl fumarate
and diroximel fumarate flagged as Michael acceptors (that *is* their Nrf2
mechanism), fingolimod's primary-amine pKa found, minocycline scored low on
tractability as a complex polycyclic natural product.

### 2. The Target Product Profile

`core/design/target_profile.py` ranks candidate proteins by

```
priority = leverage x tractability x (1 - 0.5 x liability) x (1 + 0.35 x axis_gap)
```

with leverage split 60/40 between weighted signature membership and
distance-decayed influence through STRING.

**The tractability term is the single most important guard in the pipeline.**
A new annotation, `data/targets/druggability.json`, assigns each of the 93
panel targets a class and a small-molecule prior: 22 at or above 0.6, 32 below
0.2. CD20 has among the highest leverage in MS and a tractability of **0.05**.
Without the floor the engine would confidently specify a small molecule to bind
the target of ocrelizumab. Sub-floor targets are reported as **readouts** —
the profile still wants their expression to move, but no arm is pointed at
them.

**Axis-gap analysis, the campaign's most defensible output:**

| therapeutic axis | unmet fraction among approved agents |
| --- | --- |
| remyelination | **1.00** |
| neuroprotection | 0.96 |
| cns_innate | 0.92 |
| metabolic_repair | 0.88 |
| immunomodulation | **0.00** |

This is Phase 2's stratum finding restated as a specification, and derived
independently.

### 3. Fragment selection as a Hamiltonian

Choosing which pharmacophores to fuse is a constrained binary optimisation —
coverage saturates under Bliss, chemotypes duplicate, and a shared mass budget
couples every pair. That is a QUBO, solved on three backends: exhaustive
enumeration (ground truth), an exact Ising eigensolver, and QAOA on Aer.

Two approximations, both documented rather than hidden:

1. **Coverage truncation** — exact at `k=2`, second-order beyond.
2. **Budget linearisation** — whole-molecule budgets cannot be written as a
   quadratic, so each fragment is charged against a `1/k` share. Subadditivity
   makes this **conservative**: the QUBO may reject a design that would have
   fitted, and never admits one that does not. Test-enforced.

### 4. Assembly

A design is a recipe — scaffold, an arm per attachment point, a linker each,
caps on the remainder — assembled into a real molecular graph, capped, and
measured. Search is exhaustive enumeration over assemblies followed by seeded
hill-climbing, then novelty assessment and MaxMin diversity selection.

---

## 2026-08-26 — Session B: What went wrong, and what it changed

Three failures are recorded here because they are the useful part of the
session.

### Failure 1 — Coverage-only optimisation designs molecules that cannot reach the brain

The first full run selected a **statin acid plus a sulfonylurea**: excellent
profile coverage, TPSA 170–250 against an 40–80 window, MW up to 554 against a
420 ceiling. **Every** top design failed the CNS gate.

The Hamiltonian was optimising coverage while blind to the delivery constraint.
Adding a polarity penalty fixed TPSA and **just moved the failure**: the
optimiser swapped to an adamantane and a chloroarene and landed at cLogP 5.8–7.2,
failing the same gate from the other side.

The fix required the envelope to enter the Hamiltonian on **polarity, donors,
and lipophilicity at once**, plus a mass budget that reserves atoms for the
scaffold and linkers assembly adds afterwards. Handing the arms the whole
molecular-weight ceiling was why 550 Da molecules kept appearing against a
420 Da envelope.

### Failure 2 — A silent bug in the Hamiltonian, caught by a test

Quadratic couplings were stored keyed on the fragments' **ranked** order and
looked up in **sorted** order. Every coupling whose ranking disagreed with the
alphabet read back as zero. The optimiser was solving a different Hamiltonian
than the one being reported, and nothing about the output looked wrong.

It surfaced only because a test asserted the truncated objective agrees with an
independently-recomputed exact objective at `k=2`. The lookup now raises on a
missing pair rather than defaulting to zero. **All results below are post-fix.**

### Failure 3 — QAOA returned an infeasible state that scored better than the optimum

At `reps=2`, QAOA selected **four** arms where three were allowed, scoring
*above* the true optimum precisely because the extra arm was never paid for.

| reps | ansatz depth | feasible | matches optimum |
| --- | --- | --- | --- |
| 1 | 55 | yes | no |
| 2 | 87 | **no** | no |
| 3 | 119 | yes | **yes** |
| 4 | 151 | yes | no |
| 6 | 215 | yes | no |

Textbook behaviour: `reps=3` finds the optimum, and deeper circuits do *worse*
at a fixed iteration budget as the variational landscape gets harder. Every
result now carries a `feasible` flag, because reading the objective without it
inverts the conclusion. Default depth set to 3.

### A blind spot in CNS MPO, found by a failing test

A test asserting memantine should out-score a statin acid on CNS MPO failed —
and the test premise was wrong, not the implementation. **Wager's score takes
the *most basic* pKa as its ionisation term**, correctly for the basic and
neutral compounds it was derived on, and is therefore blind to acids: a small
carboxylic acid scores above 5 of 6 while being anionic at pH 7.4 and
effectively barred from the brain.

The implementation stayed faithful to Wager. A separate `acidic_centres()`
check was added and the delivery gate now consults both.

---

## 2026-08-26 — Session C: Results

> **Superseded.** The numbers in this session predate the Phase 4 bug fixes
> (three chemistry errors found by RDKit cross-validation, and per-claim
> evidence weighting). They are left as the record of what was believed at the
> time; the current results are in the Phase 4 summary.

### Backend agreement

```
Hamiltonian : 10 binary variables, 45 couplings, choose 2 (45 feasible of 1024 states)
  enumeration  match  objective=+0.0901  feasible=True  0.024s
  eigensolver  match  objective=+0.0901  feasible=True  0.022s
  qaoa         match  objective=+0.0901  feasible=True  1.313s
```

### Target profile (top 8 of 14)

| gene | class | direction | priority |
| --- | --- | --- | --- |
| MMP9 | enzyme | down | 0.7883 |
| BTK | kinase | down | 0.7665 |
| NLRP3 | enzyme | down | 0.7651 |
| NOS2 | enzyme | down | 0.7442 |
| CSF1R | kinase | down | 0.7016 |
| RORC | nuclear receptor | down | 0.6994 |
| S1PR1 | GPCR | down | 0.6505 |
| MTOR | kinase | down | 0.5450 |

32 further targets classified as readouts and excluded from binding.

### Designs

7,350 recipes attempted, 225 rejected as chemically invalid, 7,125 distinct
structures.

| # | formula | MW | cLogP | TPSA | CNS-MPO | Tanimoto to nearest known |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | C16H17ClN4O | 316.8 | 1.45 | 48.5 | 5.83 | 0.41 |
| 2 | C18H18ClFN4O | 360.8 | 2.07 | 59.0 | 5.99 | 0.34 |
| 3 | C20H20FN3O3 | 369.4 | 3.14 | 71.5 | 5.13 | 0.30 |

Top design: `N1(CCN(CC1)c2ccc(Cl)cc2)C(=O)Nc3cccnc3`

All above the CNS-MPO gate of 4.0, all inside the property envelope, all novel
at a 0.6 Tanimoto threshold.

### The negative result, stated plainly

**Under a CNS envelope, three arms do not fit.** With the mass budget correctly
reserving atoms for the scaffold and linkers, every `k=3` selection scores
negative, while `k=2` lands at TPSA 62, cLogP 2.5, 19 heavy atoms. A
CNS-penetrant three-mechanism single molecule is not reachable from this
fragment library. The tool reports that rather than emitting 550 Da molecules
that fail the gate.

### The two stages disagree, and that is informative

The arm set the Hamiltonian ranks first is **not** the arm set that assembles
into the best molecule:

| arm set | Hamiltonian objective | best assembled fitness |
| --- | --- | --- |
| BTK + CSF1R | **+0.0901** (rank 1) | 1.5439 |
| BTK + TLR4 | +0.0822 (rank 2) | 1.4919 |
| CSF1R + TLR4 | +0.0783 (rank 3) | **1.6027** |

The QUBO scores profile coverage under a linearised envelope; design fitness
adds developability, tractability, and structural alerts to a molecule that
actually exists. The covalent BTK arm carries a Michael-acceptor alert and more
mass, and pays for both only at the second stage — so the last-ranked arm set
produced the winning structure.

Caught only because the campaign carries several arm sets forward instead of
the optimum alone. Reporting "the Hamiltonian's optimum" as though it were "the
best design" would have been wrong, and the two are now reported separately
throughout.

### Interpretation

The Hamiltonian's optimum, **BTK + CSF1R**, is the pairing of the two leading
CNS-penetrant MS mechanisms — B-cell and microglial on one arm, microglial on
the other. That it was reached from a directional signature, an interactome,
and a druggability annotation, without ever being told, is the strongest
internal-consistency check available at this stage.

The top *assembled* molecule, `N1(CCN(CC1)c2ccc(Cl)cc2)C(=O)Nc3cccnc3`, is a
**dual CSF1R / TLR4 antagonist aimed at microglial activation in progressive
MS** — CNS-MPO 5.83 of 6, coordinated suppression of the microglial activation
programme (CSF1R −0.75, AIF1 −0.69, TLR4 −0.60, CD68 −0.50, downstream
TNF/IL-1β/IL-6/CCL2), and no counter-therapeutic movement at all.

Its efficacy is **unknown and untested even for binding.** Scored on the panel's
own scale it reaches a signature reversal of 0.0519 — second of 75, above the
approved-agent median of 0.0301 — but that comparison is close to circular,
because its arms were selected to maximise alignment with this very signature
while the panel agents were not. It is an internal-consistency check, not
evidence of effect.

It also serves **one axis only** (`cns_innate`), and not the one with the
largest gap: remyelination sits at an unmet fraction of 1.00 and this molecule
does nothing for it. The envelope forced that — two arms of budget, both spent
on the same axis. The molecule the profile argues for is exactly the one that
does not fit.

The defensible outputs therefore remain structural rather than chemical: the
axis-gap asymmetry, and the demonstration that a three-mechanism CNS molecule
does not fit the envelope. **Individual molecules are the least reliable thing
the campaign produces**, for the same reason pair ranks were in Phase 2 — they
sit downstream of the most uncertain inputs.

### Artefacts

```
core/chemistry/molecule.py           SMILES parser + writer, valence, assembly
core/chemistry/descriptors.py        Ertl TPSA, Crippen logP, shape, flexibility
core/chemistry/druglikeness.py       Lipinski, Veber, CNS MPO, alerts, acids
core/chemistry/fingerprint.py        circular fingerprints, novelty, diversity
core/models/disease.py               disease model + registry loader
core/design/target_profile.py        Target Product Profile derivation
core/design/pharmacophores.py        fragment library loader + validation
core/design/quantum_assembly.py      QUBO, enumeration / eigensolver / QAOA
core/design/denovo.py                assembly, scoring, search
experiments/design/run_denovo_design.py
data/diseases/multiple_sclerosis.json
data/targets/druggability.json       93 targets annotated
data/chemistry/pharmacophore_library.json   48 fragments
data/chemistry/ms_known_structures.json     42 verified structures
tests/                               142 tests, all passing
docs/denovo_design_protocol.md
```

### Reflection

> Extended the platform from selecting existing drugs to specifying and
> assembling new ones, and made the whole stack disease-agnostic in the
> process. The most valuable outputs were again negative or structural: that
> coverage-blind optimisation produces molecules that cannot reach the tissue,
> that penalising one property just moves the failure to another, and that a
> three-mechanism CNS molecule does not fit the envelope this library can
> reach. A silent Hamiltonian bug and a documented blind spot in a published
> CNS score were both caught by tests rather than by inspection, which is the
> same argument Phase 2 made for building the validation machinery.

### Next steps

1. **Dock every proposal, then assay it.** Fragment-inherited engagement is an
   assumption, not a prediction, and it is the first thing that will break.
2. Have a chemist assess synthesisability — the tractability score is a
   declared proxy, not a route.
3. Profile `kinase_aminopyrimidine_hinge` selectivity across the kinome. A
   designed multi-target ligand and an uncontrolled polypharmacology liability
   are the same molecule seen from two sides.
4. Expand the fragment library along the remyelination axis, where the profile
   reports the largest unmet gap and the library is thinnest.
5. Add ADMET: metabolic stability, hERG, CYP, and especially P-glycoprotein
   efflux, which decides many CNS programmes and is unmodelled here.
6. Populate a second disease registry entry end to end to exercise the
   generality claim rather than asserting it.

---

---

# Phase 4 — Hardening: Generality, Validation, and Provenance

*Carried out on **2026-08-26**, in response to a review that identified five
gaps. Each is recorded below with what was actually done about it, including
where the answer was "the claim was too strong" rather than "the code was
wrong".*

---

## 2026-08-26 — Session D: The review

Five gaps, all correct:

1. **"Disease-configurable", not proven disease-agnostic.** One registry entry,
   and the new layers still imported `core.biology.ms_scoring` directly.
2. **Hypothesis generation, not drug-design prediction.** No binding,
   selectivity, assay loop, or validated ADMET.
3. **Reimplementing chemistry instead of using RDKit is the highest technical
   risk.** 42 formula checks do not validate perception, logP, or fingerprints.
4. **Fragment engagement is manually inherited** with no per-claim evidence,
   confidence, or feedback.
5. **Result files treated as durable truth** rather than reproducible run
   artifacts.

Gap 1 was the sharpest: the import graph contradicted the README. A layer that
calls itself disease-agnostic while importing MS scoring is not one, whatever
its registry says.

---

## 2026-08-26 — Session E: Chemistry validated against RDKit (gap 3)

RDKit has no wheel for the pinned interpreter, but the project's second
environment is Python 3.11, where it installs cleanly. That turned gap 3 from
unanswerable into measurable.

### Architecture change

RDKit is now the **preferred backend**, used automatically wherever importable;
`core/chemistry/` is explicitly fallback infrastructure.
`core.chemistry.backend` resolves the choice, records it on every descriptor
vector, refuses to silently downgrade when RDKit is requested and absent, and
is pinned to `local` for the test suite so assertions do not depend on the
environment.

### What cross-validation measured

`tools/validate_chemistry.py` compares both implementations over curated drugs,
library fragments, and **300 generated structures** — the molecules the
pipeline actually emits, which exercise ring systems no curated set contains.

| quantity | agreement |
| --- | --- |
| **structural perception** | **388 / 388 canonical SMILES match** |
| HBD, HBA | exact everywhere |
| TPSA | exact on designs; worst case 3.53 A^2 on curated drugs |
| rotatable bonds | 0.91–1.00 exact |
| ring count | 0.91–0.98 exact |
| **cLogP** | **MAE 0.65–0.86, worst case 2.23** |
| fingerprints | Spearman 0.982; 4/565 pairs disagree at the novelty threshold |

The perception result is the one that matters and the one that was previously
unevidenced: parser and writer together preserve the molecule as RDKit
understands it, on every structure tested.

### Three real bugs, none visible to composition checks

* **Rotatable bonds were wrong on essentially every generated structure** —
  0% exact, bias **+2.75**. Only amide C–N bonds were excluded where the strict
  definition excludes any conjugated carbonyl–heteroatom bond, and assembly
  produces carbamates and anhydrides freely. Now 0.91–1.00.
* **Aromatic rings bearing exocyclic carbonyls were mistyped in cLogP by
  +3 log units** (caffeine, uracil): a purinone ring carbon typed as an
  aromatic ether rather than a carbonyl. Any purinone design would have had its
  CNS MPO badly distorted.
* **Ring perception invented a macrocycle in adamantane.**
  `networkx.cycle_basis` is a spanning-tree basis, not a minimum one, and
  returned an eight-membered ring where three six-membered rings exist. The
  tractability proxy penalises rings above seven atoms, so every
  adamantane-containing design was charged 0.35 for a macrocycle it does not
  have. Switched to `minimum_cycle_basis`.

One over-correction is recorded too: matching RDKit's *atom-centric* strict
rotatable-bond rule exactly swung the bias negative, so the bond-level rule was
kept and the residual disagreement documented rather than chased. Chasing exact
parity with a reference implementation is a poor use of effort once the
reference itself is available — the right answer was to prefer RDKit.

---

## 2026-08-26 — Session F: A second disease, end to end (gap 1)

### Contracts extracted

`Signature`, `load_signature`, `bliss_combine`, `combine_effects`,
`alignment_metrics`, and `EVIDENCE_RANK` moved to `core/biology/signature.py`.
`ms_scoring` re-exports them, so no caller broke. Risk weighting moved out of
`ScoringConfig` and into the registry, because how strongly infection or
dyskinesia constrains use is a property of the disease and its population, not
of a scoring function. A test now asserts no module under `core/design/`
mentions `ms_scoring`.

### Parkinson's, built from scratch

| | MS | Parkinson's |
| --- | --- | --- |
| Signature | 112 genes | **90 genes, 13 pathways** |
| Panel | 74 agents | **35 agents, 29 mechanism classes** |
| Interactome | STRING v12, 261 nodes | **STRING v12, 240 nodes** |
| Druggability | 93 targets | **90 targets** |

The panel is weighted deliberately toward failures — creatine, CoQ10,
isradipine, inosine, cinpanemab — because Parkinson's has an unusually
well-documented record of neuroprotection trials that did not work, which makes
strong negative controls easy to declare in advance.

Two things the second disease forced into the design:

* **Gene aliasing.** STRING v12 still calls glucocerebrosidase `GBA`; HGNC says
  `GBA1`. Unmapped, the most common genetic risk factor in PD had **zero
  network leverage** and dropped out of the profile silently. Aliases are now
  registry data, applied at load, so the cached interactome stays a faithful
  record of what the source returned.
* **A sharper test of the tractability floor.** Alpha-synuclein — the central
  protein in the disease — is intrinsically disordered, prior 0.25. The profile
  routes around its most important target through lysosomal and autophagic
  mechanisms, which is what the clinical field actually does.

The PD profile ranks **LRRK2 first**, then SLC6A3, NLRP3, MAOB, GBA1. Axis
gaps: `symptomatic_dopaminergic` 0.00, `synuclein_proteostasis` and
`trophic_support` both 1.00.

### The library learned what it was missing

Run against PD, `unreachable_requirements()` reported **10 of 14 targets with
no chemical matter** — the library was built for MS. Eleven PD pharmacophores
were added (LRRK2 aminopyrimidine, propargylamine, nitrocatechol, aminotetralin,
iminosugar, dihydropyridine, hydroxypyridinone, ketoamide, and others), taking
the library to 59 fragments and the gap to two, both symptomatic-dopaminergic
targets a disease-modifying design would not pursue anyway.

That diagnostic driving library growth is the intended workflow: fragments are
disease-agnostic chemical matter, and the tool names what it cannot reach
rather than quietly producing something worse.

---

## 2026-08-26 — Session G: Evidence, provenance, and honest status (gaps 2, 4, 5)

### Per-claim evidence (gap 4)

All **154 fragment–target claims** now carry provenance. Each pharmacophore
declares an `evidence_tier` for its chemotype and the `primary_targets` it
claims to *bind*; everything else in its map is a downstream transcriptional
consequence and is discounted by half. Confidence weights profile coverage
directly, so the optimiser prefers well-evidenced arms instead of treating a
speculative inference like approved-drug pharmacology.

Library composition: 10 approved-drug, 8 clinical-candidate, 9
published-chemotype, 4 speculative.

This visibly changed the MS result. The top design is now a **dual BTK / CSF1R
inhibitor serving two therapeutic axes** (`cns_innate` + `immunomodulation`),
where the previous top design served one.

### Evidence status (gap 2)

The binding, selectivity, and ADMET gaps cannot be closed without models and
data this project does not have. What *was* fixed is that the caveats now
travel with the data: every candidate carries a machine-readable
`evidence_status` block enumerating **13 unassessed axes** — target binding,
selectivity, cell activity, in-vivo efficacy, permeability, BBB transport, P-gp
efflux, metabolic stability, CYP, hERG, solubility, synthetic route, freedom to
operate — each with what would satisfy it and why its absence matters.
`readiness` reads `hypothesis_only`. A consumer that never opens the protocol
still receives the caveat.

### Run artifacts (gap 5)

`core/provenance.py` records the SHA-256 digest of every input, the git
revision **and whether the tree was dirty**, interpreter, platform, tracked
package versions, and the active chemistry backend. Results moved to
`experiments/design/results/<disease>/`, are git-ignored, and carry a README
distinguishing disposable snapshots from curated `data/`.

The digest is the part that matters: a version string in a metadata block does
not survive someone editing a CSV.

---

## 2026-08-26 — Session H: QAOA's default was luck

Re-running the quantum benchmark across two diseases after the library grew
exposed something the single-disease campaign had hidden. QAOA success is
**non-monotonic in circuit depth and specific to the problem instance**: depth
3, adopted as the default because it solved the original MS instance, then
failed on *both* campaigns, while depths 2, 4, and 5 each succeeded on some
instances and not others.

| depth | MS | Parkinson's |
| --- | --- | --- |
| 1 | +0.0443 ✓ | +0.0732 |
| 2 | +0.0443 ✓ | +0.0830 ✓ |
| 3 | +0.0355 | +0.0830 ✓ |
| 4 | +0.0443 ✓ | +0.0830 ✓ |
| 5 | +0.0443 ✓ | +0.0830 ✓ |

The runner now sweeps depths 1–5, keeps the best *feasible* selection, and
prints what every depth returned. Reporting one lucky depth would have
presented a heuristic as though it were a solver.

### The stage inversion held for both diseases

| disease | Hamiltonian rank 1 | best molecule came from |
| --- | --- | --- |
| MS | BTK + RORgt (+0.0443) | BTK + CSF1R, **rank 3** |
| PD | GLUT + MAO-B (+0.0830) | caspase-1 + MAO-B, **rank 2** |

Confirming that carrying several arm sets forward, rather than the optimum
alone, is load-bearing rather than defensive.

---

## Phase 4 results

**Multiple sclerosis** — `N1(CCN(CC1)C(=O)Nc2cccnc2)C3CCN(CC3)C(=O)C=C`
C18H25N5O2, MW 343.4, cLogP 1.41, TPSA 68.8, CNS-MPO 5.33, Tanimoto 0.36.
Dual BTK/CSF1R, two axes, mean engagement confidence 0.40.

**Parkinson's** — `c1(ccc(cc1)CC(C)N(C)CC#C)OC(=O)C(=O)NCC2CC2`
C19H24N2O3, MW 328.4, cLogP 1.61, TPSA 58.6, CNS-MPO 4.93, Tanimoto 0.32.
MAO-B inhibitor fused to a caspase-1 warhead, mean engagement confidence 0.47.

Both rank first against their disease panels on signature reversal — a
**near-circular** comparison, since their arms were selected to maximise
alignment with the signature they are then scored against, and reported only
with that caveat attached.

Test suite: **174 passing** without RDKit, **198 with it**.

### Reflection

> Four of the five gaps were closed and the fifth was made honest. The most
> valuable work was again negative: cross-validation against RDKit found three
> chemistry bugs that a year of formula-and-mass tests would never have caught,
> and the second disease found a gene-alias failure that silently zeroed the
> most important target in Parkinson's. Both were invisible in the output.
> The pattern across all four phases is consistent — the machinery built to
> doubt the results is what produces the results worth keeping. Where a claim
> could not be supported, the claim moved rather than the evidence: RDKit
> became the preferred backend instead of the reimplementation being defended,
> and QAOA's depth became a swept parameter instead of a lucky constant.

### Next steps

1. **Dock the proposals.** Everything else is downstream of the
   fragment-transplantation assumption, and it remains untested.
2. Close the evidence loop: `core/design/evidence.py` has the hook for assay
   results to update fragment confidences, and nothing populates it.
3. Add ADMET where it can be done honestly — P-gp efflux first, since it
   decides CNS programmes on its own.
4. Replace both signatures with cohort-derived data; the PD one is weaker than
   the MS one and both are illustrative.
5. A third disease outside the CNS, to exercise the peripheral branch of the
   delivery constraint, which neither current entry does.

---

## Phase 5 — The second disease gets its own screen (2026-08-26)

### Objective

Parkinson's was registered in Phase 4 and exercised end to end by the *design*
campaign, but it had no combination screen of its own. The MS screen could not
be pointed at it: `core/biology/ms_scoring.py` holds MS's pathways, therapeutic
axes, risk domains, and risk weights as module constants, so running it against
Parkinson's would have scored PD combinations against MS's vocabulary and
reported `infection` as a worst risk domain for a disease where nothing is
constrained by infection risk.

Two things had to change, and one thing deliberately did not.

### 1. Vocabulary moved from the code to the registry

`core/biology/combination_scoring.py` reads every vocabulary from a
`DiseaseContext`. Adding a third disease now touches no scoring code. A test
asserts the property directly rather than by inspection: scoring a PD pair must
produce a `worst_risk_domain` that is one of PD's domains and **is not** one of
MS's.

### 2. The screen learned to count past two

The MS screen is pairwise. Parkinson's is treated with genuine polypharmacy —
levodopa plus a decarboxylase inhibitor plus a COMT inhibitor plus a MAO-B
inhibitor is an ordinary regimen — so a pairwise-only screen would have been
answering the wrong question.

`combination_metrics` scores a combination of any order. Every k-ary term is
defined as the **mean of its pairwise form**, which was the design constraint
worth holding: it makes k = 2 numerically identical to the pairwise
definitions, so nothing about the MS results is being quietly restated in a new
arithmetic, and it gives k ≥ 3 a definition that privileges no member. Tests
pin the reduction for target complementarity, safety union, and Bliss folding,
and assert that scoring is invariant to the order the members are listed in.

### 3. The comparison that needed a rule, not a number

The headline ask was monotherapy versus combination. The obvious
implementation — rank everything by `priority_score` and see what wins — is
wrong, and quietly so.

`priority_score` includes target complementarity, compartment complementarity,
and network separation. A single agent scores **zero on all three by
construction**. Ranking singles against pairs by composite score therefore
concludes that combinations win, when what they actually did was collect
bonuses that are structurally unavailable to a monotherapy. The number would
have looked like a finding.

So the rule is stated instead: **`priority_score` is comparable only within a
fixed order**, and across orders only the efficacy block — defined identically
at every k — may be compared. `monotherapy_comparison` reports that block and
nothing else, and a test asserts `priority_score` never enters the comparable
set. Three fields carry the actual comparison:
`reversal_gain_over_best_single`, `score_gain_over_best_subset` (negative means
a third agent does not earn its place), and `additivity_ratio`.

The last of those needed fixing after first sight. Its median is 1.0 at every
order, because most agents in the panel have disjoint target sets and Bliss is
exactly additive on disjoint genes — the median is the uninformative half of
the statistic. The tail is the finding: 17% of pairs and 36% of triples are
sub-additive, meaning their members are covering the same signal. The report
now carries `fraction_subadditive`, the 5th percentile, and the minimum.

### The bug that a smoke run would not have caught

The first full run returned `nan` for every order-3 stratum bootstrap interval
and a top-K Jaccard of exactly 0.0 — not a warning, a structural zero.

The order loop reassigns `prefilter` at the end of each iteration to seed the
next order. The robustness analyses ran *after* that reassignment, so at order 3
they re-enumerated the triple space filtered by a set of **3-tuples** rather
than the 2-faces that built it. Nothing matched, every re-ranking returned an
empty list, and the statistics dutifully reported `nan` over nothing.

It is the failure mode this project keeps producing: not a crash, an
answer-shaped absence of an answer. The Jaccard of 0.0 was the tell — a real
instability produces a small number, not a perfectly round one. The prefilter
in force for the current order is now held in its own name.

### The result that did not need the leaderboard

The pattern from the MS screen reproduced exactly, at both orders: individual
combination ranks are **not** stable under the curated target-effect
uncertainty, while **stratum medians are**. Same conclusion, different disease,
different vocabulary, different panel. The unit of inference is the mechanism
stratum.

Consistent with this, both screens agree on the shape of the disease they were
given: the strata that rank highest are the ones pairing a symptomatic or
metabolic agent with the axis that has **no approved agent at all** — trophic
support in Parkinson's, remyelination in MS. A screen bounded by its panel can
identify that gap and cannot, in principle, fill it. That is what the de novo
design campaign is for.

### What I chose not to build

Parkinson's regimens are chronic oral polypharmacy in an elderly population, so
the binding constraint on co-administration is pharmacokinetic — CYP and COMT
interactions, additive hypotension, serotonin syndrome with MAO-B inhibitors.
No such model exists here. Building one would have made the PD results
non-comparable with the MS v3 numbers until MS was re-run against it, so the
screen ships with route burden and half-life spread (`regimen_burden`) as the
only co-administration cost, and the gap declared in the caveats block of every
result file rather than only in prose.

It is the top-priority next model.

### Reflection

> The interesting work was again in refusing a number rather than producing
> one. `priority_score` across orders would have given a clean,
> publishable-looking answer to the question I was actually asked, and it would
> have been an artefact of which terms a single agent can structurally earn.
> Writing down the comparison rule, then testing that the rule cannot be
> violated, was worth more than the ranking it constrains. The `nan` bug makes
> the same point from the other direction: the analysis did not fail, it
> succeeded on an empty set, and only an implausibly round robustness figure
> gave it away.

### Next steps

1. **PK/DDI model** — the largest declared gap, and the one that binds hardest
   on the disease just screened.
2. Migrate the MS screen onto `combination_scoring` so both diseases run one
   contract, and re-run MS at k = 3.
3. Dock the design proposals; still the assumption everything else sits on.
4. Replace both signatures with cohort-derived data.
5. A non-CNS disease, to exercise the peripheral branch of the delivery
   constraint that neither current entry does.

---

## Phase 6 — Parkinson's drug discovery, and a silent novelty failure (2026-08-28)

### Objective

Bring the Parkinson's *design* campaign to the same standard as the MS one.
Phase 4 had already run it end to end and Phase 5 gave PD its own combination
screen, so on paper this was a re-run at parity. It was not.

### The bug: an optional field that degraded a result without saying so

Novelty in the design campaign is assessed against

```python
reference = [*disease.known_structures(), *library.parent_structures()]
```

`known_structures()` reads `data.structures` from the registry entry and
returns `[]` when the key is absent. **Parkinson's had no `structures` key.**

So every Parkinson's design was having its novelty measured against nothing but
the pharmacophore fragments it had just been assembled from. The consequence is
worse than a missing number would have been:

- the `nearest_known_compound` field was **populated**, with
  `Propargylamine MAO-B warhead`;
- the Tanimoto was **plausible** at 0.32;
- `is_novel` was **True**;
- and nothing anywhere in the JSON, the CSV, or the console output indicated
  that the comparison set had been cut in half.

A missing file would have been caught in a minute. A silently reduced
comparison survived two phases and got written into a protocol document.

With `data/chemistry/pd_known_structures.json` registered — 26 structures
covering the panel's small molecules — the same molecule's nearest neighbour is
**selegiline** at **Tanimoto 0.47**. It still clears the declared 0.6 threshold,
so the headline conclusion holds, but the margin is much narrower and it is now
measured against an approved Parkinson's drug rather than against the pipeline's
own parts. In hindsight the old answer was obviously wrong: the design *is* a
propargylamine, so of course selegiline is its nearest neighbour. The number
only looked defensible because selegiline was not in the comparison.

Curating the 26 structures went the way this repo's chemistry work usually
goes — three of my first 27 transcriptions were wrong (safinamide's amide was
on the wrong carbon, minocycline had three carbons too many, and my *literature*
mass for deferiprone was wrong rather than the SMILES). All three were caught by
re-deriving formula and mass from the structure and comparing against the
published values, which is exactly what that guard exists for. Every entry was
additionally cross-checked against RDKit before admission.

### The guards, because "remember to add the file" is not a guard

- `test_each_disease_declares_a_novelty_reference_set` — fails if a registered
  disease has no reference set, or one that barely overlaps its own panel.
- `test_novelty_is_measured_against_real_drugs_not_only_fragments` — fails if
  the set adds nothing the fragment library already contained, which would make
  registering it a no-op.
- `tests/test_chemistry.py` now iterates **every** registered disease's
  structure file rather than naming the MS one, so a third disease's set is
  formula- and mass-validated the moment it is added, without new test code.

The campaign protocol's "adding a disease" checklist gained the reference set as
step 5, with the reason attached rather than left as an instruction.

### A second reporting defect, found while writing it up

The design runner prints `design["candidates"]` — the diversity-selected
deliverable — and wrote `design["all_ranked"]` to the CSV. Those are different
lists drawn from different pools, and neither contains the other: **three of the
four designs in the report were absent from the CSV entirely**, which is the
file anyone would actually open.

The CSV now carries the union, deduplicated and ordered by fitness, with a
`selected_for_diversity` column distinguishing the deliverable from the raw
ranking. Both diseases were regenerated so their artifacts stay comparable.

### Results

Parkinson's, at parity with MS (`--k 2 --arms 3 --top 4 --quantum-benchmark`,
RDKit backend):

- 14-requirement target profile led by LRRK2 (0.965), SLC6A3, NLRP3, MAOB, GBA1.
- `DDC` and `SLC18A2` reported **unreachable** — no chemical matter in the
  library. Still the most useful line a run produces.
- 10-variable QUBO, 45 couplings, 45 feasible states of 1,024. Enumeration,
  exact eigensolver, and QAOA all reach +0.0830; QAOA needs depth ≥ 2, with
  d1 falling short at +0.0732.
- Axis gaps: `synuclein_proteostasis` and `trophic_support` both **1.00**
  unmet, `symptomatic_dopaminergic` 0.00.
- Leading design `C19H24N2O3`, MW 328.4, CNS-MPO 4.93/6 — a MAO-B inhibitor
  fused to a caspase-1 warhead.

The **stage inversion held again**: the Hamiltonian's rank-1 arm set
(GLUT + MAO-B, +0.0830) did not build the best molecule; rank 2
(caspase-1 + MAO-B, +0.0732) did. Third disease, same lesson — carrying several
arm sets forward is load-bearing.

The axis-gap result is also the *screen's* Phase 5 finding arrived at
independently: the combination screen's top mechanism strata all pair something
with `trophic_support`, and the design campaign — reading only the panel's
approved agents, never the combination scores — reports that axis at 1.00 unmet.

Test suite: **236 passing** without RDKit, **260 with it**.

### Reflection

> Both defects this phase were the same shape as the `nan` bug in Phase 5:
> nothing failed. An optional registry key degraded a scientific claim while
> leaving a well-formed, plausible number in place, and a CSV confidently
> presented the wrong four molecules. Neither would have been caught by reading
> the output, because the output looked right — they were caught by asking what
> the number was actually measured against, and by noticing that two views of
> the same result disagreed. The fix in both cases was to make the quiet failure
> loud: a test that fails when a reference set is missing, a column that names
> which rows are the deliverable. I am increasingly convinced that the useful
> question for this project is not "is this number right?" but "what would this
> number look like if it were wrong?" — and that where the answer is "exactly
> the same", that is the bug to go hunting.

### Next steps

1. **PK/DDI model** — still the largest declared gap, unchanged from Phase 5.
2. Migrate the MS screen onto `combination_scoring` and re-run MS at k = 3.
3. **Dock the proposals.** Three phases running, still the assumption
   everything else rests on, and now the only untested claim standing between
   these structures and a reason to make one.
4. Replace both signatures with cohort-derived data; PD's remains the weaker.
5. A non-CNS disease, to exercise the peripheral branch of the delivery
   constraint — and, now, to be the first disease that gets its reference set
   right on the first pass.
