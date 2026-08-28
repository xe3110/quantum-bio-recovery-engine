# Quantum Bio Recovery Engine

*A hybrid quantum–classical framework for modeling disease networks, simulating drug perturbations, and optimizing high-order therapeutic combinations using Hamiltonian-based optimization and scaling benchmarks.*

---

## Overview

Quantum Bio Recovery Engine is a research platform for discovering optimal multi-drug therapeutic combinations in complex disease systems. The framework integrates:

- Protein–protein interaction (PPI) network modeling
- Probabilistic drug response estimation
- Synergy discovery across therapeutic combinations
- Hamiltonian (QUBO) formulation of drug selection
- Exact, heuristic, and quantum-ready optimization backends
- Scaling benchmarks for computational hardness analysis
- **De novo molecular design** — deriving a target product profile from a
  disease model and assembling new chemical structures against it
- Disease-agnostic contracts, exercised end to end for **multiple sclerosis
  and Parkinson's**

The system is designed for **reproducibility, extensibility, and future quantum hardware deployment**.

---

## Conceptual Pipeline

The platform runs two campaigns over one disease model. The **selection**
campaign ranks combinations of agents that already exist; the **design**
campaign specifies and builds one that does not.

```
                    Disease Model (registry entry)
             signature · interactome · panel · druggability
                                |
             +------------------+------------------+
             |                                     |
      SELECTION campaign                    DESIGN campaign
             |                                     |
   Drug Perturbation Simulation          Target Product Profile
             |                          (leverage x tractability
   Synergy & Probability                 x liability x axis gap)
             |                                     |
   Hamiltonian (QUBO) over drugs        Hamiltonian (QUBO) over
             |                             pharmacophores
   Exact / Heuristic / QAOA                        |
             |                          Exact / Eigensolver / QAOA
   Ranked combinations                             |
                                          Fragment Assembly
                                                   |
                                     Developability · CNS gate ·
                                       Novelty · Diversity
                                                   |
                                        Proposed structures
                                          (SMILES + formula)
```

## Repository Structure
```
quantum-bio-recovery-engine/
├── core/ # Biology, chemistry, design, and Hamiltonian modules
│ ├── biology/ # Signature scoring, k-ary combination scoring, network proximity, statistics
│ ├── chemistry/ # SMILES parser/writer, descriptors, developability, fingerprints
│ ├── design/ # Target profiles, fragment library, QUBO assembly, de novo engine
│ ├── models/ # Disease model and registry loader
│ ├── quantum/ # Quantum optimizer for the drug-selection Hamiltonian
│ └── provenance.py # Input digests, environment capture, run artifacts
├── experiments/
│ ├── ms/ # Multiple sclerosis: pairwise combination screen
│ ├── parkinsons/ # Parkinson's: monotherapy + pair + triple screen
│ ├── design/ # De novo design campaigns (disease-agnostic)
│ └── benchmarks/ # Scaling and hardness benchmarks
├── tools/ # Input curation and external-database fetchers
├── tests/ # Behavioural tests for the screen and the design stack
├── data/
│ ├── diseases/ # Disease registry entries (multiple_sclerosis, parkinsons)
│ ├── chemistry/ # Fragment library and known-structure reference set
│ └── targets/ # Per-target druggability annotations
├── figures/ # Generated plots and benchmark figures
├── docs/ # Protocols and the dated lab journal
├── requirements.txt # Reproducible environment specification
├── README.md
├── LICENSE
└── CITATION.cff
```

---

## Installation

### 1. Create Environment
```bash
python -m venv qbio-env
source qbio-env/bin/activate
```

### 2. Install Dependencies
```bash
pip install -r requirements.txt
```

### One-Command Setup

For a clean, reproducible environment:

```bash
./setup.sh
```

## Reproducibility

### Scaling Benchmark (k = 4)
```bash
python -m experiments.benchmarks.scale_benchmark
```

### Generate Figures
```bash
python experiments/benchmarks/plot_scaling.py
```

Figures will be saved to:
```
figures/
├── fig_runtime.png
└── fig_quality.png
```

### MS combination-candidate screen

A deterministic, multi-parameter, discovery-stage screen over **74 MS candidates**
(2,701 pairs) against a **112-gene directional disease signature** and a cached
**STRING v12 interactome**. It scores 20 parameters per pair spanning efficacy,
pharmacology, network topology, evidence, and a 7-domain safety burden, then
attaches permutation significance (BH-FDR), a weight-free Pareto front, bootstrap
rank stability, and pre-declared positive/negative/safety controls.

```bash
python -m experiments.ms.run_extended_screen --top 30 --seed 7
python -m experiments.ms.plot_screen
pytest -q tests/
```

Outputs land in `experiments/ms/results/`: `ms_extended_screen.json` (full report
with provenance) and `ms_pairs_full.csv` (every pair, every parameter, including
the redundant pairs excluded from the primary ranking).

**Read the mechanism strata, not the leaderboard.** Individual pair ranks are not
stable under the curated target-effect uncertainty (top-30 bootstrap Jaccard
≈ 0.25), and the runner warns when that is the case. Stratum-level ordering *is*
stable, and it is the screen's defensible unit of inference.

Regenerate the inputs with `python -m tools.curate_ms_panel_v3` (validates the
panel against the signature) and `python -m tools.fetch_string_network`.

See [the MS publication protocol](docs/ms_publication_protocol.md) for data
provenance, every scoring assumption, the control design, and the validation work
required for a manuscript. This is a hypothesis generator, not a clinical
combination recommendation.

### Parkinson's monotherapy, pair, and triple screen

The same idea one order higher. A **35-agent** Parkinson's panel is scored
against a **90-gene directional signature** and a cached STRING interactome at
**k = 1, 2, and 3** — monotherapy, pairs, and triples — through a single
disease-agnostic scorer that reads its pathways, therapeutic axes, risk domains
and risk weights from the registry entry rather than from any constant in the
code.

```bash
python -m experiments.parkinsons.run_combination_screen --top 25 --seed 7
python -m experiments.parkinsons.run_combination_screen --quick      # smoke run
```

Outputs land in `experiments/parkinsons/results/`:
`pd_combination_screen.json` (full report with provenance) and
`pd_combinations_full.csv` (every combination at every order, including the
redundant ones excluded from the primary ranking).

Monotherapy runs through the **same function** as a triple, so the
single-versus-combination comparison is a comparison of biology rather than of
two implementations. Three fields carry it: `reversal_gain_over_best_single`
(what a combination adds over its own best member),
`score_gain_over_best_subset` (whether a third agent earns its place — negative
means it does not), and `additivity_ratio` (below 1 the members are covering the
same signal).

> **`priority_score` is comparable only within a fixed order.** It includes
> terms a single agent cannot earn — target complementarity, compartment
> complementarity, network separation. Across orders, compare the efficacy
> block, which is defined identically at every *k*. A test asserts the composite
> never enters the cross-order comparison.

See [the Parkinson's screen protocol](docs/parkinsons_screen.md) for the
results and the disease's own caveats, and [the disease campaign
protocol](docs/disease_campaign_protocol.md) for the standing procedure every
campaign follows.

### De novo molecular design

Where the screen above ranks **existing** agents, this campaign designs a new
one. It derives a Target Product Profile from the disease model — which
proteins are worth engaging, in which direction, inside which physicochemical
envelope — formulates pharmacophore selection as a **QUBO**, solves it by
exhaustive enumeration, an exact Ising eigensolver, and **QAOA**, then
assembles, scores, and novelty-filters real molecular structures.

```bash
python -m experiments.design.run_denovo_design --disease multiple_sclerosis
python -m experiments.design.run_denovo_design --quantum-benchmark --k 2
```

Outputs land in `experiments/design/results/<disease>/`: `denovo_designs.json`
(target profile, Hamiltonian, backend agreement, full provenance) and
`denovo_candidates.csv` (every design, every property).

These are **run artifacts, not curated data** — git-ignored, regenerated on
demand, and safe to delete. Each carries the SHA-256 digest of every input, the
git revision and whether the tree was dirty, the tracked package versions, and
the active chemistry backend. A number from one is not interpretable without
the digests beside it.

A run reports the profile it derived, the Hamiltonian it solved, whether the
backends agreed, and the structures it built:

```
Hamiltonian  : 10 binary variables, 45 couplings, choose 2 (45 feasible of 1024 states)
  enumeration  match    objective=+0.0443 feasible=True 0.061s
  eigensolver  match    objective=+0.0443 feasible=True 0.042s
  qaoa         match    objective=+0.0443 feasible=True 5.535s
               depth sweep: d1=+0.0443  d2=+0.0443  d3=+0.0355  d4=+0.0443  d5=+0.0443

  1. C18H25N5O2  fitness=1.4374
     N1(CCN(CC1)C(=O)Nc2cccnc2)C3CCN(CC3)C(=O)C=C
     MW=343.4 cLogP=1.41 TPSA=68.8 CNS-MPO=5.33 coverage=0.099
     arms: btk_acrylamide_piperidine + csf1r_picolinamide
     nearest known: Piperidinyl acrylamide covalent warhead (Tanimoto 0.3636, novel=True)
```

That structure is a **dual BTK / CSF1R inhibitor directed at compartmentalised
CNS inflammation**; the Parkinson's campaign produces a MAO-B inhibitor fused
to a caspase-1 warhead. Both are walked through in full — including what their
efficacy is not — in
[§8 of the design protocol](docs/denovo_design_protocol.md#8-reading-a-design--worked-examples).

**Read the target profile and the mechanism combinations, not the ranking of
individual molecules.** For MS the findings that survive scrutiny are
structural: approved therapy covers the immunomodulation axis completely
(gap 0.00) and remyelination not at all (gap 1.00), and under a CNS envelope a
three-mechanism single molecule does not fit the mass and polarity budget —
every `k=3` selection scores negative.

There is **no binding-affinity model anywhere in this pipeline.** Target
engagement is inherited from curated fragment annotations, not predicted from
structure. Nothing has been docked, simulated, synthesised, or assayed. Every
candidate carries a machine-readable `evidence_status` block naming the **13
assessment axes on which nothing has been done** — binding, selectivity,
permeability, P-gp efflux, hERG, CYP, metabolic stability, solubility,
synthetic route, and more — so the caveat travels with the data rather than
living only in a protocol.

Each of the library's **154 fragment–target claims** carries its own evidence
tier and confidence, distinguishing a binding claim backed by an approved drug
from a downstream transcriptional inference. Confidence weights coverage
scoring directly, so the optimiser has a reason to prefer well-evidenced arms.

See [the de novo design protocol](docs/denovo_design_protocol.md) for the full
assumption set and the validation work that would have to come first.

### Chemistry backend

**RDKit is used whenever it can be imported.** `core/chemistry/` is a
self-contained fallback — SMILES parser and writer, Ertl TPSA, Crippen logP,
circular fingerprints — for interpreters that have no RDKit wheel, including
the one this project pins. The choice is explicit, recorded on every
descriptor vector, and forceable with `QBRE_CHEMISTRY_BACKEND=local`.

The fallback's divergence is measured, not asserted:

```bash
python -m tools.validate_chemistry --designs 300   # needs RDKit
```

Current results (`docs/chemistry_validation.json`): **388/388 canonical SMILES
match** across curated drugs, library fragments, and generated designs; HBD and
HBA exact; TPSA and rotatable bonds 0.83–1.00 exact; fingerprint Spearman
0.982. **cLogP is the weak one at MAE 0.65–0.86** — wide enough to move a
molecule across a CNS MPO boundary, and the main reason RDKit is preferred.

Building that harness found three real bugs invisible to formula-and-mass
checks: rotatable bonds wrong on nearly every generated structure, aromatic
carbonyls mistyped in cLogP by three log units, and a phantom macrocycle in
adamantane inflating its synthetic-tractability penalty.

### Registered diseases

Two, so the disease-agnostic claim is checkable rather than asserted:

| | multiple sclerosis | Parkinson's |
| --- | --- | --- |
| Signature | 112 genes, 10 pathways | 90 genes, 13 pathways |
| Panel | 74 agents | 35 agents |
| Interactome | STRING v12, 261 nodes | STRING v12, 240 nodes |
| Druggability | 93 targets annotated | 90 targets annotated |
| Known structures | 42 agents | 26 agents |
| Unserved axis | remyelination (1.00) | synuclein proteostasis, trophic support (1.00) |
| Combination screen | pairs | monotherapy + pairs + triples |

They share **no therapeutic axis** and overlap on only two safety domains
(cardiac, hepatic — generic organ toxicity). Nothing in Parkinson's care is
constrained by infection risk; it is constrained by dyskinesia and psychosis.
Tests assert that separation, and assert that no module under `core/design/`
imports disease-specific scoring.

```bash
python -m experiments.design.run_denovo_design --disease parkinsons --quantum-benchmark
```

Each disease declares a **known-structure reference set** that the design
campaign measures novelty against. It is optional in the loader and that turned
out to be a trap: Parkinson's ran its first design campaign without one, so
novelty fell back to the pharmacophore library's own fragments and a design was
reported as more novel than it is. Two tests now fail if a registered disease
omits the set or fills it with chemical matter the library already had. See
[the de novo protocol](docs/denovo_design_protocol.md#the-novelty-figure-this-design-used-to-carry-was-wrong-2026-08-28).

### Adding a disease

Pathways, therapeutic axes, risk domains, risk weights, gene aliases, and the
delivery constraint all live in the registry rather than in scoring code.
Adding a disease is a data task. The full checklist — including the controls a
panel must declare and the two claims the test suite enforces — is in [the
disease campaign protocol](docs/disease_campaign_protocol.md).

1. Drop a signature CSV (`gene, logFC, direction_class, desired_direction,
   pathway, confidence`), a cached interactome TSV, and a candidate panel JSON
   into `data/`.
2. Annotate each panel target in a druggability file — `target_class` and a
   `small_molecule_tractability` prior. This is what stops the engine
   specifying a small molecule against an antibody-only antigen.
3. Write `data/diseases/<name>.json` pointing at those files and declaring the
   vocabularies, risk weights, any gene aliases (interaction databases lag
   HGNC — STRING still calls `GBA1` by its old name `GBA`), and
   `delivery.requires_cns_exposure`.

```bash
python -m experiments.design.run_denovo_design --disease <name>
cp experiments/parkinsons/run_combination_screen.py experiments/<name>/   # change DISEASE
```

The loader validates every referenced path on read and names what is missing,
so a half-populated registry entry fails at load rather than mid-campaign.
Fragments in `data/chemistry/pharmacophore_library.json` are shared across
diseases; `unreachable_requirements()` reports profile targets that no fragment
covers, which is how the library learns what it is missing.

## Documentation

| Document | What it covers |
| --- | --- |
| [Disease campaign protocol](docs/disease_campaign_protocol.md) | **Start here for a new disease.** What a registry entry must carry, the k-ary scoring contract, the cross-order comparison rule, the statistical treatment, and the checklist for adding the next disease |
| [MS publication protocol](docs/ms_publication_protocol.md) | Screen inputs, all 20 scoring parameters, statistical treatment, controls, and the validation required for a manuscript |
| [Parkinson's screen protocol](docs/parkinsons_screen.md) | The k = 1/2/3 screen: monotherapy vs combination results, mechanism strata, controls, and the caveats specific to chronic dopaminergic polypharmacy |
| [De novo design protocol](docs/denovo_design_protocol.md) | Target profile derivation, the Hamiltonian and its two approximations, the chemistry model's validation state and blind spots, and the central transplantation assumption |
| [Lab journal](docs/lab_journal.md) | Dated research log across all four phases, including the failures and what they changed |
| [Chemistry validation](docs/chemistry_validation.json) | Machine-readable cross-validation of the local chemistry stack against RDKit, regenerated by `tools/validate_chemistry.py` |

## Scientific Motivation

Therapeutic combination discovery is inherently combinatorial. As the number of candidate drugs increases, exact classical optimization becomes intractable, while heuristic methods sacrifice solution quality. This platform provides a benchmarked, Hamiltonian-based formulation suitable for evaluating hybrid and quantum-assisted optimization strategies in biological and pharmacological research.

The same combinatorial structure appears one level down. Selecting which
pharmacophores to fuse into a single multi-target molecule is a constrained
binary optimization over fragments whose benefits interact — coverage
saturates, chemotypes duplicate, and a shared mass budget couples every pair.
That is a QUBO, and the design campaign solves it on the same backends as the
combination problem. At the sizes reported here the quantum backends
demonstrate that the formulation transfers to a quantum algorithm; they
demonstrate no advantage, and the repository says so wherever the claim could
be misread.

## Disclaimer

This project is a research and educational platform only.
It does not provide medical advice, clinical recommendations, or validated treatment guidance.

The combination rankings are **not** co-administration recommendations:
concurrent immunomodulation carries additive infection and malignancy risk that
no in-silico score can quantify. The designed structures are **not** drug
candidates: they carry no binding data, no synthesis route, and no safety
assessment, and target engagement is an inherited assumption rather than a
prediction. Both outputs are hypotheses for experimental validation.
