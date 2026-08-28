# De novo molecular design — protocol

A disease-agnostic pipeline that derives a **Target Product Profile** from a
disease model, formulates **pharmacophore selection as a QUBO**, solves it
exactly and on QAOA, and **assembles novel molecular structures** against the
resulting specification. Registered and exercised end to end for two diseases:
multiple sclerosis and Parkinson's.

It is **not** a drug-discovery result. There is no binding-affinity model
anywhere in this pipeline. Nothing has been docked, simulated, synthesised, or
assayed. Every structure it emits is a hypothesis for a chemist to argue with.

This campaign shares its disease model with the
[MS combination screen](ms_publication_protocol.md) — the same signature,
interactome, and candidate panel — and inherits every input limitation
documented there. The dated development log, including the calibration failures
summarised in §7, is in the [lab journal](lab_journal.md).

---

## 1. What changed from the combination screen

The extended screen answers *"which two existing agents pair best?"* over a
fixed panel. That question has a hard ceiling: the answer can only ever be a
molecule someone already made, and in MS every approved agent is a peripheral
immunomodulator. This pipeline asks the prior question — *what should a
molecule do at all?* — and then builds one.

| | Combination screen | De novo design |
|---|---|---|
| Search space | 2,701 pairs of 74 known agents | assembled structures from 48 fragments |
| Output | a ranking of existing drugs | SMILES, molecular formula, property vector |
| Target selection | implicit in the panel | explicit, tractability-filtered profile |
| Modality realism | none — antibodies and small molecules ranked together | small-molecule tractability gates every target |
| Optimisation | weighted score + Pareto front | QUBO / Ising Hamiltonian, three backends |
| Delivery | `cns_penetration` as a scored term | CNS envelope as a **gate** on the design |
| Novelty | not applicable | Tanimoto against the disease's own known drugs (MS 42, PD 26) + 31 parent chemotypes |

## 2. The disease model

A disease is a registry entry under `data/diseases/`, not code. It names the
signature, interactome, panel, and druggability annotation, its controlled
vocabularies, and its **delivery constraint**. Adding a disease is a data task.

`data/diseases/multiple_sclerosis.json` sets
`delivery.requires_cns_exposure = true`. That flag is the whole argument of the
campaign: every approved MS therapy acts largely in the periphery and every one
plateaus in progressive disease, so a molecule designed for progressive MS has
to reach CNS parenchyma. The flag turns the CNS multi-parameter score from a
preference into a gate.

### Two diseases, and why the second one matters

An abstraction exercised by one instance is a claim, not a design. The
registry holds two:

| | multiple sclerosis | Parkinson's |
| --- | --- | --- |
| Signature | 112 genes, 10 pathways | 90 genes, 13 pathways |
| Panel | 74 agents | 35 agents |
| Interactome | STRING v12, 261 nodes | STRING v12, 240 nodes |
| Therapeutic axes | immunomodulation, cns_innate, remyelination, neuroprotection, metabolic_repair | symptomatic_dopaminergic, synuclein_proteostasis, mitochondrial_rescue, neuroinflammation_control, trophic_support |
| Safety domains | infection, malignancy, autoimmunity, teratogenicity, ocular, cardiac, hepatic | dyskinesia, impulse_control, psychiatric, orthostatic_hypotension, somnolence, gastrointestinal, cardiac, hepatic |

The two share **no therapeutic axis** and overlap on only two of their safety
domains -- cardiac and hepatic, generic organ toxicity that genuinely applies
to both. Nothing in Parkinson's care is constrained by infection risk; it is
constrained by dyskinesia and psychosis. Those vocabularies, and the weighting
over them, are registry data. Tests assert the separation, and assert that no
module under `core/design/` imports disease-specific scoring, because a layer
that calls itself disease-agnostic while importing `ms_scoring` is not one
whatever its registry contains.

Two details the second disease forced into the design:

* **Gene aliasing.** STRING v12 still calls glucocerebrosidase `GBA` where the
  current HGNC symbol is `GBA1`. Unmapped, the most common genetic risk factor
  in Parkinson's had zero network leverage and dropped out of the profile
  silently. Aliases are declared per disease in the registry rather than
  patched into the cached interactome, so the cache stays a faithful record of
  what the source returned.
* **A sharper test of the tractability floor.** Alpha-synuclein is the central
  protein in the disease and is intrinsically disordered, with a tractability
  prior of 0.25. The profile has to route around its most important target
  through lysosomal and autophagic mechanisms instead -- which is what the
  clinical field actually does.

## 3. Target Product Profile

`core/design/target_profile.py` ranks candidate proteins by

```
priority = leverage x tractability x (1 - 0.5 x safety_liability)
                     x (1 + 0.35 x unmet_axis_fraction)
```

- **leverage** — 60% direct (weighted signature membership), 40% topological
  (distance-decayed influence over signature genes through STRING, decay 0.5,
  horizon 3 steps).
- **tractability** — the per-target prior from `data/targets/druggability.json`
  (93 targets annotated: 22 at or above 0.6, 32 below 0.2).
- **safety_liability** — risk-weighted burden inherited from the panel agents
  already known to engage the target.
- **unmet_axis_fraction** — how poorly the target's therapeutic axes are served
  by *approved* agents.

**The tractability term is the most important guard in the pipeline.** CD20 has
among the highest leverage in MS and a tractability of 0.05. Without the floor,
the engine would confidently specify a small molecule to bind the target of
ocrelizumab. Targets below the floor are reported as **readouts**: the profile
still wants their expression to move, but no arm is pointed at them.

For MS the axis-gap analysis reports `immunomodulation = 0.00` (fully served by
approved drugs) and `remyelination = 1.00` (nothing approved). That asymmetry,
not any individual molecule, is the campaign's most defensible output.

## 4. The Hamiltonian

Choosing `k` arms from `n` fragments is a constrained binary optimisation:

```
maximise    sum_f b_f x_f  +  sum_{f<g} c_fg x_f x_g
subject to  sum_f x_f = k,   x_f in {0,1}
```

`b_f` is a fragment's solo profile coverage minus its costs; `c_fg` is the
interaction `cov(f,g) - cov(f) - cov(g)` minus redundancy. Three backends solve
it: exhaustive **enumeration** (ground truth), an exact **eigensolver** on the
Ising Hamiltonian via Qiskit, and **QAOA** on the Aer simulator.

### Two approximations, both deliberate

1. **Coverage truncation.** Exact at `k = 2`; a second-order approximation
   beyond it, slightly over-crediting three-way combinations of a saturating
   function.
2. **Budget linearisation.** Mass, polar surface, donors, and lipophilicity are
   whole-molecule budgets that a quadratic form cannot express, so each
   fragment is charged against a `1/k` share. Since
   `sum_f max(0, x_f - s) >= max(0, sum_f x_f - k*s)`, this **never
   under-charges**: the QUBO may reject a design that would have fitted, and
   never admits one that does not.

`exact_objective()` recomputes the true value with joint budgets, and every
solver result is reported against it.

### What the quantum backends do and do not show

QAOA reproducing the enumeration optimum shows the **formulation is faithful
and transfers to a quantum algorithm**. It does **not** show a quantum
advantage, and at ten binary variables there is none to show -- enumeration
finishes in 60 ms and the QAOA depth sweep takes about 6 s.

QAOA's behaviour is worth reading honestly, and it got less flattering the
more it was tested. An early run at depth 2 returned an **infeasible** state
selecting four arms where three were allowed, scoring *above* the true optimum
precisely because the extra arm was never paid for. Depth 3 was then adopted
as the default because it solved that instance.

That default was luck. Across the two registered diseases, success is
**non-monotonic in depth and specific to the problem instance**: after the
fragment library grew, depth 3 failed on *both* campaigns while depths 2, 4,
and 5 each succeeded on some instances and not others. Deeper circuits carry
more variational parameters and a harder classical landscape, so they do not
reliably improve on shallower ones at a fixed iteration budget.

Reporting one lucky depth would present a heuristic as though it were a
solver. The runner therefore **sweeps depths 1-5, keeps the best feasible
selection, and prints what every depth returned**:

```
qaoa   match  objective=+0.0443 feasible=True 5.5s
       depth sweep: d1=+0.0443  d2=+0.0443  d3=+0.0355  d4=+0.0443  d5=+0.0443
```

Every result carries a `feasible` flag, and reading the objective without it
inverts the conclusion.

## 5. The chemistry model

RDKit is the **preferred backend and is used whenever it can be imported.**
`core/chemistry/` is a self-contained fallback for interpreters that cannot
install it -- which includes the one this project pins, where no wheel exists.
`core.chemistry.backend` resolves the choice, `QBRE_CHEMISTRY_BACKEND=local`
forces the fallback, requesting `rdkit` where it is absent raises rather than
silently downgrading, and **every descriptor vector records which backend
produced it**.

The consequence is worth stating plainly: a campaign run with RDKit present
and one run without it can produce different numbers. Both are recorded in the
provenance block of every result file.

### Measured divergence, not asserted equivalence

`tools/validate_chemistry.py` runs both implementations over the curated
drugs, the fragment library, and 300 **generated** structures -- the molecules
the pipeline actually emits, which exercise ring fusions and substitution
patterns no curated set contains. Current numbers, regenerated into
`docs/chemistry_validation.json`:

| quantity | agreement with RDKit |
| --- | --- |
| **structural perception** | **388 / 388 canonical SMILES match** |
| HBD, HBA | exact everywhere |
| TPSA | exact on generated designs; 0.83 exact on curated drugs, worst case 3.53 A^2 |
| rotatable bonds | 0.91-1.00 exact, worst case 1 |
| ring count | 0.91-0.98 exact, worst case 1 |
| molecular weight | within 0.014 Da |
| **cLogP** | **MAE 0.65-0.86, worst case 2.23** |
| fingerprint similarity | Spearman 0.982; 4 of 565 pairs disagree at the 0.6 novelty threshold |

Structural perception is the result that matters most and the one that was
least evidenced before: the local writer's output is re-read by RDKit and
canonicalised, then compared against RDKit's canonicalisation of the input.
Agreement on all 388 means the parser and writer together preserve the
molecule as RDKit understands it.

**cLogP is the weak one.** A mean absolute error near 0.7 with a worst case
above 2 is wide enough to move a molecule across the CNS MPO lipophilicity
boundary. The local logP should not be used to make a lipophilicity call on
its own, and it is the main reason RDKit is preferred.

### What cross-validation found

Building the harness was worth more than the numbers it produced, because it
found bugs that formula-and-mass checks cannot see:

* **Rotatable bonds were wrong on essentially every generated structure**
  (0% exact, bias +2.75). Only amide C-N bonds were excluded, where the strict
  definition excludes any conjugated carbonyl-heteroatom bond -- and fragment
  assembly produces carbamates and anhydrides freely. Now 0.91-1.00 exact.
* **Aromatic rings bearing exocyclic carbonyls were mistyped in cLogP**, by
  **+3 log units** on caffeine and uracil. A purinone ring carbon was being
  typed as an aromatic ether rather than a carbonyl. Any purinone or
  pyrimidinone design would have had its CNS MPO badly distorted.
* **Ring perception reported a phantom macrocycle in adamantane.**
  `networkx.cycle_basis` is a spanning-tree basis and is not guaranteed to
  return the smallest cycles; it returned an eight-membered ring where three
  six-membered rings exist. Since the tractability proxy penalises any ring
  larger than seven atoms, every adamantane-containing design was being
  charged for a macrocycle it does not have. Now uses `minimum_cycle_basis`.

A residual difference from RDKit remains by construction: RDKit's ring count
uses a *symmetrised* SSSR, which reports extra rings on symmetric bridged
cages. Neither is wrong; they answer different questions.

### Two documented blind spots

- **CNS MPO cannot see acids.** Wager's score takes the *most basic* pKa as
  its ionisation term, correctly for the basic and neutral compounds it was
  derived on. A small carboxylic acid can score above 5 of 6 while being
  anionic at pH 7.4 and effectively barred from the brain. The delivery gate
  therefore consults `acidic_centres()` separately. This was found by a test
  whose premise was wrong -- the implementation was faithful to Wager, and
  the gap is Wager's.
- **pKa and logD are estimated** from functional-group class, not measured. A
  CNS MPO computed this way ranks within one campaign and should not be
  compared against published MPO values.

Stereochemistry is absent throughout: parsed and discarded, so every
descriptor is stereo-blind and no claim requiring a defined stereocentre is in
scope.

## 6. The central assumption, and how it is now qualified

`data/chemistry/pharmacophore_library.json` pairs each fragment with the
target-engagement hypothesis its chemotype carries. Assembly gives a molecule
the union of its arms' engagement vectors, combined under Bliss independence.

**This assumes a motif keeps its activity when transplanted onto a scaffold it
did not evolve on.** That is how medicinal chemists reason, and it is wrong
often enough that discovering where it breaks is most of what a real programme
does. No value in that file is a measured Ki, IC50, or EC50.

### Per-claim evidence and confidence

Left unqualified, that assumption is uniform: a scaffold assembled entirely
from speculative downstream inference scores exactly like one built on
approved-drug pharmacology, and looks far more convincing than its evidence
warrants. So every one of the library's **154 fragment-target claims** now
carries provenance.

Each pharmacophore declares an `evidence_tier` for its chemotype and the
`primary_targets` it is claimed to *bind*. Everything else in its engagement
map is a downstream transcriptional consequence -- a materially weaker claim.
Per-claim confidence follows:

| chemotype evidence | binding claim | downstream claim |
| --- | --- | --- |
| `approved_drug` | 0.90 | 0.45 |
| `clinical_candidate` | 0.75 | 0.38 |
| `published_chemotype` | 0.60 | 0.30 |
| `speculative` | 0.30 | 0.15 |

The library holds 10 approved-drug, 8 clinical-candidate, 9
published-chemotype, and 4 speculative pharmacophores. **Confidence weights
profile coverage directly**, so the optimiser has a reason to prefer
well-evidenced arms rather than treating all claims alike. Each design reports
its mean engagement confidence and names its weakest claim.

The tiers are coarse on purpose: the distance between an approved drug and a
speculative chemotype is what matters, not a second decimal place. They are
still curation judgements, and the weakest entries say so in their notes --
`gpr17_indole_acid`, `sirt_polyphenol`, `glut_pyridine_metabolic`, and
`trophic_aminopyrazole` are all marked speculative.

What this does **not** do is close the loop. There is no experimental feedback:
nothing in the pipeline updates a confidence from an assay result, because no
assay has been run. `core/design/evidence.py` carries the hook for it (§8).

## 7. What the pipeline found about itself

These are recorded because they are the useful part. Each was found by
machinery -- a test, a cross-validation harness, a sweep -- rather than by
inspection, which is the argument for building it.

**Coverage-only optimisation designs molecules that cannot reach the brain.**
The first full run selected a statin acid plus a sulfonylurea — excellent
coverage, TPSA 170–250 against an 40–80 window — and *every* top design failed
the CNS gate. **Penalising polarity alone just moves the failure**: the
optimiser swapped to an adamantane and a chloroarene and landed at cLogP above
6, failing the same gate from the other side. The envelope has to enter the
Hamiltonian on all of polarity, donors, and lipophilicity at once.

**Under a CNS envelope, three arms do not fit.** With the budget correctly
reserving mass for the scaffold and linkers, every `k=3` selection scores
negative while `k=2` lands at TPSA 62, cLogP 2.5, 19 heavy atoms. A
CNS-penetrant three-mechanism single molecule is not obviously reachable from
this library, and the tool says so rather than emitting 550 Da molecules.

**The Hamiltonian silently solved a different problem than the one reported.**
Quadratic couplings were stored keyed on the fragments' *ranked* order and
looked up in *sorted* order, so every coupling whose ranking disagreed with the
alphabet read back as zero. Nothing about the output looked wrong. It surfaced
only because a test asserted the truncated objective agrees with an
independently recomputed exact objective at `k = 2`. The lookup now raises on a
missing pair rather than defaulting to zero.

**Three chemistry bugs that composition checks could not see** -- rotatable
bonds wrong on essentially every generated structure, aromatic carbonyls
mistyped in cLogP by three log units, and a phantom macrocycle in adamantane.
All three were found by cross-validating against RDKit (§5), and none would
have been caught by the formula-and-mass tests that preceded it.

**A fixed QAOA depth was luck, not a default** (§4).

## 8. Reading a design — worked examples

### Multiple sclerosis

```
N1(CCN(CC1)C(=O)Nc2cccnc2)C3CCN(CC3)C(=O)C=C
C18H25N5O2   MW 343.4   cLogP 1.41   TPSA 68.8   CNS-MPO 5.33/6
arms: btk_acrylamide_piperidine (clinical_candidate) + csf1r_picolinamide (published_chemotype)
axes: cns_innate + immunomodulation
mean engagement confidence 0.40, weakest claim AIF1
nearest known compound: Tanimoto 0.36
```

A **dual BTK / CSF1R inhibitor directed at compartmentalised CNS
inflammation.** Every approved MS therapy acts largely in the periphery and
plateaus in progressive disease, so the registry sets CNS exposure as a gate;
this molecule clears it at 5.33 of 6. Both arms reach resident CNS immunity
from different directions -- BTK covers B-cell receptor signalling *and*
microglia, CSF1R covers microglial survival and proliferation -- which is why
it serves two therapeutic axes rather than one.

Predicted effects: BTK −0.85, CSF1R −0.75, TYROBP −0.58, AIF1 −0.55,
CD68 −0.50, CD79A/B −0.45, CXCL13 −0.40. No counter-therapeutic movement.

### Parkinson's disease

```
c1(ccc(cc1)CC(C)N(C)CC#C)OC(=O)C(=O)NCC2CC2
C19H24N2O3   MW 328.4   cLogP 1.61   TPSA 58.6   CNS-MPO 4.93/6
arms: caspase1_ketoamide (published_chemotype) + maob_propargylamine (approved_drug)
axes: symptomatic_dopaminergic + neuroinflammation_control + mitochondrial_rescue
mean engagement confidence 0.47, weakest claim IL1B
nearest known compound: Selegiline, Tanimoto 0.47
```

A **MAO-B inhibitor fused to a caspase-1 warhead**: symptomatic dopaminergic
benefit on one arm, inflammasome-driven neuroinflammation on the other.
Predicted effects: MAOB −0.90, CASP1 −0.80, IL1B −0.60, CASP3 −0.35,
AIF1 −0.30, with BCL2 +0.25 and SOD2 +0.20.

### The novelty figure this design used to carry was wrong (2026-08-28)

The line above previously read `nearest known compound: Tanimoto 0.32`, and the
nearest neighbour was `Propargylamine MAO-B warhead` — a fragment from this
pipeline's own library, not a drug.

The cause: novelty is assessed against
`disease.known_structures() + library.parent_structures()`, and the Parkinson's
registry entry declared **no `structures` file at all**. `known_structures()`
returns an empty list when the path is absent, so the whole first term silently
vanished and every Parkinson's design was compared only against the fragments it
was assembled from. Nothing in the output indicated a reduced comparison — the
field was populated, the number was plausible, and it was measured against the
wrong set.

With `data/chemistry/pd_known_structures.json` registered (26 structures,
covering the small-molecule members of the panel), the same molecule's nearest
neighbour is **selegiline**, an approved Parkinson's drug, at **Tanimoto 0.47**.
It still clears the declared novelty threshold of 0.6, so the conclusion does
not reverse — but the margin is materially narrower than reported, and the
comparison is now against real chemical matter. A propargylamine fused to a
ketoamide *should* look like selegiline; the previous number said otherwise only
because selegiline was not in the room.

Two guards now exist so a third disease cannot repeat this:
`test_each_disease_declares_a_novelty_reference_set` fails if a registered
disease has no reference set or one that barely overlaps its own panel, and
`test_novelty_is_measured_against_real_drugs_not_only_fragments` fails if the
set adds nothing the fragment library did not already contain. Structure
validation in `tests/test_chemistry.py` now iterates every registered disease's
file rather than the MS one by name, so a new set is checked against its
literature formula and mass automatically.

**MS was unaffected** — its reference set has been registered since the campaign
was built. Its designs' nearest neighbour is still a library fragment, which is
the honest answer there: no approved MS agent resembles a dual BTK/CSF1R
covalent.

### What their efficacy is

**Unknown, and untested even in silico for binding.** No docking score, no
assay, no animal model. Neither molecule has ever existed. Every number above
is a specification match or a computed property, and none is a measurement of
effect.

Both designs rank **first of 75 and first of 36** against their disease panels
on signed signature reversal (0.0573 and 0.0550, against approved-agent
medians of 0.0301 and 0.0104). **That comparison is close to circular** and is
reported only to explain why it should not be leaned on: the arms were
selected by the optimiser to maximise alignment with the very signature they
are then scored against, while the panel agents were not. It is evidence the
optimiser worked, not evidence the molecule does.

Since neither is a drug candidate, saying so once in a protocol is not enough.
Every candidate in every result file carries a machine-readable
`evidence_status` block enumerating the **13 assessment axes on which nothing
has been done** -- target binding, selectivity, cell activity, in-vivo
efficacy, permeability, BBB transport, P-glycoprotein efflux, metabolic
stability, CYP liability, hERG, solubility, synthetic route, and freedom to
operate -- each with what would satisfy it and why its absence matters. The
readiness field reads `hypothesis_only`. A consumer that never opens this
document still receives the caveat with the data.

### The arm set that wins is not the arm set that builds the best molecule

This held for both diseases, and it is easy to quote the wrong number:

| disease | Hamiltonian rank 1 | best assembled fitness came from |
| --- | --- | --- |
| multiple sclerosis | BTK + RORgt (+0.0443) | BTK + CSF1R, **rank 3** (+0.0355) |
| Parkinson's | GLUT + MAO-B (+0.0830) | caspase-1 + MAO-B, **rank 2** (+0.0732) |

The two stages optimise different things. The QUBO scores confidence-weighted
coverage under a linearised property envelope; design fitness adds
developability, tractability, and structural alerts to a molecule that
actually exists. An arm that is excellent on coverage can carry mass or a
reactive warhead that only shows up once assembled.

This is why the campaign carries several arm sets forward rather than the
optimum alone, and why "the Hamiltonian's optimum" and "the best design" are
reported as separate results throughout.

## 9. Result files are run artifacts, not curated data

Everything under `experiments/*/results/` is a **snapshot of one run**. It is
regenerated by re-running the recorded command, is safe to delete, and is
git-ignored. Everything under `data/` is curated, reviewed, and the thing a
reviewer is expected to argue with. Confusing the two is how a number outlives
the conditions that produced it.

Every result file carries a provenance block recording:

* the **SHA-256 digest of every input** -- signature, interactome, panel,
  druggability annotation, fragment library, registry entry. A version string
  in a metadata block does not survive someone editing a CSV; a digest does.
* the **git revision and whether the working tree was dirty**. A revision hash
  recorded from a modified tree describes code that was never committed, which
  is worse than recording nothing, so the flag travels with the hash.
* the **interpreter, platform, and tracked package versions**.
* the **active chemistry backend and RDKit version**, because that choice
  changes descriptor values and therefore rankings (§5).

A number from one of these files is not interpretable without the digests that
accompany it.

## 10. What would have to be true before any of this matters

In rough order of how much each would change the conclusions:

1. **Binding.** Dock every proposal, then assay it. Fragment-inherited
   engagement is an assumption, not a prediction, and it is the first thing
   that will break.
2. **Synthesis.** The tractability score is a declared proxy, not a route. A
   chemist has to say whether these can be made.
3. **Selectivity.** `kinase_aminopyrimidine_hinge` is deliberately promiscuous
   across the kinome. A designed multi-target ligand and an uncontrolled
   polypharmacology liability are the same molecule seen from two sides.
4. **The signatures.** The MS signature is illustrative; the Parkinson's one
   is weaker still, curated from literature knowledge rather than derived from
   a cohort. Both need replacing with versioned, batch-corrected,
   phenotype-stratified human datasets.
5. **Selection bias in the library.** Assembly can only reach chemical matter
   the library contains. `unreachable_requirements()` reports profile targets
   with no chemistry behind them, which is the most useful thing a run outputs.
6. **ADMET.** Metabolic stability, transporter efflux, hERG, CYP liability,
   and solubility are entirely unmodelled. P-glycoprotein efflux alone decides
   many CNS programmes and appears nowhere here. The `evidence_status` block
   on every candidate enumerates all thirteen unassessed axes.
7. **An experimental feedback loop.** Fragment confidences are curated priors
   that nothing updates. The value of per-claim provenance is only realised
   when assay results flow back into it; `core/design/evidence.py` carries the
   hook and nothing populates it.

## 11. Reproducibility

```bash
python -m experiments.design.run_denovo_design --disease multiple_sclerosis
python -m experiments.design.run_denovo_design --disease parkinsons --quantum-benchmark
pytest -q tests/
```

Every run is seeded and writes a full provenance block (§9). The quantum
backends need the Qiskit stack; without it the enumeration backend solves the
same problem and the run reports the skip.

To reproduce the chemistry cross-validation, in an environment with RDKit:

```bash
python -m tools.validate_chemistry --designs 300
```

The test suite pins the chemistry backend to `local` so assertions do not
depend on whether RDKit happens to be installed;
`tests/test_chemistry_backend.py` covers the RDKit path and skips itself
where RDKit is absent. Both environments are expected to pass.

To add a third disease, see the registry recipe in the README. It is a data
task: no scoring code changes.
