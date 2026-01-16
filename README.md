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

The system is designed for **reproducibility, extensibility, and future quantum hardware deployment**.

---

## Conceptual Pipeline
```
Disease / Protein Set
↓
Protein Interaction Network
↓
Drug Perturbation Simulation
↓
Synergy & Probability Estimation
↓
Hamiltonian (QUBO) Formulation
↓
Optimization (Exact / Heuristic / Quantum-Ready)
↓
Scaling & Hardness Benchmarking
```

## Repository Structure
```
quantum-bio-recovery-engine/
├── core/ # Biology, optimization, and Hamiltonian modules
├── experiments/
│ ├── ms/ # Disease-specific experiments (e.g., Multiple Sclerosis)
│ └── benchmarks/ # Scaling and hardness benchmarks
├── data/ # Synthetic and experimental datasets
├── figures/ # Generated plots and benchmark figures
├── lab_journal/ # Experimental logs and daily research notes
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

## Scientific Motivation

Therapeutic combination discovery is inherently combinatorial. As the number of candidate drugs increases, exact classical optimization becomes intractable, while heuristic methods sacrifice solution quality. This platform provides a benchmarked, Hamiltonian-based formulation suitable for evaluating hybrid and quantum-assisted optimization strategies in biological and pharmacological research.

## Disclaimer

This project is a research and educational platform only.
It does not provide medical advice, clinical recommendations, or validated treatment guidance.
