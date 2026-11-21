#  NQPCA: Quantum Noise Signature Project

A MATLAB simulation and benchmarking framework for **Noisy Quantum Principal Component Analysis (NQPCA)**. This project investigates how quantum noise affects density matrix exponentiation and Quantum Phase Estimation (QPE), enabling the extraction and validation of **quantum noise signatures** across physical noise channels.

---

## 📖 Table of Contents
- [Overview & Motivation](#-overview--motivation)
- [How the Project Works (Pipeline)](#-how-the-project-works-pipeline)
- [Repository Architecture](#-repository-architecture)
- [Getting Started](#-getting-started)
- [Validation & Reference Benchmarks](#-validation--reference-benchmarks)
- [Documentation Guide](#-documentation-guide)

---

## 🔬 Overview & Motivation

**Quantum Principal Component Analysis (qPCA)** is a foundational Quantum Machine Learning (QML) algorithm used to extract the dominant eigenvalues and eigenvectors of an unknown quantum state (represented as a density matrix \(\rho\)) exponentially faster than classical methods.

However, near-term quantum hardware (NISQ devices) suffers from environmental decoherence and gate errors. Instead of treating noise purely as an obstacle, **NQPCA (Noisy Quantum PCA)** analyzes how physical noise alters the spectral components of a density matrix:

1. **Noise Signature Extraction:** Identifies distinct eigenvalue shifts and phase peaks introduced by specific quantum noise channels.
2. **Quantum vs. Classical Benchmarking:** Evaluates how closely noisy QPE approximates ideal classical eigendecomposition.
3. **Algorithmic Validation:** Tests density matrix exponentiation (\(e^{-i\rho t}\)) and Inverse Quantum Fourier Transform (\(QFT^\dagger\)) precision under realistic noise conditions.

---

## ⚙️ How the Project Works (Pipeline)

The codebase follows a modular 5-stage scientific computing workflow:

```text
[noise_dataset.csv] 
        │
        ▼
┌───────────────────────────────┐
│ 1. Physical Noise Channels    │ ──► Simulates decoherence & noise models
└───────────────────────────────┘
        │
        ▼
┌───────────────────────────────┐
│ 2. Density Matrix Exp. Demo   │ ──► Approximates unitary evolution e^(-iρt)
└───────────────────────────────┘
        │
        ▼
┌───────────────────────────────┐
│ 3. NQPCA & QPE Execution      │ ──► Applies Inverse QFT to estimate phase peaks
└───────────────────────────────┘
        │
        ▼
┌───────────────────────────────┐
│ 4. Peak Selection             │ ──► Extracts dominant principal components
└───────────────────────────────┘
        │
        ▼
┌───────────────────────────────┐
│ 5. Classical Validation       │ ──► Compares against reference_classical_qpca.csv
└───────────────────────────────┘
```

1. **Data Ingestion & Channel Simulation (`run_physical_noise_channels.m`):** Loads quantum noise features from `data/noise_dataset.csv` and applies physical noise channel transformations.
2. **Density Matrix Exponentiation (`density_matrix_exponentiation_demo.m`):** Converts the density matrix \(\rho\) into a unitary operator \(e^{-i\rho t}\) required for phase estimation.
3. **Quantum Phase Estimation (`n_qpca_qpe.m` & `inverse_qft_matrix.m`):** Runs QPE using an explicit Inverse Quantum Fourier Transform matrix to map eigenvalues onto measurement register probabilities.
4. **Spectral Peak Detection (`select_qpca_peaks.m`):** Filters the resulting probability distribution to isolate dominant eigenvalue peaks corresponding to principal noise signatures.
5. **Benchmarking & Validation (`classical_reference.m` & `validate_nqpca.m`):** Computes exact classical eigenvalues and compares NQPCA accuracy against baseline metrics in `results/`.

---

## 📁 Repository Architecture

```text
NQPCA_Quantum_Noise_Signature_Project/
├── data/
│   └── noise_dataset.csv                     # Input dataset containing quantum noise parameters/states
├── src/
│   ├── main_NQPCA.m                          # 🚀 MAIN ENTRY POINT: Orchestrates the full simulation
│   ├── classical_reference.m                 # Computes exact classical PCA/eigendecomposition baselines
│   ├── density_matrix_exponentiation_demo.m  # Simulates e^(-iρt) density matrix exponentiation
│   ├── inverse_qft_matrix.m                  # Constructs the Inverse QFT unitary transformation matrix
│   ├── n_qpca_qpe.m                          # Core Noisy qPCA Quantum Phase Estimation engine
│   ├── run_physical_noise_channels.m         # Models physical quantum noise channels
│   ├── select_qpca_peaks.m                   # Extracts principal eigenvalue peaks from QPE output
│   └── validate_nqpca.m                      # Computes error metrics & validates against reference CSVs
├── results/
│   ├── reference_classical_qpca.csv          # Ground-truth classical qPCA eigenvalues/vectors
│   └── reference_project_metrics.csv         # Expected accuracy, fidelity, and error benchmarks
├── docs/
│   ├── REPORT.md                             # Full technical report & experimental findings
│   ├── CODE_GUIDE.md                         # Function-by-function developer walkthrough
│   ├── EXPECTED_REFERENCE.md                 # Expected console outputs and plots
│   ├── FIX_AFTER_FIRST_RUN.md                # Known post-run adjustments & troubleshooting
│   ├── TOOLBOX_NOTE.md                       # Required MATLAB toolboxes and version notes
│   ├── REFERENCES.md                         # Research papers and theoretical background
│   └── VIVA.md                               # Viva-voce / presentation Q&A study guide
├── PROJECT_MANIFEST.json                     # Structured metadata of project files
├── PROJECT_PLAN.md                           # Project phases and implementation roadmap
├── RUN_PROJECT.txt                           # Quick-start execution commands
└── README.md                                 # Project documentation
```

---

## 🚀 Getting Started

### Prerequisites
* **MATLAB** (R2021a or newer recommended)
* Review `docs/TOOLBOX_NOTE.md` before running to verify if any add-on toolboxes (such as *Statistics and Machine Learning Toolbox* or *MATLAB Support Package for Quantum Computing*) are utilized.

### Execution Steps
1. **Clone the repository:**
   ```bash
   git clone https://github.com<your-username>/NQPCA_Quantum_Noise_Signature_Project.git
   cd NQPCA_Quantum_Noise_Signature_Project
   ```
2. **Open MATLAB** and set the project root (`NQPCA_Quantum_Noise_Signature_Project/`) as your **Current Folder**.
3. **Add `src` to your MATLAB path and run the main script:**
   ```matlab
   addpath('src');
   main_NQPCA
   ```
4. *(Optional)* Check `docs/FIX_AFTER_FIRST_RUN.md` if you need to adjust parameters or file paths after your initial test run.

---

## 📊 Validation & Reference Benchmarks

To verify that your simulation reproduces the expected scientific results, `validate_nqpca.m` automatically evaluates your run against two ground-truth files in the `results/` directory:

* **`results/reference_classical_qpca.csv`:** Verifies that the estimated QPE phase peaks align with ideal classical eigenvalues.
* **`results/reference_project_metrics.csv`:** Checks reconstruction fidelity, peak estimation error, and noise channel classification metrics.

Refer to `docs/EXPECTED_REFERENCE.md` for a detailed breakdown of expected outputs.

---

## 📚 Documentation Guide

This repository comes pre-packaged with extensive documentation inside the `docs/` folder:
* **Presenting or defending this project?** Read `docs/REPORT.md` and `docs/VIVA.md`.
* **Modifying the MATLAB functions?** Read `docs/CODE_GUIDE.md`.
* **Citing academic sources?** Check `docs/REFERENCES.md`.
