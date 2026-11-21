# N-QPCA Quantum Noise Signature Analysis

**Dhanush Naidu · Amrita Vishwa Vidyapeetham**

| Release | Contents |
|---|---|
| [Version 1.0.0](https://github.com/rongalidhanush/QPCA-Noise-Signature-Analysis/releases/tag/v1.0.0) | Original published simulation, spectral estimates and noise-signature analysis. |
| [Version 2.0.0](https://github.com/rongalidhanush/QPCA-Noise-Signature-Analysis/releases/tag/v2.0.0) | Final analysis with covariance and correlation heatmaps, a top-ten component plot and reconstruction errors by component count. |

Version 2 uses classical eigenvectors for the reconstruction benchmark and simulated QPE for spectral estimation. With three components the reconstruction retains 95.91% of variance and has 20.22% relative Frobenius error. With ten components those values are 98.50% and 12.25%. [See the measured comparison](results/nqpca_reconstruction_summary_pc3_pc10.csv).

## Project overview

This version makes **N-dimensional QPCA the central algorithm**, not a side demonstration.

The supplied dataset contains:

- 3000 observations
- 15 measured response features
- 5 controlled noise conditions
- a noise-strength field

The complete pipeline is:

```text
3000 × 15 noise-response data
          ↓
standardization
          ↓
15 × 15 covariance matrix
          ↓
ρ15 = C / Tr(C)
          ↓
zero-pad to 16 × 16
          ↓
4-qubit quantum system
          ↓
canonical purification of ρ
          ↓
Quantum Phase Estimation
          ↓
N-QPCA eigenvalue spectrum
          ↓
principal-component spectrum
          ↓
noise-type and noise-strength evaluation
```

The **15-dimensional problem is retained**. It is not reduced to four hand-crafted features.

Because a quantum register represents dimensions `2^n`, the 15-dimensional density matrix is embedded into a 16-dimensional Hilbert space. The sixteenth basis state is a zero-padding state and is not a measured feature.

## The central QPCA step

The covariance-derived density matrix is:

\[
\rho = \frac{C}{\operatorname{Tr}(C)}
\]

The simulator constructs the phase unitary

\[
U = e^{2\pi i\rho_{16}}
\]

and applies QPE to a purification of the density matrix. The purification is constructed from `sqrt(rho)` so that the probability of observing an eigencomponent follows its density-matrix eigenvalue.

QPE estimates phases

\[
\phi_j \approx \lambda_j
\]

where the eigenvalues `lambda_j` are the principal spectral components of the normalized covariance matrix.

This is a **state-vector simulation of QPCA/QPE**, not a claim of execution on quantum hardware.

## Run

From the project root:

```matlab
addpath("src")
main_NQPCA
```

The script creates:

- `results/nqpca_main_results.mat`
- `results/nqpca_global_spectrum.csv`
- `results/nqpca_noise_type_results.csv`
- `results/nqpca_noise_strength_results.csv`
- `results/nqpca_qpe_distribution.csv`
- plots in `results/figures/`

## Main evaluation

The project evaluates:

1. QPCA eigenvalue agreement with the classical reference.
2. Top-3 explained variance.
3. QPCA reconstruction error.
4. Noise-type separation using QPCA scores.
5. Noise-strength sensitivity using QPCA PC1.
6. QPE probability normalization.
7. Phase-resolution error.
8. Per-noise QPCA signatures.

## Original implementation

The MATLAB implementation in `src/` was written specifically for this project. The referenced papers are used for the mathematical concepts and scientific background, not as copied source code.

## Limitation

The simulation assumes ideal state preparation and directly constructs the phase unitary from the density matrix. This is appropriate for validating the QPCA/QPE algorithmic pipeline, but it does not demonstrate hardware speedup or the full physical cost of preparing/evolving a large density matrix.
