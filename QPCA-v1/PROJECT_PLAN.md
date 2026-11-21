# Project Plan

## Title

**N-QPCA Based Quantum Noise Signature Analysis**

## Problem statement

Quantum processors and quantum-inspired sensing systems produce multi-dimensional response data under different noise mechanisms. Direct inspection of many correlated response measurements makes it difficult to identify dominant variation and compare noise conditions.

A compact spectral representation is therefore needed to identify dominant directions in the measured noise-response space.

## Aim

To implement a genuinely N-dimensional QPCA-based analysis pipeline for a 15-dimensional quantum-noise response dataset and evaluate whether the resulting quantum principal components capture dominant variation and distinguish controlled noise conditions.

## Objectives

1. Load and validate the 3000 × 15 dataset.
2. Standardize the 15 response variables.
3. Construct the 15 × 15 covariance matrix.
4. Normalize the covariance matrix into a density matrix.
5. Embed the 15-dimensional matrix into a 16-dimensional Hilbert space.
6. Use four system qubits to represent the 16-dimensional state space.
7. Construct a canonical purification of the density matrix.
8. Apply quantum phase estimation to obtain QPCA eigenvalue estimates.
9. Extract the dominant QPCA spectral components.
10. Compare QPCA estimates with a classical numerical reference.
11. Reconstruct the data using the dominant components.
12. Evaluate reconstruction error.
13. Analyze noise-type signatures.
14. Analyze noise-strength sensitivity.
15. Produce figures, CSV results, and a saved MATLAB workspace.

## Research contribution

The contribution is an engineering framework that combines:

- controlled quantum-noise data generation/analysis,
- density-matrix construction,
- N-dimensional QPCA,
- QPE-based spectral estimation,
- noise signature extraction,
- and quantitative evaluation.

The project does not claim to invent QPCA. The contribution is the application and implementation framework for quantum-noise signature analysis.

## Evaluation criteria

### QPCA correctness

Compare:

- classical eigenvalues
- QPCA/QPE estimated eigenvalues

using absolute error.

### Dimensionality reduction

Measure cumulative explained variance of the leading QPCA components.

### Reconstruction

Project the standardized data into the leading classical reference eigenvectors and reconstruct it. Report normalized reconstruction error.

### Noise discrimination

Use the first few principal scores to calculate nearest-centroid classification accuracy as an evaluation of the extracted signature space.

### Noise-strength sensitivity

Calculate the correlation between noise strength and the magnitude of the first QPCA component.

## Key limitation

The implementation is a state-vector simulation. It demonstrates the QPCA/QPE computation but not a physical quantum-hardware speedup.
