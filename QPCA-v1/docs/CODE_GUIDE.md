# Code Guide

## main_NQPCA.m

Main orchestration script. It loads the dataset, creates the density matrix, runs N-QPCA, evaluates it, analyzes noise types/strengths, and saves results.

## n_qpca_qpe.m

Central quantum algorithm.

1. Normalize rho.
2. Construct U = exp(2πi rho).
3. Build sqrt(rho).
4. Prepare the canonical purification.
5. Apply controlled U powers.
6. Apply inverse QFT.
7. Measure the ancilla distribution.
8. Extract dominant eigenvalue phases.

## inverse_qft_matrix.m

Builds the inverse quantum Fourier transform matrix.

## select_qpca_peaks.m

Finds separated peaks in the QPE distribution so adjacent finite-resolution bins are not counted as separate principal components.

## classical_reference.m

Provides a classical numerical reference only.

## validate_nqpca.m

Calculates QPCA/classical eigenvalue errors.

## density_matrix_exponentiation_demo.m

Supplementary demonstration of the density-matrix-exponentiation idea using a SWAP interaction.

## run_physical_noise_channels.m

Supplementary channel sanity check. It is not the main QPCA algorithm.

## data/noise_dataset.csv

The supplied 3000 × 17 dataset used by the project.
