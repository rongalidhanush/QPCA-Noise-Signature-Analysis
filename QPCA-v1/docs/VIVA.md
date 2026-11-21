# Viva Preparation

## 1. What is your project?

My project is **N-QPCA Based Quantum Noise Signature Analysis**. We use a 15-dimensional quantum-noise response dataset, convert its covariance structure into a density matrix, and use quantum phase estimation as the spectral core of QPCA to estimate the dominant eigenvalues.

## 2. Where exactly is QPCA used?

The main QPCA implementation is:

`src/n_qpca_qpe.m`

It receives the full 16×16 density matrix, constructs its phase unitary, prepares a density-matrix purification, applies controlled powers, performs inverse QFT, and extracts eigenvalue estimates from the phase distribution.

## 3. Why is the matrix 16×16 if the dataset has 15 features?

Because four qubits represent:

2^4 = 16

dimensions. The original 15×15 matrix is embedded into the first 15 basis states and the sixteenth state is zero padded.

## 4. Are you reducing the original data to four features?

No. The main QPCA calculation uses the complete 15-feature covariance matrix. Four refers to the number of system qubits, not the number of features.

## 5. Why convert covariance to a density matrix?

QPCA operates on a density-matrix representation. Dividing the positive-semidefinite covariance matrix by its trace gives a trace-one matrix.

## 6. Why do you use eig()?

`eig()` is not being presented as QPCA. It is a classical numerical reference used to evaluate how accurately the QPCA/QPE simulation estimates the spectrum.

## 7. What is the QPCA output?

The main output is the estimated eigenvalue spectrum. The largest eigenvalues correspond to the dominant principal spectral components.

## 8. What does QPE do?

Quantum phase estimation estimates the phase of an eigenvalue of a unitary. We construct the unitary from rho so that its eigenphases encode the eigenvalues of rho.

## 9. Why use a purification?

A purification converts the mixed density matrix into a pure bipartite state whose Schmidt weights are the square roots of the density-matrix eigenvalues. It lets the QPE simulation access the spectrum with the correct density-matrix weighting.

## 10. How many qubits?

Four system qubits represent the 16-dimensional feature-space embedding. Eight ancilla qubits provide phase resolution of 1/256.

## 11. Is this run on a real quantum computer?

No. It is a state-vector simulation. This was chosen so the complete QPCA algorithm can be inspected and evaluated without claiming hardware execution.

## 12. What is your contribution?

The contribution is the engineering pipeline that applies N-dimensional QPCA/QPE to controlled quantum-noise response data and evaluates the resulting signatures quantitatively.

## 13. Do you claim quantum speedup?

No. The project does not claim quantum speedup because the state preparation and matrix construction are performed classically in the simulator.

## 14. What are your evaluation metrics?

- QPCA/classical eigenvalue error
- phase-resolution error
- top-component explained variance
- reconstruction error
- noise-type nearest-centroid accuracy
- noise-strength correlation
- per-noise QPCA signatures

## 15. What is the limitation?

The main limitation is that the quantum circuit is simulated classically and ideal state preparation is assumed. A hardware implementation would require efficient state preparation and physically realizing the density-matrix evolution/phase operation.
