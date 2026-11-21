# First-run fixes

Two issues found during the first MATLAB run were corrected.

1. The 16-dimensional QPCA eigenvectors include a zero-padded 16th basis
   state. The original data has 15 features, so score projection now uses
   `V(1:15,1:nPC)`.

2. Finite-resolution QPE peaks do not necessarily have heights ordered by
   eigenvalue because of phase leakage. The selected QPE peak locations are
   therefore explicitly sorted by estimated eigenvalue before reporting PC1,
   PC2, PC3.

The first run already demonstrated that the QPCA/QPE stage executed:
- 3000 rows
- 15 features
- 16-dimensional quantum embedding
- 4 system qubits
- 8 QPE ancillas
- QPE probability sum = 1
- PC1 estimate 0.39453125 versus classical 0.39504677

Run `main_NQPCA` again after these fixes.
