# QPCA Noise Signature Analysis

**Dhanush Naidu · Amrita Vishwa Vidyapeetham**

This MATLAB project studies noise signatures across five conditions using a simulation of quantum phase estimation. Each version below contains its own source code, dataset, figures and results.

| Project folder | What to find |
|---|---|
| [QPCA-v1](QPCA-v1/) | Original published project, including the QPE spectrum and per-condition noise signatures. |
| [QPCA-v2-final](QPCA-v2-final/) | Final project with covariance and correlation heatmaps, top-ten component analysis and reconstruction error by component count. |

## Run the final version

Download or clone the repository. In MATLAB, set the Current Folder to the repository root and run:

```matlab
cd('QPCA-v2-final')
addpath('src')
main_NQPCA
```

To run the original, choose the QPCA-v1 folder instead. Each version writes its outputs to its own results folder.

[Download Version 1](https://github.com/rongalidhanush/QPCA-Noise-Signature-Analysis/releases/tag/v1.0.0) · [Download Version 2](https://github.com/rongalidhanush/QPCA-Noise-Signature-Analysis/releases/tag/v2.0.0)
