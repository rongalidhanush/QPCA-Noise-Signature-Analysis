# MATLAB toolbox note

The main project intentionally uses base MATLAB constructs so that the QPCA
simulation does not depend on the Statistics and Machine Learning Toolbox.

The plotting uses `scatter` rather than `gscatter`, and correlation is
computed with a small local function.
