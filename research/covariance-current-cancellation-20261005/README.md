# Exact covariance-current cancellation

Read EXACT-COVARIANCE-CURRENT-REDUCTION.md.

Finished analytical result: the existing moving-Hermite/root and covariance sampler currents combine exactly into t J Sigma J^T. The previously paid L^2 h3 term disappears. The same raw positive reference has rate L^min(m,5)e/v^(min(m,5)/2), with explicit finite-rank constants and the original normalized small-radius guard.

A concrete bounded-Hessian gradient example proves the next skew-square variance obstruction is real. The first connected correction is explicitly identified (rank two and rank four); no tensor-oracle or native executability claim is made.

No source files are modified. Native rank-five and complete-old-bank feedback admission remain open. The two-cycle return is not assumed closed. Independent audit and exact polynomial checks are included in this directory.

Reproduce finite diagnostics with `python check_covariance_reduction.py`: 92 exact assertions PASS, including noncommuting covariance/root transport. These are algebra checks, not executions of a native compiler. FINAL-STATUS.md states the precise analytical/native boundary.
