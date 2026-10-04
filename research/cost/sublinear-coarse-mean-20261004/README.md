# Coarse-mean cost: barrier and scalar law-only escape

2026-10-04. Local research package. No all-order posterior source or eventual-sublinear outer theorem is claimed.

## Results

1. STRONG-POSTERIOR-MEAN-ORACLE-BARRIER.md proves a uniform C2 original-gradient information bound, independently accepted: RMS error for the canonical posterior mean is at least c A/(N+1)^(3/2). Thus retaining the baseline strong O(A^j) mean port requires c(P)>=(P-1)/3, P=2j-1, even with arbitrary control variates, smoothing and adaptive first-action queries. It is not a fixed-potential or law-only lower bound.

2. SCALAR-POSITIVE-HERMITE-MEAN-LAW-FACTORY.md supplies an independently accepted positive scalar buffered-law construction from K complete bounded source records. The conditional density has a known exact normalization and lies between one-half and three-halves of a standard Gaussian. No unknown expectation is queried. A logarithmic rejection cap is explicit.

3. POLYLOG-SCALAR-MEAN-LAW-REFINEMENT.md extends that factory using an executed pilot anchor, fresh independent averages, clipping, and rank-uniform Hermite bounds. Its complete sample census is O((1+(S/sigma)^2) log^3((e+S/sigma)/delta)) for W2 error O(sigma delta). This permits a scalar buffer at the provider's fluctuation scale without an inverse-heat count, provided S/sigma is at most public-logarithmic. Its independent continuation audit in independent-polylog-refinement-audit/ accepts this scoped result.

## Important boundaries

- The lower bound uses a uniform Hessian-sandwich class, including potential-uniform constants and small-heat thresholds. The hard family may vary with A and N. All descendants and queried side information are charged.
- The scalar factories return marginal or explicitly conditioned buffered laws. They do not return strong mean estimates, smooth first/adjoint interfaces, actual adjacent pairs, or dimension-free vector factories.
- Fresh factor banks cannot be retained after the Hermite averaging. The pilot can be retained in the refined conditional construction, because every later factor batch is independent of it.
- Clipping, Gaussian generation and numerical thresholds need their source-specific numerical budgets. Known-operation counts are not unqualified bit-complexity theorems.
- Conditional scalar application to the new Stein packet yields Z-v_Q(Z)+sigma N. Its added sigma^2 covariance is a real outstanding debt. The first-order Stein target is not automatically the posterior.

## Diagnostics

- check_strong_mean_barrier.py: 21,315 author assertions. Exact normalized posterior sensitivities, finite conditional Walsh projections, scalar quadratic cancellation, genuine noncommuting Hessians and exponent arithmetic.
- independent-audit/: 54 independent normalized posterior cases and 81 partial sign transcripts, plus proof review.
- independent-scalar-factory-audit/: 20 initial factory diagnostics and exact proof review; continuation checks/audit are separate.
- independent-polylog-refinement-audit/: 30 parameter-bound checks and 12 quadratic identities, plus the complete continuation proof review. No enormous provider batches were simulated.
- check_polylog_scalar_refinement.py: 577 author assertions across 56 parameter fixtures, plus finite bounded-source clipping checks.

Diagnostics corroborate finite formulas; the proofs, their explicit hypotheses and their independent reviews determine the scope.
