# Integration audit and exact claim boundary

2026-10-05. Analytical checks only. No numerical OU path or history is evaluated.

## Independent work and checks

The finite six-VALUE source and its Gaussianization were independently audited in independent-audit/INDEPENDENT-PAIR-STEIN-AUDIT.md. That audit separately rederives the factor 1/2 in the double difference, the PSD-interval difference bound, the small conditional-inner derivative, the OU resolvent Stein matrix, the K-copy HS variance, and the anisotropic W2 constant. It specifically flags that every source root and the smoothing Gaussian are consumed.

The finite-order F2 theorem was derived independently of the six-VALUE source. The integration author then checked its entropy/conditioning argument, mollification stability, exact pairwise first-order representation, conditional-L2 rather than L4 projection, physical block summation, crossover exponent and all-scale summation. The full coherent-block extension is proved separately, rather than treating a small mark energy as a small derivative.

The cluster-translation and gap-obstruction note was independently derived from the 2-by-2 OU correlation matrix, with full analytic proofs. Its strong claim is a lower bound for a GENERIC source-independent gap rule over arbitrary bounded positive fields. It is not a lower bound specific to the true history field.

## Noncommuting and normalization checks

1. The source's derivative is (H11-H21)H_z in exactly that order. Neither symmetry of its product nor positivity of its signed derivative is claimed.
2. Four orthogonal copies of the pair-Hoeffding component appear in the rectangular difference; division by 2, rather than 4 or sqrt(2), gives its exact covariance.
3. Covariance control is ell^2 I, ell=cA^2 tau/sqrt(2). The complete first is sqrt(A^2 sigma^2+c^2 A^4 tau^2)/sqrt(2), not ell.
4. The Gaussianized sum is K^(-1/2), preserving the target covariance. Shared retained-caller derivatives therefore grow by sqrt(K); complete private derivatives do not.
5. The Stein matrix can be nonsymmetric. Its operator norm, mean identity and HS variance suffice for the stated W2 estimate; it is not used as a PSD executing matrix.
6. Nuisance averaging shares the four midpoint roots, uses R independent complete nuisance banks, and has a PSD 1/R excess covariance. Treating R replicas as independent midpoint banks would change that target.
7. Positive sums of interacting modes must share the midpoint source before squaring/Gaussianizing. Separately squaring their modes would delete cross-mode terms.
8. The exact ordered leading history pair is A_I B_J, with covariance E[A_I Cov(B_J) A_I*]. No matrix factor is commuted.
9. The block readout [I I] costs 2 on covariance error and sqrt(2) on LAW error; these are different factors.

## Tail theorem checks

For a block of size n, the paid proof-only smoothing gives a comparison with L4 squared norm C A^5 n^2 J_k. Its high-order conditional covariance is controlled by an L2 projection in every physical direction. No L4 operator bound for high-order Hoeffding projections is used.

The fine-scale Poincare cap is C A^2 d_k^2 sqrt(n). Physical block summation gives the displayed b^(3/2) sqrt(D) coarse term. Balancing A^2 d^2 with A^5 b^(3/2) d^(-1/2) gives d_*=A^(6/5)b^(3/5), and both dyadic series have scale A^(22/5)b^(6/5). The coarse d>1 scales cost the separate log_+T term.

The full coherent mark has L4 squared energy C A^6 n, absorbed by C A^5 n^2 J_k for T>=1. Adding it to the STACKED comparison controls every omitted high-order cross block at once. The fine-scale block cap uses the complete O(A) sensitivity kernel, so no false O(A^3) derivative enters.

The result improves A^4 sqrt(D) uniformly only for fixed bounded physical block size, or other explicit small-dimension regimes satisfying the displayed factors. For an arbitrary source, b=D is legitimate but the resulting dimension factor remains in the estimate.

## Finite source, fixed shape, and actual-history source are distinct

The packet proves an executing finite pair interaction and its buffered positive covariance law. It also proves that a fixed analytical order cutoff is adequate at a stated block-dependent error for the true history covariance. These two theorems are NOT silently joined: a finite original-VALUE representation of the exact singleton/pair conditional coefficients, including their entire old skeleton and nonlinear coherent ancestors, is still required.

Likewise, the common-translation theorem accepts one fixed complete cluster field. Its endpoint projection is analytical, not a unit-cost original-VALUE oracle. Translating only the selected pair while leaving correlated old background fixed changes its source/conditioning. Varying the relative gap changes the physical stationary pair law, and the independent gap obstruction proves why the generic dimension-free theorem does not cover that operation.

Accordingly the finite source has the literal 6K (or 6KR with nuisance averaging) bill and the explicit growing caller first. The true-history coefficient count remains unsupplied; literal order-two enumeration is order 4^K. No new logarithmic m4 compiler, marked cubic/quartic cut theorem, retained-mark Gaussianization, same-readable-keep LAW comparison, or unrestricted-dimension accuracy-raising induction is claimed.

## Sealed input check

The four input MANIFEST.json hashes exactly match the pinned requests. The new packet does not modify those inputs. Its mathematical deliverables are text proofs; hashing and file assembly perform no simulation or force-history evaluation. Nothing is uploaded or published externally.
