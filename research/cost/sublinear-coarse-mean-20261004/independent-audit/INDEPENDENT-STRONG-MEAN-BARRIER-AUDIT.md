# Independent audit: strong posterior-mean oracle barrier

2026-10-04. Verdict: **accepted within the stated uniform black-box, strong-mean scope. No fatal mathematical gap found.** The exponent consequence is `c(P) >= (P-1)/3`. It is not a lower bound for a single fixed potential, arbitrary posterior sampling, or a final buffered-law comparison.

Audited barrier SHA256: `518c1eb8a623acb0268873cd8579643b60b13dd4d98d37f22f2bfe233911ddc8`.

Revision verification, 2026-10-04 19:36 UTC: accepted both scoped Section 5 edits, making the endpoint-node realization and the potential-uniform small-heat threshold explicit. Re-read the complete revised argument and found no change to the mathematical verdict; the numerical checks remain applicable. The original audit pin was `f8f4fef16c07a90796dc86aa7908dfcfa6eca6e11354146736f0fd934a93e472`. A reverse-edit byte-hash reconstruction did not reproduce that pin, so this verification rests on reviewing the current complete text rather than claiming a verified two-edit-only byte diff.

Audited baseline SHA256: `f7514fbaf777a573d262416759f2776dcce7bbc909741ceeb08a6276671785f4`.

Also inspected the complete-family numerical qualifiers, the earlier coarse-bank gate, the linear-cost map, and LOW30's literal collocation rows and hidden-source contract. The baseline is `../growing-order-outer-20261004/direct-mean/ANY-ORDER-DIRECT-MEAN-OUTER.md` relative to the coarse-mean directory.

## 1. Exact target and normalization

With `U=X/sqrt(A)`, the standardized density is proportional to `exp[-u^2/2-V(sqrt(A)u)]`. Integrating its derivative gives

    E U + sqrt(A) E V'(X) = 0.

Hence `mu_A=-E U` exactly. The tail boundary term vanishes because the perturbation is bounded and the quadratic precision is in `[1,3/2]`.

Write `R=eta A n^-2 sum_i theta_i Phi_i`. Its absolute value is at most `eta A B/n`; both the unnormalized weight and the normalizer change. Their ratio yields

    exp(-2 eta A B/n) <= p_theta/q_A <= exp(2 eta A B/n).

Differentiation of the normalized expectation therefore gives

    d mu_A/d theta_i = eta A n^-2 Cov_theta(U,Phi_i(U)).

This includes the normalizing constant. Differentiation under the integral is justified uniformly by a Gaussian integrable envelope, since every `Phi_i` is bounded by `B`.

The independent-copy covariance proof is correct. Both cross-events, `U in [2,5/2], U' in [-1,0]` and its reversal, contribute; together they give at least `2 B P_+ P_-`. All remaining contributions are nonnegative.

A completely explicit uniform scalar bound is available. For `s in [1,3/2]`, the centered Gaussian density is at least `exp(-75/16)/sqrt(2 pi)` on `[2,5/2]` and `exp(-3/4)/sqrt(2 pi)` on `[-1,0]`. The density ratio is at least `exp(-1/120)`, since `eta=1/4`, `B=1/30`, and `n>=2`. Thus one may take

    c0 = exp(-87/16 - 1/60)/(60 pi) > 0,
    k = eta c0,
    d mu_A/d theta_i >= k A/n^2.

The constants are intentionally loose. Crucially, they are independent of `A,n,i`, and the interpolated sign vector.

## 2. Regularity, side information, and oracle locality

`phi=t^2(1-t)^2` on `[0,1]`, zero elsewhere, is `C1`; both endpoint values and endpoint first derivatives vanish. Consequently `Phi` is `C2`, and the displayed potential is globally `C2`. The exact maximum `|phi'|` is `sqrt(3)/9 < 1`, so the claimed Hessian sandwich follows with room to spare. Only one open cell can contribute to a gradient or Hessian at a point.

Outside the cells, the original gradient is exactly `lambda x`. The cumulative plateaus in the potential itself do not appear in that gradient. Inside a cell, a gradient value or an original Hessian/first/adjoint action depends on at most that cell's sign. Allowing the stronger oracle to identify this sign perfectly cannot make the lower bound harder to prove.

At the caller zero the potential, gradient and Hessian are the same for every sign vector. Both the potential mode and the posterior mode at caller zero are exactly zero by strict convexity. Free knowledge of these objects leaks no hidden sign. Exact private random-number generation also leaks none. Any additional potential-dependent preprocessing must be included in the charged queries or independently shown sign-independent.

The exclusion of a potential-value oracle is material: potential values contain accumulated preceding-cell information. Likewise an expectation, posterior-sample, normalization, or general smoothed-score oracle is not covered. This is an oracle-information theorem, so arbitrary computation means computation independent of the unknown potential except through permitted oracle responses and the declared side information.

## 3. Adaptive randomness and stopping

Condition first on the entire private random seed. An adaptive point-query program is then a decision tree whose next point and stopping decision use only previously revealed signs. Each leaf determines the signs revealed along that branch and imposes no restriction on unrevealed signs. Thus the latter remain independent uniform Rademachers, including when stopping is random or depends on previously observed signs.

Given a transcript and unknown index set `I`, each first Walsh coefficient satisfies

    b_i = E[f(theta_I) theta_i | transcript]
        = (1/2) E[f(theta_i=+1)-f(theta_i=-1) | transcript]
        >= k A/n^2.

Conditional orthogonality, not an unjustified sum of arbitrary influences, supplies

    Var(mu_A | transcript) >= sum_i b_i^2
                              >= (n-m) k^2 A^2/n^4.

The algorithm's conditional squared loss is no smaller than that variance, even when its output is the optimal posterior mean. Averaging over signs and private randomness then lower-bounds its worst-case risk.

For a deterministic cap, `n=2 max(1,N)` gives the displayed `N^-3` rate (with the usual harmless treatment of `N=0`). For a finite uniform worst-case expected query budget, choose `n=2 max(1,ceil(N))`; then `E m <= E Q <= N <= n/2`, so the same bound follows after expectation. A deterministic cap is not needed. The expectation must include all algorithm randomness and apply to every potential in the class, or at least to the finite prior used in this proof. Merely knowing the work on a favorable potential would not suffice.

## 4. Link to the actual baseline mean port

The baseline requires each complete executable force statistic to be strongly close to its deterministic conditional exact posterior force. Its finite-source bias is `O(A^(j+1))`; if its centered/finite-source strong error is `O(A^j)`, triangle inequality gives the exact-target strong error `O(A^j)`. No access to that finite-source expectation is assumed.

There is a useful explicit way to realize caller zero inside the actual hidden-node graph. LOW30 defines

    A_iℓ = integral_0^(c_i) (c_i-u) ell_ℓ(u) du.

At the Chebyshev–Lobatto endpoint `c_i=0`, every `A_iℓ=0`. Therefore `y_i^[d]=x` at every Picard depth, independently of the captured momentum refresh. Setting the original position `x=0` realizes exactly the hard posterior used here. One must not assert that all other hidden Picard targets are zero; they generally are not. The verified Section 5 revision now includes this endpoint observation and makes the baseline application explicit.

All original-gradient work used to prepare that force statistic, including other-node banks if consulted, belongs in the total transcript budget. Shared or adaptive work does not defeat the bound: the complete output component still has only the total number of queried sign cells available to it.

The separately selected numerical allowance must actually be chosen at `O(A^j)` or below to claim this strong target. A larger fixed numerical tolerance would no longer instantiate the target being lower-bounded. At caller zero the relevant original caller profile is uniformly bounded, so unbounded-center numerical envelopes are not a loophole.

## 5. Uniform scope and exponent

The finite hard family depends on both `A` and `n`. This is legitimate for a supremum over the Hessian-sandwich class at each heat, but does not establish a single potential that remains hard at every heat. Uniform error/work constants and a uniform small-heat threshold over this class are essential. Potential-specific thresholds or uncontrolled potential-specific constants are outside the argument.

From the mean-square lower bound,

    C_j A^j >= c A/(N(A)+1)^(3/2).

Thus for `j>1`, `N(A)+1 >= c_j A^[-2(j-1)/3]`. A uniform upper bound `N(A) <= C_j log(1/A)^d_j A^-c_j` is incompatible with `c_j < 2(j-1)/3`, because a fixed logarithmic power cannot compensate for a nonzero power of `A`. With `P=2j-1`, this is exactly

    c(P) >= (P-1)/3.

Constants and logarithmic degrees may grow arbitrarily with the fixed requested order. The conclusion still prevents `c(P)/P -> 0` for a family retaining this same uniform strong-mean contract. It leaves open different contracts, especially a comparison of the fully buffered output law. It also does not prove that the necessary slope is attainable.

## 6. Noncommuting two-dimensional fixture

`T` has eigenvalues `7/16` and `9/16`. Adding `eta theta_i phi' E11` with `eta=1/8` places the Hessian spectrum inside `[5/16,11/16]`, hence inside the requested sandwich. In any one sign configuration, two points with distinct bump derivative values give distinct coefficients of `E11`, and

    [T+a E11,T+b E11] = (b-a)[T,E11] != 0.

This is noncommutation within one potential, not merely between different members of the family.

The Schur complement obtained by integrating out `U2` is exactly

    s_A = 1+A/2 - (A/16)^2/(1+A/2).

It lies in `[1,3/2]` for `0<A<=1`. The remaining one-dimensional perturbation is unchanged. Vector integration by parts gives the first force component `-E U1`, so the same sensitivity and information bound apply. A full gradient or HVP at any point depends on only one sign through its first coordinate; the off-diagonal coupling cannot reveal another cell's sign.

## 7. Independent finite checks

`check_strong_mean.py` uses normalized, split-interval quadrature; it checks scalar and two-dimensional cases, the exact integration-by-parts identity, normalization-sensitive finite differences, every partial sign transcript for `n=4`, and a same-potential nonzero Hessian commutator. Results are in `numerical-check-results.json`.

- 54 scalar/2D posterior cases and 81 partial transcripts passed.
- Maximum integration-by-parts absolute error: `1.36e-16`.
- Maximum relative finite-difference derivative discrepancy: `1.23e-8`.
- Minimum observed covariance: `7.81e-4`, above the explicit uniform bound `2.27e-5`.
- Minimum projection slack was `-2.07e-25`, floating-point roundoff.
- Checked noncommuting Hessian eigenvalues were between `0.4243` and `0.5757`; commutator norm was `0.00425`.

These checks corroborate formulas and implementation only. The lower bound is established by the analytical arguments above, not by those numerical tests.
