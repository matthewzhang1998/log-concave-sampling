# Audit: local 14-VALUE pair producer

Date: 2026-10-05. Reviewed candidate: `TRUE-PAIR-HISTORY-RESUMMATION.md`, as present at 10:50 UTC. This supplements the transport audit and the pair-history covariance comparison check.

## Verdict

The local Markov join, two-replica nuisance-sharing pattern, seven VALUES per replica, uniform inner Lipschitz estimate, level weights, and stated fixed-block covariance exponents check out. No substantive mathematical gap was found in those parts. The result remains a marginal pair-component producer with all the scope limits stated in the candidate.

## 1. Locality, Markov order, and weights

For ordered distinct cells, cancellation of the time weights leaves the smooth response integrand \(e^{-s}Dh(X_t)h(X_s)\). The triangular restriction on times disappears across ordered cells, so the exact pair Hoeffding component factors into its separately centered local matrix and vector factors. All other cell midpoints are absent from that component.

Starting from \(Y\), generating only the four selected cell endpoints via the OU transitions over the actual gaps has exactly the required marginal distribution. Adjacent cells must literally identify their common endpoint. The semigroup order

\[
 e^{-2b}P_a\Gamma_dP_{b-a-d}B_d
\]

is correct. In particular, the matrix sandwich is not commuted through an intervening transition.

The finite geometry normalization is also correct:

\[
 M_{N,d}=\sum_{i+r\le N-2}e^{-2(i+r+1)d}
 ={q\,[1-Nq^{N-1}+(N-1)q^N]\over(1-q)^2},\quad q=e^{-2d}.
\]

The infinite resolvent formula follows by separately summing the nonnegative \(i\) and \(r\) geometric series, preserving the operator order.

## 2. Rectangular secant and conditional integration

In one replica, the three inner VALUES are the shared anchor and two inner evaluations; the four outer VALUES are the rectangular secant. Hence there are exactly seven original source evaluations.

For the two replicas, geometry, selected endpoints, and all four standardized midpoint roots must be shared. Local integration times, bridge residuals, and smoothing roots must be independently resampled. This is exactly the candidate's prescription. It produces the square of the conditional integral. Sharing the integration times instead would generally replace that square by an average of squared integrands and would be incorrect.

The midpoint-pair rectangular difference is correctly normalized by \(1/2\): averaging its Gram over an independent copy of each midpoint doubles the inner covariance and then doubles the outer covariance, cancelling the factor \(1/4\).

There are two harmless but worthwhile notational corrections to the description of the secant limit:

- The displayed secant tends to the negative rectangular derivative response, owing to the shifts \(-\delta u\). This does not change its Gram.
- Since the earlier definition of \(B_J\) uses the absolute-time weight, the local-time sampler targets \(e^b B_J\). The geometry probability separately supplies \(e^{-2b}\). The text should state this explicitly rather than identify the unweighted local limit directly with the globally weighted expression.

## 3. Uniform inner Lipschitz bound

The Hessian interval is useful here: for any \(x,x'\),

\[
 \|Dg(x)-Dg(x')\|_{\rm op}\le A.
\]

Differentiating the secant with respect to \(\zeta\) therefore gives operator norm at most \(A^2\tau_s/2\); the same holds for \(\zeta'\). The joint derivative bound is \(A^2\tau_s/\sqrt2\). It is independent of the finite-difference scale and all nuisance roots. Swapping \(\zeta,\zeta'\) reverses the secant, giving conditional centering.

With \(L_d=d(1-e^{-d})\), \(p_k=w_k/S\), and \(w_k=M_{N,d}L_d^2\tanh(d/2)/2\), the scaled Lipschitz constant is exactly bounded by

\[
 a_kA^2\sqrt{\tanh(d/2)/2}=A^2\sqrt S.
\]

The dyadic sum \(S\) is uniformly bounded: its small-cell summands are \(O(d^3)\), and its large-cell summands are \(O(d^2e^{-2d})\).

## 4. Paid secant bias

After integrating the outer smoothing root, the derivative of the smoothed source has Lipschitz constant \(CA/e\). The inner anchor gives the pointwise bound \(|u|\le A|z|\), independent of the smoothing-root magnitude. Thus the unscaled, conditionally integrated secant bias has integrated \(L^4\) bound \(C\delta A^3n/e\).

The chosen finite-difference scales satisfy the exact cancellation

\[
 {a_k\delta_k\over e_k}=A^2\sqrt{2S}.
\]

Therefore the scaled conditional mean-field error is \(CA^5n\) in integrated \(L^4\). The target field has \(L^4\) norm \(CA^2\sqrt n\), from the uniform centered Gaussian inner-root bound. The covariance error is consequently \(C(A^7n^{3/2}+A^{10}n^2)\), as claimed.

This argument is integrated over stationary \(Y\) and the selected endpoint distribution. It does not claim a uniform-in-endpoint Taylor-error moment bound. The cap/transport concentration, in contrast, really is uniform in the fixed non-inner roots.

## 5. Block execution and numerical bookkeeping

For the stated independent-block Wasserstein aggregation, all primitive roots of different physical block draws, including their clock choices, should be independent conditional on \(Y\). The blocks can still be packed into fourteen full-vector calls because the physical splitting basis is supplied. Shared clock choices across blocks would need a different joint-law argument.

The Gaussian count of at most \(17n\) per active block is correct: four endpoint innovations, four midpoint roots, four nuisance vectors in each of two replicas, and one output Gaussian. The inactive-slot buffer is additional, as declared.

Clock-weight setup is polylogarithmic and can be shared. Clock selection is polylogarithmic per physical block; total selection work has an additional factor equal to the number of blocks. Computing \(HG\) uses the two rank-one products directly and requires only \(O(n)\) arithmetic per block.

The total \(O(D\,\mathrm{polylog})\) arithmetic assertion also presumes the physical splitting is already in the source oracle's coordinates, or that its basis transforms are fast or separately charged. An arbitrary supplied dense orthogonal change of basis generally adds \(O(D^2)\) arithmetic per full-vector transform. It does not change the fourteen-call original-VALUE count.

For the final scaled law, the VALUE-error budget must include \(a_k\), not just the unscaled secant factor \(1/\delta_k\). Outer VALUE error has amplification \(O(a_k/\delta_k)\); inner VALUE errors have amplification \(O(a_kA)\). These remain polynomial in the retained cutoff parameters, so the correction does not undermine the claimed polylogarithmic precision count. Stable evaluation of the positive geometric normalizations is advisable at small \(d\).

The total error decomposition is consistent: the true-history pair covariance replacement, third-increment pair omission, fine-level cutoff, finite secant bias, cap bias, and marginal two-replica transport error are distinct paid terms.

## Final revision check, 10:54 UTC

The main candidate now explicitly requires \(R_j^2/v\le1/8\) in every physical block, which matches the proved transport bound. The sign/local-time normalization, final scaled VALUE-error amplification, and dimension-dependent clock work were clarified.

The newly supplied frozen-clock complete/caller first bound also checks: \(Ca_k(A/\delta_k+A^2)\) is at most \(CA^{-3/2}d_K^{-5/4}\), up to fixed coarse-scale constants. The transported source derivative has the claimed \(CRL_{\rm source}|G|/\sqrt v\) moment profile, and hence is not a globally bounded native first.

The initially supplied \(A^{22/5}\) higher-order tail theorem controlled \(F_2\) alone. The parent subsequently identified `/workspace/shared/dyadic-pair-secant-port-20261005/COHERENT-BLOCK-ORDER-TWO-EXTENSION.md`, which I read at 10:55 UTC. That extension does prove the stated higher-order tail for the full coherent \((F_2,\Delta_3)\) block, including its cross block. Explicitly citing this input resolves the coherent interaction-order-at-least-two corollary.

The extension's displayed summed theorem assumes \(A^3b^{3/2}\le1\). Retain that hypothesis, or cover the complementary regime by the trivial full-covariance bound \(CA^2\sqrt D\); when \(b\ge A^{-2}\), the displayed leading bound \(A^{22/5}b^{6/5}\sqrt D\) is already at least \(A^2\sqrt D\).

An independent alternative justification for the extra coherent slots is to apply the same third-increment omission proof to \(P_{\ge2}\) on retained levels, using \(P_{\ge2}F_1=0\), and pay the cross and diagonal covariances by \(CA^5K\sqrt{bD}\). This is absorbed into the existing term by changing its constant. Below the cutoff, use the full coherent complete-history fine-scale cap.
