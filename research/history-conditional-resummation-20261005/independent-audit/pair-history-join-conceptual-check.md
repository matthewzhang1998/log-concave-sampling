# Conceptual check: fixed-block true-history pair replacement

Date: 2026-10-05. Inputs: the stated pair-replacement proposal and the existing `FINITE-ORDER-TWO-C2-TAIL.md` theorem. This check does not inspect or certify the proposed 14-VALUE implementation graph.

## 1. Projection argument is sound

At one dyadic level, let \(S\) be the coarse skeleton, let \(P_2\) denote the conditional orthogonal projection onto the sum of exact two-midpoint Hoeffding slots, and write

\[
 r=P_2\mathbb E[R_g\mid\mathrm{fine}],\qquad
 q=-P_2\mathbb E[B_h\mid\mathrm{fine}],\qquad
 w=r-q=P_2V.
\]

All these slots have conditional mean zero given \(S\). Conditional \(L^2\) contraction gives

\[
 \mathbb E[|r|^2\mid S]\le\mathbb E[|R_g|^2\mid S],\qquad
 \mathbb E[|w|^2\mid S]\le\mathbb E[|V|^2\mid S].
\]

Expand \(rr^T-qq^T=rw^T+wr^T-ww^T\). The conditional Cauchy–Schwarz inequality, conditioning further down to \(Y\), and then Jensen and Hölder give

\[
 \left\|\mathbb E[rr^T-qq^T\mid Y]\right\|_{L^2(Y;\mathrm{HS})}
 \le2\|R_g\|_4\|V\|_4+\|V\|_4^2.
\]

No \(L^4\)-boundedness or \(L^4\)-contraction of \(P_2\) has been assumed. Inserting the supplied blockwise bounds gives

\[
 C\left[A^{9/2}n^{3/2}\sqrt{J_d}+A^5n^2J_d\right],
 \qquad J_d\lesssim1+d^{-1/2}.
\]

For a known physical block decomposition of maximum block size \(b\), squaring and summing gives global factors \(b\sqrt D\) and \(b^{3/2}\sqrt D\), respectively. Fixed-block uniformity means these are allowed constants in the final fixed-\(b\) assertion.

## 2. Claimed cutoff exponent is correct

For \(d\le1\), the retained-level dyadic sums are controlled by

\[
 A^{9/2}d_{\min}^{-1/4}+A^5d_{\min}^{-1/2}.
\]

The omitted true pair refinements below the cutoff have their own positive-covariance bound \(CA^2d_{\min}^2\), obtained by bounding the pair covariance by the complete true increment covariance and summing the fine-scale geometric series.

At \(d_{\min}=A^{10/9}\), the first and third powers are both \(A^{38/9}\), and the second is \(A^{40/9}\). Coarse levels \(d>1\) additionally cost \(C_b\sqrt D(A^{9/2}+A^5)(1+\log_+T)\). Discretizing the cutoff to a neighboring dyadic level changes constants only. Thus the claimed fixed-block leading gain is consistent.

This argument is a covariance comparison. Any final Gaussian-law comparison must separately charge its covariance-to-transport conversion and any baseline spectral scale.

## 3. Omitting the pair slots of the third increment

If \(\Delta_3=F_3-F_2\), \(\|\Delta_3\|_4\lesssim A^3\sqrt n\), and \(P_2F_1=0\), then the pair slot of \(F_2\) is controlled by \(R_g=F_2-F_1\), whose \(L^4\) norm is \(O(A^2\sqrt n)\). The same conditional \(L^2\) argument bounds the pair cross covariance by \(CA^5n\); the omitted \(\Delta_3\)-diagonal covariance is \(CA^6n\). Over \(K\) retained levels and physically independent blocks this is \(CA^5K\sqrt{bD}\), up to a constant when \(A\le1\).

The assertion concerns those pair slots only; it does not dispose of the coarse, singleton, or higher-order pieces of the coherent block.

## 4. Exact two-cell locality of the smooth comparison

Writing the ordered term as

\[
 B_h^T=\int_0^T\int_t^\infty e^{-s}Dh(X_t)h(X_s)\,ds\,dt,
\]

consider distinct coarse cells \(I_i<I_j\). Their order makes the condition \(s\ge t\) automatic. Their contribution is exactly

\[
 \left(\int_{I_i}Dh(X_t)\,dt\right)
 \left(\int_{I_j}e^{-s}h(X_s)\,ds\right).
\]

Given the coarse skeleton, the cell bridges are independent. Conditional expectation onto fine midpoints therefore factors as a matrix-valued local function of midpoint \(i\), times a vector-valued local function of midpoint \(j\). The exact pair projection is the product of their separately centered functions. Each function reads only that cell's endpoints and midpoint. Same-cell terms and terms with the second time beyond \(T\) produce only singleton or constant midpoint dependence and hence no pair contribution.

After averaging the local midpoint functions, the pair covariance is a completely positive matrix sandwich. Integrating the finitely many local endpoints uses the OU Markov transition over the actual gap, so this comparison covariance need not read a full coarse skeleton. The weight \(e^{-2b}\), when \(b\) denotes the inner cell's start time, is consistent with squaring the inner factor's \(e^{-b}\) time weight. One should use a different symbol for this start time and for the maximum physical block size in the final statement.

This locality belongs to the explicitly factored smooth comparison. It is not a claim that the original true-history pair coefficient itself already has a local readset. The preceding covariance comparison is what pays for replacing the true coefficient by the local one.

## Follow-on construction audit

The actual 14-VALUE graph, nuisance-sharing pattern, centering symmetry, uniform inner Lipschitz constant, finite secant bias, finite geometry normalization, root counts, and covariance-to-law conversion were subsequently checked in `local-fourteen-value-graph-audit.md` and `conditional-gram-transport-audit.md` in this directory. Those items are no longer open. The reports retain the stated physical-block, unrestricted-dimension, native-first, and marginal-law scope limitations.
