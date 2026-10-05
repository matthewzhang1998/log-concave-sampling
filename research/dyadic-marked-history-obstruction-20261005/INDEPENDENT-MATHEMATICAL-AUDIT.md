# Independent mathematical audit

Date: 2026-10-05. Audited `DYADIC-MARKED-HISTORY-NONCLOSURE.md` directly. This review uses analytic calculations and deterministic one-dimensional scalar quadrature only. No OU path grid, stochastic path simulation, nested-history sampler, external research, or external upload was used.

Audited final snapshot SHA-256: `2007b8a22fdc6ee1bbd954902bdb2d97d04ad41bcfaf48e253799420a486c42f`.

## Verdict

The mixed-midpoint obstruction, the critical A⁴ covariance refinement, and the one-Hilbert-Schmidt orthogonal-residual bound are mathematically sound under the stated stationary OU normalization and target assumptions. The document correctly confines the obstruction to direct reuse of the additive-prefix, exact-singleton architecture. It does not establish a general impossibility for covariance representations, interacting marked histories, or other finite-source constructions.

## Independent checks

1. **Admissible original gradient.** The witness satisfies `g(0)=0`, is smooth, and obeys `A/4 <= g' <= 3A/4 < A`. Consequently `U'=g` has the stipulated convexity and Hessian bound. Every derivative of `g` of order at least two is bounded. The witness is unbounded and is not odd, but neither property conflicts with the target assumptions stated in Section 1. This checks mathematical admissibility against those stated assumptions; it does not independently authenticate the pinned source manifests.

2. **Exact signs and mixed derivative.** Write `L_h=R(g'(X)h)` and `L_k=R(g'(X)k)`. Differentiating `F2=R g(X-F1)` gives the integrand `g''(X-F1)(h-L_h)(k-L_k)`, because `D²F1[h,k]=R(g''(X)hk)=0`. Ordered disjoint supports imply both `hk=0` and `k L_h=0`, leaving exactly `g''(X-F1)(-h L_k+L_h L_k)`. Thus Equation (7) has all terms and correct signs.

3. **Conditional lower bound.** At the zero finite coarse skeleton and the two zero midpoint values, the remaining conditional Gaussian process is centered and has point variances at most one. Since `g(0)=0` and its Lipschitz constant is `u`, Minkowski gives `||F1,t||₂ <= u` and `||X_t-F1,t||₂ <= 1+u`. Markov's inequality supplies the stated probability `3/4`, hence Equation (8). An 80-decimal-digit mpmath evaluation gives

   `(3/4) sech²(2.015) = 0.05147679585244062351344992936763936673623211989027285938091875905076523`,

   exceeding `1/20` by `0.001476795852440623513449929367639366736232119890272859380918759050765228`.

4. **Resolvent overlap and final constant.** If `h` is earlier than `k`, Tonelli gives

   `J0 = Rk(0) ∫h(s)ds`,

   `K0 = Rk(0) ∫(1-e^(-s))h(s)ds`.

   Thus `0 <= K0 <= J0`, with `J0>0`. The resulting upper bound is

   `A² J0[-1/320 + 9A/64] <= -(11/6400) A² J0`

   for `0<A<=1/100`. Correlations between `g''(X-F1)` and `L_h,L_k` do not invalidate this step: deterministic lower and upper envelopes for the latter quantities are used before expectation.

5. **Positive pair variance.** The exact conditional projection is smooth in the two midpoint values. Its negative mixed derivative excludes an additive decomposition almost everywhere under the product of the two nondegenerate conditional Gaussian laws. Therefore its conditional pair Hoeffding component has positive variance. The document's continuity argument in the finitely many coarse endpoints correctly extends this beyond the probability-zero zero skeleton to a positive-probability open set. This is an interaction within the F2 marginal, so it is distinct from merely omitting the F2/Delta cross block.

6. **Hat identities.** Direct integration independently yields `∫_I h = 2 tanh(d/4)` and `∫_I e^(-s)h(s)ds = (d/2)e^(-(a+d/2))`. Splitting deterministic quadrature at the midpoint, at 80-digit precision, checked `(a,d)=(0,0.2),(1.3,0.7),(0,2),(2,0.1)`; absolute discrepancies were zero at displayed precision or at most `8.24×10^(-84)`. These confirm the displayed closed-form `J0` and `K0/J0`; they are scalar checks, not path discretization.

7. **Full-history sensitivity and residual.** The derivative recursion gives `beta_j(s)=A e^(-s)∑_(r=0)^(j-1)(As)^r/r!` by convolution and operator-norm bounds without commuting Jacobians. Stacking F2 with F3-F2 is bounded by `2 beta_2+beta_3`, exactly the displayed beta. Given the endpoint skeleton, the bridge pieces and future tail are independent Gaussian sources. Their row bounds produce `C d∑_I(∫_I beta)² + C(∫_T^∞ beta)²`; Cauchy-Schwarz bounds the first term by `C d²∫beta²`, and exponential-polynomial integration bounds the second by `C A²(1+T⁴)e^(-2T)`. Applying this in every physical unit direction proves the uniform residual operator bound, followed by one HS conversion with factor `sqrt(2D)`.

8. **Exact scope of the remainder.** The bound concerns `E[Cov(W|S_K)|Y]`, equivalently the covariance loss under the exact orthogonal projection `E[W|S_K]`. It does not, by sensitivity alone, establish the same bound for evaluating W on the mean/interpolated path. The document uses the former throughout. Including time zero and T among conditioned endpoints is essential to this derivation and is explicit in the document.

9. **Added explicit constants.** For `A<=1/2`, direct integration gives `A^(-2)∫beta²=9/2+(9/2)A+3A²+(9/8)A³+(3/16)A⁴ <=1959/256<8`. Also `A^(-1)e^T∫_T^∞beta <=19/4+(7/4)T+(1/8)T²`. Squaring with Cauchy-Schwarz and using `T²<=(1+T⁴)/2` proves the stated bound with constant 80. An OU bridge's maximum point variance is `tanh(d/2)<=d/2`, so its displayed row norm `sqrt(d/2)` is valid.

## Minor caveats and interpretation

- Gaussian smoothing generally changes `g(0)`. The revised note correctly restores it by subtracting the smoothed value at zero. This preserves the derivative bounds; sufficiently small smoothing preserves the strict interaction by continuity. The initial normalization caveat is therefore resolved in the audited snapshot.
- The displayed positive identities are exact analytic identities, not a finite original-VALUE implementation. The note correctly leaves the conditional-mean supplier and its complexity unsupplied.
- The nonzero pair term prevents equality with the sum of the exact singleton midpoint covariances. It cannot prohibit reallocating interaction variance into different effective fields, nor does it prohibit sequential martingale representations. The note explicitly retains these possibilities.
- This audit verifies a distinct, scoped contribution within the supplied development. It does not claim priority or novelty relative to external mathematical literature.

**No substantive correction is required for Equations (1)–(10).** The one minor normalization qualification identified during review is resolved in the audited snapshot.

## Additional audit: the critical A⁴ covariance refinement

Equations (11)–(13) are valid, with the fixed horizon and chosen cells independent of A as explicitly assumed.

1. **Taylor remainder, including moment scope.** Since `f''<=1/4`, writing `Q=Rf(X)` gives the pointwise bound

   `|F2,0-A Rf(X)(0)+A² B2(X)| <= (A³/8) R(Q²)(0)`.

   The linear-growth bound `|f(x)|<=3|x|/4`, Minkowski, and stationary Gaussian moments imply `||R(Q²)(0)||_p<∞` for every fixed finite `p>=1`. Hence the remainder is uniformly `O_Lp(A³)` as A tends to zero. Under fixed finite coarse values the same reasoning applies to a Gaussian process with a bounded conditional mean and point variances at most one; the conditional constants may depend on those coarse values. No skeleton-uniform Taylor-moment constant is asserted or needed. Integrating under the original stationary Gaussian caller law restores finite global constants.

2. **Mixed derivative of B2.** Differentiating `B2=R[f'(X)Q]` twice produces

   `R[f'''(X)hk Q + f''(X)h R(f'(X)k) + f''(X)k R(f'(X)h) + f'(X)R(f''(X)hk)]`.

   The first and last terms vanish because `hk=0`; the third vanishes by support ordering. The remaining term is exactly Equation (12) and is strictly positive on every continuous path. Differentiation under the conditional Gaussian expectation is justified by the bounded derivatives and resolvent envelopes. Thus the pair projection Psi is nonzero for every fixed finite coarse skeleton, strengthening the earlier zero-skeleton witness.

3. **Projection and conditional covariance remainder.** The pair projection is the sum of four conditional expectations, so its operator norm on every `L^p`, `p>=1`, is at most four. The order-A term has zero pair projection by disjoint bridge-cell additivity. Write the result as `d_A=-A² Psi+A³ e_A`, where `sup_A||e_A||_4<∞`. For any retained sigma field G contained in the coarse skeleton,

   `E[d_A²|G]=A⁴ q_G+A⁵ r_A,G`, where `q_G=E[Psi²|G]`,

   `r_A,G=-2E[Psi e_A|G]+A E[e_A²|G]`.

   Conditional-expectation contraction and Hölder give

   `||r_A,G||_2 <= 2||Psi||_4||e_A||_4 + A||e_A||_4²`.

   This proves the integrated L² conditional-covariance claim using p=4. In particular one may take G to be the full coarse skeleton or just Y. The conditional means of the pair component vanish, so these second moments are its actual covariance contributions.

4. **Dimension dependence and exact norm.** For independent coordinatewise copies of the scalar witness, the selected pair's omitted F2 covariance contribution conditional on `Y=(Y_1,...,Y_D)` is diagonal, with entries `v_A(Y_j)=E[d_A²|Y_j]`. Set `q(y)=E[Psi²|Y=y]` and `c=||q||_(L²(Y))>0`. Then

   `||diag(v_A(Y_1),...,v_A(Y_D))||_(L²(Y;HS)) = sqrt(D)||v_A||_2 = c A⁴ sqrt(D)+O(A⁵ sqrt(D))`.

   Constants are independent of A and D; they may depend on the fixed horizon and chosen cells. This is an equality at the indicated asymptotic order for the selected pair contribution. For the entire omitted PSD sum, it provides a lower bound unless the other contributions are also included in the definition of c. Allowing the chosen cell geometry to vary with A would require further control of c and is not covered by this fixed-cell witness.

The refinement therefore locates a genuine omitted covariance contribution at the claimed A⁴ one-HS scale. It remains a witness against deleting actual Hoeffding interactions, not an obstruction to a construction that supplies or reorganizes them. The final main-note snapshot explicitly defines `q_pair(Y)` and the integrated `L²(Y;HS)` norm, resolving the norm-interpretation ambiguity. Equations (11)–(13) require no substantive correction.
