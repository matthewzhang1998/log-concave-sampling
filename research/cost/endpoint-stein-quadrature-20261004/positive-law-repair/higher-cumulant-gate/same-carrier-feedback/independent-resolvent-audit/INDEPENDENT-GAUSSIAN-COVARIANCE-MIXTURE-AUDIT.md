# Independent audit: one-energy Gaussian covariance-mixture law

Date: 2026-10-04.

## Verdict and pin

**PASS for the stated Gaussian-input covariance-mixture theorem, its constant, and its finite-B_Q specialization.** This is an analytical law consumer, not a finite coefficient/mean producer by itself.

Reviewed `../GAUSSIAN-COVARIANCE-MIXTURE-ONE-ENERGY-LAW.md`, SHA256
`7476219bacabfacc44030a14710408735d3e31b78103376d57f525a3518726dc`.
This corrected pin has the proper Stein-kernel orientation. The earlier superseded pin was not admitted by this audit.

The theorem is

    W2(Law(Y),N(0,nu I+E C))
       <= Lip_(X -> HS)(C) ||C-E C||_(L2;HS)
                                      /(sqrt(2)nu^(3/2)),

where X is a standard Gaussian tape, `C(X)>=0`, `Y|X~N(0,nu I+C(X))`, and `nu>0` is fixed. The hypothesis is a Lipschitz bound into the physical Frobenius Hilbert space, not merely an operator-norm Lipschitz bound.

## 1. Exact Stein indices

Choose an orthonormal basis `E_a` of Sym_D under the Frobenius inner product, write `F_a=<C-E C,E_a>`, and let `h_a=(-L)^(-1)F_a`. Gaussian integration by parts gives

    E[F_a phi(F)] = E[sum_k D_k h_a D_k(phi(F))]
                  = E[sum_b sum_k D_k h_a D_k F_b partial_b phi(F)].

Therefore the kernel entries are

    tau_ab = E[sum_k D_k(-L)^(-1)F_a D_k F_b | F].

The frozen corrected source has precisely this order. Transposing it does not in general preserve the Stein identity. In particular, for the diagnostic `F=(X,X^2-1)`, the unconditioned matrix is

    tau = [[1,2X],[X,2X^2]],

and F already determines X. For the test `phi(f)=f_1 f_2`, its correct first-row expectation is 2, whereas the reversed-order matrix gives 1. The independent checker detects this difference. This polynomial diagnostic is an index check, not a globally Lipschitz covariance fixture.

The conditional Jensen and matrix product estimate are

    ||tau||_(L2;HS)
       <= ||D(-L)^(-1)F (DF)*||_(L2;HS)
       <= Lip_HS(C) ||D(-L)^(-1)F||_(L2;HS).

Centered Hermite degree n contributes `1/n` to the squared norm in the last expression, so the spectral gap implies that it is at most `||F||2`. This proves the one-energy kernel bound with constant one, independently of both input tape dimension and dim(Sym_D).

## 2. Posterior covariance fluctuation and Wick contraction

For `S=nu I+C`, differentiation in a Frobenius symmetric direction E gives

    D_C phi_S(y)[E]
       = (1/2)phi_S(y)<E,S^(-1)yy*S^(-1)-S^(-1)>.

The basis must be Frobenius-orthonormal, including the `1/sqrt(2)` normalization for off-diagonal symmetric coordinates. Applying the above Stein identity to this density proves, coordinate by coordinate,

    E[F_a|Y=y]
       = (1/2)E[sum_b tau_ab
          <E_b,S^(-1)yy*S^(-1)-S^(-1)> |Y=y].

For a fixed output row a of tau define the physical symmetric matrix `T_a=sum_b tau_ab E_b`. Conditional on C, tau is fixed, and `Y=S^(1/2)N` with N an independent standard Gaussian. Its contraction is exactly

    <S^(-1/2)T_a S^(-1/2), NN*-I>.

Wick's formula gives squared L2 norm `2||S^(-1/2)T_a S^(-1/2)||HS^2`. Since `S>=nu I`, this is at most `2 nu^(-2)||T_a||HS^2`. Summation over a costs `||tau||HS^2`, not a second physical dimension factor. Conditional Jensen and the factor 1/2 in the density derivative therefore give

    ||E[C-E C|Y]||_(L2;HS) <= ||tau||_(L2;HS)/(sqrt(2)nu).

The source's general smooth/Sobolev approximation is legitimate. One can first truncate by metric projection onto a bounded convex subset of PSD matrices, then Gaussian-smooth the matrix-valued map; positivity and Lipschitz control are preserved. The displayed L2 majorants permit the limit. The actual finite-B_Q covariance target is already uniformly bounded in operator norm.

## 3. The Stein-to-W2 interpolation has the stated orientation and constant

Conditional Gaussian integration by parts gives `T_Y(y)=E[nu I+C|Y=y]` as a Stein matrix of Y. With `Sigma=nu I+E C`, take independent `G~N(0,Sigma)` and

    Y_t=tY+sqrt(1-t^2)G.

For a smooth test phi, differentiation and the two Stein identities yield

    d/dt E phi(Y_t)=t E[(T_Y-Sigma):D2 phi(Y_t)].

Integrating once in G produces the continuity velocity

    v_t(Y_t)=t/sqrt(1-t^2)
              E[(T_Y-Sigma)Sigma^(-1)G |Y_t].

The original G is independent of Y and therefore of `T_Y(Y)`. Its covariance computes the exact squared matrix factor

    E_G |(T_Y-Sigma)Sigma^(-1)G|^2
       = ||(T_Y-Sigma)Sigma^(-1/2)||HS^2.

The dynamic W2 bound, conditional Jensen and
`integral_0^1 t/sqrt(1-t^2)dt=1` prove

    W2(Law(Y),N(0,Sigma))
       <= ||(T_Y-Sigma)Sigma^(-1/2)||_(L2;HS).

Since `Sigma>=nu I`, this contributes exactly one additional `nu^(-1/2)`. Together with Section 2 it proves the claimed `1/(sqrt(2)nu^(3/2))` constant. The interpolation Gaussian is analytical only; it is not an actual program root to retain afterward.

## 4. The finite B_Q specialization is dimension-free in the required norm

For the already audited Lipschitz field j_Q,

    B_Q(X)=sum_i w_i r_i P_(r_i)Dj_Q(X),
    ||B_Q||op<=L_j/2.

Although j_Q need not be a gradient, the two input derivative slots of each scalar output commute after Gaussian smoothing. Thus, for fixed input direction h,

    D_X B_Q(X)[h]
      =sum_i (w_i r_i^2/c_i)
          E[(Dj_Q(r_iX+c_iG)h)G*].

The expectation's HS norm is at most `L_j |h|` by Gaussian Bessel. The admitted dyadic coefficient sum bounds `sum w_i r_i^2/c_i`, proving a Lipschitz constant `C A` from X into HS. No second original derivative is required.

For `Q(X)=B_Q(X)B_Q(X)*`, the product rule gives

    ||DQ(X)[h]||HS <= 2||B_Q(X)||op ||DB_Q(X)[h]||HS
                   <= C A^2 |h|.

For the complete independent owned tape `G=(G_1,...,G_N)` and positive weights of total one,

    C_mix(G;Z)=sum_k a_k Q(q_k Z+sqrt(1-q_k^2)G_k),

its directional derivative is bounded by

    C A^2 sum_k a_k |h_k|
      <= C A^2 (sum_k a_k^2)^(1/2)||h|| <= C A^2 ||h||.

This bound is uniform in Z. Its pointwise operator bound `C_mix<=C A^2 I` gives centered HS energy at most `C A^2 sqrt(D)`. The theorem therefore yields `C_nu A^4 sqrt(D)` without a clock-count or extra sqrt(D) loss.

Taking `a_k=2v_kq_k` makes the averaged target exactly the transpose-retaining C_Q from the coefficient source. The specialization replaces a random covariance by its mean in a full law estimate; covariance matching alone was not used as a Gaussianity argument.

## 5. Evidence and scope

`check_mixture_and_gram_independent.py` supplies 277 new finite assertions in addition to the rerun 1,930 coefficient-source assertions. The mixture checks include the asymmetric Stein-index separator, physical symmetric-direction Gaussian density derivatives, exact finite-Gaussian Wick calculations, and scalar quantile-mixture diagnostics with the proved bound. They supplement the analytical proof and do not replace it.

This theorem introduces no original g queries. Its use in a finite algorithm remains dependent on an independently admitted completed conditional mean service, actual owned roots, independent complete banks, fixed positive variance shares, and propagated absolute floors. The accompanying rectangular-Gram audit checks one such bounded service separately. K, m3 and the full order-four endpoint remain outside this result.
