# Independent audit: resolvent covariance stability and finite j coefficient source

Date: 2026-10-04.

## Verdict and exact scope

**PASS for the dimension-free resolvent covariance theorem, its application to the pinned nested mean, the positive finite covariance target, and the literal original-gradient VALUE coefficient source.** No mathematical correction to the frozen source is required for these claims.

This is **not** a PASS for a covariance action, a positive order-four reserve, or the full order-four endpoint. An admitted VALUE action for the nonsymmetric product `B_Q B_Q*`, the localized mixed current `K`, and the same-carrier mean `m3` remain open exactly as stated in the source. No count is assigned to the missing action.

Frozen source reviewed:

- `../RESOLVENT-COVARIANCE-STABILITY-AND-FINITE-J-KERNEL.md`
- SHA256 `359cdfffa1ce4b7467d24ed24d6b8efdedb85a16c5de0b650d3b4ea37aa52cfd`

Imported nested-mean source read, including its pre-outer comparison and independent audit:

- `../../THIRD-ORDER-NESTED-FORCE-MEAN-WITH-C2-DECOUPLING.md`
- SHA256 `b056e213ac28c60869727589d3ee6d1b2c53d6666908cdbe1961fa3121161ced`

The independent checker `check_resolvent_kernel_independent.py` passes **1,930 assertions**, seed `202610042337`. It imports no author checker. Its JSON output is `resolvent_kernel_independent_checks.json`. These finite tests supplement the proofs below; they do not establish a universal estimate by experiment.

## 1. Row-Hessian smoothing really is dimension-free

Use the convention `Df[output,input]`. If `Lip(f)<=L`, `v=R1 f`, and `m` is a fixed output row, then

    Dv = integral_0^1 r P_r(Df) dr,
    ||Dv||op <= L/2.

For `q=sqrt(1-r^2)`, differentiating one further time by Gaussian integration by parts, rather than differentiating `Df`, gives

    D2(m dot v)(x)
      = integral_0^1 (r^2/q) E[(Df(rx+qG)*m)G*] dr.

Depending on the displayed Hessian index convention, the matrix inside the expectation can be transposed; its HS norm is unchanged. No second derivative of `f` is used.

For an arbitrary vector-valued square-integrable `a(G)`, apply scalar Bessel to each component against the orthonormal coordinates `G_1,...,G_D` and add:

    sum_(i,k) (E[a_i G_k])^2 <= sum_i E[a_i^2].

Here `a=Df(rx+qG)*m`, so the right side is at most `L^2 |m|^2`. Minkowski and the exact integral

    integral_0^1 r^2/sqrt(1-r^2) dr = pi/4

prove

    ||D2(m dot v)(x)||HS <= (pi/4)L|m|.

This is precisely the row-contracted tensor norm required later. It is not an operator-to-HS estimate with a hidden dimension multiplier. Squaring and summing over the rows of a frozen matrix `N(x)` gives

    sum_i ||D2(N_i(x) dot v)(x)||HS^2
       <= (pi L/4)^2 ||N(x)||HS^2.

One may choose these rows separately at each x because the assertion is pointwise and linear in m. In the weak formulation, a countable dense collection of rows gives a common full-measure set; continuity in m extends it to all rows. Derivatives of the varying N are a separate product-rule term, not omitted from this argument.

The calculation holds for Lipschitz f by Gaussian smoothing/Sobolev closure. In the application, `g=grad U` is C1 with bounded first under U's original C2 hypothesis, and the source uses only this first. The new proof introduces no uniform continuity modulus or bounded higher derivative of g.

## 2. The divergence bound preserves Gaussian cancellation

For each output row of T, Gaussian divergence is the L2 adjoint of differentiation. Its exact finite-dimensional Sobolev isometry is

    ||delta T||2^2
      = ||T||_(L2;HS)^2
        + sum_(i,j,k) E[(partial_j T_ik)(partial_k T_ij)].

The second term can have either sign. Its upper bound by `||DT||_(L2;HS)^2` gives

    ||delta T||2 <= sqrt(||T||2^2+||DT||2^2)
                   <= ||T||2+||DT||2.

This estimate is dimension-free. Separately bounding `xT` and the trace derivative would lose the cancellation and is not what the source does.

For `T=N Dv`, the true product rule has two terms:

- `(partial_l N) Dv`, whose tensor norm is bounded by `(L/2)||DN||HS`.
- `N partial_l Dv`, whose output-row tensor is the Hessian of the scalar frozen-row function `N_i(x) dot v` and is bounded by `(pi L/4)||N||HS` after summing rows.

Consequently

    ||delta(N Dv)||2
      <= (L/2+pi L/4)||N||2+(L/2)||DN||2.

The finite Hermite tests independently verify the full isometry, including its crossed derivative indices, in several input dimensions. They do not replace its Sobolev proof.

## 3. Paraproduct estimate, constant, and covariance identity

Let `P_t=P_(exp(-t))`. For a matrix test M, Hilbert-valued Gaussian Bessel and invariance give

    ||D P_t M||_(L2;HS)
      <= exp(-t)/sqrt(1-exp(-2t)) ||M||_(L2;HS).

This norm includes all matrix output and derivative indices. It has no extra matrix-size factor.

For `u in W^(1,2)` and `v=R1 f`, with first smooth truncations if necessary, the matrix orientation is

    (Du Dv*):P_t M = Du:((P_t M)Dv).

Thus self-adjointness and divergence duality give

    <R2(Du Dv*),M>
      = integral_0^infinity exp(-2t)
          E[u dot delta((P_t M)Dv)]dt.

Use the preceding divergence estimate and `||M||2=1`. The two scalar integrals are

    integral_0^infinity exp(-2t)dt = 1/2,
    integral_0^infinity exp(-3t)/sqrt(1-exp(-2t))dt = pi/4.

The second is obtained by `r=exp(-t)`. Therefore the resulting coefficient is exactly

    (1/2+pi/4)(1/2)+(1/2)(pi/4) = (1+pi)/4,

and

    ||R2(Du D(R1 f)*)||_(L2;HS)
       <= [(1+pi)/4] Lip(f) ||u||2.

The heat singularity is integrable at t=0. Smooth approximations, then the displayed uniform bound and Sobolev closure, justify Fubini, integration by parts and duality for the actual Lipschitz sources. Here `u=R1(f1-f0)` is itself in the necessary first-order Sobolev class. No estimate on the derivative of the small error is asserted.

For completeness, the covariance functional is the true Markov-path covariance. If

    H_f=integral_0^infinity exp(-t)f(X_t)dt,
    v=R1 f,

then the conditional martingale is

    E[H_f|F_t] = integral_0^t exp(-s)f(X_s)ds + exp(-t)v(X_t).

Using `(1-L_OU)v=f`, its Brownian coefficient is `sqrt(2)exp(-t)Dv(X_t)`. Conditional Ito isometry yields

    Cov(H_f|Z)=2R2[Dv Dv*](Z).

This identity holds for vector-valued, non-gradient f. The transpose is required.

With `u=R1(f1-f0)` and `v=R1(f1+f0)`, exact noncommutative polarization is

    C(f1)-C(f0)=R2[Du Dv*+Dv Du*].

The two terms are transposes of each other after R2. Applying the proved estimate to each gives exactly

    ||C(f1)-C(f0)||_(L2;HS)
      <=2[(1+pi)/4][Lip(f1)+Lip(f0)]||R1(f1-f0)||2.

There is no missing factor two, reversed product, symmetry assumption or second dimension-sized energy.

## 4. The pinned nested proof supplies precisely the imported error

The pinned prior Sections 2–3 compare the two inner conditional laws before the outer quadrature is introduced. The true conditional Markov inner law at the ancestor x defines j(x); the cheap common-H inner law defines j_Q(x). Both are independent of the outer r once x is exposed. Hence the two integrated conditional force means are exactly `R1 j` and `R1 j_Q`.

The prior decoupling estimate is used twice, once for each private inner law. Multiplication by the A-Lipschitz outer g gives `C A^3(|z|+sqrt(D))/sqrt(1-r^2)`, whose r integral is finite. The middle deterministic-shift comparison uses

    ||E_H H_Q - R1 g||2 <= delta_in A sqrt(D)

and therefore costs `delta_in A^2 sqrt(D)`. This proves the imported statement

    ||R1(j-j_Q)||2 <= C[A^3+delta_in A^2]sqrt(D).

The prior's additional outer quadrature error `delta_out A` is absent here because that operation has not yet occurred. The present source does not import the full finite force-mean bound and mistakenly drop such a term.

Both actual conditional fields have first bounded by `L_j=A(1+A/2)`: the conditional inner x first is a positive average of `tau Dg`, lying in `[0,A I/2]`, and the outer chain rule is `Dg(x-inner)(I-D_x inner)`. Thus the new theorem gives

    ||Cov(H_j|Z)-Cov(H_(j_Q)|Z)||_(L2;HS)
       <=C[A^4+delta_in A^3]sqrt(D).

The original prior proof integrates its private roots in a marginal conditional-law comparison. Nothing here preserves or reuses those old roots. This application does not identify the covariance of the cheap raw common-root chord with the covariance of H_j. It first obtains a deterministic conditional-mean field and then inserts that field into the true Markov covariance functional. The earlier shared-root separator is therefore unaffected.

## 5. Finite r and q clocks retain the transpose and the Gaussian caller law

Write `B=D R1 j_Q=R2(Dj_Q)`. On Hermite degree n, the r approximation to B uses the moment `sum_i w_i r_i^(n+1)`, whereas the exact multiplier is `1/(n+2)`. Thus a uniform moment error `delta_r` really implies

    ||B_Q-B||_(L2;HS) <= delta_r ||Dj_Q||_(L2;HS)
                      <= delta_r L_j sqrt(D).

Positivity and the exact first moment give `||B_Q||op,||B||op<=L_j/2`. They do not imply symmetry.

The appropriate noncommutative product identity is

    B_Q B_Q* - B B* = (B_Q-B)B_Q* + B(B_Q*-B*).

Hence its L2 HS norm is at most `delta_r L_j^2 sqrt(D)`. The finite positive q operator `2 sum_k a_k q_k P_(q_k)` has mass one and is an L2 contraction. Its shifted-degree multiplier error against `2R2` is at most `2 delta_q`. Finally

    ||B B*||_(L2;HS) <= L_j^2 sqrt(D)/4.

An explicit valid bound is therefore

    ||C_Q-C(j_Q)||_(L2;HS)
       <= L_j^2(delta_r+delta_q/2)sqrt(D).

Also `C_Q` is symmetric PSD and `C_Q<=L_j^2 I/4`, pointwise. Combining this with Section 4 and tolerances at most A^2 gives the source's order-four L2 HS coefficient-target accuracy.

All these clock replacement statements use a fresh standard Gaussian endpoint Z. They neither assert pointwise endpoint accuracy nor authorize applying the same L2 estimate at an arbitrary approximate caller law without the stated host join/restoration.

## 6. Literal original-g VALUE graph and mean Jacobian

For each fixed captured X, the source is

    V_X(u,H)=sum_i alpha_i[J_Q(r_iX+c_i u,H)-J_Q(r_iX,H)],
    c_i=sqrt(1-r_i^2), alpha_i=w_i r_i/c_i.

Each of the two branches at each r_i has exactly the following graph:

1. Its own exposed x, either `r_iX+c_i u` or `r_iX`.
2. N_in original VALUES at `tau_l x+sqrt(1-tau_l^2)H`.
3. One terminal original VALUE at x minus that branch's computed positive inner sum.

The same H is deliberately retained across all inner nodes and both branches within the occurrence. The terminal branch cannot substitute another branch's inner ancestor.

The safe raw count is therefore exactly the advertised upper bound

    Q_V <= 2 N_r(N_in+1).

It does not assume a global cache of the anchors. The Gaussian raw input is precisely two D-roots, `(u,H)`. Clock nodes add deterministic sites rather than new private roots.

Differentiation in u yields

    D_u V_X=sum_i w_i r_i D_xJ_Q(r_iX+c_i u,H).

Since u and H are independent standard roots when this coefficient is evaluated, its expectation is exactly B_Q(X). No HVP-valued producer is being introduced by this target identity.

The independent finite checker executes original VALUES, separately checks every inner and terminal site, and verifies firsts by central differences. A separate tensor Gauss-Hermite integration verifies the coefficient through Stein's identity

    E[V_X(u,H)u*]=E[D_u V_X]=B_Q(X),

with maximum entry error `2.7149e-11` in its nonlinear two-dimensional fixture. Its quadratic calibration gives the exact matrix coefficient `B_Q=(1/2)S(I-S/2)` for `g(x)=Sx`.

## 7. Actual firsts, coisometry, curl, and private-H anchors

Let `beta_in=sum_l v_l sqrt(1-tau_l^2)`, `S1=sum_i alpha_i`, and `S2=sum_i alpha_i r_i`. The original C2 chain rule gives

    ||D_x J_Q||op <= L_j,
    ||D_H J_Q||op <= A^2 beta_in,
    ||D_x J_Q-(D_x J_Q)*||op <= A^2.

The last expression is the commutator of two symmetric matrices; they cannot be commuted. The complete source bounds follow directly:

    ||D_u V_X||op <= L_j/2,
    ||D_H V_X||op <= 2 A^2 beta_in S1,
    ||D_X V_X||op <= 2 L_j S2.

There is no reason to replace the caller bound by O(A^2). The first difference of outer Hessians need only be O(A), consistent with the original C2 class.

At u=0, the two complete branch graphs coincide for every X,H. With coherent same-key reuse their difference is literally zero, including in the numerical program. Integrating the u first at fixed H gives

    |V_X(u,H)| <= L_j |u|/2,
    ||V_X||_(Lp(u,H)) <= C_p A sqrt(D),

uniformly in X. This is a one-energy bound depending on the D-dimensional u, even though the raw input dimension is 2D. The latter dimension still belongs in any later consumer guard.

For `P=(I,0)`, the square lift has derivative and antisymmetric part

    D(P*V) = [[D_u V,D_H V],[0,0]],
    Curl(P*V) = [[D_uV-(D_uV)*,D_HV],[-(D_HV)*,0]].

The u-u block is bounded by `A^2 sum_i w_i r_i=A^2/2`; the off-diagonal skew block has norm `||D_HV||op`. A valid complete bound is

    ||Curl(P*V)||op <= A^2/2+2A^2 beta_in S1.

Thus the source is near-gradient with O(A^2) curl and O(A) full/private first. It is not licensed for a genuine-gradient action. The nonlinear fixtures explicitly exhibit nonzero curl.

The term `J_Q(r_iX,H)` is independent of u but generally depends on the private H. It is not a caller-only capture. Its H derivative is `-Dg(r_iX-H_Q)D_H H_Q`, and the test fixtures exhibit nonzero anchor changes whenever H changes. Changed H/source/complete-bank records must replay its inner values. Safe full replay also covers changed X or u banks; any exact reuse must be proved at the complete semantic key. The checker makes no cross-bank reuse.

If X is subsequently owned as `q_kZ+sqrt(1-q_k^2)G_k`, the actual D_X port gives its G_k and Z chain-rule paths. Neither root can be left frozen during a comparison that promises to integrate the complete bank.

## 8. Dyadic coefficients, absolute precision, and complete counts

The bounded coefficient sums use the **stated admitted positive dyadic Gaussian rule family**, not merely positivity, mass one and first moment one half. This distinction is important: arbitrary positive moment-calibrated rules can put too much mass extremely close to r=1 and have unbounded `sum w_i/c_i`.

For the admitted family set `s=1-r`. A dyadic panel `[a,2a]` has weight mass a, and `c=sqrt(s(2-s))>=sqrt(s)>=sqrt(a)`. Its contribution to `sum w_i/c_i` is at most sqrt(a), independently of within-panel Gauss order. The terminal midpoint on `[0,h]` contributes at most sqrt(2h). Summation gives

    sum_i w_i/c_i <= 1+sqrt(2)-sqrt(h) < 1+sqrt(2).

Both S1 and S2 are bounded by this sum. Hence there is no inverse endpoint heat factor in the aggregate VALUE or caller bounds. Individual alpha_i remain in the numerical leaf ledger. The inherited fixed-grade dyadic rule has `O(log^2(1/A))` nodes; the source has not asserted a new uniform moment theorem for an arbitrary rule.

An explicit absolute VALUE-error calculation is available. If the inner leaf error in branch sign s is `eps_(i,s,l)` and the terminal leaf error is `eps_(i,s,term)`, the output error is at most

    sum_(i,s) alpha_i [eps_(i,s,term)
                         + A sum_l v_l eps_(i,s,l)].

The A term comes from moving the terminal query by the error in its computed inner packet. A uniform leaf tolerance eps therefore gives `2 S1(1+A)eps`. This is absolute, never divided by a small actual energy. Coherent same-key evaluation preserves the exact u=0 cancellation. The numerical checker tests this bound and the literal numerical zero.

This VALUE-error calculation does **not** by itself prove derivative accuracy of a numerical implementation. The actual exact-source firsts are the chain-rule bounds in Section 7; numerical first/adjoint, mode, original-caller and finite-version restoration remain the imported source/oracle obligations. This is compatible with the source's explicit separation of those floors and is not an admission of a missing consumer.

All original first/adjoint requests can be made as HVP sweeps at the already recorded inner and terminal VALUE points. Reverse or forward sweeps include the feedback derivative through the inner packet. No derivative of an HVP is needed, and no HVP is a new primal producer.

With N_action complete occurrences per q bank, the transparent safe recurrence is

    2 N_q N_r(N_in+1)N_action
      + actual caller/known/numerical/replay work.

This remains conditional on an admitted action and its complete tape. The raw 2D dimension is not the final action dimension. Original mode/caller restoration, discarded-primal replay, actual first/adjoint work, coisometry and frozen tolerances are not silently removed from the bill.

## 9. Checks and remaining gates

The 1,930 passing assertions cover:

- Exact frozen source and pinned-prior hashes.
- Both pi/4 integrals and the displayed paraproduct constant.
- Exact finite-Hermite Gaussian divergence isometry, including crossed indices.
- Dimension-free row-Hessian tests for non-gradient vector fields up to D=41.
- Hermite-resolved scalar paraproduct tests at multiple chaos degrees/frequencies.
- Exact nonsymmetric covariance polarization and an explicit ordinary-square separator.
- Positive dyadic clocks, mass/first moment, shifted Hermite moments and coefficient sums.
- Correct original VALUE sites/counts, u/H/X chain-rule checks, complete lifted curl, one-energy amplitude, private-H anchor dependence, coherent zero and absolute precision.
- Matrix-quadratic calibration and nonlinear Gaussian coefficient/Stein verification.

The test suite's finite numerical checks are not a proof of uniform clock error, a substitute for the pinned decoupling theorem, or a covariance action. No such inference is made.

The following are still substantive unfilled gates:

1. A finite admitted original-VALUE action realizing `B_Q B_Q*`, at the actual near-gradient radius, curl, raw/completed dimensions, padding and caller guards.
2. A positive covariance reserve whose error meets the order-four allowance. The previous cubic reserve error does not automatically improve.
3. A native same-I service for the mixed K-current.
4. The same-carrier Gaussian-history mean m3, with its actual source-zero/caller/floor/count ports.
5. The full host join and independent reviews of these services.

In particular, replacing `B_Q B_Q*` by `B_Q^2` is invalid for the actual nonsymmetric source, and no completed Gaussian mean-law coupling is a strong value at a saved old root. The frozen source respects both restrictions. The present audit closes the resolvent theorem and finite coefficient-source admission only.
