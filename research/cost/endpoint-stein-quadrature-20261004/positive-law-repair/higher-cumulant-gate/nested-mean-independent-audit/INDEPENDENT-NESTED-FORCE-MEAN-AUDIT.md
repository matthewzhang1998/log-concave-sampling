# Independent audit: nested-force mean with uniform C2 decoupling

Date: 2026-10-04.

## Verdict and exact scope

**PASS for the finite force-mean component and its guarded conditional Gaussian mean-law consumer.** The source's actual positive finite VALUE graph satisfies

    ||E_(G,H) F_Q(Z,G,H) - m_cont(Z)||_(L2(gamma_Z))
       <= C [A^3 + delta A + delta A^2] sqrt(D).

The comparison is uniform over anchored gradients with `0<=Dg<=A I`, `0<A<=1/2`. It requires no modulus of continuity of Dg, no derivative of Dg, no pathwise accuracy of the cheap inner packet, and no root-count-dependent constant. Setting `delta<=A^2` gives the stated third-order mean accuracy. The outer carrier Z must have its actual standard Gaussian law for this integrated estimate.

The actual finite source, split into a D-dimensional genuine-gradient baseline and a 2D-dimensional near-gradient chord, meets the imported fixed-order mean-consumer ports after its recorded private origins are subtracted. Under the actual gradient threshold, near-gradient covariance-gap bound, padding and all imported finite-clock/filter/numerical guards, the complete mean program has conditional target `N(E_(G,H)F_Q(Z,G,H),I)` and error `Lambda A^4(|Z|+sqrt(D))` plus absolute floors. The two completed branches use independent banks at the same retained Z. This gives the source's integrated conditional third-order Gaussian mean-law return.

This audit does not convert the raw F_Q into the full continuous force variable, identify their conditional covariances, retain the integrated G,H bank, supply a strong mean oracle, or establish an unbuffered third-order posterior sampler. The separate covariance and buffered-law join must be audited on their own terms.

Inspected source:

- `../THIRD-ORDER-NESTED-FORCE-MEAN-WITH-C2-DECOUPLING.md`
- Final inspected SHA256: `b056e213ac28c60869727589d3ee6d1b2c53d6666908cdbe1961fa3121161ced`
- The earlier inspected pin `5c41f1619125e23f8ac16f613aab0ebae7bddbe4b52fd1cc1c2f09c4f0710b07` was strengthened with the explicit sufficient normalized near-gradient first/curl guard derived in Section 6 below; the final bytes have been checked.

The independent program `check_nested_force_mean.py` passes **541 assertions**. It does not import an author checker. It tests finite Gaussian Markov ancestry, complete firsts and square-lift curl on a genuinely noncommuting Hessian fixture, captured origins, exact matrix-quadratic force means and translated Euler-decoupling laws, and uniform Hermite multipliers. Its output is `nested_force_mean_checks.json`. Those checks supplement, rather than replace, the continuum and consumer proofs below.

## 1. The conditional past genealogy is the correct one

Fix `0<r<1`, let `q=sqrt(1-r^2)`, and retain Z. For every `tau in (0,1]`, define the residual

    R_tau = X_(r tau) - tau X_r.

The reverse-OU covariance gives

    Cov(R_tau,X_r)=0,
    Cov(R_tau,Z)=r tau - tau r=0,
    Cov(R_tau,R_sigma)
      =[min(tau,sigma)/max(tau,sigma)-tau sigma] I.

The finite vector of all such residuals is jointly Gaussian with `(X_r,Z)`. Therefore it is independent of that entire pair. Because `G=(X_r-rZ)/q`, it is also independent of `(Z,G)`. This is stronger than an informal conditional independence assertion and proves the source's finite representation using a new standard Gaussian bank V:

    X_(r tau)=tau(rZ+qG)+L_tau V,
    L_tau L_tau*= (1-tau^2) I.

The rows L_tau share the residual covariance above. They cannot be replaced by independent node rows. The cheap common-H packet deliberately has another covariance; the proof never identifies it with this past law.

For finite positive weights a_j of mass one, write

    J(x,V)=sum_j a_j g(tau_j x+L_j V).

Then

    0<=D_x J<=A(sum_j a_j tau_j) I,
    ||D_V J||<=A sum_j a_j sqrt(1-tau_j^2)<=A.

The private derivative estimate uses the operator norm of each literal Gaussian row and the positive sum, without an extra square root of the number of rows. For `x=rz+qG`, each force argument has mean `r tau_j z` and covariance `(1-r^2 tau_j^2)I`, so

    ||J(rz+qG,V)||_p<=C_p A(|z|+sqrt(D)).

Midpoint Riemann rules on [0,1] have strictly positive nodes, mass one and exact first moment 1/2. Such rules may therefore be used for the analytical continuous approximation, giving `D_xJ<=A I/2` uniformly. The source's weaker bound A also suffices.

On every closed interval away from zero, the Gaussian OU process is continuous in L2 as a function of tau. The force is A-Lipschitz. Its positive Riemann sums converge in L2 on that interval; the interval `[0,eta]` costs at most `C_p eta A(|z|+sqrt(D))`. This establishes the continuous integral without inventing an X_0. The finite bounds survive the limit. No path approximation rate or growing-root implementation is being billed as an executable algorithm.

## 2. The imported decoupling theorem really reaches Euler

The independently audited conditional-decoupling lemma applies to `Y=G-K(G,V)`, where

    K(G,V)=J(rz+qG,V)/q.

At each exposed z and fixed `r<1`, its parameters may be chosen as

    a=||D_G K||<=A,
    b=Lip_V K<=A/q,
    E_K<=C A(|z|+sqrt(D))/q,
    J_K=(E||D_G K||_HS^2)^(1/2)<=A sqrt(D).

The q in `D_G J=q D_xJ` cancels the displayed division by q. Thus `a<=1/2<1` uniformly even when the other-private first b is large. The lemma only requires a<1; it places no smallness requirement on b or E_K.

The lemma first compares the random endpoint to the flow of `v_K(G)=E_V K(G,V)`. It separately proves the flow-to-Euler error `a exp(a) E_K/2`. Keeping that term yields a direct marginal comparison to `G-v_K(G)`. A convenient combined bound is

    W2(Law(G-K),Law(G-v_K(G)))
       <=C [a E_K + b(E_K+J_K)].

Its source hypotheses are C1 with bounded firsts and linear growth, which hold for each finite approximation. The finite-dimensional theorem is therefore sufficient; no unproved infinite-dimensional Gaussian transport statement is needed.

Multiplying by q and translating by rz gives

    W2(Law(X_r-J(X_r,V)|z),Law(X_r-E_VJ(X_r,V)|z))
       <= C A^2(|z|+sqrt(D))/q.

For clarity, substituting the parameters produces terms of orders

    A^2(|z|+sqrt(D)),
    A^2(|z|+sqrt(D))/q,
    A^2 sqrt(D).

All are dominated by the stated right side because `q<=1`. There is no omitted translation term: the lemma is applied to standard G, with rz kept inside the function K.

For the continuous Markov integral, use the preceding finite approximations, all with the same bound. Their random shifted endpoints converge in W2 by the literal common-process L2 coupling. Their deterministic conditional means converge by conditional Jensen. Passing to the limit in the triangle inequality proves the continuous version. The private-bank dimension does not enter the constant at any stage.

This is a marginal law comparison after V has been integrated. Keeping the old V or H after using the comparison would invalidate its use. The proof does not do so.

## 3. The force comparison gains the required additional A

The two inner objects have analytical means

    E_private I_r | X_r=x = v(x)=integral_0^1 P_tau g(x)d tau,
    E_H H_Q(x,H)=v_Q(x).

Apply the preceding physical decoupling bound separately to I_r and H_Q. Since g itself is A-Lipschitz, a vector expectation of g changes by at most A times the W2 displacement. This gives, at fixed r and z, two errors bounded by

    C A^3(|z|+sqrt(D))/sqrt(1-r^2).

For the middle deterministic-shift comparison, the source's one-clock Hermite theorem gives `||v_Q-v||_L2(gamma)<=delta A sqrt(D)`. The unconditional X_r is standard Gaussian whenever Z and G are independent standard Gaussians. Conditional Jensen therefore yields

    || E_G[g(X_r-v_Q(X_r))-g(X_r-v(X_r)) | Z] ||_2
       <= delta A^2 sqrt(D).

This is an L2 estimate under the actual carrier law, rather than a pointwise estimate uniform over arbitrary z.

Minkowski in r and `||( |Z|+sqrt(D) )||_2<=2sqrt(D)` give an integrated decoupling cost bounded by

    C A^3 sqrt(D) integral_0^1 (1-r^2)^(-1/2)dr
      =(C pi/2) A^3 sqrt(D).

The blow-up is integrable. The endpoint r=1 has measure zero; no executable normalized source is constructed there. In particular the proof does not use the decoupling formula at outer quadrature nodes and then hide a quadrature sum of inverse q factors. It compares the exact r integrals first.

## 4. The outer quadrature is a legitimate single Mehler quadrature

Set, in analysis only,

    Psi_Q(x)=E_H g(x-H_Q(x,H)).

The inner quadrature is the same frozen Q_in for all outer r, and H is independent of the outer Gaussian pair. Thus Psi_Q has no hidden r-dependence. Even though the actual finite program uses the same G,H for all outer terms, linearity of expectation gives exactly

    E_(G,H)F_Q(Z,G,H)=sum_i w_i P_(r_i) Psi_Q(Z).

No node independence is needed for this mean identity. It would not establish an identity of covariance operators.

For independent standard X,H, each `tau_j X+sqrt(1-tau_j^2)H` is standard Gaussian. Positivity and mass one yield `||H_Q||_2<=A sqrt(D)`. Since `|g(y)|<=A|y|`,

    ||Psi_Q||_L2(gamma)<=A(1+A)sqrt(D).

The uniform Hermite-multiplier error bound applies to every vector-valued L2 function, so it applies to Psi_Q without any high-derivative or pointwise clock-analyticity hypothesis. It gives outer error `delta A(1+A)sqrt(D)`. Together with Section 3 this proves the advertised force-mean estimate.

Choosing the inner and outer tolerances at a fixed heat power makes each node count logarithmic squared. The number of raw original VALUES is their product plus the outer terminals. No estimate has changed an L2 error into a pointwise bound, and no inverse power of heat appears in this quadrature count.

## 5. Complete actual firsts, curl, energy and origins

Let

    beta_out=sum_i w_i q_i, q_i=sqrt(1-r_i^2),
    beta_in=sum_j v_j sqrt(1-tau_j^2).

Both beta values lie in `[1/2,sqrt(3)/2]`, by their exact first moments, the inequality `sqrt(1-t^2)>=1-t`, and Jensen. The following slightly weaker bounds would also suffice.

The baseline is

    B(Z,G)=sum_i w_i g(r_i Z+q_i G).

Its G-Jacobian is symmetric positive semidefinite, with norm at most `A beta_out`. Because every q_i is positive, an analytical potential is `sum_i (w_i/q_i) U(r_i Z+q_i G)`; that potential is never queried. Subtracting the recorded B(Z,0) preserves gradient status. The anchored baseline has energy at most `A beta_out sqrt(D)` and private dimension D.

For a single chord E_i, denote the two outer Hessians by B1,B0 and the inner derivatives by Jx,JH as in the source. The literal chain rule is

    D_G E_i=q_i[(B1-B0)-B1 Jx],
    D_H E_i=-B1 JH.

The two B matrices both lie in `[0,A I]`, so `||B1-B0||<=A`. The inner Jx is symmetric, lies in `[0,A I/2]`, and `||JH||<=A beta_in`. Moreover the full private derivative of the inner packet is at most A: each inner row on (G,H) is `(tau_j q_i I,sqrt(1-tau_j^2)I)` and has norm at most one. Consequently the aggregate chord obeys the useful explicit bound

    ||D_(G,H) E||<=A beta_out + A^2.

The same upper bound applies to F_Q. In the square lift `(E,0)` on R^(2D), the skew G-G block is

    sum_i w_i q_i(Jx B1-B1 Jx),

and the G-H block is the actual D_HE. Hence

    ||D(P_G*E)-D(P_G*E)*||
       <=A^2(beta_out+beta_in).

This bound does not commute any local Hessians. The leading Hessian difference contributes no skew because it is symmetric. No claim that this difference is O(A^2) is needed.

At a fixed z, positive weights and the chord inequality give

    ||E(z,G,H)||_p<=C_p A^2(|z|+sqrt(D)).

The private origin is generally nonzero:

    E_origin(z)=sum_i w_i [g(r_i z-H_Q(r_i z,0))-g(r_i z)].

Using `|H_Q(r_i z,0)|<=A r_i |z|/2` gives `|E_origin(z)|<=A^2|z|/4`. The executable anchored source `E-E_origin` has the same energy profile and exactly the same private derivative and curl. It has dimension 2D under the square lift.

For the exposed carrier,

    D_z E_i=r_i[(B1-B0)-B1 Jx],

so `||D_z E||<=A/2+A^2/4`. Its origin has the same O(A) caller bound; subtracting it preserves the required O(A) scale. The baseline and full force have the analogous O(A) caller first. Small VALUE energy does not imply a quadratic caller derivative and is not used that way.

## 6. Guarded positive mean-consumer join

The source states that the imported consumer guards are required. They are substantive restrictions; `A<=1/2` alone is not enough. For the fixed variance shares `v_B=v_E=1/2`, explicit sufficient source radii are

    rho_B=sqrt(2) A beta_out,
    ell_E=sqrt(2) A(beta_out+A).

Require `rho_B<=r_*(k,...)` at the actual fixed-order gradient compiler, `ell_E<=1/4` for the admitted near-gradient covariance gap, and all its finite-clock/filter/precision guards. For example the latter inequality is a fully executable sufficient check based only on frozen A and the clock rule. The normalized full curl is at most `sqrt(2) A^2(beta_out+beta_in)`, so its relative curl is O(A). No measured source energy is used to choose the program.

Use padding `mu=A`. The previously audited near-gradient mean error has the form

    Lambda e_E { A(A+mu)+A^3(1+mu^(-1/2)) }
       + absolute floors.

At the actual anchored energy `e_E<=C A^2(|z|+sqrt(D))`, this is bounded by `Lambda A^4(|z|+sqrt(D))`. The powers and fixed constants coming from share normalization are absorbed in the declared Lambda, not in a hidden inverse-energy choice. The gradient bank with order k>=4 has error at most `Lambda A^4 sqrt(D)` plus its own floors.

Each complete output uses an independent full bank conditional on z and every earlier external caller. Restoring its literal recorded origin changes neither covariance nor W2 error. Independent convolution adds their analytical target means and the strictly positive variances `I/2+I/2=I`. This proves the scoped conditional Gaussian mean-law return. Integrating the conditional error over fresh Gaussian Z, and then using the force-mean estimate, proves the total third-order mean-law grade.

The conditional W2 comparison preserves Z as a caller; it integrates G,H and all complete response/filter/clock banks. It cannot be reused as if it retained an old force value, a private Gaussian observer, or a common-root covariance with another branch.

## 7. Complete work and numerical/restoration ledger

A raw baseline occurrence uses n_out original g VALUES. A raw chord occurrence uses `n_out(n_in+2)` VALUES: every inner ancestor, its shifted outer terminal, and its baseline. The outer and inner roots are shared exactly as stated inside that occurrence. A changed consumer argument requires a new whole occurrence, not reuse of an old partial inner packet. Origins may be reused only as identical captured caller-only graphs with complete semantic keys.

Thus the source's bill

    Q_captured+n_out N_B+n_out(n_in+2)N_E
       + actual known/numerical work

is a safe complete original-VALUE bill, conditional on paying the imported completed-consumer N_B,N_E. Their larger replay/filter/clock Gaussian dimensions are not 2D. A requested first or adjoint sweep applies original HVPs only at recorded VALUE sites, follows all feedback paths, and pays discarded-primal replay. No HVP is a producer VALUE and no saved HVP is differentiated.

For a uniform absolute anchored g-VALUE error epsilon at exact sites, the raw F_Q VALUE error is at most `(1+A)epsilon`; the raw chord error is at most `(2+A)epsilon`. Positive weights of total mass one, rather than a node count, control this propagation. Centering the chord includes the separately propagated origin error. Errors in nodes, known Gaussian rows, weights, anchors, modes and consumer filters retain their separate actual allowances. Approximate quadrature moment identities must not be treated as exact if their numerical debt has not been paid.

At Z=G=H=0, every source value is literally zero under saved-anchor reuse. At fixed nonzero z the privately centered maps have their literal zero because the same recorded origin is subtracted. Independently perturbed repeated evaluations would instead leave a numerical-zero debt. VALUE convergence gives no convergence of Hessians; requested numerical firsts keep the original-first primitive contract and the actual finite graph's first bounds.

The theorem is stated for anchored g. Original finite-mode linear tilts, conditional-mode offsets when this source is rescaled to another conditional problem, physical coordinate factors and numerical caller profiles remain separate inherited floors. They are not absorbed into the intrinsic `A^3 sqrt(D)` force-mean estimate. All versions and numerical parameters must be frozen before caller differentiation. The analytical factor 1/q appears only inside the decoupling proof and has no executed source-width or replay cost.

## 8. Diagnostic limits and conclusion

The checker independently verifies all finite Gaussian covariance identities for several inner bank sizes and outer clocks, including clocks extremely close to one. A smooth globally bounded-Hessian gradient in dimension three has genuinely noncommuting local Hessians; its literal chain-rule checks have maximum directional finite-difference error below `4e-12`, and its full square-lift curl meets the stated bound. Exact matrix quadratics give mean `(B/2-B^2/4)z` for both the continuous target and the finite graph, independently of the cheap inner common-root covariance. Their Euler-decoupling covariance distance also obeys the physical inverse-q estimate, including near-endpoint clocks.

These finite tests do not instantiate the full fixed-order mean compiler and do not prove its imported theorems. The written genealogy, finite-to-continuous argument, force-Lipschitz gain and Hermite operator estimate establish the new analytical step. The prior audited gradient/near-gradient consumer is applied only after its actual source ports, dimensions, origins, error profile and guards have been checked above.

**No mathematical blocker remains for this bounded mean component at the inspected pin.** The separate covariance and cubic Gaussianization components are necessary for any full buffered bridge claim.
