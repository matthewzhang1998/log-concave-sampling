# Common auxiliary terminal heat: exact native source return and coherent-drift screen

This is a NEW bounded audit of the literal candidate proposed by P on 2026-10-04. It uses the native K3 source c8b20bf4, the six-VALUE ambient lift 64ed4369, and the held coherent-drift family f8363e7d. It is neither a whole-source OU substitution nor a stationary-law theorem on a nonlinear constraint graph.

## 1. The exact modified program

Keep the entire original record W=(S,U,Z), including

    x0=cS+sU, Delta=g_t(S)-g_t(S+epsilon Z),
    d=a Delta, x1=x0+M* d,
    t_i=S+a M g_0(x_i), lambda=0.

Add ONE independent standard V and change only the two constraint coordinates to

    p_i=g_0(x_i)+sigma V,  i=0,1.

The SAME V is used in both terminal queries. The physical S-output of the ambient gradient is exactly

    E_sigma(W,V)=r[g_t(t0+a sigma M V)-g_t(t1+a sigma M V)].       (1)

All original feedback and query aliases are unchanged. In particular the added terminal covariance is a² sigma² M M*, not an arbitrary isotropic primitive heat. The deterministic conditional mean is E_V E_sigma(W,V); the random source in (1) is also considered below, with V included in its complete tape.

## 2. Actual mark, firsts and replay

Assume the terminal primitive has the supplied uniform strong-gradient bounds

    m I <= Dg_t <= L I,  m>0.

For every W,V, strong monotonicity and Lipschitzness give

    m|t0-t1| <= |g_t(t0)-g_t(t1)|,
    |E_sigma(W,V)| <= r L |t0-t1| <= (L/m)|E(W)|.                 (2)

This is a pointwise preservation of the actual UNcentered original E mark, not its larger nominal sqrt(d) bound. Consequently every finite Lp energy is preserved by the same numerical ratio. A centered-only energy statement needs a mean bound; the held coherent-drift source supplies that separately in section 4.

The old native pointwise bound also survives:

    |E_sigma| <= r L² a² epsilon ||M||² |Z|.

The full old-W Jacobian has the same two terminal Jacobians as E, evaluated at shifted queries, and hence the same ordinary O(r) bound under the original native query-path certificate. The new V row is exactly

    D_V E_sigma = r a sigma [H_t(t0+a sigma MV)-H_t(t1+a sigma MV)]M,

and has norm at most 2r a sigma L||M||. With geometry/r/sigma fixed in theta, old caller paths have their original ordinary bounds; source-dependent moving geometry would require its separately supplied rows.

There are still two old terminal evaluations, the original two feedback evaluations, and the original two ancestor evaluations. The new Gaussian row and known M application are charged. Thus the collapsed K3 source uses six original VALUES, plus the actual known-row arithmetic and a new d-dimensional standard block; there is no inverse-sigma replica or callback count. A general finite K uses its actual full old replay plus the two shifted final terminal calls. Numerical errors are passed through the literal readout. No derivative-of-error guarantee is inferred.

## 3. What the ambient companions do

On this modified constraint graph, the other ambient gradient blocks are

    G_Z=-r a epsilon g_t(S+epsilon Z),
    G_x0=-r sigma V,       G_x1=r sigma V,
    G_p0=r[a M* g_t(t0+a sigma MV)-x0],
    G_p1=r[-a M* g_t(t1+a sigma MV)+x1],
    G_d=-rS,              G_lambda=0.                              (3)

Thus the complete ambient gradient is not an E-energy source: its d block alone has L2 norm r sqrt(d), and the new x blocks have norm r sigma sqrt(2d). These are actual outputs of the supplied genuine-gradient lift. Retaining or canceling them requires an explicit program/current; the physical estimate (2) does not assign them the E mark.

The constrained input remains a nonlinear, generally singular function of (W,V), not a standard ambient Gaussian. The zero constraint-output blocks do not change that fact.

There is one exact auxiliary-coordinate gradient available without differentiating an HVP:

    h_W(V) = M* E_sigma(W,V),  sigma>0.                   (4)

It is the gradient of

    r/(a sigma) [V_t(t0+a sigma MV)-V_t(t1+a sigma MV)].

Its active first is at most 2r a sigma L||M||², its captured-W first is O(r||M||), and its VALUE energy is bounded by ||M|| times the physical energy. The displayed inverse width occurs only in the analytical potential, which is not queried; evaluating (4) uses the original gradient VALUES and a known M* action. Recovering the whole physical E_sigma from (4) requires invertible M and the actual known map M^{-*}; for singular M, (4) provides only M*E_sigma and loses transverse output components. A positive projected mean using this inverse must charge its norm, carrier covariance and fill. Even in the invertible case, conditional stationarization of (4) only supplies a mean field depending on W; it has not averaged the original roots or repaired their orientation. If a consumer instead normalizes (4) by 1/(a sigma), then its active first becomes O(r) and its captured-W first O(r/(a sigma)); those readouts must not be mixed.

## 4. The held coherent-drift family passes this heat screen

For f8363e7d, M=e0e0*, r=a, epsilon=a^.9, and

    V_d(x)=|x|²/4 + delta x0(sqrt(1+|x|²)-1) - beta cos x0,
    delta=beta=1/40.

The auxiliary displacement is only a sigma V0 e0. It does not heat the large orthogonal radial bath responsible for the original coherent ancestor drift.

In fact this particular potential has a dimension-uniform Hessian Lipschitz bound. With R=sqrt(1+|x|²),

    ||D²R||op <= 1/R,      ||D³R||op <= 6/R²,
    ||D³(x0 R)||op <= 3/R + 6|x0|/R² <= 9.

Therefore ||D³V_d||op <=9 delta+beta=1/4. Applying the exact mixed finite-difference integral to (1), then the strong monotonicity in (2), gives

    |E_sigma-E| <= (a sigma/(4m)) |V0| |E(W)|,
    ||E_sigma-E||_2 <= (a sigma/(4m)) ||E||_2.                    (5)

The same bound holds for the conditional mean by Jensen. The held proof gives both ||E||_2<=C ra and centered e>=c ra, so the right-hand side is at most C a sigma e with numerical dimension-independent C.

For precision, use the Gaussian output/input-swap orientation identity on the complete tape, padding the physical row by zero on V. If A=I-E-T is the centered swap-defect operator, ||A||_(2->2)<=2, and

    O(F)=Sym E[F_centered tensor A(F)].

Consequently, if delta_F=||F-G||_2,

    ||O(F)-O(G)||HS <= 2(delta_F)(||F_centered||_2+||G_centered||_2).  (6)

The analogous mixed bound with the ORIGINAL intact mark is

    ||Sym E[E_centered tensor A(F-E)]||HS <= 2 e delta_F.          (7)

These are full Gaussian Hilbert bounds; no Hessian convergence or unproved root derivative is used. Equations (5)-(7), together with e comparable to kappa0=ra for this family, give

    ||O(E_sigma)-O(E)||HS <= C a sigma kappa0 e,                  (8)

and the same conclusion for O(E_V E_sigma) and for the mixed original-mark coefficient. Thus EVERY sigma=A^p with fixed p>0 restores this example. The original orientation remains of order kappa0 e, but this auxiliary heat does not change its leading value.

The earlier independent isotropic primitive-heat counterexample a025c546 is not applicable to (1): its dimension-amplifying bath is not in range(M) for this source. Applying it without the actual M M* covariance would be a genealogy error.

## 5. Exact boundary of the conclusion

The native physical signal has a genuine one-energy return under the stated strong terminal hypothesis, ordinary complete first/caller, and finite six-VALUE cost. The held coherent-drift fixture produces no restoration obstruction for this source-aligned heat. The uniform Hessian modulus used in section 4 belongs to that explicit test potential; the native class supplies only bounded Hessians, so (8) is not asserted for every native source.

To obtain a positive stationary repair from the ambient lift, one still needs a law/current theorem for its actual nonlinear constraint graph, or a separately certified source completing (4)'s output/old-root interface. The large companion blocks in (3), any inverse-M readout, and the retained original W-dependent mean must be owned. Independent marginal Gaussian substitutions, or a full-gradient pullback through the nonlinear graph, do not provide this admission.
