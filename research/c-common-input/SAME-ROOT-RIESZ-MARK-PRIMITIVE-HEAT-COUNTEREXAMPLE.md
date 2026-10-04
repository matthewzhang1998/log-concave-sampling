# A same-root, one-energy test of independent primitive heat

## 0. Precise shortcut tested

Freezing a common root Y before invoking a primitive gradient pair preserves its saved centers, but the pair with source g(q(Y)+sZ) averages an EXTRA primitive Gaussian heat. This is not the whole-input OU heat of the original finite graph. The following two-terminal source shows that this distinction can cost a full kappa e, even when the marked factor is the ACTUAL same-root Gaussian Riesz Jacobian of the complete source.

The result rejects an unqualified replacement of a terminal Hessian by independent additive primitive heat. It does not reject a known variance-preserving pullback/adapter, a source-specific correlated extraction, or the proved whole-source OU restoration port. The displayed state is affine and its exact known Gaussian geometry is available, so this is deliberately a test of the proposed extraction step, not an impossibility for this source.

## 1. One original uniformly convex potential and a small complete VALUE source

Let n=m+2 and write x=(u1,u2,z), z in R^m. Define

    T_m(z)=(|z|²-m)/sqrt(m),
    psi_m(x)=sin(u1)sin(u2) exp(-T_m(z)²),
    V_m(x)=|x|²/4+epsilon0 psi_m(x),
    g_m=grad V_m,
    epsilon0=10^(-8),
    D=e1 e1*.

There is a uniform bound ||D²psi_m||op<=100 for all m>=1 and all x. One direct check uses h(t)=exp(-t²), |grad T_m|²<=4(1+|T_m|), ||D²T_m||<=2, and the2-by-m Hessian blocks. The u block is at most2, twice the cross-block norm is at most16sqrt(2), and the z block is at most68. Thus the total is less than100.

Consequently

    (.5-100epsilon0)I <= Dg_m <= (.5+100epsilon0)I,

so every original gradient query uses a uniformly convex, unit-first potential, with constants independent of dimension.

Let 0<a<=1, r>0, B=I+aD, and execute the finite source

    E(Y)=r[g_m(BY)-g_m(Y)],          Y standard,
    kappa0=r a.                                         (1)

Both terminal calls use the SAME original potential and the SAME Gaussian root. The state perturbation is the known small affine map aDY. The source is odd and centered. Its actual energy obeys

    (.5-100epsilon0)kappa0 <= e=||E||2
                                   <=(.5+100epsilon0)kappa0.           (2)

It has full first at most2r, and curl

    C_E(Y)=kappa0[H_m(BY)D-DH_m(BY)],
    H_m=Dg_m,                                           (3)

with operator norm at most200epsilon0 kappa0, in particular at most kappa0. The second terminal in(1) is a genuine gradient and contributes no curl. E(0)=0 exactly.

This is a literal terminal-gradient VALUE graph with known affine state, complete first and adjoint. It is not asserted to be the named CW7 weak-proxy batch or to satisfy a more restrictive scalar-block-only state grammar without its own embedding certificate.

## 2. Actual same-root Riesz mark and the replacement under test

Let N be the Gaussian number operator and

    R_E(Y)=D N^(-1) E(Y).

The exact true orientation is

    O_E=-Sym E_Y[R_E(Y) C_E(Y)].                          (4)

This identity is only an analytical coefficient identity; R_E is not an emitted differentiable oracle.

The freeze-Y primitive-pair step replaces H_m(BY) in(3) by

    H_(m,s)(BY)=E_Z H_m(BY+sZ),

with a NEW independent standard Z and nonzero primitive width s. Define

    C_(E,s)(Y)=kappa0[H_(m,s)(BY)D-DH_(m,s)(BY)].

Keep the ACTUAL same-root R_E(Y), without changing or independently resampling it. The coefficient calibration error is exactly

    Delta_s=-Sym E_Y[R_E(Y)(C_(E,s)(Y)-C_E(Y))].          (5)

For a->0, s->0, m->infinity with s²sqrt(m)->1, this error has

    ||Delta_s||HS >= c_* kappa0 e                       (6)

for one fixed positive numerical c_*, once the parameters are sufficiently small/large. In particular no uniform bound o(1)kappa0 e follows just because s shrinks.

## 3. The nonzero leading block

Write h_m=E exp(-T_m(Z)²), and

    h_(m,s)=E exp(-[((1+s²)|Z|²-m)/sqrt(m)]²),
    beta_a=exp(-[(1+a)²+1]/2).

A direct Gaussian calculation gives the COMPLETE mean of the curl replacement error: it is supported on the first two coordinates and equals

    E[C_(E,s)-C_E]
      =kappa0 epsilon0 beta_a [exp(-s²)h_(m,s)-h_m]
                            [[0,-1],[1,0]].             (7)

The quadratic Hessian commutes with D. Every other expected entry vanishes by a sine oddness in u1 or u2. The added primitive heat inflates the m bulk-coordinate variance from1 to1+s²; keeping Y shared does not undo this.

The central limit theorem and bounded convergence give

    h_m -> 1/sqrt(5),
    h_(m,s) -> exp(-1/5)/sqrt(5)

along s²sqrt(m)->1. Put

    delta_h=(exp(-1/5)-1)/sqrt(5) != 0.

Therefore the bracket in(7) tends to delta_h, while beta_a tends to exp(-1).

## 4. Why the WHOLE actual mark cannot cancel this effect

The complete source has the exact decomposition

    E(Y)=kappa0[.5 D Y+epsilon0 v_(a,m)(Y)],
    v_(a,m)(Y)=[grad psi_m(BY)-grad psi_m(Y)]/a.

The uniform Hessian bound gives ||v_(a,m)||2<=100. Gaussian Riesz contraction gives

    R_E=kappa0[.5 D+epsilon0 R_v],
    ||R_v||_(L2;HS)<=||v_(a,m)||2<=100.                  (8)

Also

    ||C_(E,s)-C_E||op<=400epsilon0 kappa0.                (9)

Thus the contribution of epsilon0 R_v to(5) is bounded in HS by

    40000epsilon0² kappa0².                             (10)

This is a bound on the actual same-root product; no independence is used.

The exact .5D part of(8), combined with(7), has Hilbert--Schmidt norm

    kappa0² epsilon0 |beta_a[exp(-s²)h_(m,s)-h_m]|
                                     /(2sqrt(2)).       (11)

Its limiting coefficient divided by epsilon0 kappa0² is

    exp(-1)|delta_h|/(2sqrt(2)) > .0105.

For epsilon0=1e-8, the remainder(10) is less than .0004epsilon0 kappa0². Equations(10)--(11) therefore give, for example,

    ||Delta_s||HS >= .004 epsilon0 kappa0²

eventually. Together with(2) this proves(6), with a fixed positive c_*.

The mark in this proof is the full R_E, not its constant approximation. The constant component is used to prove a lower bound, and the entire remaining marked component is controlled at the SAME root. The source energy is O(kappa0), independent of m, so an ambient sqrt(m) energy cannot absorb the defect.

## 5. What can escape the counterexample

Whole-input OU heat substitutes cY+sZ into BOTH complete terminal queries in(1), retaining the state B and all source aliases. It preserves the original Gaussian input law and is covered by the separate bound

    ||O_E-O_(P_t E)||HS<=kappa e min(1,sqrt(2t)).

The independent primitive heat in section2 is a different operation.

Since B is KNOWN and affine in this fixture, an exact variance-preserving Gaussian pullback can also be used: heat the common original input before applying B, or use the genuine-gradient pullback B* g_m(Bx) with its actual known rows and any required known inverse/readout. Such a construction must be priced and its joint input law retained, but is not excluded here. The fully affine-query case can be easier than a genuinely nonlinear q(Y)=PY+a g0(CY).

The lesson for that nonlinear case is precise. Conditioning on Y before adding primitive-pair heat is insufficient by itself. One needs a common-input genealogy identity or a genuinely variance-preserving joint adapter, together with whole-source first/energy control. A root Bessel estimate for a standalone terminal Hessian, or the fact that both source labels were saved from the same Y, does not supply that identity.
