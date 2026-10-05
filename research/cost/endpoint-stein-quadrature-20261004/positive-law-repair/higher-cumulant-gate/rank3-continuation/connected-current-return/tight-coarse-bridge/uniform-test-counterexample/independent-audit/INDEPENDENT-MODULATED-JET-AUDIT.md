# Independent audit of the modulated original-gradient covariance witness

2026-10-05. Mathematical audit of the proposed uniform scalar-test counterexample. No sealed source was changed.

## Verdict

PASS for a source-qualified counterexample to a dimension-free or fixed-public-log scalar-test covariance bound for the normalized tight center-center six-tree. The construction also survives adding the other derivative-hit pairs, physical symmetrization, and a positive normalized common-bank average over the complete short Price interval.

Two precision points matter:

1. At a fixed retained original Z=0, the common modulation is asymptotically cos²(sqrt(v) Z0), where v is the actual nonzero center coarse variance. It is not literally cos²(Z0) unless the remaining retained root is also averaged or an equivalent unit-variance convention is used.
2. Retain either the full center injection or a fixed nonvanishing original ancestry injection such as G1. The W2-only injection has squared coefficient sigma² and does not give the claimed growing variance on its own.

This refutes the proposed uniform covariance scalar-test port, not all possible integrated/grouped current methods. It does not construct a positive current producer or settle the recurrence.

## 1. Literal original gradient and Hessian guard

Let n>=3, D=n+1, sigma=n^(-1/2), kappa=1/2, and fix 0<eta<=1/8. Define

    U(x)=A[kappa |x|²/2 + eta sigma² cos(x0) sum_(i=1)^n cos(xi/sigma)],
    g=grad U.

The perturbation K in Dg/A=kappa I+K has entries

    K00=-eta sigma² cos(x0) sum_i cos(xi/sigma),
    Kii=-eta cos(x0) cos(xi/sigma),
    K0i=Ki0=eta sigma sin(x0) sin(xi/sigma),
    Kij=0 for distinct positive i,j.

Its diagonal part has operator norm at most eta; its off-diagonal block has norm at most eta sigma sqrt(n)=eta. Hence

    (1/4)I <= Dg/A <= (3/4)I.

In particular this is a globally bounded PSD original Hessian, and g(0)=0. The potential need not have a bounded gradient; the admitted original source is Lipschitz, and its anchored VALUE source has the required radius.

The same argument works for n=floor(sigma^(-2)), because n sigma²<=1. This variant allows sigma to be taken from the actual finite clock rule rather than requiring an exact reciprocal-integer clock.

## 2. Exact original Gaussian genealogy

Use the source adapter's original q,s,r1,r2,r3 clocks. A transparent representative choice is

    q=s=1/2, r1=r3=1/sqrt(2), r2=sqrt(1-2/n).

Thus the leaf shields are 1/2 and the center shield is sigma. Fix retained Z=0. Write the complete coarse bank as (G0,G1,W1,W2,W3), exactly as in the source adapter. The center query is the scalar row L_n applied coordinatewise to that bank, where

    L_n=(r2 s sqrt(1-q²), r2 sqrt(1-s²), 0, sigma, 0).

Consequently

    v_n=|L_n|²=r2²(1-s²q²)+sigma²
       =15/16-7/(8n) -> v=15/16.

Take the whole-bank Price endpoints G and G_t=tG+sqrt(1-t²)H, with G,H independent standard 5D-dimensional Gaussians. An outer bridge scale equal to any fixed positive constant only changes v and a fixed prefactor; take scale one here. The two center queries X=L_n G and Y=L_n G_t therefore satisfy:

- The pairs (Xi,Yi), i=0,...,n, are independent across physical coordinates.
- Each pair is centered Gaussian with marginal variance v_n and correlation t.
- The low-coordinate pair (X0,Y0) is independent of all n high-coordinate pairs.

This is a consequence of the literal original shared ancestry. It is not a substitution of independent centers or enlarged private shields.

The summed center-center derivative injection is lambda_n=|L_n|²=v_n. If injection histories are kept separate, the G1/G1 history has

    lambda_n=r2²(1-s²) -> 3/4,

and is already sufficient. The W2/W2 history instead has lambda_n=sigma² and should not be used as the growing-variance witness.

The exact numerical interior clocks are not essential: any positive nodes with q,s and both leaf shields in fixed compact interior intervals give v_n bounded away from zero and fixed positive leaf heat. Actual endpoint nodes tending to sigma=0 work with n=floor(sigma^(-2)).

## 3. Exact heat jets, including all cross derivatives

For physical heat P_b f(x)=E f(x+bZ), every trigonometric perturbation term has the same multiplier

    d_b=exp[-b²(1+sigma^(-2))/2].

In particular the exact leaf Hessian is

    P_(1/2)Dg/A=kappa I+d_(1/2)K,
    ||P_(1/2)Dg/A-kappa I||op<=2eta exp[-(n+1)/8].

Let the normalized center C2 analytic tensor be

    T(x)=sigma² P_sigma D³g(x)/A,
    d_n=exp[-(1+sigma²)/2].

Because g is a gradient, T is fully symmetric in its four physical indices. Its complete list of possible nonzero components, up to index permutation, is

    Tiiii= eta d_n cos(x0) cos(xi/sigma),
    T0iii=-eta d_n sigma sin(x0) sin(xi/sigma),
    T00ii= eta d_n sigma² cos(x0) cos(xi/sigma),
    T000i=-eta d_n sigma³ sin(x0) sin(xi/sigma),
    T0000=eta d_n sigma⁴ cos(x0) sum_i cos(xi/sigma).

Components with two distinct positive indices vanish. The quadratic potential contributes nothing to this jet.

These tensors have uniformly bounded proper cuts. One direct verification expands the displayed terms as diagonal tensors with 0-slots inserted. The worst newly degenerate cut is the three-positive-slots versus one-0-slot cut of the sigma-weighted term, whose norm is at most sigma sqrt(n)=1. All other cross terms have at least as much suppression. The diagonal four-tensor has every proper cut at most one. There are only finitely many slot placements. The same argument applies to the normalized C1 center tensor sigma P_sigma D²g/A. Their HS norms are O(sqrt(n)).

Any fixed nonzero native calibration convention multiplies this analytic target by a fixed numerical constant and does not change the conclusion. This audit concerns the analytical coefficient defined by the source adapter; finite calibrated VALUE realizations keep their separate absolute error budgets.

## 4. Exact diagonal six-slot contraction

Replace the four leaves temporarily by kappa I. Let J_n^(0) be the resulting center-center tree, including the chosen injection lambda_n. For an external all-i physical index tuple, only contracted center indices i and 0 can contribute. Exactly,

    (J_n^(0))_(iiiiii)
      =lambda_n kappa⁴ eta² d_n² [
          cos(X0)cos(Y0)cos(Xi/sigma)cos(Yi/sigma)
          +sigma² sin(X0)sin(Y0)sin(Xi/sigma)sin(Yi/sigma)].

Thus none of the unlisted mixed center derivatives can change this diagonal entry. In particular there is no cancellation of the leading term by a cross derivative.

Let

    H_n=n^(-1/2) sum_(i=1)^n e_i^(tensor 6).

Then ||H_n||HS=1, and H_n is symmetric. Full physical symmetrization leaves each all-i component unchanged, so it leaves this test unchanged exactly.

The true leaf matrices differ from kappa I by epsilon_n=O(exp[-(n+1)/8]) in operator norm, uniformly in every query center. The 1|3 cut of each center tensor is O(1). A telescoping replacement of the four leaf matrices therefore changes each displayed scalar entry by O(epsilon_n), uniformly in i and the entire coarse bank. As a result,

    |<H_n,J_n-J_n^(0)>/sqrt(n)|<=C epsilon_n.

The fact that the leaf centers have their actual shared ancestry causes no problem: this is a pointwise estimate.

## 5. Uniform short-Price-interval limit

Let delta²=A_n, with A_n/sigma² ->0; for example A_n=n^(-2) or n^(-3). Take any

    1-2A_n <= t <=1.

Define

    S_n=(1/n) sum_i cos(Xi/sigma)cos(Yi/sigma).

The independent Gaussian coordinate pairs give the exact mean

    m_n(t)=E S_n
      =[exp(-v_n(1-t)/sigma²)+exp(-v_n(1+t)/sigma²)]/2.

It tends to 1/2 uniformly over the entire declared Price interval, since A_n/sigma² ->0 and v_n stays positive. Also

    E|S_n-m_n(t)|²<=1/n.

This empirical fluctuation is independent of (X0,Y0). The second sine average in Section 4 is bounded in absolute value by sigma², so it tends uniformly to zero.

Using a common representation X0=sqrt(v_n)Z0 and Y0=sqrt(v_n)[tZ0+sqrt(1-t²)Z1],

    cos(X0)cos(Y0) -> cos²(sqrt(v)Z0)

in L², uniformly over the Price interval. Hence, for either lambda_n=v_n or a fixed nonvanishing ancestry injection,

    <H_n,J_n>/sqrt(n) -> c cos²(sqrt(v)Z0) in L²,
    c=lambda_* kappa⁴ eta² exp(-1)/2 >0.

A quantitative bound is O(n^(-1/2)+A_n/sigma²+sqrt(A_n)+1/n+epsilon_n), with fixed constants. This uniform convergence also permits an arbitrary positive normalized Price-clock quadrature or an integral over the interval when the same whole-bank G,H family is used: the same limit remains. A randomly selected independent Price clock is allowed as well. If different quadrature histories are instead executed with independent banks and subsequently summed, the resulting variance must retain their squared weights; that is a distinct coefficient and is not identified with this common-bank average. No endpoint atom t=1 is required.

Finally,

    Var(cos²(sqrt(v)Z0))=(1-exp(-4v))²/8>0.

Therefore

    Var(<H_n,J_n>)/n -> c²(1-exp(-4v))²/8>0.

Since ||H_n||HS=1, the proposed dimension-free variance bound fails. Taking A as a fixed negative power of n makes n dominate every fixed polynomial in the public logarithms.

## 6. Pooling derivative-hit histories cannot rescue this witness

Let F_n denote the complete normalized cubic node. Its center C1 tensor has fixed proper cuts. At a leaf, the normalized derivative sigma D(P_(1/2)Dg/A) has fixed proper cuts times exp[-(n+1)/8]; this follows either from the explicit third derivatives or from the same diagonal tensor expansion.

Thus sigma D_QF_n equals its center-hit contribution plus an exponentially small tensor, including all actual ancestry rows. Contracting the two derivatives and testing the all-i six-slot diagonal shows that the fully pooled normalized coefficient

    sigma² sum_a (partial_a F_n)(G) tensor (partial_a F_n)(G_t)

has the same limit. Leaf-hit terms cannot cancel an order-one normalized common modulation. Symmetrization before differentiating or after forming the six-slot coefficient likewise preserves the all-i test.

This statement is about the sum of derivative-hit histories at the chosen original clock node. It does not claim a lower bound for an arbitrarily weighted sum across different original clock nodes or for a new grouped common-heat coefficient that changes the target itself.

## 7. Clock/admissibility and scope checks

- Choose a center endpoint panel above the inherited cutoff with sigma ->0 and A/sigma² ->0. For a cutoff sigma_min=A^b, any fixed exponent 0<a<min(b,1/2) gives the relevant scale sigma comparable to A^a. Endpoint dyadic positive quadrature nodes supply such scales. Use n=floor(sigma^(-2)).
- Every fixed original interior clock can be replaced by an actual quadrature node in a fixed compact interior interval. The conclusion only needs fixed positive bounds, not the illustrative exact fractions.
- All private original shields remain unchanged. Every original root and both entire Price endpoint banks are retained.
- The counterexample is an analytic original-gradient coefficient and requires no derivative oracle or external tensor source. The source guard follows from the explicit Hessian sandwich.
- Fixed-order/public-log native calibration constants and separately chosen arbitrarily small absolute floors cannot create a uniform scalar-test theorem for an analytic target whose variance grows linearly with dimension.
- The result rules out the uniform proposed port. It does not prove that the previously generic sigma^(-2) estimate is always sharp under all other normalizations, nor does it rule out a nonuniform integrated clock-energy substitute with additional structure.

## Sources inspected

- ../../COVARIANCE-TEST-PORT-REDUCTION.md
- ../../TIGHT-BRIDGE-LEDGER-AND-REROOT-GATE.md
- ../../../../EXACT-FIVE-CLOCK-NATIVE-SOURCE-ADAPTER.md
- ../../../../grouped-kernel-response/PAIRED-ROTATION-COVARIANCE-ALL-CUT-LEMMA.md
- ../../../../grouped-kernel-response/independent-audit/INDEPENDENT-FRACTIONAL-B-PRICE-AUDIT.md
