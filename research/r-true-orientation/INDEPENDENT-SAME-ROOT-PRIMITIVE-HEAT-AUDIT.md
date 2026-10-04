# Independent audit of the complete same-root primitive-heat counterexample

Publication copy: nonmathematical context and/or local paths were sanitized. Original and public SHA-256 values are recorded in `INVENTORY.json`; source/audit pins below identify their historical versions, not these edited bytes.

## Source and scoped verdict

Scoped PASS against `research/c-common-input/SAME-ROOT-RIESZ-MARK-PRIMITIVE-HEAT-COUNTEREXAMPLE.md`, SHA-256 `a025c546a6cb8010fb714bac4c6dcdaa97c28175dfd954873fef20387a5ecf29`.

This is a new audit, not a reconstruction of an earlier review. It verifies the FULL actual same-root Riesz mark, a single uniformly convex original potential, dimension-independent source energy, and a nonvanishing one-energy error from independent additive primitive heat. The source has known affine query geometry and explicitly permits a variance-preserving known-row repair. It does not establish membership in the named CW7 weak-proxy batch or in a more restrictive isotropic/scalar-block state grammar.

## 1. Original source and uniform first bound

Put x=(u1,u2,z), T=(|z|^2-m)/sqrt(m), h(T)=exp(-T^2), psi=sin(u1)sin(u2)h(T), and V=|x|^2/4+epsilon0 psi with epsilon0=1e-8. The gradient g=grad V is the single original oracle used by both terminal VALUES.

The displayed Hessian bound 100 is safe uniformly in m. Indeed |grad T|^2<=4(1+|T|), ||D^2T||<=2, |h'|=2|T|exp(-T^2), and |h''|<=(4T^2+2)exp(-T^2). The radial Hessian block is bounded by exp(-t^2)(16t^3+16t^2+12t+8), t=|T|. Bounding each monomial by its elementary Gaussian maximum already gives less than 26. The u block is at most 2; the two off-diagonal blocks have total norm less than the source's bound 16sqrt(2). Thus 100 leaves a large uniform margin. V is uniformly convex and its gradient has first less than one.

With B=I+aD, D=e1e1*, E(Y)=r[g(BY)-g(Y)], and kappa0=ra, the complete source is a literal two-VALUE graph on one shared Gaussian input. Its exact decomposition is

    E=kappa0[.5 D Y+epsilon0 v],
    v=[grad psi(BY)-grad psi(Y)]/a.

The Hessian bound gives ||v||2<=100||DY||2=100. Hence (.5-100epsilon0)kappa0<=e<= (.5+100epsilon0)kappa0 independently of m. The source is odd, centered and exactly zero at the origin. No sqrt(m) source energy is available to absorb a bad estimate.

Differentiating the complete source gives

    Curl E = kappa0[H(BY)D-DH(BY)],

since the second terminal Hessian and the identity part of B are symmetric. Its operator norm is at most 200epsilon0 kappa0. In particular kappa0 is a valid declared curl bound. Both complete source first and actual adjoint are finite.

## 2. The actual marked coefficient and sign

For R_E=D N^-1 E, Gaussian covariance/Riesz integration gives

    O_E = -Sym E[R_E Curl E].

Replacing only the first terminal Hessian by the independently heated H_s(BY)=E_Z H(BY+sZ), while keeping the original R_E(Y) on the SAME Y, changes this target by

    Delta_s=-Sym E[R_E(C_s-C)].

This is the actual marked coefficient being tested. It does not replace R_E by an independent variable or a standalone test vector.

Gaussian Riesz contraction gives the exact decomposition

    R_E=kappa0[.5D+epsilon0 R_v],   ||R_v||L2,HS<=100.

Also ||C_s-C||op<=400epsilon0 kappa0. The contribution of the COMPLETE remaining R_v mark is therefore bounded by 40000epsilon0^2 kappa0^2, by HS times operator norm and Cauchy-Schwarz on the same root. No independence is used.

## 3. Nonzero leading block and dimension limit

The mean curl difference is supported on the first two coordinates:

    E[C_s-C]=kappa0 epsilon0 beta_a
      [exp(-s^2) h_(m,s)-h_m] [[0,-1],[1,0]],
    beta_a=exp(-((1+a)^2+1)/2).

This follows by integrating the off-diagonal terminal Hessian cos(u1)cos(u2)h(T); every other candidate entry has an odd sine in one of the first two coordinates. The independent primitive shield raises each radial bulk variance from 1 to 1+s^2.

Under s^2 sqrt(m)->1, the radial statistic converges to sqrt(2)Z+1 instead of sqrt(2)Z. Since h is bounded, the expectations converge to exp(-1/5)/sqrt(5) and 1/sqrt(5), respectively. This is a genuine nonzero constant change despite s->0.

Multiplication by the exact .5D component of R_E, followed by symmetrization, has HS norm

    epsilon0 kappa0^2 |beta_a[exp(-s^2)h_(m,s)-h_m]|/(2sqrt(2)).

The limiting coefficient is 0.0105438605374851... times epsilon0 kappa0^2. The bound on the full remainder is only 0.0004 times epsilon0 kappa0^2. The eventual lower bound .004epsilon0 kappa0^2 is therefore conservative. Combining it with e comparable to kappa0 proves a fixed positive multiple of kappa0 e; the full same-root mark cannot cancel the effect.

## 4. What the example rules out and what remains legal

It rules out a uniform o(1)kappa0 e calibration based only on shrinking an independent additive primitive heat while freezing the saved root. Sharing that root before the new heat does not preserve its Gaussian covariance.

It does not rule out common-input OU heat applied to the whole graph, or a correctly variance-preserving known B adapter. In this affine-query example, B is known and invertible; B* g(Bx) is itself a genuine-gradient pullback with known readout. A joint known-row construction can preserve the query law. A genuinely nonlinear q(Y) requires its own such construction rather than the same label with an extra independent Gaussian.

The exact common-input orientation-stability lemma remains consistent with the example: it heats the COMPLETE E under one substituted Gaussian input, not the isolated terminal Hessian in the curl factor.

## 5. Diagnostic checks

I reran the NEW checker `c-common-input/check_primitive_heat_counterexample.py`. All 1,000 exact 2-by-2 sign/normalization cases passed. The five radial chi-square integrations through m=10^6 approach the stated limit, and the largest matrix discrepancy was below 7e-18. These numerical checks supplement the CLT/bounded-convergence and full-mark remainder proof. They do not establish named native batch eligibility.
