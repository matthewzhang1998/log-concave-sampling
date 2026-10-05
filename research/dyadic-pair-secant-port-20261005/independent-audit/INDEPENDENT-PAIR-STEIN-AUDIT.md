# Audit of the six-VALUE mixed-secant Gaussianization

2026-10-05. Analytical source and law certificate. This proves the source's exact pair-Hoeffding Gram and its smoothed iid Gaussian approximation. It does not identify a finite sum of these sources with the complete true-history covariance or assert logarithmic copy count.

## 1. Literal source and covariance

Let g: R^D -> R^D be a gradient with symmetric Jacobian H satisfying 0<=H<=A I. Fix captured a,b,c,sigma,tau, with c>=0 and sigma,tau>=0. Take four independent standard D-dimensional Gaussian roots xi,xi',zeta,zeta'. Put

    x=a+sigma xi,      x'=a+sigma xi',
    z=b+tau zeta,      z'=b+tau zeta',
    h(x,z)=g(x-c g(z)),
    Psi=(1/2)[h(x,z)-h(x',z)-h(x,z')+h(x',z')].         (1)

Two inner g VALUES and four outer g VALUES suffice, with same-caller reuse of each inner value. All displayed parameters and numerical scalar versions are frozen first.

Let d(x,z) be the two-variable Hoeffding component of h with respect to the independent Gaussian laws of x and z. The singleton and constant terms cancel in (1). The four remaining d terms are pairwise orthogonal because d is centered in each argument. Hence

    E Psi=0,
    Cov(Psi)=E[d(x,z)d(x,z)^*].                         (2)

This is exact equality of covariance, not equality of the laws of Psi and d.

Conditional on the two outer roots (xi,xi'), swapping zeta and zeta' negates Psi, so its conditional mean is zero. Differentiating only the inner roots gives, for example,

    D_zeta Psi
       =(c tau/2)[H(x'-c g(z))-H(x-c g(z))] H(z).

Since two matrices in [0,A I] differ in operator norm by at most A,

    ||D_(zeta,zeta') Psi||op <= L := c A^2 tau/sqrt(2).

Conditional Gaussian Poincare and the zero conditional mean imply the dimension-safe bound

    0<=Sigma:=Cov(Psi)<=L^2 I = c^2 A^4 tau^2 I/2.       (3)

The full first derivative remains separately paid. For instance,

    ||D_(xi,xi') Psi||op <= A sigma/sqrt(2).

No small complete-first conclusion follows from the small inner Stein radius.

## 2. A conditional-inner Stein kernel

For fixed outer roots o=(xi,xi'), write f_o(Z)=Psi, where Z=(zeta,zeta') is standard Gaussian in R^(2D). It is centered and has ||D f_o||op<=L. Let P_t be the standard Gaussian OU semigroup on Z and set

    U_o(Z)=integral_0^infinity e^(-t) P_t(D f_o)(Z) dt.

Then ||U_o||op<=L. The Gaussian resolvent integration-by-parts identity gives, for smooth scalar phi,

    E[f_(o,i)(Z) phi(f_o(Z))]
       =sum_j E[(U_o D f_o^*)_(ij) partial_j phi(f_o(Z))].

Consequently

    T(Psi)=E[U_o(Z) D f_o(Z)^* | Psi]

is an exact Stein kernel after also averaging over the random outer roots. It may be nonsymmetric; symmetry is unnecessary. Uniformly,

    ||T||op<=L^2,       E T=Sigma,
    E||T-Sigma||HS^2 <= D L^4.                         (4)

The last inequality is the Hilbert-space variance identity combined with ||T||HS^2<=D L^4. The construction differentiates only the inner Gaussian roots. Centering conditional on all outer roots is what permits their derivative costs to be absent from (4).

## 3. iid copies and a Gaussian buffer

Take K independent complete copies Psi_1,...,Psi_K at the same captured parameters, and an independent standard D-dimensional Gaussian G. Define

    S_K=K^(-1/2) sum_i Psi_i,
    Y_K=S_K+sqrt(v) G,       v>0,
    Q=Sigma+v I.

A Stein kernel of S_K is the conditional expectation of K^(-1) sum_i T_i given S_K. Independence and (4) yield

    E||T_(S_K)-Sigma||HS^2 <= D L^4/K.

A Stein kernel of Y_K is v I+E[T_(S_K)|Y_K], so

    E||T_(Y_K)-Q||HS^2 <= D L^4/K.                     (5)

For any centered Y with target covariance Q>0 and Stein kernel T_Y, the anisotropic Gaussian Stein bound is

    W2(Law(Y),N(0,Q))
       <= [E||(T_Y-Q)Q^(-1/2)||HS^2]^(1/2).

For completeness, interpolate Y_t=e^(-t)Y+sqrt(1-e^(-2t))Q^(1/2)Z with independent standard Z. Integration by parts produces a continuity-equation velocity with L2 norm bounded by

    e^(-2t)/sqrt(1-e^(-2t))
       * [E||(T_Y-Q)Q^(-1/2)||HS^2]^(1/2).

Its scalar time factor integrates to one. The Wasserstein path-length bound proves the assertion. Applying Q>=v I to (5) therefore gives

    W2(Law(Y_K),N(0,Sigma+v I))
       <= L^2 sqrt(D/(v K))
       = (c^2 A^4 tau^2/2) sqrt(D/(v K)).               (6)

The source is a finite original-VALUE graph with 6K VALUES before any additional caller replay and (4K+1)D complete private Gaussian coordinates. Equation (6) has algebraic K^(-1/2) convergence. It does not make K logarithmic in the requested error.

## 4. Caller and source qualifications

- Equations (2)-(6) hold pointwise in fixed captured parameters, hence uniformly when those parameters are later integrated under an actual caller law, provided the stated uniform envelopes apply.
- All four root banks per copy are consumed for the marginal-law comparison. An external caller cannot inspect xi,xi',zeta,zeta' afterward unless a separate joint/conditional comparison has been proved.
- The variance-v Gaussian is actual independent smoothing noise in the output. The v-dependent gain does not establish a comparison conditional on that G or a coupling that must preserve the same caller-readable G exactly.
- If one conditions on the outer roots instead, one obtains their conditional covariance as the Gaussian target, rather than the unconditional Sigma in (6).
- Every changed captured parameter or source root replays every affected inner VALUE, outer shifted VALUE, and same-caller record. Small covariance does not erase these complete-first and replay paths.
- This source has the pair-Hoeffding Gram of the literal two-site shifted composition h. It becomes a source for the true coherent history interaction only after a separate exact representation or approximation theorem, with all temporal cross terms and coherent ancestors retained.

No simulated path, history grid, external message, or upload was used.
