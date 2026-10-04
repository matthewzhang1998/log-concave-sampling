# Exact two-clock current and a uniform-C2 analytical reference

2026-10-04. A reusable full-current identity, an exact remainder reference, and the finite native source it calls for. The stationary continuous process below is used only in analysis. It is not an executed infinite oracle, an unpriced path quadrature, or a restart of the paused unsmoothed path-discretization route.

## 1. Result and precise finite gate

The true second-order tilt current has a positive two-clock Gaussian representation which preserves every bilinear Hermite interaction:

    X_r=rZ+sqrt(1-r^2)G,
    X_(r,tau)=tau X_r+sqrt(1-tau^2)H,

    C2(phi)=integral_0^1 dr integral_0^1 d tau
       E[ Dg0(X_r)g0(X_(r,tau)) dot grad phi(Z)
          +r Sym(g0(X_r) tensor g0(X_(r,tau))):Hess phi(Z)]. (1)

Z,G,H are independent standard Gaussian vectors. The inner point uses the same actual X_r. The formula is analytical; Dg0 is not a permitted VALUE producer leaf. There is no missing tau weight.

The rank-two term is exactly the rank-two current of a Gaussian **Markov OU** force integral, with covariance min(r,s)/max(r,s). The rank-one term is its nested VALUE feedback. An exact stationary Langevin resolvent supplies a uniform O(A^3 sqrt(D)) law remainder for two continuous substitutions, using only Lipschitz/monotonicity and no Hessian modulus.

What is not yet proved is a finite polylogarithmic positive law/current realization of those Markov two- and three-root correlations with the same remainder. A finite program and all its raw ports are explicit below. Its quadrature error must be proved in the relevant one-energy operator/current norm. Matching finitely many scalar moments, or assuming that the clocks see a smooth force, does not do this.

## 2. Exact second-order current under the C2 sandwich

Let U=A U0, g=A g0=grad U, g(0)=0, and 0<=Dg<=A I, with 0<A<=1/2. For a smooth compactly supported test phi, write mu_A proportional to exp(-|x|^2/2-AU0(x)). Its amplitude expansion is

    E_muA phi=E_gamma phi-A Cov_gamma(U0,phi)
                      +A^2 C2(phi)+o(A^2),
    C2(phi)=(1/2)Cov_gamma(phi,(U0-EU0)^2).              (2)

Here U0 has at most quadratic growth. Gaussian integrability and differentiation under the integral establish (2); the fixed-potential coefficient identity is not itself a uniform third-order remainder claim.

The Gaussian covariance identity gives

    Cov(F,G)=integral_0^1 E[grad F(X) dot
                  grad G(tau X+sqrt(1-tau^2)N)]d tau.   (3)

Apply it first to (U0-EU0)^2/2 and phi. This yields

    C2(phi)=integral_0^1 E[(U0(X_r)-EU0)
                              g0(X_r) dot grad phi(Z)]dr.

At fixed r, view (Z,G) as one Gaussian vector W and X_r=R_r W, R_r=(rI,sqrt(1-r^2)I). Apply (3) to U0(R_r W) and g0(R_r W) dot grad phi(PW), P=(I,0). Its derivative has two terms:

    R_r* Dg0(X_r) grad phi(Z),
    P* Hess phi(Z) g0(X_r).

Since R_r R_r*=I and R_r P*=rI, their contractions give exactly (1). One may equivalently condition W on its correlated copy; the old projected root then becomes tau X_r+sqrt(1-tau^2)H with a fresh standard H. This proves the displayed literal correlation, not an independent force replacement.

The proof uses only Dg0, which is continuous and bounded under C2, plus derivatives of the test. Approximation by smooth functions and dominated Gaussian moments extend the identity to the Sobolev/compact-test closure needed here. It neither differentiates Dg0 nor requires U0 to be C3.

Polarizing U0=U_i+U_j gives every second-order bilinear interaction. The independent symbolic checker verifies 52 exact Wick identities for potential Hermite pairs (1,2),(2,3),(2,5),(3,4),(3,3),(4,4),(4,5),(5,5), with carrier ranks through eight. For example pair (2,3), carrier rank three, gives rank-one kernel 36 r^2 tau^2+72 r^2 tau and rank-two kernel 72 r^2 tau^2+72 r^2 tau; their double integral is 36, the exact target coefficient. These polynomial tests diagnose the two-variable geometry; they are not an operator-norm quadrature proof for arbitrary C2 sources.

## 3. Why the common-two-root packet has the wrong next genealogy

The original packet uses X_r=rZ+sqrt(1-r^2)G with the **same G at all r**. For r,s its covariance is

    r s+sqrt(1-r^2)sqrt(1-s^2).

In the true two-clock current, if s=r tau<=r then

    Cov(X_r,X_(r,tau))=tau=s/r,
    Cov(Z,X_r)=r, Cov(Z,X_(r,tau))=s.                   (4)

These two root geometries agree in the conditional first mean but disagree in the joint current. The same-root chord and two-stage cubic repair expose some consequences of that discrepancy; neither changes the whole covariance kernel to (4).

Let {X_r:0<r<=1} be a stationary Gaussian OU process in logarithmic time, with

    Cov(X_r,X_s)=[min(r,s)/max(r,s)]I, X_1=Z.            (5)

For each fixed ordered pair s=r tau, its triple (Z,X_r,X_s) has precisely the representation in (1). Put h=integral_0^1 g0(X_r)dr as an L2 integral. Then

    (1/2)E[h tensor h:Hess phi(Z)]
       =integral_0^1 dr integral_0^1 d tau
          r E[Sym(g0(X_r) tensor g0(X_(r,tau)))
                                                     :Hess phi(Z)]. (6)

This is just the ordered triangular domain s<=r with ds=r d tau, including the shared root law (5). The exact rank-one part is the amplitude coefficient of the nested feedback

    Z-integral_0^1 g(X_r-integral_0^1 g(X_(r,tau))d tau)dr. (7)

Thus the continuous expression (7) matches all of (1), not only one covariance or cumulant rank. It is still an analytical object until its finite realization is supplied.

## 4. Uniform C2 remainder by an exact stationary resolvent

This section supplies a comparison reference only. It executes no continuum force evaluations.

On a two-sided Brownian probability space, let X_t be stationary Gaussian OU,

    dX_t=-X_t dt+sqrt(2)dB_t,

and couple the stationary Langevin process Y_t for the target by the same noise,

    dY_t=-(Y_t+g(Y_t))dt+sqrt(2)dB_t.

The drift is globally Lipschitz and strongly monotone. Starting finite-past solutions with the common noise and taking the stationary L2 limit supplies a stationary joint coupling. The target marginal is mu_U. Its anchored monotonicity gives E|Y_t|^2<=D by target integration by parts.

Variation of constants and the vanishing stationary L2 boundary term imply the exact identity

    Y_t=X_t-integral_0^infinity e^(-s)g(Y_(t-s))ds.      (8)

Consequently ||X_t-Y_t||2<=A sqrt(D). Define comparison processes only in analysis,

    P0_t=X_t,
    P_(j+1),t=X_t-integral_0^infinity e^(-s)g(P_j,(t-s))ds.

Stationarity, Minkowski, and Lip(g)<=A give the exact contraction

    sup_t ||P_(j+1),t-Y_t||2
       <=A sup_t ||P_j,t-Y_t||2,
    ||P_j,t-Y_t||2<=A^(j+1)sqrt(D).                    (9)

For j=2, set r=e^(-s), tau=e^(-u). The expression P2_0 is exactly (7), with Gaussian joint law (5). Therefore

    W2(Law(P2_0),mu_U)<=A^3 sqrt(D).                   (10)

This is uniform over the C2 Hessian sandwich (indeed the contraction uses less). No expansion of the Hessian, high derivative, or dimension-dependent displacement guard enters (9). More analytical substitutions have the displayed contraction, but this is not a finite-cost recurrence: each continuous source has an unimplemented integral and Gaussian path.

For the anchored conditional source f(y)=s[g(x+s y)-g(x)], alpha=As^2, the same reference gives normalized error alpha^3 sqrt(D), physical conditional error A^3 s^7 sqrt(D), and after exact affine reverse-OU gluing

    C (Delta/t) A^3 s^5 sqrt(D)

plus the existing finite-mode/tilt/numerical floors. These are analytical target rates for a finite compiler to realize. They are not attributed to an unproved finite quadrature.

## 5. Explicit finite native graph, with the unresolved law comparison exposed

Freeze finite positive quadratures (w_i,r_i) and (v_j,tau_j), with both weight sums one and 0<r_i,tau_j<=1. Zero clock nodes are excluded from this finite Gram because they have no finite logarithmic OU time; the admitted interior dyadic rules satisfy the condition. Form the finite set of times

    T={1} union {r_i} union {r_i tau_j}.

Identify equal times under exact complete keys. Draw a Gaussian vector {X_t:t in T} with the known PSD scalar Gram

    Gamma_(t,s)=min(t,s)/max(t,s),

tensor I_D. It can be drawn exactly in the Gaussian real-arithmetic convention by successive known OU transitions at sorted times; a finite numerical version pays its explicit covariance-factor precision. No potential data are used to choose Gamma. These nodes cannot be sampled independently or replaced by a common two-root cosine bridge.

Execute

    I_i=sum_j v_j g(X_(r_i tau_j)),
    Y_Q=X_1-sum_i w_i g(X_(r_i)-I_i).                  (11)

The VALUE count is at most n_r n_tau+n_r, plus captured original/conditional mode and anchor sites. Equal original sites may be reused only under their full semantic keys. A first/adjoint sweep uses those original HVP sites and the actual feedback through I_i; replaying a discarded primal pays the full graph. No HVP is a producer VALUE.

The root dimension is D times the number of distinct times in T, not 2D. If n_r,n_tau are public-log polynomials this is also a public-log root count. Known Gaussian covariance arithmetic is separately paid. The actual Gaussian row at each time has norm one. For any positive weights a_l summing to one, the stack of rows weighted by sqrt(a_l) has operator norm at most one, because its scalar Gram has trace one. This yields dimension-free raw first bounds

    Lip(Y_Q-X_1)<=A(1+A),
    ||Y_Q-X_1||p<=C_p A(1+A)sqrt(D),
    Y_Q(owned roots zero)=0.                           (12)

For the first estimate, take the outer weighted L2 stack and the inner stack weighted by w_i v_j; both have scalar Gram trace one. The chain through the inner force contributes A^2. No sqrt(number of nodes) is silently absorbed. At original conditional scaling, the same positive weights preserve the finite-mode caller bound and give private residual first CAs^3(1+As^2); captured caller origins and original physical coordinates are retained exactly as in the packet service.

Every point is a literal original-gradient VALUE site. Absolute numerical errors propagate with total weight one and feedback factor 1+A. All clocks, Gram versions, precision allocations, and source versions are frozen before differentiation. The zero is literal with recorded-anchor reuse. The actual original caller and every external label remain outside this fresh Gaussian draw; no joint law with discarded private nodes is assumed after a comparison.

The missing theorem is now specific:

    W2(Law(Y_Q|caller),Law(P2_0|caller))
       <=C A^3 sqrt(D)+allocated floors

with n_r,n_tau only public-logarithmic and all native ports above preserved, possibly after a finite covariance/moment/current calibration. Raw conditional-mean quadrature accuracy is insufficient for this statement.

## 6. Quadratic covariance screen for the finite graph

For g(x)=Bx, assume the rules have first moments sum w_i r_i=sum v_j tau_j=1/2. Then

    Y_Q=Z-B H_Q+B^2 J_Q,
    H_Q=sum_i w_i X_(r_i),
    J_Q=sum_(i,j)w_i v_j X_(r_i tau_j),
    Cov(Z,H_Q)=I/2, Cov(Z,J_Q)=I/4.

Writing

    S_Q=sum_(i,l)w_i w_l min(r_i,r_l)/max(r_i,r_l),

the actual order-B^2 covariance coefficient is

    S_Q+1/2.                                          (13)

The continuous coefficient is one because S_cont=1/2. Thus an uncalibrated finite rule must explicitly pay S_Q-1/2 even on quadratics. This is not the shared-two-root beta_Q^2 coefficient. A positive finite covariance or VALUE-polynomial repair may remove (13), but its higher-Hermite/root descendants remain to be priced. A complementary finite-clock analysis is investigating calibrated representations and the high-chaos diagonal constraint; no broad obstruction to all positive current constructions follows from (13).

This screen prevents the analytical contraction (9) from being misapplied to an unrelated discretized law. It is not an invitation to resume an ordinary unsmoothed path grid: the active finite goal is law-level Markov current compression/calibration with an actual one-energy C2 norm proof.

## 7. Exact local sources requested by a current-based alternative

Instead of realizing the whole Markov path, (1) offers a two-clock local target at each (r,tau), conditional on the carrier Z. Its native VALUE record is

    X=rZ+sqrt(1-r^2)G,
    B=tau X+sqrt(1-tau^2)H,
    F=g(X), G1=g(B), F_shift=g(X-G1).

The chord F_shift-F is an actual local source with one marked energy O(A^2 sqrt(D)) under Gaussian averaging, full first O(A), captured-carrier first O(A), and exact anchored zero. Its analytical leading current is -Dg(X)g(B); a uniform third-order *finite* restoration cannot be concluded by pointwise Taylor under C2. The rank-two target is the actual full product r Sym(F tensor G1), and the old base packet contributes its whole common-root half-square. Neither is replaced by a covariance trace.

In particular, conditional covariance alone is not the whole product:

    E[Sym(F tensor G1)|Z]
       =Sym Cov(F,G1|Z)
          +Sym(E[F|Z] tensor E[G1|Z]).                  (14)

A covariance-action consumer must retain the second term through a legitimate joint mean/pair construction. The means in (14) are analytical, not strong oracles. Every mixed source must use the same X in B. Independently resampling that ancestor changes the exact target. A finite positive pair/current consumer that realizes the whole difference of (1) and the base packet's rank-two current, with a C2-safe restoration tied to (8)-(10), would close this specific next-order gate.

The all-rank Hermite consumer is not invoked on a formal tensor here. Its missing input is precisely this finite original-VALUE rank-one/rank-two producer with its actual root, caller, origin, first, and one-energy restoration return. Equations (1), (8), (11), and (14) identify that input completely enough to test proposed positive constructions.
