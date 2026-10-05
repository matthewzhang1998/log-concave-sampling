# Same-carrier P3 mean and full covariance: exact contract and failed finite gate

2026-10-04. Bounded target-identity result and a uniform-C2 counterexample to dropping the self-covariance. Independent review requested. This is not a fourth-order endpoint theorem.

## Result first

For the analytical third stationary substitution P3=Z−F3, the correct fourth-order buffered-law target uses

- m3(Z)=E[F3|X1=Z] on the original standard Gaussian endpoint carrier;
- C2(Z)=Cov(F2|X1=Z), INCLUDING Cov(E2|Z), where F2=H+E2;
- the connected third cumulant of H, subject to the separately proved one-marked-energy third-cumulant comparison and native rank-three packet.

The posterior conditional-force re-entry E_nu(r,z) g is a different mean target. It cannot replace m3 in this join. Moreover, a claimed uniform bound Cov(E2|Z)=O(A4 sqrt(D)) is false in the stated Hessian class. Section 5 gives a C-infinity family, A=D^(-1/2), for which its L2(Z;HS) norm is Omega(A2), while A4 sqrt(D)=A3. Thus the existing O(A3 sqrt(D)) one-energy self-covariance bound is genuinely sharp at this scale. Any fixed public-log multiplier leaves the separation intact.

The analytical covariance replacement Cov(F3|Z) -> Cov(F2|Z) DOES cost O(A4 sqrt(D)). The missing finite target is the full latter covariance and the same-carrier m3, not merely a mixed term. No finite polylogarithmic producer for those two targets has been proved here. No paused unsmoothed grid is executed or revived.

## 1. Source pins and exact scope

The following inputs were read and hashed locally:

1. `../THIRD-ORDER-FULL-REVERSE-OU-LAW-AND-ENDPOINT.md`, SHA256 e83d737da13d05054ff0dda129a05e72d350c6e3e33b019d548cc517b816e702.
2. `../THIRD-ORDER-NESTED-FORCE-MEAN-WITH-C2-DECOUPLING.md`, SHA256 b056e213ac28c60869727589d3ee6d1b2c53d6666908cdbe1961fa3121161ced.
3. `../EXACT-TWO-CLOCK-CURRENT-AND-UNIFORM-C2-REFERENCE.md`, SHA256 0a966fbd6076dff94d2b928ec1a4a05aa6d33c9bc8222eff051b221ee9881477.
4. `../clock-compression/POSITIVE-SQUARE-CLOCK-COVARIANCE-SERVICE.md`, SHA256 108af2124e3da0e2a04af92e9f0403bfce2f82f7e2912ba6b7e1975ef1c8b18b.
5. `../order-reentry/FOURTH-ORDER-FORCE-MEAN-REENTRY.md`, SHA256 e1231f72136283b978895dc8763c08eddd35f73fc1aa3955a1e1b7fbb9888d85.

Input 1 compiles the conditional mean m2(Z)=E[F2|Z] and the leading covariance Cov(H|Z) to an integrated order-three full buffered law. Input 2 only proves a conditional-MEAN comparison for its cheap shared G,H source, with bias O(A3 sqrt(D)); it does not preserve the joint law of H and F2. Input 4 compiles the full continuous Cov(H|Z), not Cov(F2|Z). Input 5 compiles the posterior force mean after a full endpoint-law re-entry and explicitly does not retain the old Gaussian endpoint through its law coupling.

## 2. Literal Gaussian-carrier hierarchy

Let g=grad U, g(0)=0 and 0<=Dg<=A I. Let {X_r:0<r<=1} be the analytical stationary Gaussian OU process with Cov(X_r,X_s)=min(r,s)/max(r,s), and X1=Z. All coordinates use the same Markov history in every descendant. Define, for j>=1,

    F1,r = integral_0^1 g(X_(r tau)) d tau,
    F_(j+1),r = integral_0^1 g(X_(r tau)−F_j,(r tau)) d tau.

Put H=F1,1, F2=F2,1, F3=F3,1, and Pj,r=X_r−Fj,r. This indexing agrees with the source P2 in input 3. Then

    P3,1=Z−F3,
    W2(Law(P3,1),mu_U)<=A4 sqrt(D).

The contraction is only an analytical comparison, as in input 3. It does not execute an infinite Gaussian history.

By Gaussian Markov conditioning, the pre-r history, conditioned on X_r=x, is independent of Z. Consequently, with the SAME Gaussian-history conditional target,

    chi2(x)=Law(x−F2,1 | X1=x),
    psi2(x)=integral g(y) chi2(x)(dy),
    m3(Z)=E[F3|Z]=integral_0^1 P_r psi2(Z) d r.       (2.1)

Here P_r is the ordinary Gaussian Mehler operator. There is no posterior density claim about chi2. The forward Gaussian endpoint x is the captured conditional datum, and the private record integrated by chi2 is its past Markov Gaussian history.

A fixed positive Hermite-multiplier quadrature can approximate the OUTER resolvent in (2.1) to arbitrary fixed order if a qualifying finite psi2 source has already been supplied. It does not supply that source. In particular, the input-2 source gives m2 at order three, not a conditional chi2 sampler at order three.

## 3. Exact quadratic target-identity check

For one-dimensional g(x)=a x, write

    H0 = integral_0^infinity e^(−s) X_(−s) ds,
    J0 = integral_0^infinity s e^(−s) X_(−s) ds,
    K0 = integral_0^infinity (s2/2)e^(−s) X_(−s) ds.

Then F2=a H0−a2 J0 and F3=a H0−a2 J0+a3 K0. Conditioning on X0=z gives

    m2(z)=(a/2−a2/4)z,
    m3(z)=(a/2−a2/4+a3/8)z,                         (3.1)
    E chi2(z)=(1−a/2+a2/4)z,
    Var chi2(z)=a2/4−a3/2+5a4/16.                   (3.2)

Indeed Var(H0|Z)=1/4, Cov(H0,J0|Z)=1/4 and Var(J0|Z)=5/16. The posterior conditional-force target from input 5 is instead

    E_nu(r,z) g = a r z/[1+a(1−r2)],
    Var nu(r,z)=(1−r2)/[1+a(1−r2)].                 (3.3)

For example r=1/2 gives mean coefficient a/2−3a2/8+9a3/32+..., already differing from m3 at second order. More fundamentally, chi2 has conditional variance O(a2), whereas a fixed-heat posterior has variance of order one. An A-dependent r fitted to one quadratic identity would not establish a general target identity.

## 4. Orders that do follow, and the covariance species that remain

Conditioned on Z=z, positive weights and actual Gaussian row norms give, by finite approximations followed by L2 limits,

    Lip_private Fj <= A(1+A+...+A^(j−1)),
    ||F2−H||_(L2|z) <= C A2 (|z|+sqrt(D)),
    ||F3−F2||_(L2|z) <= C A3 (|z|+sqrt(D)).          (4.1)

The small VALUE energy does not improve the first of E2=F2−H: its terminal Jacobian difference is Dg(X−F1)−Dg(X), which can be O(A) under C2.

For Gaussian-root U,V, a one-marked-energy covariance bound is

    ||Cov(U,V)||HS <= Lip(U) ||V−EV||2.              (4.2)

For completeness, dualize against an HS-unit matrix M and use Gaussian integration by parts with the Ornstein–Uhlenbeck resolvent of V. The factor DU* M has HS norm at most Lip(U), and the gradient of the resolvent of centered V has L2 HS norm at most ||V−EV||2 by the Gaussian spectral gap. This proves (4.2) without a second dimension-sized energy. Finite Gaussian approximations suffice for the path targets.

Set DeltaF=F3−F2. The identity

    Cov(F3)−Cov(F2)=Cov(DeltaF,F3)+Cov(F2,DeltaF)

and (4.1)-(4.2) imply

    ||Cov(F3|Z)−Cov(F2|Z)||_(L2(Z);HS)
        <= C A4 sqrt(D).                            (4.3)

In contrast, the exact covariance decomposition is

    Cov(F2|Z)=Cov(H|Z)+Cov(H,E2|Z)
                    +Cov(E2,H|Z)+Cov(E2|Z).          (4.4)

Both the mixed pair and the self-covariance have the currently justified allowance C A3 sqrt(D), since Lip(E2)=O(A) and its energy is O(A2 sqrt(D)). The self-covariance cannot uniformly be assigned to the A4 remainder: Section 5 is a counterexample. Thus any new service must retain all three correction species in (4.4), or directly realize Cov(F2|Z).

The conditional mean difference m3−m2 has norm O(A3 sqrt(D)); it too must be corrected, not dropped at order four. Its exact VALUE expression is

    E[integral_0^1 {g(X_r−F2,r)−g(X_r−F1,r)}d r | Z].

Again, its small amplitude does not give an O(A3) first, and it is not the posterior-force mean.

## 5. C-infinity counterexample: the E2 self-covariance can saturate A3 sqrt(D)

Set epsilon=1/2, c=1/4, v_D=(1,...,1)/sqrt(D), A_D=1/sqrt(D), and

    b(t)=t+epsilon log cosh(t),
    g_D(x)=c A_D (b(x1),...,b(xD))
                         +c A_D v_D tanh(v_D dot x). (5.1)

This is the gradient of

    c A_D sum_i integral_0^xi b(t)dt
                            +c A_D log cosh(v_D dot x).

It is anchored at zero, smooth, and

    Dg_D=c A_D diag(1+epsilon tanh(xi))
              +c A_D sech2(v_D dot x) v_D v_D*,
    0<=Dg_D<=c(2+epsilon)A_D I < A_D I.               (5.2)

Thus it lies in the admitted class with no high-derivative assumption or violation of the small radius.

Use independent coordinate copies of the stationary Gaussian OU path. Set U_r=v_D dot X_r, so U is a standard one-dimensional OU path for every D. Write

    K_r^i=integral_0^1 b(X_(r tau)^i)d tau,
    L_r=integral_0^1 tanh(U_(r tau))d tau,
    I_r^i=F1,r^i=c A_D [K_r^i+D^(−1/2)L_r].

With ell(t)=log cosh(t), c0=E ell(N)>0 for N standard normal, and d=c epsilon c0>0, the coordinate law of large numbers gives

    v_D dot I_r = c[D^(−1)sum_i K_r^i+D^(−1/2)L_r]
                           -> d                       (5.3)

in L2, uniformly in r by stationarity, hence also after integrating r. Bounded tanh and the finite fourth moments of K justify every limit below.

The exact projection of E2=integral[g_D(X_r−I_r)−g_D(X_r)]dr, divided by c A_D, is

    integral_0^1 [-v_D dot I_r
      +epsilon D^(−1/2)sum_i {ell(X_r^i−I_r^i)−ell(X_r^i)}
      +tanh(U_r−v_D dot I_r)−tanh(U_r)]d r.            (5.4)

Because ||ell''||infinity<=1,

    ell(x−u)−ell(x)=−u tanh(x)+R, |R|<=u2/2.

The L2 norm of the total remainder after D^(−1/2) summation is O(D^(−1/2)): use Minkowski and E|I_r^i|4=O(D^(−2)). The leading coordinate sum in (5.4) converges by the ordinary L2 law of large numbers to the constant −epsilon c k, where

    k=E[K_r^1 tanh(X_r^1)],

independent of r by stationarity. The rank-one contribution inside that sum is O(D^(−1/2)). Combining this with (5.3),

    (v_D dot E2)/A_D − c[C0+J_D] -> 0 in L2,
    C0=−d−epsilon c k,
    J_D=integral_0^1 h(U_r)d r,
    h(u)=tanh(u−d)−tanh(u).                           (5.5)

The law of J_D does not depend on D. Its conditional variance given the complete endpoint Z_D is strictly positive. Here is a quantitative argument that also avoids any interchange of conditional limits. Fix 0<a<1 and put V_D=U_a−a U1. This Gaussian is independent of ALL coordinates of Z_D and has variance 1−a2. Gaussian integration by parts gives

    E[J_D V_D] = E[h'(N)] integral_0^1
                 [min(r,a)/max(r,a)−r a]d r
               = E[h'(N)] (−a log a).                (5.6)

The derivative expectation is strictly negative:

    E[h'(N)] = E sech2(N−d)−E sech2(N)<0.

For example, strict decrease of the convolution follows by expressing the positive even strictly decreasing function sech2 as a positive layer-cake integral of interval indicators; each Gaussian interval mass P(|N−d|<t) strictly decreases for d>0,t>0. Therefore the right side of (5.6) is a nonzero fixed constant.

Since V_D is independent of Z_D, conditional centering of E2 does not change its inner product with V_D. Equations (5.5)-(5.6) and Cauchy–Schwarz imply, for all sufficiently large D,

    E Var(v_D dot E2 | Z_D) >= c_* A_D2              (5.7)

for a constant c_*>0 independent of D. Consequently

    ||Cov(E2|Z_D)||_(L2(Z_D);HS)
      >= E[v_D* Cov(E2|Z_D)v_D]
      >= c_* A_D2.                                  (5.8)

But A_D4 sqrt(D)=A_D3. The ratio is at least c_*/A_D and diverges faster than any fixed polynomial in log(1/A_D), log(D). This disproves the suggested A4 remainder even with such public-log losses. At the same time A_D3 sqrt(D)=A_D2, matching the existing one-marked-energy bound.

This counterexample rejects ONLY the dropped self-covariance estimate. It is not an impossibility theorem for a positive fourth-order law construction that actually realizes this covariance.

## 6. Precise constructive admission contract, currently unfilled

A sufficient bounded route must supply two actual finite original-gradient VALUE sources/services conditional on exactly the same fresh standard Z and earlier captured callers:

M3. A complete positive mean-law source M3(Z) with integrated conditional target N(m3(Z),I), error Lambda A4 sqrt(D), and separately restored absolute floors.

C2. A complete positive covariance reserve with target N(0,v0 I+h2 Cov(F2|Z)), for fixed v0 bounded below and h<=1/2, error Lambda A4 sqrt(D), plus absolute floors. Directly compiling the full covariance or realizing every term of (4.4) is acceptable. A signed algebraic mixed action is permissible inside an ultimately positive Gaussian-root pushforward; a signed probability or tensor oracle is not.

Their complete private banks must be independent after Z and all original caller-only anchors are captured. A rank-three packet must target the raw-state sign −h3 kappa3(H|Z), with its induced mean/covariance feedback paid separately. The analytical error (4.3) then allows replacing Cov(F3) by Cov(F2); it does not erase any finite service error or feedback term.

One possible source qualification would be a finite JOINT VALUE record (Hhat,F2hat) with a proved same-Z conditional coupling to (H,F2) strong enough to give the full covariance difference at O(A4 sqrt(D)), together with actual one-energy and first bounds. Marginal mean matching is insufficient. A law-level compiler can instead furnish C2 directly, but must name and prove that exact target. No such finite joint/covariance service is supplied by inputs 1-5. In particular, their independent output couplings may not be joined while retaining an integrated old force/root.

For M3, a conditional chi2 sampler at integrated order-three state error would suffice after outer Lip(g)=A ONLY if its estimate is actually conditional on the captured Gaussian endpoint used in (2.1), all roots/callers are retained as stated, and outer quadrature is separately controlled. The posterior endpoint theorem has a different target and supplies no such conditional sampler.

## 7. Native records, counts, guards and stopping condition

No new producer is claimed, so there is no invented finite count. The exact missing count is the original-VALUE cost of an admitted M3 source and an admitted C2 source. The following are mandatory before either receives a PASS:

- Every changed outer/inner argument replays the entire old source graph, including all mode, anchor, filter, stage, action and nested ancestor calls. Cache only identical genuinely caller-only records under complete finite keys.
- Root dimension is the literal sum of all Gaussian banks. Z is exposed first. No old force, private clock, source root or law-coupling Gaussian becomes an observer after that bank has been integrated.
- Any derivative/adjoint uses original HVPs only at recorded original-VALUE sites, through all ancestors. No HVP is a producer and no saved HVP is differentiated.
- The actual first of a chord/correction is retained, generally O(A), despite an O(A2) or O(A3) VALUE energy. All coarse-caller and saved-origin firsts must be propagated explicitly. The fixed nonzero-Z origins are executed and restored.
- Floors are absolute with their actual coefficients and caller profiles, including global/conditional finite-mode tilt. They are not normalized by a realized correction energy.
- All clock orders, source versions, shares, padding, response widths, selector rules and numerical tolerances are frozen before differentiation. Every inherited radius, covariance-gap, clock, caller and precision guard is imposed at its actual active dimension and public-log majorant.
- Original VALUE query counts must be a fixed public-log polynomial for this bounded grade, including all complete replays; the old analytical OU path is not a free source occurrence.

The independent review should verify (2.1), (4.3)-(4.4), the full counterexample proof, and the no-join scope. Until both M3 and C2 are constructively filled and independently reviewed, the same-carrier fourth-order endpoint join is FAILED/OPEN, even if the separate native rank-three packet passes.
