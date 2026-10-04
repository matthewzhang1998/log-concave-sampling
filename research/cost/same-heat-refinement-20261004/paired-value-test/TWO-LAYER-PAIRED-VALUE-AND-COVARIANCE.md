# Two-layer paired VALUE packet and its covariance companion

Publication copy: nonmathematical provenance wording and/or local paths were sanitized. Original/public SHA-256 values are recorded in `INVENTORY.json`; historical source/audit pins refer to the original versions.

2026-10-04. Independent bounded test. Original potential regularity is C2, so the force b is C1. This is a scalar terminal/source-identity analysis, not a posterior-law lower bound or a new native family.

## Result first

The proposed finite packet

    C(q1,q2) = 2 b((q1+q2)/2) - [b(q1)+b(q2)]/2

cancels the quadratic Taylor contribution to its conditional force mean. It does not cancel the exact nonlinear-feedback rank-one current under the available C1 force contract. Its leading covariance is the original propagated covariance divided by two, including all cross-time terms. Thus it does not cancel both currents.

A concrete finite positive, bounded-first scalar covariance companion does exist for an endpoint-only observer: duplicate the complete packet endpoint, estimate its conditional covariance by half its squared difference, and use a known odd capped Gaussian fill with one untouched Gaussian keep. This preserves the packet's own mean and exactly restores its scalar variance before separately charged known-moment/numerical errors. It cannot change the remaining nonlinear mean.

At R=4, the available conditional certificate still charges the nearest previous layer N_near >= a^(-1). This remains true even if every covariance and higher centered current is granted ideal zero cost. The paired force program costs 2 N_near + 3 N_final original force queries plus the cached prefix, not their product. A duplicated complete endpoint for the explicit covariance companion costs 4 N_near + 6 N_final plus the same cached prefix. No incumbent or earlier cached prefix is rebuilt.

This is not a universal obstruction. Averaging the actual Gaussian root together with its actual later observers might improve the force-mean current. That would require a new coupled identity/estimate stated in Section 8, not conditional centering alone.

## 1. Exact records and actual observers

Freeze the entering caller, the earlier cached Picard graph and its labels, and the actual next-layer evaluation times t_l and weights w_l. The complete records of two independently drawn current-layer clock banks generate path vectors

    q1 = q + D1,   q2 = q + D2,
    E[Di | R] = 0,   D1 independent of D2 conditional on R.

Here q is the analytical conditional path mean. It is never queried. Each Di is a whole finite vector of the path values used by the actual next force calls. Its coordinates at different t_l are correlated. Independence is between complete banks, not between times within a bank. The actual later clock bank is independent of these two banks conditional on the earlier records; it may be included in R.

Put A=(D1+D2)/2 and

    C_l(t) = 2 b(q_l+t A_l) - 1/2 sum_i b(q_l+t Di_l),
    Y_t = B - sum_l w_l C_l(t).

Only the actually executed future force observations are included. An exterior read of a raw fresh clock or a paired child must be included explicitly and cannot silently be frozen into R while preserving its centering claim. This distinction is also needed for the endpoint-only covariance companion below.

## 2. The exact paired currents

Let

    J_t = -sum_l w_l [2 Db(q_l+t A_l) A_l
                       - 1/2 sum_i Db(q_l+t Di_l) Di_l],
    J_0 = -sum_l w_l Db(q_l) A_l,
    V_t = J_t - J_0,
    Q_s = J_0 tensor J_s.

The all-layer identity gives, exactly,

    E[phi(Y_1)-phi(Y_0) | R]
      = E integral_0^1 V_t dot grad phi(Y_t) dt
        + E integral_0^1 (1-s) Q_s : Hess phi(Y_s) ds.

In particular the force-mean error at coordinate l is

    m_l(R) = E[C_l(1)-b(q_l) | R]
      = E integral_0^1 {
          2 [Db(q_l+t A_l)-Db(q_l)] A_l
          - 1/2 sum_i [Db(q_l+t Di_l)-Db(q_l)] Di_l
        } dt.                                                (M)

This is a difference of original firsts evaluated at literal original-VALUE query paths. It uses no derivative of Db. The integral/current is an analytical identity, not an executed new differentiable HVP source.

The leading covariance source is

    1/2 E[J_0 tensor J_0 | R].                              (Q)

For scalar b this contains every term

    w_l w_k b'(q_l)b'(q_k) E[A_l A_k | R],
    E[A_l A_k | R] = 1/2 E[D_l D_k | R].

Thus pair averaging divides the leading covariance by two; it does not annihilate it. For affine b, V_t=0, J_t=J_0, and this covariance is exact. The quadratic-force Gaussian scale-mixture issue therefore survives as a constant-factor reduction unless a separate covariance action is supplied.

## 3. What Richardson cancellation actually establishes

For a formal quadratic force b with second derivative B,

    C(q+D1,q+D2)
      = b(q) + Db(q) A + 1/2 B[D1,D2].

Consequently its quadratic conditional mean bias vanishes. This algebra does not supply a Taylor remainder rate when b is merely C1. In particular, the two evaluations of a difference of Db in (M) cannot be assigned an extra power without a modulus or a coupled marginal theorem.

For the strict-sandwich smooth test

    b(x) = a [x/2 + c (sqrt(x^2+eta^2)-eta)],  0<c<1/4,
    q=0,  D1,D2 iid uniform[-epsilon,epsilon],

as eta/epsilon tends to zero,

    E[C(D1,D2)-b(0)]/(a c epsilon) -> 1/6.

Indeed E|D_i|=epsilon/2 and E|(D1+D2)/2|=epsilon/3. The original one-copy bias is a c epsilon/2, so the packet improves this fixture only by a factor of three. For iid centered Gaussian Di of standard deviation epsilon the corresponding factor is sqrt(2)-1 times E|Di|, again nonzero.

This is a conditional caller/source fixture, not a realization of every full posterior genealogy and not an unconditional W2 lower bound for the proposed complete sampler. In particular it does not preclude small-region cancellation from an unobserved nondegenerate Gaussian caller.

## 4. Which imported modules match

The basic two-vertex rectangular covariance reserve in `exact-slack/exact32-proof-v2.tex`, Proposition `prop:k2`, accepts an actual finite rectangular C1 VALUE source f with complete fresh first l and conditional centered energy e, and a positive allocated numerical variance v with Lambda l^2 << v. It is a possible conditional terminal covariance consumer for the literal paired source, after all its hypotheses and rescaling losses are checked. It emits a centered reserve targeting N(0,v I-Cov(f|R)). Its remaining lower/channel/root rows and the main source's el^2 rank-three carry remain charged. It returns no asserted RAW first. It does not change E[f|R].

The `SCALAR-K3-NO-COPY-MEAN-COROLLARY.md` module is restricted to its actual signed-twin K3 source with the proved true-Gram identities. No identity equating the present paired Picard source with that native K3 source has been demonstrated. Matching only a centered-energy/first power table does not authorize that substitution.

## 5. A direct positive scalar covariance companion

The following narrow construction uses original VALUES, their ordinary firsts, known scalar operations and capped Gaussian fills. It is independent of an imported K3 claim.

Let H be one actual complete paired endpoint before its final Gaussian keep, and let H' be an independent complete copy conditional on common records R. They use the same old caller, original heat, common momentum Z and earlier cached prefix. Their fresh paired banks are independent. For a full two-layer endpoint they also use independent final quadrature banks. Put

    S = (H-H')^2/2,   mu = E[H|R],   v= sigma/sqrt(2).

Only S is executed; mu is analytical. Take a known odd bounded C1 cap c_B, with bounded first, and independent standard K,L, unread by every force. Define the known scalar moments

    r=E c_B(K)^2,   s=E K c_B(K)>0,
    d(S)=(-v s + sqrt(v^2 s^2-r S))/r,
    Y=H+v K+d(S)c_B(K)+v L.                              (F)

Allocate the actual positive gap

    S <= v^2 s^2/(2r).                                    (G)

Known moments r,s are not expectations of an original source. Their finite numerical calculation/calibration and every square-root evaluation have their assigned numerical allowances. The following exact statements refer to the ideal known-moment graph before those allowances.

Because both K and c_B(K) are odd, the fill is conditionally centered. The defining quadratic identity is

    2 v s d + r d^2 = -S,
    Var(v K+d c_B(K)+v L | H,H',R) = sigma^2-S.

Independence gives E[S|R]=Var(H|R), and therefore

    E[Y|R] = mu,    Var(Y|R) = sigma^2.                   (V)

The fill can be correlated with H through S; (V) remains exact because its conditional mean is zero. It is unnecessary and generally false to regard it as an independent covariance-estimator reserve. In particular this construction is not an invocation of an independent-reserve source contract.

The source-zero coefficient is exact: d(0)=0. At S=0 the fill is the known Gaussian baseline v(K+L). The L coordinate has constant coefficient v and remains a genuine independent unread keep, with normalized variance sigma^2/2 and physical variance a sigma^2/2.

The first bound is explicit:

    d'(S) = -1/[2 sqrt(v^2 s^2-r S)],
    |d'(S)| <= 1/(sqrt(2) v s),
    DS = (H-H')(DH-DH').

Thus the potentially problematic root/clock term is d'(S) DS c_B(K), and it is bounded by the executed cap. No term proportional to an unbounded K multiplies a caller derivative. The K-first is v+d(S)c_B'(K); the L-first is exactly v. H and H' are differentiated only once, through their original gradient queries. No original HVP child is differentiated as a VALUE oracle. Uniform caller/clock firsts still include the displayed 1/v and the cap size, and a complete native hidden-host/retained-domain return is not asserted.

Only endpoint plus R is compared. If an exterior observer reads H, H', S, or the fresh paired banks, this is a different joint law and needs its own contract. Retaining records internally for correct first/adjoint execution is not permission to erase exterior observers.

### Remaining higher centered currents

With the uncapped Gaussian fill Y_g=H+sqrt(sigma^2-S) G, centered about mu, exact moment algebra gives

    E[(Y_g-mu)^3|R] = -1/2 E[(H-mu)^3|R],
    fourth_cumulant(Y_g|R) = -1/2 fourth_cumulant(H|R).

So exact covariance restoration does not remove all higher currents. The bounded cap version adds explicitly controlled cap-dependent higher moments; it does not inherit zero third or fourth current.

For completeness, a narrow scalar Gaussian-density bound is available if |H-mu| <= c sigma for a sufficiently small fixed c:

    W2(Law(Y_g|R), N(mu,sigma^2))
        <= C E[|H-mu|^3 |R]/sigma^2.                    (W)

To see this, standardize u=(H-mu)/sigma and t=S/sigma^2. The Gaussian likelihood ratio of N(u,1-t) relative to N(0,1) has its L2 Taylor expansion

    1+u He_1 + (u^2-t)He_2/2
      + O_L2(|u|^3+|u|t+t^2).

The first two nonconstant coefficients average to zero. For t=(u-u')^2/2 and iid bounded u,u', the averaged remainder is at most C E|u|^3. Gaussian transport-entropy then proves (W). This is scalar, uses a small bounded-displacement hypothesis, and does not establish dimension-uniform one-energy scaling.

The actual capped fill (F) differs from this Gaussian fill by a coupling cost

    C ||S||_2 ||c_B(K)-K||_2 / sigma,

provided r,s stay in fixed neighborhoods of 1. Bounds on |r-1| and |s-1| follow from the same cap-tail norm. The keep is split before this comparison. The numerical approximation of r,s has its additional charged error.

## 6. Heat and query ledger at R=4

Use the deterministic M>=3 Picard reference and the original same-heat incumbent once. For M=3 the prior scalar ledger assigns

    sigma = a^(5/4),
    deterministic Picard tail = O(a^(9/2)),
    keep law price = O(a^(3/2) sigma^2)=O(a^4),
    physical untouched keep variance = a sigma^2/2 = a^(7/2)/2.

For a capped previous path, its stratified clock error has scale

    e_D = O(a B N_near^(-3/2)),
    sup error = O(a B N_near^(-1)).

The paired force feedback mean has the still-valid bound from (M)

    |m| <= C a e_D,
    physical propagated endpoint price
       <= C a^(5/2) N_near^(-3/2)                      (P)

with scalar moment constants/public cap logs included. Consequently the available conditional proof requires N_near >= a^(-1) for a^4 error. The conditional fixture in Section 3 explains why cancellation of its quadratic Taylor coefficient does not by itself improve (P). It does not prove (P) sharp for the full marginal posterior.

With the existing final-bank buffered estimate one has

    physical final-clock price = a^(5/2)/(sigma N_final^3),
    N_final >= a^(-11/12).

For M=3 the remote first layer can use N_1>=a^(-1/3), since its propagated price is a^(7/2)N_1^(-3/2). Thus the nearest paired layer dominates, and the certified total original-query exponent remains 1.

The direct covariance companion's centered endpoint error may be budgeted separately using (W). A full paired endpoint has scalar centered Lp scale for every fixed p used here (in particular p=3,4)

    e_H <= C [a N_final^(-3/2) + a^2 N_near^(-3/2)],

and sup centered displacement at most

    C B [a/N_final + a^2/N_near].

The latter, divided by sigma, must be small to certify the actual fill gap. The existing N_near=a^(-1), N_final=a^(-11/12) more than suffice. Formally applying (W), including its explicit hypotheses, gives at sigma=a^(5/4)

    physical centered law error
       <= C [a N_final^(-9/2) + a^4 N_near^(-9/2)]

plus cap-tail and numerical rows. This can improve the scalar final-only allocation to N_final=a^(-2/3), which also meets the gap. It does not lower the overall exponent 1 because (P) remains. No optimized general covariance/native source bill is inferred from this optional scalar terminal estimate.

The actual variable covariance branch is not as small as its terminal weak-law remainder. Its normalized L2 energy is

    ||d(S)c_B(K)||_2 <= C ||S||_2/sigma <= C e_H^2/sigma.

At N_final=a^(-2/3), one has e_H=O(a^2), so this is O(a^(11/4)) normalized and O(a^(13/4)) physically. Deleting that branch is therefore certified only at retained mark 11/4=2.75 in the convention physical error sqrt(a) a^K; it does not return a mark 3.9. Keeping its cap/square-root/product graph instead requires its own native retained admission. Neither terminal law control nor its ordinary firsts imply that admission. Likewise the final protected variance a^(7/2)/2 is real but is not automatically compatible with the heat/width required by a later protected proxy.

Literal queries without the direct duplicate companion:

    Q_old(a) + Q_cached_prefix + 2 N_near + 3 N_final
      + Q_mode + Q_zero + Q_guard + Q_numerical.

For (F), duplicating the complete paired endpoint gives

    Q_old(a) + Q_cached_prefix + 4 N_near + 6 N_final
      + Q_mode + Q_zero + Q_guard + Q_numerical.

Each copied packet owns its fresh clocks, but all copies reuse the old caller and cached earlier forces. The dense known-kernel work may be O(N_near N_final) and is real arithmetic; it is not an original-force query product. A first/adjoint sweep uses one original HVP at each executed gradient point. The fresh tape contains 4 N_near+2 N_final clocks for the full duplicate construction, plus K,L; old and earlier clocks remain charged. Any later source consuming that tape pays its actual dimension.

## 7. Why this does not consume the retained-path zero-Gram obstruction for free

The covariance companion acts after the actual two-layer force readouts and compares only their scalar endpoint and common conditioning records. Its new L keep supplies positive variance at that endpoint. It does not create variance in earlier path coordinates. If the future caller retains or reconstructs the intermediate path as an additional observer, its transverse zero-Gram directions remain and require a separate joint-path theorem/reserve. Conversely, one should not require unused intermediate path coordinates as artificial exterior observers when testing this narrow endpoint construction.

## 8. The concrete coupled identity still needed

The missing improvement is an observer-qualified treatment of

    E_R E_D integral_0^1 V_t dot grad phi(Y_t+fill) dt,

together with the same-record Q_s rank-two current and every current generated when moving its source point. The expectation must retain the actual Gaussian q,D genealogy, the real future force evaluations and labels, the relevant caller/root firsts, and the unread keep chronology. It must provide a quantitative extra power under original C2, without differentiating Db, and either cancel the mean (M) or control its actual marginal transport effect at a better grade.

The finite paired VALUE identity and positive covariance fill provide concrete source programs on which to test such a result. They do not currently supply it. Neither the fixed-caller near-kink test nor the conditional exponent ledger rules out a different valid marginal or coupled source theorem.

## Inputs and provenance

- `randomized-weak/ALL-LAYER-MARTINGALE-CURRENTS-AND-PATH-GAP.md`: SHA256 `dbf2b9b04c72cd83970273cb487b2ac800ff9358dec945e1ef489f054b829d6b`.
- `randomized-weak/SCALAR-RANDOMIZED-QUARTER-FLOW-AND-ALL-LAYER-DEBTS.md`: SHA256 `14ff69c491ab5cbcc3dfcfa196156cefced1c685f4905e5df663a6f8e48f8f27`.
- Imported contract scope checked against `research/exact-slack/exact32-proof-v2.tex`, Section `sec:k2`, and `research/r-native-twin/SCALAR-K3-NO-COPY-MEAN-COROLLARY.md`.
- The exact capped-fill parameterization in Section 5 was independently checked here.
