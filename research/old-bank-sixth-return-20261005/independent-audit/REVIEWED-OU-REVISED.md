# Exact symmetric OU covariance grouping with the literal old clock weights

2026-10-05. Normalization-clarified revision; the earlier frozen artifact is preserved at SHA256 404c7829fb172cf1ccd95a8ead875f75537b7920259a465c949b222d901bc722.

Convention: g is always the literal original gradient, with ||Dg||op<=A. The service amplitude alpha=qA/sqrt(u) is distinct from A. Every source normalizes original g by A; alpha appears only in service/native amplitudes and their ledgers. The bridge numerical-error certificate in this revision conservatively permits M_old^2 D and tightens the scalar operator tolerance. The substantive grouped one-HS estimates are unchanged.

A new bounded analytical/source candidate for the old C0-C1-C0 cubic node. The exact covariance and Gaussian redistribution identities, endpoint estimates, and source normalization below are proved. They do not by themselves certify the full old physical-readout return. In particular, Section 8 is a real boundary condition, not an omitted routine step.

## 1. What changes relative to the unsuccessful reheating estimate

Let L(Q) be the exact tensor-valued old cubic coefficient, including its original scalar clock weight and source amplitude. Q is the complete standard old structural bank. Define P_r on this WHOLE bank, without changing any original private shield.

The exact identity is

    Cov(L) = 2 integral_0^1 r E[(P_r D_Q L) contracted_Q (P_r D_Q L)] dr.       (1)

Both gradient factors receive their own independent Mehler noise, conditional on the same outer Q. Their two copies have cross covariance r^2 I in the old bank. This differs materially from leaving one derivative factor unsmoothed.

The useful operation is to redistribute some of this already present joint Gaussian covariance into independent per-force shields. This preserves every original force-query covariance. It is not the replacement of L by a differently heated target.

For a cutoff a<1,

    Cov(L) = Cov(P_a L) + R_a,
    R_a = 2 integral_a^1 r E[(P_r DL) contracted (P_r DL)] dr.               (2)

R_a is positive semidefinite as an operator on the old tensor space. The literal six-tree estimate gives its proper cuts and HS/sqrt(D) at most

    Lambda alpha^6 (1-a^2) sum_old_nodes w^2/sigma_c^4,                     (3)

with the corresponding mixed-hit shield denominators for the other eight hit pairs. This is a valid one-HS near-endpoint remainder, using actual squared old weights. In particular 1-a^2=tau^2 gives Lambda alpha^6 tau^2, with no proposed normalized-six-tree variance theorem. A tensor coefficient allowance is not automatically a complete law-error allowance; the consumer still has to supply that implication.

Identity (1) follows from d/dr E[(P_r L) tensor (P_r L)] = 2r E[(P_r DL) contracted (P_r DL)], using Gaussian integration by parts. It also follows term by term in Hermite chaos. It applies to a sum of old nodes sharing a bank only after that sum is formed; cross-node terms must then be retained. Bounds below can absorb their finite public-log multiplicity by triangle/Cauchy-Schwarz, but cannot replace their correlations by independence.

## 2. A center-independent spectral lower bound for the exact five-clock bank

Use the exact old five-clock genealogy, conditional on retained Z:

    Y = q Z + sqrt(1-q^2) G0,
    X = s Y + sqrt(1-s^2) G1,
    q1 = r1 Y + sigma1 W1,
    q2 = r2 X + sigma2 W2,
    q3 = r3 X + sigma3 W3,
    sigma_i^2 = (1-r_i^2)/2.

Each force then has its separate original sigma_i-private Gaussian shield. Let Gamma be the 3 by 3 scalar covariance matrix of (q1,q2,q3) over Q=(G0,G1,W1,W2,W3), before these private shields are added. Put

    m = min(sigma1^2, sigma3^2, 1-s^2).

Then, at every interior clock tuple,

    Gamma >= (m/16) I_3.                                                   (4)

Crucially, m does not contain sigma2.

Proof. Write vY=1-q^2, vX=1-s^2 q^2. Regress the random part of r1Y on the random part of X:

    r1Y = a X + independent residual,  a=r1 s vY/vX,  0<=a<=1.

Discarding the independent residual and W2 gives the covariance lower bound

    sigma1^2 x1^2 + sigma3^2 x3^2 + vX(a x1+r2 x2+r3 x3)^2.

If r2>=1/2, write z=a x1+r2 x2+r3 x3. Then x2^2<=12(z^2+x1^2+x3^2), so |x|^2<=13(x1^2+x3^2+z^2). Since vX>=1-s^2, the lower bound is at least m|x|^2/13. If r2<1/2, sigma2^2>=3/8 and the independent W1,W2,W3 covariance alone is at least (3m/4)I, since m<=1/2. Both cases imply (4).

A uniform heat floor independent of ALL clocks is false. For example, as s,r1,r3 approach one, the query covariance approaches a low-rank common-shift covariance. The extra endpoint weights are therefore essential, not optional bookkeeping.

## 3. Exact independent-shield extraction inside the symmetric bridge

At bridge r let h^2=1-r^2. One factor P_r DL has joint original query-noise covariance, conditional on the outer Q,

    h^2 Gamma + diag(sigma1^2,sigma2^2,sigma3^2).

Define

    t_i^2 = sigma_i^2 + h^2 m/32,
    Gamma_res = h^2 (Gamma - m I_3/32).                                    (5)

Gamma_res is positive semidefinite by (4), and

    Gamma_res + diag(t_i^2)
       = h^2 Gamma + diag(sigma_i^2).                                     (6)

Thus draw a known correlated 3D-dimensional residual bank with covariance Gamma_res, and add independent t_i-private shields at the three original forces. The affine mean includes the unchanged deterministic retained-Z row and r times the original Q-to-query row. Do this with two independent residual/shield banks conditional on the same full Q. The cross-copy covariance remains exactly r^2 Gamma. All original query means, individual variances, within-copy cross covariances, and cross-copy cross covariances are unchanged from (1).

The original derivative injection in DL is the original Q-to-query row, not r times that row: P_r is acting on an already differentiated coefficient. The bridge contraction still sums EVERY coordinate of the original Q. Expanding the nine hit pairs only after applying (5) retains the literal source genealogy.

This is a Gaussian-law identity for the full coefficient. It does not state a conditional identity after exposing the new shield variables individually, and no already integrated private tape may be appended as an observer.

## 4. The actual positive old weights pay the spectral degeneracy

On a center endpoint panel of size Delta, sigma_c^2 is comparable to Delta and the original positive panel masses obey sum_panel w_c^p <= C_p Delta^p for p>=1, modulo the already declared public-log multiplicity. Let b^2=h^2 m/32. Dyadic summation gives

    sum_c w_c^2/(sigma_c^2+b^2)^2       <= Lambda,
    sum_c w_c^2/(sigma_c^2+b^2)^(5/2) <= C/b,
    sum_c w_c^4/(sigma_c^2+b^2)^5     <= C/b^2.                             (7)

The first bound can have the usual logarithm. The second follows by splitting Delta above/below b^2: above, terms are Delta^(-1/2), and below they are Delta^2/b^5. The third is analogous, with Delta^(-1) above and Delta^4/b^10 below. These are original w_c^2 and w_c^4, not independent center clocks substituted for a repeated clock.

The outer/leaf old weights absorb the m factors:

    m^(-1/2) <= sigma1^(-1)+sigma3^(-1)+(1-s^2)^(-1/2),
    m^(-1)   <= sigma1^(-2)+sigma3^(-2)+(1-s^2)^(-1).

Their squared or fourth positive panel masses make each corresponding sum finite, uniformly in inherited endpoint cutoffs. The q clock needs no inverse factor in this estimate. Original mixed-hit histories have t_c^2 t_v t_v' in place of t_c^4 and can be handled by the same product-weight bounds; this statement should be checked against each exact original history normalization when importing it.

There is a second important weight: the NEW symmetric bridge measure is 2r dr. Writing h=sqrt(1-r^2),

    integral_0^1 2r/h dr = 2.                                             (8)

A positive dyadic bridge rule with mass u_j comparable to its endpoint panel width has the uniform discrete analogue

    sum_j u_j/h_j <= C.                                                   (9)

Consequently the summed coefficient derivative allowance in (7) is finite without a fixed alpha heat floor. Dropping this bridge weight would incorrectly reintroduce an inverse-cutoff loss. Any finite bridge rule still needs its own positive quadrature/calibration error proof; (9) is the needed mass envelope, not that entire proof.

## 5. Literal six-vertex VALUE candidate and its first ledger

For the center-center hit, the two centers are C2 and four leaves are C0. Use the original gradient VALUES and same-center anchors at the new exact query representation (5):

    F_i(u) = [g(qbar_i+t_i u)-g(qbar_i)]/(A t_i).

These are genuine square-gradient sources with private radius at most one. Their C0/C1/C2 adapters are the usual finite bounded scalar filters with freshly frozen calibration. Analytical derivative tensors are used only to identify coefficients.

At bridge node u and old node w, the normalized six-tree scalar factor is

    beta = u w^2/t_c^4                                                       (10)

up to the fixed original derivative-injection and physical-readout constants. More generally it is u w^2/(t_c^2 t_v t_v'). Set five nonroot amplitudes to alpha and a genuine physical-leaf root amplitude to minus alpha beta divided by the positive six-mark readout product. This realizes the negative leading coefficient alpha^6 beta with positive probabilities and positive readout variance shares.

The source radius sum is Lambda alpha. The root's own coarse derivative is paid by its leaf t1>=sigma1 and the squared leaf-clock mass. The nearest center's actual derivative has one strict alpha ancestor and is bounded in the complete sum by

    Lambda alpha^2 sum_bridge u/h <= Lambda alpha^2.                        (11)

Other center paths have at least as many ancestors. Every changed affine caller rebuilds its full original gradient source and anchor. This is a path bound for the candidate actual program, not an inference from a small ideal-law cross. No external probe is promoted to a Gaussian channel terminal.

The six-mark root-spine has four forces, retaining the two center physical publics and both side substitutions. Its ordinary conditional second cumulant begins with eight forces, not twelve, at scale alpha^8 beta^2; all physical tilts and side covariance completions must be retained. The previous four-probe Wick expansion remains the appropriate conditional coefficient generator.

## 6. Bounded descendant ledgers that the new weights suggest

For a fixed r, the generic derivative-cut certificate for the complete six-tree coefficient has a center contribution sum w^2/t_c^5 <= Lambda/h after the other old weights are summed. Including (9), the full derivative allowance is Lambda alpha^6. The connected Price contraction of two such derivatives therefore yields the candidate bound

    new main's whole-coarse covariance:
        proper cuts <= Lambda alpha^12,
        HS <= Lambda alpha^12 sqrt(D).                                    (12)

Mixed with the old cubic coefficient, whose summed derivative allowance is Lambda alpha^3, the corresponding candidate bound is Lambda alpha^9 (and one sqrt(D) in HS). Correlated node banks can be embedded in their complete common Gaussian record and bounded by the same finite triangle sum; independence of nodes is not assumed.

These are bounds on the correctly grouped leading coefficients, conditional on the native derivative-cut facts and exact readout boundary being available. They are not yet a law theorem for the finished sampler. In particular the old main and correction can have additional public-tilt mixed coefficients, finite bounded filters need uniform frozen-label Sobolev/Lp floors at their actual t_i, and the same-endpoint consumer must account for every returned term.

There is no contradiction with the modulated C2 six-tree counterexample. That counterexample concerns a fixed normalized six-tree and a near-identity Price panel. Here both old cubic gradients are separately smoothed before contraction, their fresh common-heat covariance is redistributed exactly, and the NEW positive bridge weights are used. Nothing claims a shield-uniform scalar-test inequality for an arbitrary isolated J.

## 7. Finite work and numerical scope

For each old node, bridge node, and derivative-hit pair, the leading candidate uses six COMPLETE native pair programs, with the correct C0/C1/C2 pattern; the tight pair is four C0 plus two C2. Include all original g VALUES and anchors, finite filters, pair/calibration replicas, old Q, two residual bridge banks, readout publics, untouched keep, and discarded-primal replay. Known small scalar covariance square roots in (5) are not source coefficient oracles. Original HVPs are used only in requested first/adjoint sweeps at recorded original VALUE sites.

The new bridge count and filter precision may be public-log only if their finite positive quadrature/calibration guarantees are supplied. A chosen near-r=1 cutoff has the certified coefficient remainder (3). For a fixed target grade P one may take tau a fixed power of alpha; this does not blow the first-path sums (7)-(9). Higher derivative or arbitrary-order recurrence bounds are not asserted. The current result is a bounded sixth-current candidate and its first descendant ledger.

## 8. Physical readout tilts are an actual remaining gate

If the original physical Gaussian readout is independent of Q, the construction above has the appropriate unshifted bank target. If it has a row C Q, its exponential physical tilt shifts Q by mu=C^T theta. The exact identity (1) must then be applied with the shifted semigroup

    P_r^mu f(x)=E f(mu+r(x-mu)+hZ).

Naively using rQ+hZ under the old readout gives mean r mu, not mu. This is wrong even at the formal coefficient level for its public-tilt descendants.

A useful exact Gaussian feasibility test is as follows. Suppose W is the intended unit-covariance physical carrier and one wants two standard bank copies Q+,Q- with cross r^2 I and Cov(W,Q+)=Cov(W,Q-)=C. Such a joint Gaussian exists iff

    (1+r^2) I - 2 C^T C >= 0,                                              (13)

besides the standard dimensions. A concrete construction uses common W and independent plus/minus noises:

    Q± = C^T W + (1/sqrt(2)) Splus^(1/2) U
                    ± sqrt((1-r^2)/2) V,
    Splus=(1+r^2)I-2C^TC.

Thus ||C||^2<=1/2 suffices uniformly in r. With ||C||>1/sqrt(2), the full r range cannot preserve those rows in this representation. This is a real PSD obstruction to the naive coupling, not a no-go for another current construction.

Even when (13) holds, the full native program must place W and the bridge noises in actual known public/keep rows, keep selected-input/private independence, and retain the old packets' unchanged physical rows. One may not simply condition on a Gaussian carrier containing a new packet's own private terminal noise. A fresh correlated-bank realization can also change mixed old/new currents; these must be generated, not presumed equal. The full old retained-bank/public-readout producer is therefore NOT certified by this note.

## 9. Checks

check_exact_ou_grouping.py tests 20,000 random endpoint-heavy clock tuples, the covariance decomposition, and a finite Hermite-chaos version of (2). exact_ou_grouping_checks.json records the results. Observed minimum lambda_min(Gamma)/m was approximately 0.38787, above the proved 1/16; covariance reconstruction error was 2.3e-16, and the Hermite identity error was 3.6e-15. These numerical checks are only sanity checks; Sections 1-4 give the arguments.

## 10. Canonical retained boundary, exact rows, and nodewise maxima

The parent has verified that the intended canonical old packet has zero physical readout row on Q. Its physical publics are independent of Q, and the retained endpoint y is a caller held fixed. Under this boundary the shifted-readout issue in Section 8 does not arise. The following identities are conditional on that unchanged y and all unchanged old physical publics.

Let A be the 3 by 5 scalar old injection matrix, so Gamma=A A^T, and let d=(r1 q,r2 s q,r3 s q). At the bridge node use one complete old bank Q and, for each copy epsilon=+,-, independent standard 3D residual roots U_epsilon and independent private Z_i,epsilon. Define the known scalar matrix

    B_r = h (Gamma-mI/32)^(1/2),
    qbar_i,epsilon = d_i y + r (AQ)_i + (B_r U_epsilon)_i,
    X_i,epsilon = qbar_i,epsilon + t_i Z_i,epsilon.                        (14)

All U+,U-,Z+,Z- are mutually independent and independent of Q and the old physical publics. Then

    E[X_i,epsilon | y] = d_i y,
    Cov(X_i,epsilon,X_j,epsilon | y)=Gamma_ij+delta_ij sigma_i^2,
    Cov(X_i,+,X_j,- | y)=r^2 Gamma_ij.                                    (15)

These are exactly the means and covariances in (1). The same full Q, not independent Q copies, occurs in both factors. Old y rows and physical-public rows are unchanged. Both copies' independent residual banks must remain separate until their conditional coefficient is formed.

Each Gamma_ii<=1. Each complete affine qbar_i row on (Q,U_epsilon) has squared norm

    r^2 Gamma_ii + h^2(Gamma_ii-m/32)
       = Gamma_ii-h^2m/32 <=1.

The y caller row is |d_i|<=1. The original DL injection is A_i and has norm <=1; its contraction A_i A_j^T is a known scalar of magnitude <=1. Hence neither source replay nor derivative injection secretly costs 1/h, 1/m, or a dimension factor. The only inverse widths are the explicitly displayed native t_i-normalizations.

For exact C2 normalization, an old unnormalized coefficient has center factor P_sigma D^2g. Its original-Q derivative has factor P_sigma D^3g contracted with A_c. After (14), the exact privately averaged factor is P_t D^3g at qbar_c. The calibrated C2 target is (t_c^2/A)P_t D^3g with its two probe slots. Thus the two centers together contribute t_c^4/A^2 in original-derivative units; four C0 leaves contribute A^-4. Equivalently, in tensors made from the normalized original gradient g/A, the center product supplies exactly t_c^4. The service amplitudes supply alpha^6, so matching the old coefficient scale alpha^6 w^2 in those normalized derivatives requires precisely w^2/t_c^4, multiplied by bridge mass u and the known old injection. There is no leftover sigma_c denominator: introducing one would double-count the old C1 normalization.

Here are useful nodewise bounds in addition to the sums. Write w=v w_c, absorbing bounded clock polynomials into the positive factors. Let x=sigma_c^2, x1=sigma1^2, x3=sigma3^2, z=1-s^2. Endpoint node bounds give w_c<=C x and v<=C x1 x3 z (with the q weight bounded). Consequently v<=C m and v^2<=C m^2. On an endpoint bridge panel, u<=C h^2. Uniformly over the whole finite clock record,

    beta = u v^2 w_c^2/(x+h^2m/32)^2 <= C u v^2 <= C,
    beta/t_c <= C (u/h) v^2/sqrt(m) <= C,
    beta/t_1 <= C u v^2/sigma1 <= C,
    beta/t_3 <= C u v^2/sigma3 <= C,
    beta^2/t_c^2 <= C (u^2/h^2) v^4/m <= C.                              (16)

The constants include the fixed number of quadrature nodes per old panel when appropriate. Away from endpoint panels all widths are bounded below and the same bounds are immediate. These pointwise bounds show that a bad individual node is not concealed by the integrated estimates. The stronger summation statements (7)-(9) supply the complete actual caller budget.

## 11. An explicit positive finite bridge rule with an operator certificate

There is a dimension-free positive rule for the bridge that does not require sampling or differentiating an unbounded coefficient. Let

    R = integral_0^1 P_t dt

on L2 of the complete old Gaussian bank, where P_t is its Mehler semigroup. Since the kth Hermite chaos has multiplier t^k, R has multiplier 1/(k+1). With t=r^2, identity (1) is

    Cov(L) = E[DL contracted (R DL)].                                    (17)

Choose eta=2^-J and partition [0,1-eta] into panels

    I_j=[1-2^-j,1-2^(-j-1)],  j=0,...,J-1.

Use an n-point Gauss-Legendre rule on each panel, with positive weights u_jl and nodes t_jl, and set r_jl=sqrt(t_jl). This gives

    R_num=sum_jl u_jl P_(t_jl),
    ||R-R_num||_(L2->L2) <= eta + C 4^-n.                                (18)

Proof: for every integer k>=0, z^k is bounded by one on the fixed Bernstein ellipse of parameter 2 around each I_j. Indeed its center is 1-3d/2 and its major semiaxis is 5d/8, where d=|I_j|, so the ellipse lies strictly within the unit disk. Uniform degree-(2n-1) polynomial approximation on that ellipse gives panel quadrature error at most C d 4^-n. Sum d<=1. The omitted endpoint interval contributes at most eta for each scalar multiplier. Taking the supremum over Hermite chaoses proves (18).

Take J=O(log(1/epsilon)) and n=O(log(1/epsilon)). Then the number of positive bridge nodes is O(log^2(1/epsilon)), and the operator error is at most epsilon. The weights satisfy

    sum_l u_jl=|I_j|=d,
    0<u_jl<=d,
    d <= h_jl^2=1-t_jl <=2d,
    sum_jl u_jl/h_jl <= sum_j sqrt(d) <=1/(sqrt(2)-1).                    (19)

Thus the same explicit positive rule has both the required covariance accuracy and the endpoint mass envelope used in the actual-source proof. Its endpoint strip is eta=1-a^2 in (3).

For the numerical bridge discrepancy only, this revision uses a conservative Hilbert estimate. Let M_old be a deterministic finite envelope satisfying

    ||D_Q L||_(L2;HS) <= M_old sqrt(D),

where L is the complete old coefficient (sum first when nodes share Q), and M_old includes all original amplitude, clock, inverse-shield, injection and normalization factors. Such an envelope is computable from the frozen old source ledger; it is not a realized energy or an oracle evaluation. If E=R-R_num has L2 operator norm at most delta_op, then the HS norm of the uncontracted coefficient error is bounded conservatively by

    || E_Q[DL contracted (E DL)] ||HS
         <= delta_op ||DL||_(L2;HS)^2
         <= delta_op M_old^2 D.                                          (20)

The inequality follows by Cauchy-Schwarz on the full derivative tensors and the Hilbert operator bound; it also bounds every physical matrix cut by the same HS quantity. This numerical-error estimate is allowed to spend two dimension-sized energies. It is separate from the substantive original-law, grouped-current and near-endpoint one-HS estimates, which are unchanged.

Given a separately allocated absolute coefficient tolerance e_bridge, choose

    delta_op <= min(1, e_bridge/max(1,M_old^2 D)),
    J = ceil(log2(2/delta_op)),
    n = max(1, ceil((1/2) log2(16/delta_op))).                             (21)

For the Bernstein ellipse of parameter 2, the best degree-(2n-1) approximation error is at most 4*4^-n. Positivity and exact panel mass of Gauss-Legendre quadrature therefore give panel error at most 8d*4^-n. Equations (18) and (21) imply ||E||<=delta_op and hence numerical discrepancy <=e_bridge. The positive node count remains O(log^2(1/delta_op)), logarithmic in D, original inverse cutoffs, amplitude bounds and 1/e_bridge. When target-grade cutoffs are fixed powers of alpha, this remains a fixed-order/public-log cost.

A stronger one-HS estimate for the bridge quadrature error would require its own accepted operator-valued cut proof; this revision does not rely on it. The previous frozen note's proposed internal-edge one-HS quadrature argument is superseded for purposes of the complete positive return by the conservative bound (20). Finite native filters/pair priors remain separate absolute errors and need the complete-program tolerance propagation from the parent construction.
