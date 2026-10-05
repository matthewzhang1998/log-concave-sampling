# A bounded-Hessian original-gradient counterexample to the uniform six-tree test

2026-10-05. Source-qualified countertest to the OPEN hypothesis in `../COVARIANCE-TEST-PORT-REDUCTION.md`, not an impossibility theorem for a different grouped construction. Independent audit is recorded separately.

## Result

For arbitrarily large physical dimension D=n+1 there is a C-infinity original gradient g_n with

    (A_n/4) I <= Dg_n <= (3 A_n/4) I

at every point, using the literal original five-clock query genealogy and its unchanged private shields, for which the normalized center-center six-tree J_n and an HS-unit deterministic symmetric rank-six tensor H_n satisfy

    Var <H_n,J_n> >= n / 2^40,      n >= 128.                 (1)

The construction uses A_n=n^-3, delta_n=sqrt(A_n), center shield sigma_n=n^-1/2, fixed interior leaf shields, and ANY Price clock

    1-2delta_n^2 <= t <= 1.

Thus the lower bound is uniform over the entire declared Price panel in this example; it does not exploit just its endpoint t=1. It also survives positive averaging over that panel. All proper cuts of J_n remain bounded by numerical constants, and its HS norm is at most a numerical constant times sqrt(D). The new covariance operator really has an eigenvalue of order D. Since sigma_n^-2=n, the previously generic inverse-square-shield scale is attained along this family.

No fixed-order/public-polylogarithmic Lambda can satisfy the proposed estimate `Var <H,J> <= Lambda ||H||HS^2` uniformly over this original-gradient class. Fixed numerical normalization factors, symmetrization, positive readout constants, and the nonvanishing selected old-root injection do not change the conclusion.

## 1. The original source is genuinely admissible

Put kappa=1/2, eta=1/8, sigma=n^-1/2, and define on R^(n+1), with coordinates 0,1,...,n,

    U_n(x)=A_n [ kappa |x|^2/2
             + eta sigma^2 cos(x_0) sum_(i=1)^n cos(x_i/sigma) ],
    g_n=grad U_n,                  A_n=n^-3.                (2)

The source is C-infinity, g_n(0)=0, and no higher derivative oracle is used. Let E=Dg_n/A_n-kappa I. Its nonzero entries are

    E_ii = -eta cos(x_0) cos(x_i/sigma),                 i>=1,
    E_0i = E_i0 = eta sigma sin(x_0) sin(x_i/sigma),
    E_00 = -eta sigma^2 cos(x_0) sum_i cos(x_i/sigma).

The diagonal block has norm at most eta because n sigma^2=1. The off-diagonal arrow has norm equal to the Euclidean norm of its vector, at most eta sigma sqrt(n)=eta. Hence ||E||op<=2eta=1/4, proving the global sandwich preceding (1). This is stronger than the required bounded original Hessian and convex-gradient qualification.

Write P_h f(q)=E f(q+hZ) for additive Gaussian smoothing. Every perturbation entry above is multiplied under P_h by the SAME factor

    d_h=exp[-h^2(1+sigma^-2)/2].                           (3)

Consequently a leaf with h=1/2 has normalized Jacobian

    B(q)=P_(1/2) Dg_n(q)/A_n = kappa I+E_B(q),
    ||E_B(q)||op <= 2eta exp[-(n+1)/8].                   (4)

This bound is uniform over every leaf query and all of its old ancestors.

## 2. Exact old genealogy and Price endpoints

Use the source record in `../../../EXACT-FIVE-CLOCK-NATIVE-SOURCE-ADAPTER.md` (the path is relative to this note). To avoid ambiguity, the five-clock formulas are written in full. Fix retained Z=0, old clocks q=s=1/2, r1=r3=1/sqrt(2), and

    r2=sqrt(1-2/n),   sigma2=sigma=n^-1/2,
    sigma1=sigma3=1/2.

For the complete standard coarse bank Q=(G0,G1,W1,W2,W3),

    Y=(sqrt(3)/2)G0,
    X=Y/2+(sqrt(3)/2)G1,
    q1=r1 Y+W1/2,
    q2=r2 X+sigma W2,
    q3=r3 X+W3/2.                                       (5)

Every original force still has its independent original private shield: 1/2 at each leaf and sigma at each center. Nothing replaces a private shield by a larger common heat.

The whole-bank Price endpoints are Q and Q'=tQ+sqrt(1-t^2)Qtilde, with Qtilde an independent complete standard bank. This is the allowed affine center zero and outer bank width one in the Price identity. The center pair (q2_i,q2'_i) is an iid Gaussian pair across i=0,...,n, with variance

    v_n = (15/16)r2^2+sigma^2
        = 15/16-7/(8n),                                (6)

and correlation t. In particular the pair with i=0 is independent of all pairs with i>=1.

Select the old G1-block derivative injection on both center hits. Its contraction is the nonvanishing scalar

    gamma_n = r2^2(1-s^2) = (3/4)(1-2/n).               (7)

G1 also changes the right leaf query, but the specified product-rule history hits the center on each old tree. Other derivative-hit histories are separate, exactly as in the sealed ledger. Retaining G1 and every other root is essential; no ancestor has been replaced by an independent copy. Alternatively, summing all center-center old-bank injections replaces gamma_n by v_n. If the chosen convention absorbs this known injection into the scalar weight, set gamma_n=1. Each convention obeys gamma_n>=1/2 for n>=16 and gamma_n<=1, which is all the proof uses.

The chosen original clock numbers make the calculation transparent. The construction is robust to any fixed interior old q,s,r1,r3, and to any endpoint center-node sequence sigma->0 with n=floor(sigma^-2): all limiting constants stay positive, n sigma^2 stays bounded and approaches one, and (3)-(6) change continuously. Thus a fixed positive endpoint quadrature need not contain these exact illustrative clock values. This counterexample applies whenever its node family contains the usual arbitrarily small center shields and fixed interior other clocks. With A=sigma^6, any original cutoff A^b with b>1/6 includes such shields. No conclusion is being claimed for a separately restricted clock family that excludes them.

## 3. Literal analytical C2 tensor and six-tree contraction

The native analytical center tensor is

    T(q)= (sigma^2/A_n) P_sigma D^3 g_n(q)
        = sigma^2 P_sigma D^4(U_n/A_n)(q).               (8)

This is the C2 heat coefficient represented by finite original-VALUE filters, not a request to evaluate D^3g. Let

    d_n=exp[-(1+sigma^2)/2].

For i>=1 its only entries of the form T_(i i i a) are

    T_(i i i i)= eta d_n cos(q_0) cos(q_i/sigma),
    T_(i i i 0)=-eta d_n sigma sin(q_0) sin(q_i/sigma).  (9)

There is no term with a different nonzero coordinate a.

For clarity, label the six physical slots a,b,c,d,e,f. An ordered center-center six-tree with four leaves is

    J_(a b c d e f)
      = gamma_n sum_(u,v,w,z,h)
           B1_(a u) B3_(c v) T_(b u v h)(q2)
           B1'_(d w) B3'_(f z) T_(e w z h)(q2').        (10)

The prime is evaluation on Q', with its own independently averaged private original shield. Gradient symmetry makes any transpose convention immaterial. Equation (10) is precisely two three-marked old paths joined at their center derivative indices. All six original physical marks are still present.

Let J0 replace its four leaves by kappa I. The deterministic tensor

    H_n=n^-1/2 sum_(i=1)^n e_i^(tensor 6)                (11)

has HS norm one and is symmetric. Thus physical symmetrization of (10) leaves its pairing with H_n unchanged. From (9), writing U=q2_0, V=q2'_0 and alpha_i=q2_i/sigma, beta_i=q2'_i/sigma,

    <H_n,J0>
      = c_n n^-1/2 sum_i [ M X_i + sigma^2 N Y_i ],
    c_n=gamma_n kappa^4 eta^2 d_n^2,
    M=cos U cos V,       N=sin U sin V,
    X_i=cos alpha_i cos beta_i,
    Y_i=sin alpha_i sin beta_i.                         (12)

No unlisted center cross term is discarded in (12).

## 4. Quantitative variance lower bound

Condition on U,V. The pairs (X_i,Y_i) are iid and independent of U,V. Their means are

    mu_n = E X_i
         = [exp(-v_n(1-t)/sigma^2)
            +exp(-v_n(1+t)/sigma^2)]/2,
    nu_n = E Y_i
         = [exp(-v_n(1-t)/sigma^2)
            -exp(-v_n(1+t)/sigma^2)]/2.                 (13)

Therefore conditional variance decomposition gives the exact lower bound

    Var <H_n,J0>
       >= c_n^2 n Var(mu_n M+sigma^2 nu_n N).            (14)

Set delta_n^2=A_n=n^-3 and allow any t in [1-2n^-3,1]. For n>=16, v_n>=7/8 and v_n<=1. We have mu_n>=1/3 and |nu_n|<=1/2. Also

    sd(cos^2 U)=(1-exp(-4v_n))/sqrt(8) >= 1/3,
    ||M-cos^2 U||L2 <= ||V-U||L2 <= 2 n^-3/2.

Since standard deviation is 1-Lipschitz under L2 perturbations,

    sd(M)>=1/4,
    sd(mu_n M+sigma^2 nu_n N)
       >= (1/3)(1/4)-1/(2n) >= 1/32.                   (15)

Here sd(N)<=1 was sufficient. Moreover gamma_n>=1/2, kappa^4 eta^2=1/1024, and d_n^2>=exp(-2)>1/8, so c_n>=2^-14. Equations (14)-(15) imply

    sd <H_n,J0> >= 2^-19 sqrt(n).                       (16)

This is a true variance estimate over the complete old/Price bank. Captured retained Z=0 is an allowed caller value; a claimed bound uniform in the retained caller must cover it.

For completeness the sharper asymptotic is

    Var <H_n,J0>/n
      -> (gamma_infty kappa^4 eta^2 e^-1)^2
           (1-exp(-4v_infty))^2/32 > 0,                (17)

where v_infty=15/16 and gamma_infty=3/4 for (7). The limit is uniform across the Price panel. It follows directly from (12)-(13), bounded iid sampling variance divided by n, and U-V->0 in L2. The coherent surviving scalar is `cos^2(U)/2`, not an unsmoothed higher-derivative artifact.

## 5. Pointwise cut bounds and the exponentially small actual-leaf error

Every proper flattening of T has operator norm at most 16eta<=2. Here is a direct proof rather than an appeal to source qualification. Expanding its four derivative indices yields 16 placements of the zero coordinate. For k zero-coordinate indices and 4-k repeated i indices, its coefficients have absolute value at most eta d_n sigma^k. When 4-k>=1, each flattening has norm at most eta d_n times either sigma^k or sigma^k sqrt(n); the latter occurs only when all repeated i indices are on one side of the cut. For k=0 both sides contain an i index, so the cut norm is at most eta d_n. For k>=1, sigma^k sqrt(n)<=1. The k=4 tensor has only the all-zero entry, of absolute value at most eta d_n n sigma^4<=eta. Summing the 16 placements proves the claim.

The standard one-edge contraction lemma from the sealed Price note now bounds every proper cut of (10) by gamma_n times 2^2 times the four leaf operator norms. Those norms are at most one. Thus J and J0 have numerical proper-cut bounds and HS norms at most a numerical constant times sqrt(D), as required by the premise of the proposed covariance-test port.

Telescoping the four leaves in (10), using (4), gives every proper cut of J-J0 at most

    4 gamma_n (2)^2 [2eta exp(-(n+1)/8)]
       <= 4 exp(-(n+1)/8).

Isolating one physical index bounds its HS norm, hence

    |<H_n,J-J0>| <= 4 exp(-(n+1)/8) sqrt(n+1)
                  <= 6 exp(-(n+1)/8) sqrt(n).           (18)

For n>=128 the last coefficient is at most 2^-20. Combining (16)-(18) by the L2 triangle inequality for centered scalars gives

    sd <H_n,J> >= 2^-20 sqrt(n),

which proves (1). Leaf smoothing has not silently altered the source: the actual B matrices in (10) come from the very same g_n and literal fixed original leaf shields, and (18) accounts for their difference from the comparison kappa I.

## 6. Positive Price averaging, finite filters, and exact scope

The asymptotic coherent term in (12) is uniform in t over the full Price panel. If a positive normalized quadrature or the normalized integral over that panel is formed using its declared common G,H bank, the same coherent term remains, and Jensen bounds on the vanishing errors give the positive limit (17). This observation is only about that source-qualified center-center analytical coefficient. It does not assert that arbitrary signed combinations of other correction histories cannot cancel it.

The native finite C2 filters use original g VALUES and can approximate their analytical coefficient targets at separately frozen absolute Sobolev/Lp floors. The present proof tests that target, which is the J in the sealed OPEN hypothesis. Arbitrarily accurate finite calibration cannot turn an order-sqrt(n) scalar standard deviation into a dimension-free one: any L2-HS coefficient error o(sqrt(n)) preserves the asymptotic counterexample. No finite filter is falsely declared exactly multilinear, and no analytical coefficient is requested as a producer input.

The source's high frequency is 1/sigma and its potential amplitude is sigma^2. This is why its original Hessian stays uniformly bounded while its TWO normalized higher-heat jets reveal a common low-frequency modulation. A single native Jacobian matrix test does not see the same coherent product. Neither the source C2 assumption nor allowing smoothing forbids this family.

What has been refuted is the proposed shield-uniform scalar-test estimate for this normalized center-center six-tree, including the source's complete old roots and marks. The result does NOT refute the general any-order positive-law program, an integrated clock cancellation theorem, a different grouped common-heat coefficient, or a constructive current correction. It does not settle positive grouped VALUE realization, same-carrier means, higher joint physical-R tilts, generic chunk/spine preservation, or eventual-sublinear complexity. Those remain distinct obligations. In particular this countertest cannot be relabeled as a positive grouped VALUE producer.

## 7. What weighting does and does not follow

For the fixed source-qualified history in this proof, restoring the literal scalar amplitude from the sealed ledger gives

    C_n=A_n^7 beta_n J_n,       beta_n=w_n^2/sigma_n^4,
    Var <H_n,C_n> >= A_n^14 beta_n^2 n / 2^40.            (19)

Fixed known normalization/readout factors can be restored on both sides. Equation (19) is an exact fixed-node consequence. For any node whose beta_n is bounded below modulo public logarithms, the desired shield-free covariance certificate fails there as well. The sealed upper bound `Delta_2 <= C sigma_2^2` by itself does NOT supply such a lower bound on the actual w_n; a quadrature-specific lower bound must be checked separately.

A positive normalized average over the declared Price t-panel is covered by Section 6 because the same complete old/Price roots can be coupled there and the coherent term converges uniformly. This is not a statement about an arbitrary aggregate over the five ORIGINAL clocks. For a sum of completed nodes with genuinely independent complete coarse banks, its scalar variance is the sum of the corresponding weighted scalar variances; then (19) gives a lower bound from any included node. For shared-node banks, signed combinations, or a regrouped coefficient with additional exact cancellations, the actual weights, cross-covariances, and history identities must first be calculated. No universal weighted integrated no-go is asserted in this package.

## 8. A surviving restricted test and the grouping obligation

The uniform scalar-test port DOES hold for an orthogonally coordinate-separable original potential, with the same coordinatewise isotropic Gaussian genealogy and a deterministic retained caller. In a separating orthonormal basis write U(x)=sum_i u_i(x_i), with the bounded Hessian assumptions. Every normalized C0 leaf and C2 center is diagonal in all its physical/source indices. The six-tree is therefore

    J(Q)=sum_i j_i(Q_i) e_i^(tensor 6),

where Q_i denotes all old/Price roots in coordinate i. These banks are independent across i, and the scalar heat-derivative bound gives |j_i|<=K for a dimension-free fixed K. Consequently

    Var <H,J> = sum_i H_(i i i i i i)^2 Var j_i
              <= K^2 ||H||HS^2.                         (20)

Equation (20) is a restricted-source fact, not a replacement of the authorized general class and not a VALUE realization theorem.

For the nonseparable family (2), conditional on its modulation-coordinate bank (U,V), the ideal diagonal test (12) has variance at most a numerical constant: its n summands are conditionally independent and bounded after the n^-1/2 normalization. The entire order-n variance comes from the conditional mean `c_n sqrt(n) [mu_n cos U cos V+sigma^2 nu_n sin U sin V]`. Thus a different grouped strategy would have to keep this coherent old-root conditional mean as an explicit current and constructively account for it. Merely centering locally or invoking one-Hilbert proper cuts does not eliminate the required observer-boundary return. The existence of such a general positive grouped current construction remains open.
