# A smoothed finite bridge prefix beyond the frozen-prefix grade

2026-10-05. Independently reviewed constructive source theorem under the explicitly imported native guards. No existing sealed file is changed.

## 1. Result and imported input

Assume g=grad U, g(0)=0, 0<=Dg<=A I in the anchored C2 class. Import the actual native guards and complete LAW-only delayed-tail services from `/workspace/shared/law-only-variance-join-20261005/LAW-ONLY-SHRINKING-BUFFER-JOIN.md` (manifest SHA256 673bdfcc6dc54c00ea380c942039add5ecda9f32686e6af0398243923dc6f3a7).

Replace its frozen short prefix by the finite Gaussian bridge source below. Let w=1-q, delta=-log q, h in (0,1), r=sqrt(1-h²), and q>=1/2. The new prefix and smoothing error is

    C sqrt(D)[A²h+A h²+A³w+A³w³/h+A²w epsilon_bridge].     (1)

A second useful refinement narrows the fresh-RAW endpoint interval independently of the delay: use LAW services for t<=1-eta, and the actual RAW source for 1-eta<t<1, where 0<eta<=w. The complete mean error is (1) plus

    C sqrt(D)[A³(sqrt(eta)+eta/sqrt(w))
                +A⁴(w^-1+log(1/eta)+1)]
      +Lambda sqrt(D)[A⁵(w^(-3/2)+eta^(-1/2))
                       +A⁶(w^(-7/6)+eta^(-1/6))]+e_abs. (2)

Thus w=A^(3/5), eta=A, h=A^(7/5), epsilon_bridge<=A², and epsilon_outer<=A³ yield

    ||E T-m3||_(L2 gamma) <= C A^(17/5)sqrt(D)
                              +Lambda A⁴sqrt(D)+e_abs.    (3)

The two leading prefix terms and bulk Gaussianization have grade 17/5. Near-endpoint RAW debt has grade 7/2, coherent prefix corrections grade 18/5, and terminal smoothing drift grade 19/5. The logarithmic A⁴ term is included in the declared public-log Lambda. With the original wider endpoint branch eta=w, the same prefix construction already gives grade10/3; the independent cutoff is what improves the final grade to17/5. This is a guarded canonical-m3 finite positive original-VALUE mean theorem. The same direct graph ports below give its optional completed own-mean LAW. It is not an all-order theorem, zero-buffer theorem, or RAW covariance producer.

## 2. A separately charged terminal smoothing

Write B=F2 for the genuine nested force and Kf(x)=E[f(x-B(x,Z))]. Its actual Gaussian source has caller first and L2 energy bounded by L=C A and L sqrt(D), respectively. These assertions follow directly from the positive nested history: Lip_x B<=A/2+A²/4 and ||B||2<=A(1+A)sqrt(D).

For OU correlation r and c=sqrt(1-r²)=h, coherent use of the same complete B bank gives

    ||(K P_r-P_r K)g||2
       <=A L[(1-r)+sqrt(2(1-r))]sqrt(D).

Also ||R1(P_r-I)||_(2->2)<=1-r and ||Kg||2<=A(1+L)sqrt(D). Therefore

    ||R1 K(P_r-I)g||2<=C[A h²+A²h]sqrt(D).              (4)

Here P_r g is only analytical notation. Every executing terminal occurrence is the ONE original VALUE

    g(r u+h G_h),

with a fresh standard D-root G_h. No smoothed-value oracle, derivative of g, or expectation is executed.

## 3. Exact short-prefix bridge ancestry

Split the genuine force at delta:

    B=P_nested+b,
    P_nested=integral_0^delta e^-s g(X_s-H_s)ds,
    b=q F2_future(Y), Y=X_delta.

The future force b is independent of the entire past bridge conditional on Y. Replacing P_nested by P=integral_0^delta e^-s g(X_s)ds costs <=A²w sqrt(D) in force, hence <=A³w sqrt(D) at the smoothed terminal VALUE. This retains the true F2 tail; there is no global F2-to-F1 replacement.

Conditional on X=X0 and Y=X_delta, the OU bridge has

    X_s=a_s X+b_s Y+sigma_s Z_s,
    a_s=sinh(delta-s)/sinh(delta),
    b_s=sinh(s)/sinh(delta),
    sigma_s²=(1-e^-2s)(1-e^-2(delta-s))/(1-e^-2delta).     (5)

The Z_s are the normalized centered bridge, not independent in s. They are jointly independent of the entire future-tail bank given (X,Y). Both a_s and b_s are nonnegative. The maximum bridge variance is w/(2-w)<=w. Its full isonormal-bank first therefore obeys

    Lip_bridge P <=A integral_0^delta e^-s sigma_s ds
                  <=A w^(3/2).                         (6)

All infinite-history statements here are analytical and pass through finite strong Gaussian approximations preserving these bounds. No path grid is executed.

## 4. Finite positive bridge-marginal quadrature

There is a deterministic rule with beta_j>=0, sum beta_j=w, s_j in [0,delta], J_bridge=O(log²(1/epsilon_bridge)), for which the actual source

    P_Q(x,y,N)=sum_j beta_j g(a_j x+b_j y+sigma_j N)     (7)

uses ONE fresh D-root N. It has

    Lip_N P_Q<=A w^(3/2),
    ||E[P|X,Y]-E_N P_Q(X,Y,N)||2
                   <=epsilon_bridge A w sqrt(D).      (8)

The same N is deliberately shared across all nodes; no true multi-time path law is claimed or needed. It is independent of (X,Y) and the entire future force bank. The exact conditional mean in (8) is an analytical comparison, not a leaf.

Here is a constructive positive polylog proof. Standardize the endpoints as Y=q X+sqrt(1-q²) Z. The conditional-mean operator for time s is the second-quantization contraction Gamma(l_s), with

    l_s=(e^-s, 2q sinh(s)/sqrt(1-q²)).

On Hermite degree n it acts by the nth symmetric tensor power of this row; this statement is dimension-independent. For complex s=x+iy,

    ||l_s||²=1-sigma_x² + 4q² sin²(y)/(1-q²).

Since 1-e^-2x>=2x e^-2x and 1-e^-2(delta-x)>=2(delta-x)e^-2(delta-x), the product on the right of (5)'s numerator is >=4q²x(delta-x). Hence Gamma(l_s) is an analytic operator contraction whenever 0<x<delta and y²<=x(delta-x). In particular it is uniformly bounded in the wedges |y|<=min(x,delta-x). The weighted operator e^-s Gamma(l_s) has norm <=1 there.

Take endpoint cutoff a=epsilon_bridge w/64. Replace each endpoint interval by its exact scalar mass at s=0 or s=delta. Each replacement costs at most twice its interval length in operator norm. Partition [a,delta/2] into doubling panels of width at most their left endpoint; reflect the panels around delta/2 for the right half. A fixed Bernstein ellipse of parameter 2 around each panel remains inside the contraction domain. Positive n-point Gauss-Legendre quadrature on each panel has operator error <=C(length)4^-n, by degree-(2n-1) polynomial approximation and positivity. Choose n=C+ceil(log_4(1/epsilon_bridge)). Sum the errors and normalize the positive weights to total w. Increasing the numerical constant C if needed gives (8). There are O(log(1/epsilon_bridge)) panels, and O(log(1/epsilon_bridge)) nodes per panel. A small fixed precision decrease in epsilon_bridge absorbs all scalar setup error. This construction does not assume a Hessian modulus.

## 5. A one-energy comparison using only the new terminal root

Fix (X,Y) and the entire future force b. The true P and finite P_Q are independent of b given (X,Y). For a complete Gaussian source p with private first L_p, mean mu, and an independent hG root, the elementary one-energy mean-only Stein identity gives

    |E g(a-rp+hG)-E g(a-rmu+hG)|
                <= A r² L_p² sqrt(D)/(2h).             (9)

Proof: center p, take its Gaussian Stein field tau=[D_G(-L_G)^(-1)(p-mu)](D_Gp)^T, where L_G is the OU generator on its complete private Gaussian bank, and use ||tau||_(2;HS)<=L_p²sqrt(D). The mean interpolation gives an integrated factor 1/2 times tau:D²g. Move one derivative to independent hG, leaving Dg of operator norm A. Gaussian isometry gives E|tau G|²=E||tau||HS². No derivative of Dp or of Dg is taken or assumed, and no extra physical dimension appears. Smooth first only for the identity and pass to the limit preserving first bounds.

Apply (9) twice with (6),(8), then compare the conditional means with the A-Lipschitz terminal. Integrated conditional Jensen and R1 contraction give

    C A³w³ sqrt(D)/h+epsilon_bridge A²w sqrt(D).        (10)

This comparison is precisely where the independent terminal smoothing is spent. It does not spend conditional X variance and it does not assert a strong bridge approximation.

## 6. Coherent commutation and lawful tail replacement

Since partial_x P_Q=sum beta_j a_j Dg(...) is symmetric PSD and bounded by Aw I, u->u-P_Q(u,Y,N) is a contraction whenever Aw<=1. Moreover,

    ||P_Q(X,Y,N)-P_Q(X-b,Y,N)||2
        <=Aw ||b||2<=C A²w sqrt(D).                   (11)

Consequently commuting b into EVERY X-dependent bridge site costs <=C A³w sqrt(D). Define the actual terminal consumer

    F_(Y,N,G_h)(u)=g(r[u-P_Q(u,Y,N)]+h G_h).            (12)

It has Lip_u F<=A. Its retained roots N,G_h are independent of (X,Y) and of the future-tail bank before the comparison. All descendants in (12) are updated coherently when u changes.

At each outer t, use exactly the sealed joint disintegration (z,Y,X), variance v=(1-t²)(1-q²)/(1-q²t²), and positive shares v/2,v/4,v/4. Extend the LAW branch to t<=1-eta. For eta<=w<=1/2 this has c eta<=v<=2w, rather than v comparable to w at every node. The analytical true b is independent of conditional-X noise given (z,Y,N,G_h). The actual mean/Gram/K banks are COMPLETE and fresh given Y, so their LAW comparisons can retain unread (z,N,G_h). Thus (12) is a legal A-Lipschitz consumer of the same complete whole-LAW variance allocation. No service's carrier is subtracted during this LAW comparison, and no terminal LAW is relabeled RAW.

For t>1-eta keep the original fresh actual F_Q tail and consume it with (12). All finite nodes stay below t=1. Their one-energy mean-only comparison uses conditional-X noise while retaining the fresh bridge and smoothing roots. The same finite outer rule first approximates the fixed new analytical source, then replaces each node; it never treats the t-dependent service program as a semigroup.

The exact identity

    1/v=q²/(1-q²)+1/(1-t²)

gives the integrated bills, with numerical constants,

    integral_0^(1-eta) v^-1 dt <=C[w^-1+log(1/eta)],
    integral_0^(1-eta) v^(-3/2)dt<=C[w^(-3/2)+eta^(-1/2)],
    integral_0^(1-eta) v^(-7/6)dt<=C[w^(-7/6)+eta^(-1/6)],
    integral_(1-eta)^1 v^(-1/2)dt
                         <=C[eta/sqrt(w)+sqrt(eta)].   (14)

The same estimates hold for the actual finite positive outer rule. Insert both t=q and t=1-eta as panel boundaries. A dyadic panel of endpoint gap d has total weight O(d), and every node gap is comparable to d; splitting such a panel preserves these inequalities. For p>1, sum d^(1-p) down to d~eta is O(eta^(1-p)); for p=1 the count is O(log(1/eta)). On the RAW side, sum sqrt(d) is O(sqrt(eta)). The early midpoint gap is positive, chosen <=c epsilon_outer<=c A³, and is charged in the same square-root sum. Thus no endpoint singularity is replaced by a hidden pointwise constant or an inverse-A node count.

The pointwise native mean/Gram/K service bills are exactly the sealed ones at their ACTUAL local buffer v. Integrating those pointwise bounds by (14), plus the unchanged A⁴ tail-mean and target restorations, gives (2). At eta=A all allocated buffers are >=cA, which retains the full O(Lambda A) actual source first needed for the final ports. Combining (4),(10),(11), the nested-prefix omission and these finite sums proves (1)-(3).

## 7. Full first, curl, energy and optional own-mean LAW ports

At each outer node, let the actual tail-service preterminal output be X_t-d_t with its known source-zero carrier X_t=t z+c_t G_t. The sealed actual VALUE-path theorem gives first(d_t)<=Lambda A and energy(d_t)<=Lambda A(|z|+sqrt(D)), including all private and retained-label paths. The new complete terminal VALUE is

    g(r[X_t-d_t-P_Q(X_t-d_t,Y,N)]+h G_h).

Its source-zero carrier inside g is rX_t+hG_h, with caller row r t and Gaussian row norm sqrt(1-r²t²)>0. Rotate each node's full known Gaussian row to one common D-root G, retaining every perpendicular coordinate. This is an orthogonal row rotation only; one-node distributions and the legal comparisons above are unchanged. Then

    B(z,G)=sum_outer omega_t g(r t z+sqrt(1-r²t²)G)

is a genuine gradient in G, with first<=A. Let E=T-B. The displacement has full first Lambda A and energy Lambda A(|z|+sqrt(D)); the prefix adds at most Aw times the actual u/Y/N profiles. Hence

    first(E)<=Lambda A,
    curl(P_G*E)<=Lambda A²,
    ||E||_p<=Lambda_p A²(|z|+sqrt(D)),
    first_z(T)<=Lambda A.                              (13)

For the curl statement, the Hessian difference multiplying the scalar common carrier is symmetric. Every remaining nonsymmetric product has a terminal Dg factor O(A) and a displacement first O(Lambda A). Off-G paths also pass through that displacement. This argument never differentiates a Hessian.

All original sites vanish at the whole source zero because g(0)=0; numerical zero aliases use identical records. Capture and restore the true nonzero-caller origins on their full caller readsets. Apply the same imported order-four raw-split own-mean compiler with actual r=rho=Lambda A, delta=a_curl=Lambda A, mu=A and actual full dimension. It adds Lambda A⁴sqrt(D) and its fully enumerated floors, giving the completed N(m3(z),I) LAW with (3).

## 8. Roots, replays, VALUE counts and guards

Use the sealed Q_M,Q_H,Q_K and complete native dimensions, evaluated at the new delay w=A^(3/5) and local buffers down to c eta=cA. A bulk node costs Q_M+Q_H+Q_K+J_bridge+1 original VALUES. A near node costs Q_FQ+J_bridge+1. Each node adds exactly two D-roots, N and G_h, before alignment. Every native reentry of T replays its full original graph and gets a fresh complete tape. A retained current value never supplies a changed-argument replay.

If d_t^old is the sealed unaligned node dimension, the new aligned whole-source dimension is

    d_T=D+sum_nodes(d_t^old+2D-D).

The final raw-split baseline uses N_outer original VALUES per invocation; its residual uses Q_T+N_outer, with fully expanded native occurrence counts and capture/restoration costs exactly as in the sealed equation (31). Dimensions include all service coarse, filter, passive, fill, keep and pair roots. Actual complete private dimension divided by D and all VALUE counts are fixed public-log polynomials.

The new small scales w,h,positive interior bridge node gaps and row norms affect scalar precision only. There is no inverse-A replication, Monte Carlo averaging, strong path grid, or hidden mean evaluation. Count every original VALUE, requested HVP at a recorded VALUE site, adjoint, origin, mode, private-bank and discarded-primal replay. No HVP is an executing producer or differentiated again.

Check actual guards, not just asymptotic orders:

- q>=1/2, Aw<=1, 0<h<1, eta=A<=w, and every LAW-branch allocated u>=c eta=cA.
- Mean/Gram actual native radius Lambda A/sqrt(u) is at most C Lambda A^(1/2) and must satisfy its imported order-four radius window. Its complete residual first Lambda[A+A^(3/2)/sqrt(u)] is therefore at most C Lambda A.
- Its raw-split self-reserve first Lambda r_native²/sqrt(A) is at most a public-log factor times A^(1/2), with its own numerical guard.
- Mixed K actual radius C A/u^(1/3) is at most C A^(2/3); impose rho<=1/2, its imported order-five native radius window, and u>=A³. All covariance gaps retain their stated fixed fraction of u.
- The final own-mean raw-split is checked at actual radius/curl, complete dimension and public logarithm; the new roots and counts enter Lambda.
- To absorb Lambda A⁴ into C A^(17/5), additionally require Lambda A^(3/5)<=1. Otherwise retain Lambda explicitly as in (3).

Freeze every clock, positive weight, variance share, smoothing scale, bridge coefficient, orthogonal rotation, origin and numerical version before first/adjoint evaluation. Enumerate the finite full graph; assign each primitive tolerance by its actual downstream norm/profile so the sum is <=e_abs. Include every bridge row and weight, positive numerical square root, h,r, native calibration floor, terminal VALUE and caller restoration. No tolerance is divided by realized source energy. Intrinsic errors in (1),(2) remain substantive.
