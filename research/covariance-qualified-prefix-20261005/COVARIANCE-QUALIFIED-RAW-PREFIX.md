# A covariance-qualified RAW bridge prefix, with explicit inverse-A work

2026-10-05. Constructive finite source result. This packet separates a complete cost-qualified construction from the still-unproved public-logarithmic covariance-LAW route. It does not claim an any-order/public-log advance.

## Result

Let g=grad U, g(0)=0, 0<=Dg<=AI, q=1-w>=1/2, delta=-log q, and let P=integral_0^delta e^-s g(X_s) ds be the true ordinary OU prefix conditioned on its original endpoint pair (X,Y). For every integer M>=1 and 0<epsilon<1/4 there is a positive finite original-g-VALUE source Q_M(X,Y,Z), with no expectation, derivative, covariance or tensor producer, satisfying

    complete private first <= L := A w^(3/2),
    each endpoint first <= Aw,
    0 <= D_X Q_M <= Aw I,
    mean defect <= epsilon Aw sqrt(D),
    covariance defect <= C A^2 w^3 M^(-1) sqrt(D).       (1)

The moment defects are respectively the L2 and L2(HS) norms over the same stationary Gaussian endpoint pair. There are

    J <= C [M+log(1/epsilon)] log(1/epsilon)

original VALUE sites and at most JD independent Gaussian coordinates. Every source partial Jacobian is symmetric. The covariance is an analytical target only.

The complete two-source terminal-h comparison, paying both third Stein defects, replaces the previous mean-only prefix debt by

    C sqrt(D) [A^2 w epsilon
       + A^3 w^3/(h M) + A^4 w^(9/2)/h^2].       (2)

The terminal smoothing, nested-prefix omission and coherent-tail commutation still cost C[A^2 h+Ah^2+A^3w]sqrt(D). The Gaussian endpoint covariance certificate is used BEFORE commuting the tail displacement; the actual Aw endpoint first makes the subsequent commutation legal.

In particular the cost-qualified parameters

    w=A^(4/5), h=A^(9/5), M=ceil(A^(-1/5)),
    eta=A^(8/5), epsilon<=A^2,
    variable-buffer raw-split padding mu=A^(1/2)          (3)

give a finite positive canonical-m3 source at

    C A^(19/5)sqrt(D) + Lambda A^(19/5)sqrt(D) + e_abs,  (4)

under the actual native guards below. The public-log Lambda stays at leading order. This improves the fixed 15/4 ledger ceiling ONLY by changing the source and paying an explicit inverse-A node count. It is not the desired public-logarithmic covariance port. The original fixed-tape ceiling was never a lower bound on all possible graphs.

## 1. Positive clocks and exact joint bridge execution

Use the existing analytic conditional-mean quadrature construction with endpoint cutoff a<=min(epsilon w/64,delta/M) and its doubling panels. Further split every interior panel until its width is at most delta/M. A panel still lies in the same fixed analytic contraction wedge. Use n=C+ceil(log_4(1/epsilon)) positive Gauss nodes per panel. Replace each endpoint panel by its exact e^-s mass at the endpoint. On every interior panel multiply the positive quadrature weights by the known ratio between its exact e^-s mass and its quadrature mass. The scalar ratio is 1+O(4^-n), preserving the analytical mean error after adjusting the numerical constant. All panel masses are now exact, all weights beta_j are positive, sum beta_j=w, and the claimed J count follows. These are public scalar setup operations, frozen before differentiation.

Sort the interior nodes 0<s_1<...<s_J<delta. Starting with x_0=X, sample recursively

    x_i = A_i x_(i-1)+B_i Y+S_i Z_i,
    A_i = sinh(delta-s_i)/sinh(delta-s_(i-1)),
    B_i = sinh(s_i-s_(i-1))/sinh(delta-s_(i-1)),
    S_i^2 = (1-exp(-2(s_i-s_(i-1))))
              (1-exp(-2(delta-s_i)))
              /(1-exp(-2(delta-s_(i-1)))),             (5)

where s_0=0 and the Z_i are independent standard D-roots. Endpoints are literal X,Y values and use no root. Every denominator and every interior innovation variance is strictly positive. Formula (5) is the ordinary scalar-block Markov bridge disintegration, so the executed x_i have the TRUE joint conditional bridge law. Return

    Q_M = sum_j beta_j g(x_j).                         (6)

Every node is an original VALUE. Unlike the earlier common-root marginal source, (6) preserves the multitime covariance of these actual nodes. No Gaussian root of an unknown matrix is requested: every row in (5) is a known scalar row.

## 2. Full bank first, caller first, energy and moment proof

The expanded x_j has endpoint coefficients

    a_j=sinh(delta-s_j)/sinh delta,
    b_j=sinh(s_j)/sinh delta

and a scalar row on (Z_1,...,Z_J) of Euclidean norm sigma_j, where

    sigma_j^2=(1-exp(-2s_j))(1-exp(-2(delta-s_j)))
                         /(1-exp(-2delta)) <= w.

Weighted triangle inequality therefore gives the complete private first A sum beta_j sigma_j<=A w^(3/2), independently of J. All coefficients a_j,b_j and all innovation-row coefficients are nonnegative. Thus each private-block partial Jacobian is symmetric PSD; D_X Q_M and D_Y Q_M are symmetric PSD and bounded by Aw I. This is a direct graph bound, not a consequence of a LAW estimate.

The conditional centered L2 energy is at most L sqrt(D) by Gaussian Poincare. Uncentered energy is at most C Aw(|X|+|Y|+sqrt(wD)), directly from the positive sum and the marginal row norms. At X=Y=Z=0 all original sites are zero; identical numerical zero records are reused. The graph's caller origin is the literal same-key graph value and is captured/restored when a native consumer needs it.

One-time marginals and the original bridge quadrature operator proof give the mean statement in (1). For the covariance proof, couple the executing nodes to the genuine continuous stationary OU path. Let

    nu(ds)=e^-s ds-sum beta_j delta_sj(ds),
    F(s)=nu([0,s]), Delta=delta/M.

Each panel has exact mass and width at most Delta. Hence F vanishes at panel boundaries, |F(s)|<=Delta, and integral F(s)^2 ds<=delta Delta^2. The vector-valued stationary OU spectral measure mu_g satisfies

    E[g(X_s) dot g(X_t)]=integral exp(-lambda|s-t|)dmu_g(lambda),
    integral lambda dmu_g(lambda)=E||Dg||_HS^2<=A^2D.

The constant-chaos term vanishes because nu has total mass zero. For each positive spectral frequency, distributional integration by parts gives

    integral integral exp(-lambda|s-t|)nu(ds)nu(dt)
       =2lambda integral F^2
          -lambda^2 integral integral exp(-lambda|s-t|)F(s)F(t)dsdt
       <=2lambda delta Delta^2,

because the exponential covariance kernel is positive semidefinite. Therefore

    ||P-Q_M||_2 <=sqrt(2) A delta^(3/2)sqrt(D)/M
                  <= C L sqrt(D)/M.                  (7)

This is a stationary-OU one-energy estimate, using only the original bounded Hessian through E||Dg||_HS^2. It is stronger than the elementary pointwise-in-time Brownian modulus estimate, which would yield only M^(-1/2). No new smoothness, temporal derivative oracle, or mean replication is assumed.

The same statement remains true for endpoint replacement panels because their chosen sites lie in the panel and their masses are exact. This strong coupling is only used to prove the covariance certificate; it does not assert public-log strong path approximation.

For any vector U and V,

    ||Cov(U,V)||_HS <= sqrt(||Cov(U)||op) ||V-EV||_2.   (8)

Indeed test against any HS-unit matrix K and apply Cauchy-Schwarz to (K^T(U-EU)) dot (V-EV). Conditional Gaussian Poincare gives ||Cov(P|X,Y)||op and ||Cov(Q_M|X,Y)||op at most L^2. Put R=P-Q_M and use

    Cov(P)-Cov(Q_M)=Cov(P,R)+Cov(R,Q_M).

Applying (8) conditionally and then integrating (X,Y), with conditional centering a contraction, proves the covariance statement in (1) from (7). This is a one-HS-energy proof with sqrt(D), not D.

## 3. Legal h-buffer comparison and retained endpoints

Retain the genuine future force b, independent of the entire past bridge given (X,Y). The new roots in (5) and the independent terminal hG_h are fresh and unread by every future LAW service. First compare P and Q_M at the original Gaussian endpoint pair. The established two-source covariance-preserving Stein estimate gives

    A delta_m + C A delta_Sigma/h + C A L^3 sqrt(D)/h^2.

Both complete third Stein defects are included. It yields (2), up to constants and r=sqrt(1-h^2)<=1. No conditional-X variance is spent in this comparison.

Only afterward move the future displacement into EVERY X-dependent site of (5)-(6). Since D_X Q_M is symmetric PSD <=Aw I, the actual consumer

    u -> g(r[u-Q_M(u,Y,Z)]+hG_h)

has Lipschitz constant <=A whenever Aw<=1. The change Q_M(X,Y,Z)-Q_M(X-b,Y,Z) costs at most Aw|b|, hence CA^3w sqrt(D) after the terminal g. All descendants and their anchors replay under the changed u. The source retains the actual (X,Y) ancestry throughout; a certificate at unshifted endpoints is never silently applied at shifted endpoints.

The independent complete bulk mean/Gram/mixed-K/cubic LAW banks and the near-endpoint original RAW F_Q branch can now be substituted exactly as in the sealed combined source. Z and G_h are unread retained labels during those comparisons. No part of a consumed LAW carrier is subtracted or reused as RAW data.

## 4. Integrated exponent ledger and native guards

The new prefix exponents at (3) are

    A^2h: 19/5,
    Ah^2: 23/5,
    A^3w: 19/5,
    A^2w epsilon: at least 24/5,
    A^3w^3/(h M): at least 19/5,
    A^4w^(9/2)/h^2: 4.

The RAW endpoint exponents are 19/5,21/5,28/5 for A^3sqrt(eta), A^3eta/sqrt(w), A^4eta. Fixed-target outer quadrature remains O(A^4 sqrt(D)).

At variable buffer u, use the actual six mean/raw-split debts AFTER the terminal A factor with mu=A^(1/2):

    A^5u^(-3/2), A^5u^-1,
    A^(9/2)u^(-1/2), A^5u^(-1/2),
    A^(23/4)u^(-3/2), A^6u^(-3/2).                 (9)

The smallest integrated exponents are respectively 19/5,21/5,41/10,23/5,91/20,24/5. The separately retained Gram covariance-mixture term A^5 sum omega u^(-3/2) remains 19/5. The existing cubic/fourth-current terms have minimum exponents 19/5 for A^5u^(-3/2), 22/5 for A^6u^-2, 23/5 for A^7u^(-5/2), and 29/5 for A^7u^(-3/2). Mixed-K terms are higher. Exact finite positive dyadic sums, not an unproved integral substitute, are the ones in the sealed combined source.

Every executing LAW share is >=c eta. The actual normalized mean/Gram radius Lambda A/sqrt(u) has power A^(1/5); the native self-reserve first Lambda A^2 mu^(-1/2)/u has power A^(3/20). The per-node physical residual first is Lambda[A+A^(7/4)/sqrt(u)], which may exceed A near the endpoint, but its positive weighted sum is <=Lambda[A+A^(27/20)]<=C Lambda A. The actual mixed-K radius A/u^(1/3) has power A^(7/15). Cubic amplitude A/sqrt(u) has power A^(1/5); include every actual inverse shield, readout share and clock coefficient in its native Lambda and guard. Its reference covariance correction relative to u has power A^(2/5), and third-order reference gap parameter A^3/u^(3/2) has power A^(3/5). All margins are positive; they must still satisfy the literal numerical native windows. q>=1/2, eta<=w<=1/2, Aw<=1 and 0<h<1 are also required.

The changed mu is used ONLY in the variable-buffer mean services receiving the later terminal A. If a final own-mean compiler is subsequently used, it must retain mu=A and be rechecked against the actual complete first/curl/energy, full dimension and source-block count. This packet exports the needed direct source ports but does not certify a uniform final-compiler cost or its dimension/block-dependent constants when J grows as an inverse power. Equation (4) is the uncompiled canonical-mean source theorem; it is not by itself a completed own-mean LAW theorem. No leading Lambda is absorbed without a separate window.

## 5. Actual first/curl and complete execution bill

At every outer node replace the old one-root prefix by (6) with its J independent bridge roots. The untouched terminal root G_h remains separate. The preterminal displacement still has complete weighted first Lambda A and energy Lambda A(|z|+sqrt(D)): the prefix adds at most Aw times actual caller/root profiles, independent of J by the row-norm proof above. Align the known Gaussian source-zero terminal row exactly as in the sealed combined source, preserving all perpendicular coordinates. The resulting common-carrier baseline is

    B(z,G)=sum_t omega_t g(rtz+sqrt(1-r^2t^2)G).

The Hessian-difference term multiplying its scalar carrier is symmetric; every remaining skew product includes terminal Dg of norm A and one displacement first. Therefore the actual whole source has baseline first <=A, residual first Lambda A, square-lift curl Lambda A^2, residual energy Lambda A^2(|z|+sqrt(D)), and caller first Lambda A. No derivative of a Hessian or of a LAW error is taken. Source norms use the actual positive marginal profiles, so no sqrt(J) energy is introduced by replacing these estimates with the loose norm of the full root tape.

Let Q_Mean,Q_Gram,Q_K,Q_3 and d_Mean,d_Gram,d_K,d_3 denote the fully expanded original native counts and dimensions at their actual local buffers and padding. They include all response/filter/mean/pair/calibration/origin/private/coarse/fill replays.

    Q_LAW-node = Q_Mean+Q_Gram+Q_K+Q_3+J+1+Q_capture,
    d_LAW-node = (J+3)D+d_Mean+d_Gram+d_K+d_3,
    Q_RAW-node = Q_FQ+J+1+Q_capture,
    d_RAW-node <= (J+5)D.

Here a LAW node has Y, keep and terminal roots plus J bridge roots; the RAW source retains its original four-D endpoint/tail bank plus terminal and J bridge roots. Endpoint quadrature values can reduce J but are conservatively counted. After source-zero carrier alignment,

    d_T = D + sum_nodes(d_node-D).

Every changed caller or private argument reexecutes (5), all original sites (6), the entire affected native graph, origin and numerical version. No path, covariance, mean, tensor or sample is cached across changed arguments. Requested first/adjoint sweeps use original HVPs only at already recorded VALUE sites; HVPs are not producer leaves.

The new J is O(A^(-1/5) log(1/A)+log^2(1/A)), before absolute-precision logarithms. This inverse-A count is disclosed and not renamed a public logarithm. All imported native counts remain their literal functions of actual buffers, dimensions, radii, orders and public logs. An optional final own-mean replay multiplies the full Q_T+Q_B costs at every native residual occurrence and allocates its complete d_T tape; no uniform public-log final-work claim follows from this packet.

Freeze every clock, weight, scalar row/root, share, rotation, origin and numerical source version. Enumerate the finite downstream absolute path profiles and choose leaf tolerances so their weighted sum is <=e_abs. Include the potentially small positive innovation gaps in (5), h,r, all native shields/filter floors and complete origin restoration. No floor is divided by realized energy. Every term in (1)-(4),(9) is substantive, not a numerical tolerance.

## Scope

This source is an actual covariance-qualified RAW producer and is reusable with the stated complete tape/caller interface. It is not an assumed conditional expectation or a Gaussian-law output relabeled RAW. It breaks the old fixed37/10-to15/4 ledger by paying a different, inverse-A source cost. A public-log covariance-qualified producer, and a proved order/cost/radius recurrence for arbitrary order, remain open. The separate martingale-square note identifies a positive original-VALUE route to that stronger port but does not certify its unresolved clock compression.
