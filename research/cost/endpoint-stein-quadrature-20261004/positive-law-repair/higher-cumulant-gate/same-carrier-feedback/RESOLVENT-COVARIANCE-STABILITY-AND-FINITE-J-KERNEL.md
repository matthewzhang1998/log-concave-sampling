# Resolvent covariance stability admits a finite j coefficient kernel

2026-10-04. Constructive finite coefficient-source admission, independent review requested. This is not yet the covariance reserve or the full fourth-order endpoint. The matrix orientation, K-current, and same-carrier m3 ports remain separate.

## Result

Let R_b f=integral_0^1 r^(b−1)P_r f dr for Gaussian Mehler P_r, and define the true conditional Markov covariance functional

    C(f)=2 R2[(D R1 f)(D R1 f)*].

For Lipschitz vector fields f0,f1, put u=R1(f1−f0). Then

    ||C(f1)−C(f0)||_(L2(gamma);HS)
       <=2 C_* [Lip(f1)+Lip(f0)] ||u||2,
    C_*=(1+pi)/4.                                    (1)

This is a resolvent covariance stability theorem. It differentiates no error estimate and pays only one dimension-sized energy. It applies to non-gradient vector fields.

Let j(x)=E[g(x−I)|Gaussian OU endpoint x] be the exact first-shift field in the audited innovation localization, and let

    H_Q(x,H)=sum_l v_l g(tau_l x+sqrt(1−tau_l2)H),
    J_Q(x,H)=g(x−H_Q(x,H)),
    j_Q(x)=E_H J_Q(x,H).                              (2)

The already audited nested-mean decoupling proof gives

    ||R1(j−j_Q)||2<=C[A3+delta_in A2]sqrt(D).          (3)

Both j and j_Q have Lipschitz constant at most A(1+A/2). Hence (1) proves

    ||Cov(H_j|Z)−Cov(H_(j_Q)|Z)||_(L2;HS)
       <=C[A4+delta_in A3]sqrt(D).                   (4)

This does NOT use Cov of the cheap common-root chord. It admits its CONDITIONAL MEAN into a freshly constructed true-Markov covariance functional. The earlier shared-root covariance separator remains fully valid.

Sections 4-5 give an actual finite original-gradient VALUE source whose conditional mean Jacobian is the needed finite B_Q. It has a small but nonzero curl. Producing B_Q B_Q* by an admitted positive VALUE action is a further orientation/action gate; producing B_Q2 instead is insufficient.

## 1. Dimension-free row-Hessian smoothing

Let f:R^D->R^D have Lip(f)<=L and v=R1 f. Ordinary Gaussian Sobolev differentiation gives

    ||Dv||op<=L/2.                                   (5)

For any fixed row m in R^D, the scalar m dot f has gradient Df* m with norm at most L|m|. For r<1, q=sqrt(1−r2), one additional Gaussian integration by parts gives

    D2[m dot v](x)
      =integral_0^1 (r2/q)
          E[(Df(rx+qG)*m) G*]dr.                    (6)

For any square-integrable vector a(G), Gaussian Bessel inequality gives

    ||E[a(G)G*]||HS2<=E|a(G)|2,

because the coordinate Gaussians are orthonormal in scalar L2, separately for each output coordinate. Therefore

    ||D2[m dot v](x)||HS <= (pi/4)L|m|               (7)

pointwise, with no sqrt(D), since integral_0^1 r2/sqrt(1−r2)dr=pi/4.

In particular, for any fixed matrix M with rows m_i,

    sum_i ||D2[m_i dot v](x)||HS2
                         <=(pi L/4)2 ||M||HS2.      (8)

The row M may be chosen anew at each x when applying the pointwise inequality; its derivatives are handled separately below. Formula (6) requires only the bounded first derivative of f. Approximation and Gaussian Sobolev closure justify it for Lipschitz f. In our application f=j or j_Q is C1 with globally bounded first under the original C2 assumption on U; no D2g is queried or assumed.

## 2. Gaussian divergence estimate with no hidden trace factor

For a matrix field T:R^D->R^(D by D), define the Gaussian divergence rowwise by

    (delta T)_i=sum_k [x_k T_ik−partial_k T_ik].

It is the L2(gamma) adjoint of the Jacobian: E[Du:T]=E[u dot delta T]. The Gaussian divergence isometry and Cauchy–Schwarz give

    ||delta T||2 <= ||T||_(L2;HS)+||DT||_(L2;HS),    (9)

where the last norm sums all three derivative/output indices. This is the ordinary matrix-valued Gaussian Sobolev divergence estimate; it does not separately bound the xT and derivative terms and thereby lose their cancellation.

Let N be a smooth matrix field and T=N Dv. Then

    ||T(x)||HS<= (L/2)||N(x)||HS,
    ||DT(x)||HS<= (L/2)||DN(x)||HS
                                  +(pi L/4)||N(x)||HS. (10)

The first product-rule term follows from ||Dv||op<=L/2. For the second, freeze each row m_i=N_i(x): its tensor entries are precisely D2[m_i dot v], and (8) applies. This is the crucial one-energy row-Hessian step; replacing it by a trace or by D times an operator bound would be invalid.

Combining (9)-(10),

    ||delta(N Dv)||2
      <=(L/2+pi L/4)||N||2+(L/2)||DN||2.             (11)

## 3. Resolvent paraproduct estimate and covariance stability

Let P_t=P_(exp(−t)) denote OU time notation. For any matrix test M of L2 HS norm one, self-adjointness of P_t and Gaussian integration by parts give

    <R2(Du Dv*),M>
      =integral_0^infinity e^(−2t)
          E[Du:(P_t M)Dv]dt
      =integral_0^infinity e^(−2t)
          E[u dot delta((P_t M)Dv)]dt.               (12)

Use the dimension-free Hilbert-valued heat gradient estimate

    ||D P_t M||2 <= e^(−t)/sqrt(1−e^(−2t)) ||M||2.  (13)

It follows either by the same Gaussian Bessel inequality and invariance, or by the Hermite spectral bound. Applying (11), then integrating,

    ||R2(Du Dv*)||_(L2;HS)
      <=L||u||2 [(1/2+pi/4) integral_0^infinity e^(−2t)dt
            +(1/2) integral_0^infinity
                     e^(−3t)/sqrt(1−e^(−2t))dt]
      =[(1+pi)/4] L ||u||2.                          (14)

The singularity in (13) is integrable; its displayed integral is pi/4. Smooth compact/spectral truncations justify all steps first. The bound extends by density and the first-order Sobolev closure to the actual sources. It is an L2(Gaussian endpoint;HS) theorem, not a pointwise endpoint statement or a uniform arbitrary-caller claim.

For (1), set v=R1(f1+f0). The exact polarization is

    C(f1)−C(f0)=R2[Du Dv*+Dv Du*].                   (15)

Apply (14) and transpose invariance, with Lip(f1+f0)<=Lip(f1)+Lip(f0). This proves (1).

## 4. Apply to the actual cheap j kernel and finite clocks

The source J_Q in (2) is a literal original-g VALUE program, with the SAME inner H at every tau node. Its x derivative is

    D_x J_Q=Dg(x−H_Q)(I−D_x H_Q),
    0<=D_x H_Q<=A I/2.

Thus ||D_x J_Q||op<=A(1+A/2)=:L_j, and the skew part is at most A2. Its H first is at most A2 beta_in, beta_in=sum_l v_l sqrt(1−tau_l2)<=1. The exact field j has the same Lipschitz bound by conditional Markov differentiation. Equation (3) is exactly the pre-outer-quadrature comparison in Sections 2-3 of the pinned `../THIRD-ORDER-NESTED-FORCE-MEAN-WITH-C2-DECOUPLING.md` (SHA b056e213ac28c60869727589d3ee6d1b2c53d6666908cdbe1961fa3121161ced). It does not assert that j_Q is pointwise close to j.

Choose a positive outer rule (w_i,r_i) with mass one, exact first moment one half, and uniform Hermite moment error delta_r. Define only as a coefficient target

    B_Q(X)=sum_i w_i r_i E_(u,H)
                      D_x J_Q(r_i X+sqrt(1−r_i2)u,H). (16)

Then

    ||B_Q||op<=L_j/2,
    ||B_Q−D R1 j_Q||_(L2;HS)<=delta_r L_j sqrt(D).    (17)

No symmetry of B_Q is asserted. Take a positive q rule (a_k,q_k) with the same mass/first-moment identities and error delta_q. Define

    C_Q(Z)=2 sum_k a_k q_k P_(q_k)[B_Q B_Q*](Z).       (18)

This is PSD as a target, and 0<=C_Q<=L_j2 I/4. The same noncommutative product estimate as the audited square-clock service, with a transpose retained, gives

    ||C_Q−C(j_Q)||_(L2;HS)
       <=C(delta_r+delta_q) A2 sqrt(D).              (19)

For delta_in,delta_r,delta_q<=A2, (4) and (19) show

    ||C_Q−Cov(H_j|Z)||_(L2;HS)<=C A4 sqrt(D).         (20)

Every quadrature count is O(log2(1/A)) at this fixed grade. These L2 bounds refer to the fresh standard endpoint Z, exactly as required by the buffered-stage join.

## 5. Literal finite coefficient VALUE source and native ports

For captured X, let c_i=sqrt(1−r_i2)>0 and execute the vector source on two D-dimensional roots u,H:

    V_X(u,H)=sum_i (w_i r_i/c_i)
        [J_Q(r_i X+c_i u,H)−J_Q(r_i X,H)].             (21)

Both copies of J_Q use the SAME H and their own correct inner ancestors. This source uses only original g VALUES. At u=0 reuse the identical original sites under their complete keys, so V_X(0,H)=0 literally for every H, including frozen numerical versions with coherent reuse. It does not query or subtract the uncomputed expectation j(0).

Its mean u-Jacobian is exactly B_Q(X). Its actual native ports are

    ||D_u V_X||op<=L_j/2,
    ||D_H V_X||op<=2 A2 beta_in sum_i w_i r_i/c_i
                               <=C A2,
    |V_X(u,H)|<=L_j |u|/2,
    ||V_X||Lp<=C_p A sqrt(D),
    ||D_X V_X||op<=2L_j sum_i w_i r_i2/c_i<=C A.     (22)

The aggregate coefficient sums are uniformly bounded for the admitted positive dyadic Gaussian rule, by the exact same dyadic square-root argument as the existing square-clock source. No inverse endpoint heat factor survives those sums; the individual coefficients still enter absolute leaf precision allocation.

Under coisometry P=(I,0), the square lift P*V_X has u-u curl at most A2/2, and its off-diagonal curl blocks are controlled by D_H V_X. Therefore its complete curl is O(A2). Its full/private first is O(A), not O(A2). It is a near-gradient source, not a genuine gradient; this distinction is retained at the next action gate.

One raw source occurrence has the safe original VALUE count

    Q_V <= 2 N_r (N_in+1).                            (23)

The putative anchors J_Q(r_i X,H) depend on the PRIVATE H. They may not be moved into a global caller-only capture or cached across changed H/compiler banks. Genuine same-key aliases inside one complete record may be reused; every changed u,H,X bank replays all its inner ancestors. The raw input dimension is 2D. An outer covariance clock additionally owns G_k with X_k=q_kZ+sqrt(1−q_k2)G_k; its derivative path uses the actual D_X port in (22).

All original HVPs are first/adjoint requests at the VALUE sites in (21), including every inner H_Q ancestor; none is a producer and no saved HVP is differentiated. Numerical VALUE errors propagate through the bounded absolute coefficient sums and the actual feedback factor 1+A. All numerical floors remain absolute; original mode/caller restoration remains separate. Quadratures, coefficients, versions, tolerances and coisometry are frozen before differentiation.

If an orientation-correct action uses N_action COMPLETE V occurrences per q bank, its safe original work is

    2 N_q N_r (N_in+1) N_action
                +actual caller/known/numerical/replay work.     (24)

No value is assigned to N_action until such an action is admitted at the actual near-gradient radius, curl, source dimension, padding and caller guards. Equation (24) is a transparent recurrence, not a claim that the missing action is already supplied.

## 6. Exact remaining gates

The first finite coefficient-source port is now explicit: (21) targets (16), and its positive covariance target (18) approximates the correct H_j covariance at order four. The missing covariance ACTION must realize B_Q B_Q* with the required positive reserve and A4 sqrt(D) error. The old genuine-gradient square action realizes a symmetric square; silently using B_Q2 would incur an uncontrolled orientation/curl contribution, and its earlier A3 reserve allowance does not automatically improve to A4.

The localized mixed current 2 Sym R2 K, K(x)=Cov(I,g(x−I)|x), remains unfilled. So does the same-carrier mean m3=R1 psi2 with psi2(x)=E[g(x−F2)|X1=x]. The new theorem does not equate that Gaussian-history kernel to a posterior. It does not preserve old private roots through completed mean-law couplings.

Only after these remaining native services, their actual source-zero/caller/floor ports, complete counts and independent reviews are supplied can the full order-four endpoint join pass. The known source restrictions and analytical-only OU reference are unchanged.
