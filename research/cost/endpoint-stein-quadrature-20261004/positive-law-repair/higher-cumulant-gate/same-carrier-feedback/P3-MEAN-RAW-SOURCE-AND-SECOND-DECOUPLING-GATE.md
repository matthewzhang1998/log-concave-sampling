# A finite P3 mean source and the precise second-decoupling gate

2026-10-05. Constructive raw-source contract and exact analytical current. This note does not assert that the raw source has the correct P3 mean, or that a fourth-order same-carrier mean service is complete. The proposed bias correction in Section 4 remains unproved.

## 1. Exact target, unchanged

On the actual Gaussian OU history conditioned on X1=x, let F2 be the second-substitution force. The third-substitution conditional mean is

    m3(Z)=R1 psi2(Z), psi2(x)=E[g(x−F2)|X1=x].

This is the Gaussian-history conditional target, not a posterior conditional-force mean. Its exact analytical P3 reference has endpoint error A4 sqrt(D), but no history is executed for free.

## 2. Literal finite three-root raw source

Freeze three positive interior clock rules: (w_i,r_i), (v_j,tau_j), and (u_k,sigma_k), each of mass one and exact first moment one half. Write their respective complementary square roots c_i,d_j,e_k. At captured standard Z, use three independent D-roots G,H,J, shared across their entire levels:

    x_i=r_i Z+c_i G,
    y_ij=tau_j x_i+d_j H,
    z_ijk=sigma_k y_ij+e_k J,
    L_ij=sum_k u_k g(z_ijk),
    I2_i=sum_j v_j g(y_ij−L_ij),
    F3_Q(Z,G,H,J)=sum_i w_i g(x_i−I2_i).             (1)

All descendants remain on the same literal record. This is a new finite VALUE source, not an execution of a continuous path. No claim that its common-root genealogy equals the Markov genealogy is made.

Let B_Q=sum_i w_i g(x_i), E_Q=F3_Q−B_Q. The baseline is a genuine gradient in G. The exact split uses the same x_i and the same complete ancestors; independent complete mean banks may be introduced only after that split identifies their target means.

The positive weights and bounded Gaussian rows give

    ||E_Q||_(Lp|Z)<=C_p A2(|Z|+sqrt(D)),
    first(E_Q)<=C A,
    curl(P_G*E_Q)<=C A2,
    captured-Z first(E_Q)<=C A,                     (2)

where P_G=(I,0,0). In detail, the inner I2 map has x first at most A(1+A/2)/2, H first O(A), and J first O(A2). For one outer term,

    D_G E_i=c_i[(Dg(x_i−I2_i)−Dg(x_i))
                                  −Dg(x_i−I2_i)D_x I2_i].

The first Hessian difference is symmetric regardless of its magnitude. The remaining G-G skew is O(A2), the H off-diagonal block is O(A2), and the J block is O(A3). This proves the stated small curl without differentiating a Hessian. The actual full first stays O(A).

At fixed nonzero Z, B_Q(Z,0) and E_Q(Z,0,0,0) are executed caller-only origins. They must be subtracted/restored under the identical saved version. They are not set to zero by a Gaussian argument. At total Z=G=H=J=0 the source zero follows from g(0)=0 and coherent same-site anchor reuse.

The raw counts are

    Q_B=N_out,
    Q_E=N_out[N_mid(N_in+1)+2].                      (3)

A completed gradient/near-gradient mean with fixed positive shares, gradient order at least four and padding mu=A therefore has the conditional target N(E F3_Q|Z,I) with error Lambda A4(|Z|+sqrt(D)), subject to the actual C A radius, C A2 curl, active dimension, finite-clock, precision and caller guards. This is the same admitted source-class argument as the earlier nested mean, now on a literal 3D private tape. Its safe complete count is

    Q_captured+N_B Q_B+N_E Q_E+known/numerical/replay work. (4)

Every changed raw argument replays every original inner g site. Requested first/adjoint sweeps use original HVPs only at these recorded VALUE sites. No stored HVP is differentiated.

This is only a mean-law compiler for the raw source's OWN mean. The target-bias gate below is not supplied by (2)-(4).

## 3. Exact innovation current for the missing bias

The independently audited innovation localization defines

    j(x)=E[g(x−I)|X0=x],
    H_j=integral_0^infinity e^(−t)j(X_t)dt,
    F2=H_j+S, E[S|x]=0, Cov(S|x)<=A4 I/8.

Put F_theta=H_j+theta S, theta in [0,1], and psi_j(x)=E g(x−H_j). The exact first-derivative identity is

    m3−R1 psi_j
       =−R1 integral_0^1 E[Dg(x−F_theta)S|x]dtheta.  (5)

This is legitimate under C2 and keeps the true same history. A mere Lipschitz bound on (5) costs O(A3 sqrt(D)), one order too large. Zero conditional mean and small innovation covariance alone do not imply the desired A4 force bound.

For use as a WHOLE analytical current, let B_s be the conditional OU Brownian motion. The earlier Clark–Ocone proof gives

    S=integral_0^infinity L_s dB_s,
    ||L_s||op<=A2 s e^(−s)/sqrt(2).                  (6)

The literal first-program derivatives obey

    ||D_s H_j||op<=Lip(j)e^(−s)/sqrt(2),
    ||D_s F2||op<=A e^(−s)/sqrt(2)
                       +A2(s+1/2)e^(−s)/sqrt(2).   (7)

The latter follows by differentiating only the VALUE graph: the direct X_t term integrates over t>=s; the I_t ancestor derivative is bounded by A e^(−|s−t|)/sqrt(2). The clock identity integral_0^infinity e^(−t)e^(−|s−t|)dt=(s+1/2)e^(−s) closes the feedback bound.

For smooth regularizations ONLY, (5) can be regrouped as

    m3−R1 psi_j
       =R1 integral_0^1 E[D2g(x−F_theta):T_theta|x]dtheta,
    T_theta=integral_0^infinity L_s(D_s F_theta)*ds.  (8)

The coefficient has operator norm O(A3) and HS energy O(A3 sqrt(D)), uniformly in the regularization. Formula (8) is a complete law/current identity; D2g is not assumed bounded, is not an executed producer, and is not asserted to converge as a separate pointwise factor under C2. The original expression (5) remains the target definition.

The potential closure is to transfer the remaining derivative through a Gaussian score that preserves every coefficient ancestor. Each positive-time coefficient vertex depends only on the history from its first time onward (and on its correctly conditioned future copies). The augmented early-bridge score used in the audited K reduction can preserve those future ancestors. However, a vector trace in (8) requires an actual one-marked-energy/Hermite contraction bound. Naive multiplication of two vector energies can introduce D rather than sqrt(D). No such missing estimate is hidden in this note.

## 4. A narrower leading-covariance correction candidate

Let the finite inner I2_Q(x;H,J) be the two-level source in (1), before the final outer g. The existing nested-mean theorem gives its conditional mean error relative to true F2 at L2(gamma) order A3. Also, one-marked-energy covariance estimates give

    Cov(F2|x)=C_H(x)+O_L2HS(A3 sqrt(D)),
    Cov(I2_Q|x)=C_cheap(x)+O_L2HS(A3 sqrt(D)),

where C_H is the true conditional Markov covariance of I, and C_cheap is the covariance of the finite one-level common-H packet sum_j v_j g(tau_j x+d_j H).

This suggests the bounded mean-current target

    m3 ?= R1 E[g(x−I2_Q)|x]
       +(1/2)R1[(C_H−C_cheap):D2g]
       +O(Lambda A4 sqrt(D)).                        (9)

The sign is TRUE covariance minus CHEAP covariance. This candidate does not require the full order-three covariance correction K for an order-four force mean: an A3 covariance allowance would already be multiplied by the outer force's A scale after a valid one-energy heat estimate.

Equation (9) is NOT proved here. It needs a C2-safe conditional second-moment expansion with coherent means retained, actual Gaussian genealogy/cuts, and a uniform one-energy remainder. A pointwise Taylor expansion of Dg, or an argument based only on moment matching, is invalid. Centering around a common conditional mean may be necessary to avoid a spurious product of dimension-sized deterministic means. Any outer/earliest-clock cut must be analytical with a priced tail and a finite positive quadrature implementation; it must not resume the paused ordinary path grid.

If this current is established, its leading coefficient has only the original affine g covariance genealogies: true C_H is the independently admitted square-clock covariance; C_cheap is a finite genuine-gradient covariance source. The D2g current still requires a proper original-VALUE response/marked-law consumer with all ancestor cuts preserved. It cannot be invoked as a tensor or Hessian producer.

## 5. Current stopping boundary

The exact target and the raw finite source/mean compiler are explicit. The unresolved step is the source-qualified conditional second-decoupling/current comparison from (1) to m3, potentially via (5) or (9). The separately closing K covariance service does not by itself establish this mean. The full order-four endpoint remains unproved until this gate and the full final feedback/restoration join receive independent review.
