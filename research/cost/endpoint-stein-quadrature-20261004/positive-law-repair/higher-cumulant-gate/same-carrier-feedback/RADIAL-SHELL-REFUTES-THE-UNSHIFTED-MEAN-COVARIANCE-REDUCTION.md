# A radial-shell same-source counterexample to the unshifted mean reduction

2026-10-05. Uniform-C2 separator. Independently derived asymptotics; exact-text review requested. Unlike the earlier port-only example, every force and covariance here comes from the SAME anchored gradient. The centered/resummed current is not refuted.

## Result

Let I_H be the true conditional Gaussian OU force integral and I_Q the actual finite positive common-root packet, with the admitted Hermite clock tolerance delta<=A2 and exact first moment 1/2. Put

    j_a(x)=E[g(x−I_a)|x], C_a(x)=Cov(I_a|x), a=H,Q.

The proposed unshifted reduction

    R1(j_H−j_Q)
       ?=(1/2)R1[(C_H−C_Q):D2g(x)]+O(Lambda A4 sqrt(D))          (1)

is FALSE uniformly in the original C2 Hessian class. There is a C-infinity, globally anchored, monotone gradient family with A=D^(−1/2) for which the L2(gamma_D) residual in (1) is at least c_* A2. The proposed allowance is A4 sqrt(D)=A3, and any fixed polynomial-log multiplier still fails.

The mechanism is a coherent radial mean displacement combined with the high-dimensional covariance trace. The correct leading force derivative is evaluated at the radius shifted by a fixed amount. Expanding it back at the original x loses an order. This does not refute a centered/resummed target or the literal finite F3_Q source.

## 1. One genuine gradient for all branches

Choose a nonnegative symmetric C-infinity bump b supported in [−2,2], strictly decreasing on (0,2), with integral one. For example, normalize exp(−1/(1−s2/4)) on |s|<2. Let

    psi(s)=integral_(−infinity)^s b(u)du,
    M_k=||psi^(k)||infinity.

Then psi is zero below −2 and one above 2. Fix c=1/4 and any fixed 0<d<=1/[8(1+M1)]. For integer D sufficiently large, set

    R=sqrt(D), A=1/R, n_x=x/|x|,
    f_D(x)=psi(|x|−R)n_x,
    g_D(x)=A[c x+d f_D(x)].                           (2)

Set f_D=0 near the origin, as its displayed formula already prescribes when R>=4. It is globally C-infinity. It is the gradient of a smooth radial potential, and g_D(0)=0.

The radial and tangential Hessian eigenvalues of g_D are

    A[c+d psi'(|x|−R)],
    A[c+d psi(|x|−R)/|x|].                            (3)

They are nonnegative and at most A for the fixed c,d and all sufficiently large D. Thus every branch uses one original gradient satisfying 0<=Dg_D<=A I, with no A sqrt(D) restriction.

## 2. Exact decomposition of the two inner sources

At captured x, the conditional OU history is X_t^x=e^(−t)x+L_tV. Its linear integral is Gaussian with conditional mean x/2 and covariance I/4. For the finite common-root rule, let

    beta_Q=sum_j v_j sqrt(1−tau_j2),
    sum_j v_j=1, sum_j v_j tau_j=1/2.

Both inner sources have the exact representation

    I_a=c A[x/2+sigma_a G_a]+d A K_a,
    sigma_H=1/2, sigma_Q=beta_Q,                       (4)
    |K_a|<=1.

G_a is a standard Gaussian conditional on x, possibly correlated with K_a. K_H is the positive integral of f_D along the true history; K_Q is its finite positive common-root sum. Their conditional means obey the admitted one-variable Hermite multiplier bound

    ||E K_H−E K_Q||_(L2(gamma_D))<=delta ||f_D||2<=delta.          (5)

No joint-history approximation is used.

The covariance gap stays away from zero uniformly over the actual finite rules. Since sqrt(1−tau2)>=1−tau2 and the second-moment error is at most delta,

    beta_Q>=1−sum_j v_j tau_j2>=2/3−delta.             (6)

For small delta, beta_Q2−1/4 has a fixed positive lower bound. No continuum common-root limit or beta_Q convergence theorem is needed.

## 3. Gaussian radial expansion with the mean shift retained

Write r=|x|, S_D=r−R, n=n_x and

    y0=(1−c A/2)x,
    y_a=y0−c A sigma_a G_a.

All following conditional-expectation errors are measured in L2(x~gamma_D), with constants uniform for sigma_a in [1/2,1]. Standard Gaussian norm tails let us work on r>=R/2, with exponentially small discarded contributions.

The displacement zeta=−c A sigma_a G_a has radial component n dot zeta of Lp size O(A), tangential squared norm of size O(1), and

    |y0+zeta|−|y0|
       =n dot zeta+|P_n zeta|2/(2|y0|)+O_Lp(A3),
    E[|y0+zeta|−|y0| | x]
       =(c2 sigma_a2/2)A+O_L2(A2).                  (7)

Here P_n=I−nn*. The radius fluctuation is O_Lp(A), despite the tangential displacement norm being O(1). Its direction changes by O_Lp(A). Rotational symmetry about n makes the conditional vector expectation parallel to n, and its directional attenuation contributes only O(A2). Bounded derivatives of psi therefore give

    E[f_D(y_a)|x]
       =f_D(y0)+(c2 sigma_a2/2)A psi'(S_D−c/2)n
                                            +O_L2(A2).          (8)

The argument S_D−c/2 is essential: |y0|−R=S_D−c/2−(cA/2)S_D.

For clarity, (7) follows by expanding the Euclidean norm at radius comparable to R. The cubic remainder is controlled by |n dot zeta||zeta|2/R2+|zeta|4/R3. Their Gaussian moments are O(A3). The exceptional event |zeta|>r/2 has super-polynomially small probability and bounded f_D; it does not change the stated estimates.

## 4. The nonlinear inner part cancels at leading order

The derivatives of f_D satisfy ||Df_D||op<=C and ||D2f_D||op<=C, uniformly in D. Moreover the radial/directional geometry above gives

    ||Df_D(y_a)−Df_D(y0)||_(Lp(x,G_a);op)<=C_p A.     (9)

The same Gaussian marginal bound applies even though K_a is correlated with G_a, since |K_a|<=1. Taylor's formula with the bounded second derivative gives

    E[f_D(y_a−dA K_a)|x]
       =E[f_D(y_a)|x]−dA Df_D(y0) E[K_a|x]
                                                +O_L2(A2).       (10)

The two conditional K means cancel up to (5). The linear c A part of g also sees only their mean difference. Combining (4)-(10),

    j_H−j_Q
       =k_D A2 psi'(S_D−c/2)n+O_L2(A3+A2 delta),
    k_D=(d c2/2)(1/4−beta_Q2).                       (11)

By (6), k_D is negative and bounded away from zero; it is also uniformly bounded in magnitude. This derivation keeps the complete original nonlinear I_a. It does not replace it by its linear part without paying the contribution: (5),(9),(10) are the required cancellation and remainder estimates.

## 5. Actual covariance and the unshifted current

Gaussian first-chaos Bessel, conditional on x, yields

    ||Cov(G_a,K_a|x)||HS<=||K_a−E(K_a|x)||2<=1.

Also ||Cov(K_a|x)||HS<=tr Cov(K_a|x)<=1. Hence the actual covariance in (4) satisfies the UNIFORM-in-x estimate

    C_a=c2 A2 sigma_a2 I+E_a,
    ||E_a||HS<=C A2.                                (12)

This is a one-energy bound and uses no false independence of G_a and K_a.

For any symmetric matrix B, differentiating the radial field gives exactly

    D2 f_D(x):B
      =psi''(r−R)n(n*Bn)
       +a(r)[n tr(P_n B)+P_n(B+B*)n],
    a(r)=psi'(r−R)/r−psi(r−R)/r2.                    (13)

Thus the HS-to-vector norm is at most |psi''|+(sqrt(D)+2)|a(r)|<=C on the nonzero support, and is zero inside the inactive ball. Consequently

    ||D2 g_D(x):(E_H−E_Q)(x)||<=C A3.               (14)

For the isotropic leading covariance,

    Delta g_D
       =A d[psi''(S_D)+(D−1)a(r)]n
       =d psi'(S_D)n+O_L2(A).                       (15)

Combining (12)-(15), the unshifted current in (1), before the outer R1, is

    (1/2)(C_H−C_Q):D2g_D
       =k_D A2 psi'(S_D)n+O_L2(A3).                 (16)

All derivatives in this counterexample are ordinary smooth functions; no weak-current existence issue is involved.

## 6. A nonzero R1 witness

Subtract (16) from (11) and put

    h(s)=psi'(s−c/2)−psi'(s).

The unshifted residual before R1 equals

    k_D A2 h(S_D)n+O_L2(A3+A2 delta).               (17)

Take T_D(x)=A x, whose L2(gamma_D) norm is one. Self-adjointness and R1 T_D=T_D/2 give

    <R1[h(S_D)n],T_D>
          =(1/2)E[h(S_D)|X|/sqrt(D)].              (18)

The Gaussian radial central limit theorem gives S_D -> N(0,1/2), with uniform moment bounds; |X|/sqrt(D)->1. Since h is bounded and continuous, (18) converges to (1/2)E h(S), S~N(0,1/2).

This limit is strictly negative. The bump psi'=b is nonzero, symmetric and strictly decreasing away from zero, so its convolution with the centered Gaussian density is strictly maximized at zero. Therefore E b(S−c/2)<E b(S).

Use (6), (17)-(18), the R1 contraction on the remainder, and delta<=A2. For all large D,

    ||R1(j_H−j_Q)
        −(1/2)R1[(C_H−C_Q):D2g_D]||2 >=c_* A2.     (19)

At A=D^(−1/2), this exceeds A4 sqrt(D)=A3 by an inverse-A factor. This proves the claimed same-source separator, including for the actual finite cheap rule.

## Exact scope and next target

Equation (19) rejects the UNshifted terminal argument in (1), which was the exploratory candidate in the earlier m3 notes. It does not reject a centered leading term at x−mu(x), a resummed finite shifted kernel, or the exact covariance interpolation/Stein current that keeps x−F_theta. In this family the missing displacement is visibly c/2 in the radial variable, and the centered term restores that leading location.

Every branch here uses the same original g_D, the true Markov I_H and the genuine finite common-root I_Q. No ordinary Markov path-grid algorithm is executed. The admitted finite full-covariance service remains valid. The same-carrier m3 mean and its native shifted-current consumer remain open.
