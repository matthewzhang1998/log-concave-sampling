# A third-order VALUE-only reserve for the exact reverse-OU curvature

2026-10-04. Constructive candidate; independent audit pending. This is a Gaussian curvature reserve, not a complete non-Gaussian reverse transition or an all-order recurrence.

## 1. Target, scope, and finite source

Capture the caller a=rz, r, t, original anchored gradient g=grad U, and every source version before sampling. Assume g(0)=0, 0<=Dg<=A I, 0<A<=1/2. Write

    s²=1-r², Delta=t²-r²<=1/4, 0<=r<t<=1,
    v0=Delta/t², alpha=A s².

Thus s>0 and Delta<=s². The exact conditional law is nu(dx) proportional to exp(-|x-a|²/(2s²)-U(x)) dx. Its curvature is

    H=E_nu Dg-Cov_nu(g)=s^(-2) Cov_nu(X,g(X)),
    0<=H<=A I.

Cov(X,Y) means E[(X-EX)(Y-EY)^T]; Sym(B)=(B+B^T)/2. Every target below includes the full true conditional covariance. Nothing replaces nu by a Gaussian without a displayed restoration price.

The claimed normalized reserve is an actual finite VALUE map R with

    W2(Law(R|a,r,t),N(0,(1-Delta)I-Delta H))
       <= Lambda sqrt(D) Delta A³ s⁴ + restored finite floors.       (R)

Multiply R by sqrt(v0) for the reserve required by the reverse-OU gluing ledger. The imported finite-gradient and near-gradient compilers require their actual small-radius/public-log guards; A<=1/2 alone does not imply those guards.

Use the exact finite conditional mode and packet of the audited conditional-velocity service: x_0=a, x_(j+1)=a-s²g(x_j), capture x=x_M and g(x). Draw W=(u,v) standard in R^(2D). For the positive dyadic quadrature (w_i,r_i), put q_i=sqrt(1-r_i²), R_i=(r_i I,q_i I), and

    gbar(y)=s[g(x+s y)-g(x)],
    K(W)=sum_i w_i gbar(R_i W),
    Q(W)=x+s[u-K(W)],
    F(W)=g(Q(W)), B(u)=g(x+s u), E=F-B,
    f(W)=[F(W)-g(x)]/s, b(W)=[B(u)-g(x)]/s.

K is the same two-root packet, with the same v at every node. The source service proves

    W2(Q,nu)<=C A²s⁵ sqrt(D)+eta_Q,
    Lip Q<=C s, Lip K<=alpha, ||K||p<=C_p alpha sqrt(D),
    ||E/s||p<=C_p A²s² sqrt(D),
    Lip f<=C A,
    ||D(P^*f)-D(P^*f)^*||op<=C A alpha,
    P=(I,0), PP^*=I.

All raw sources vanish exactly at W=0 with recorded same-site anchor reuse. No Hessian value is executed to form them. eta_Q contains the actual conditional finite-mode residual, original-mode linear-tilt restoration where applicable, and state numerical floors.

## 2. Restore the true conditional covariance with one marked energy

For random vectors A_1,B_1,A_2,B_2 on any common coupling,

    ||Cov(A_1,B_1)-Cov(A_2,B_2)||HS
      <= sqrt(||Cov(B_1)||op) ||A_1-A_2||2
          +sqrt(||Cov(A_2)||op) ||B_1-B_2||2.              (1)

Subtract means in the differences if desired; this only improves the bound. The inequality ||E[ab^T]||HS<=||a||2 sqrt(||Cov b||op), for centered b, follows by summing first-chaos/Bessel dual bounds over the rows of a. It does not multiply two Euclidean sqrt(D) energies.

Gaussian Poincare for Q,F and conditional Poincare for X,g(X) give covariance operator bounds C s² and C A²s². Couple Q to X optimally and use |g(Q)-g(X)|<=A|Q-X|. Therefore

    H_Q=s^(-2) Sym Cov(Q,F),
    ||H_Q-H||HS<=C A³s⁴ sqrt(D)+C(A/s)eta_Q.             (2)

Gaussian integration by parts in u, with the displayed covariance orientation, gives the exact identity

    H_Q=Sym E D_u f-Sym Cov(K,f).

Write f=b+E/s and define

    H_star=Sym E D_u f-Sym Cov(K,b).                     (3)

The actual omitted covariance is retained-and-paid, not asserted zero:

    ||Cov(K,E/s)||HS
       <=sqrt(||Cov K||op) ||E/s||2
       <=C A³s⁴ sqrt(D).                                (4)

Thus H_star approximates the full exact H to the allowance in (2). No derivative of the small mean error or the small E VALUE error has been used.

## 3. Exact full-gradient lift and balanced polarization

All r_i are strictly positive. Let

    L_Q=sum_i w_i/r_i,
    G_H(W)=sum_i (w_i/r_i) R_i^* gbar(R_i W),
    G_B(W)=P^* b(W).

G_H and G_B are genuine full gradients on the SAME W. They satisfy

    P G_H=K, P G_B=b,
    Lip G_H<=alpha L_Q, Lip G_B<=A.

Indeed each summand is the gradient of (w_i/r_i)[U(x+sR_iW)-s g(x) dot R_iW]. G_B is the gradient of [U(x+s u)-s g(x) dot u]/s². These potentials are analytical certificates; the producer evaluates only g.

The dyadic rule is Gauss-Legendre on t=1-r panels. The smallest r occurs on t in [1/2,1], and r_min=(1-x_max)/4 for the largest zero x_max of P_m. Markov's inequality gives max|P'_m|<=m² and hence 1-x_max>=1/m². Thus r_min>=1/(4m²), L_Q<=4m². There is no inverse-alpha query replication here; this lift has a public-log coefficient bill. The tiny endpoint panel near r=1 is harmless for 1/r.

Use the balanced full gradients

    J_+=s G_B+G_H/s,  J_-=s G_B-G_H/s.                  (5)

Their full firsts and centered/anchored Lp energies are at most Lambda A s and Lambda A s sqrt(D). Their exact zeros are zero. The balancing is essential: an unbalanced covariance action would lose the s powers needed below.

The exact same-tape covariance polarization is

    P[Cov J_+-Cov J_-]P^*=4 Sym Cov(K,b).                (6)

Independence may be imposed between COMPLETE + and - action banks only. Within each J occurrence, G_H and G_B use the same W and anchors. This retains all cross-node and mixed ancestor products in (6).

## 4. The finite positive derivative block

Set v_D=(1-Delta)/2 and v_+=v_-=(1-Delta)/4. These are fixed numerical positive shares bounded away from zero over Delta<=1/4. Let s0 in (0,1/4) be the fixed scale in the admitted near-gradient VALUE pair. Form the square source on R^(2D)

    h=(Delta/(s0 v_D)) P^* f.                            (7)

Run the actual near-gradient VALUE pair Y=Pair_v(h;p) with full incoming p standard in R^(2D), retained through the pair comparison. Return

    R_D=sqrt(v_D/2) P(p-Y).                              (8)

The pair target has conditional law N(Mp,I-MM^T), M=s0 E Dh. On integrating p, direct joint Gaussian algebra gives

    Cov R_D^ref=v_D I-Delta Sym E D_u f.                 (9)

This is an exact linear-covariance extraction from a native retained-input pair. It is not a square interpreted as a transpose and does not replace the pair's actual carrier by its comparison noise. The whole p is retained until (8) is formed.

At the actual scaling (7), the source has radius ell<=Lambda Delta A, curl relative a_c<=Lambda alpha, and centered energy <=Lambda Delta A sqrt(D). Set the finite pair padding mu=alpha. Its literal admitted error is bounded by

    Lambda sqrt(D) [Delta² A²(alpha+alpha)
                   +Delta⁴ A⁴(1+alpha^(-1/2))]
       <=Lambda sqrt(D) Delta A³s⁴.                    (10)

For the last inequality use Delta<=s² and A<=1/2 in each term. In particular the inverse-s term is Delta⁴ A^(7/2)/s; its ratio to Delta A³s⁴ is Delta³ sqrt(A)/s⁵<=sqrt(A)s. There is no inverse-s replication.

The actual source and returned residual firsts are respectively Lambda Delta A and Lambda[Delta A+Delta² A²/sqrt(alpha)]<=Lambda Delta A. The literal original caller-a first is <=Lambda Delta A/s<=Lambda A s. This is a graph bound, not a derivative of the law error.

## 5. Two signed full-gradient covariance-action blocks

Use the admitted finite forward covariance action C_cov(J;p) on a full-gradient source J, with padding mu=alpha, all its finite positive covariance clocks, and its original VALUE response graph. For a full gradient its analytical target matrix is EXACTLY Cov J, since every averaged Jacobian is symmetric. Its error/energy/first bounds are

    ||E_private C_cov(J;p)-(Cov J)p||L2_p
          <=Lambda ell e mu+e_num,
    ||C_cov||2<=Lambda ell e+e_num,
    Lip C_cov<=Lambda ell²/sqrt(mu),

where ell=Lip J and e=||J-EJ||2. These are the literal t30 forward-action port, including full response roots, finite first-chaos filters and unexpanded positive covariance clocks. They are not covariance-oracle instructions.

For sigma=+1 and -1, use J_sigma from (5), draw independent full roots p_sigma,z_sigma and all fresh action tapes, set eta_sigma=zeta_sigma=sqrt(v_sigma/2), and execute

    T_sigma=eta_sigma p_sigma
             +sigma [Delta/(8 eta_sigma)] C_cov(J_sigma;p_sigma)
             +zeta_sigma z_sigma,
    R_sigma=P T_sigma.                                  (11)

The independent-buffer current estimate used in the admitted near-gradient mean proof gives a law comparison to the Gaussian with linearized covariance

    v_sigma I+sigma (Delta/4) P Cov(J_sigma)P^*,          (12)

with error at most

    Lambda [Delta ell e mu
              +Delta² ell³ e(1+mu^(-1/2))]+e_num.        (13)

For explicit provenance, first condition on p_sigma and interpolate the private action to its conditional mean; the independent z_sigma buffer prices its fluctuation by (complete private first)*(action energy), giving the mu^(-1/2) term in (13). Calibration gives the first term. The deterministic linear Gaussian has an extra positive quadratic covariance [Delta/(8eta_sigma)]²(Cov J_sigma)²; the fixed-gap Gaussian-root inequality prices it by the remaining Delta² ell³e term. No arbitrary random covariance mixture is replaced merely by its expectation.

A short standalone form of this private-action estimate is useful. If X=h(V) is centered, V is standard Gaussian, and Lip h<=L, the Gaussian Riesz representation supplies a Stein matrix tau for X with ||tau||L2(HS)<=L||X||2. Along Y_t=zeta Z+t c X, first use this Stein identity and then Gaussian integration by parts in the independent Z. Its continuity velocity is (t c²/zeta)E[tau Z|Y_t], with L2 norm at most t c² L||X||2/zeta. Integrating t gives W2(Law(zeta Z+cX),Law(zeta Z))<=c²L||X||2/(2zeta). Apply at each fixed p_sigma to C_cov minus its private conditional mean, then square and integrate its conditional energy. Centering reduces energy and does not change the private first. This supplies the exact one-marked-energy factor in (13).

Here ell<=Lambda As, e<=Lambda As sqrt(D), mu=As². Thus (13) is at most

    Lambda sqrt(D)[Delta A³s⁴
       +Delta² A⁴s⁴+Delta² A^(7/2)s³]
      <=Lambda sqrt(D) Delta A³s⁴.                      (14)

The two signed reference covariances (12) are uniformly positive once the actual Lambda As radius guard is applied; for example require Delta ell_J²<=2v_sigma using the declared deterministic bound ell_J>=Lip J_sigma. This sufficient guard does not query Cov J_sigma. All implemented programs are positive maps regardless; the gap is also needed for their Gaussian comparison and numerical restoration. The actual action residual first in (11) is <=Lambda Delta A^(3/2)s<=Lambda Delta A, and its caller-a first <=Lambda Delta A^(3/2). No derivative of a saved first action is taken.

## 6. Join, exact centering, and target restoration

Execute (8), (11)+, (11)- on independent COMPLETE banks, conditional only on the same captured caller and parameters. Return R_raw=R_D+R_++R_-. Their independent Gaussian references have covariance

    (1-Delta)I-Delta Sym E D_u f
           +(Delta/4)P[Cov J_+-Cov J_-]P^*
       =(1-Delta)I-Delta H_star.                        (15)

The target in (15) is uniformly gapped: ||H_star||op<=A(1+alpha)+A alpha from Gaussian Poincare, so a numerical universal gap holds already for A<=1/2, Delta<=1/4. The stricter individual-source compiler guards remain mandatory.

If literal zero mean is required, return (R_raw^(1)-R_raw^(2))/sqrt(2) from two independent COMPLETE copies. It has exactly zero mean, preserves the target Gaussian, and does not enlarge the L2 coupling error (the squared difference becomes the variance of the prior coupling error). Counts double by this known constant. All firsts remain bounded up to an absolute factor.

Equations (2)-(4), the fixed-gap Gaussian-root inequality, and (10),(14) prove (R), with additional normalized covariance restoration C Delta(A/s)eta_Q. After multiplication by sqrt(v0), this becomes C sqrt(v0)Delta(A/s)eta_Q. For an original finite-mode tilt, restore its effect in eta_Q against the actual tilted conditional law before applying (1); do not silently reuse the anchored target.

This supplies a covariance reserve independent of the completed unit-buffer conditional mean service. Joining sqrt(v0)R to cz-d M_r has the exact-mean/exact-covariance Gaussian REFERENCE of the reverse transition, with error <=d epsilon_M+sqrt(v0)epsilon_R. It still does not reproduce the true conditional cumulants of order >=3. No all-order reverse transition or sublinear complete-cost recurrence is claimed.

## 7. Original-query, caller, zero, and numerical ledger

The caller-only mode and last anchor cost M+1 original g VALUES and may be cached only under the identical complete captured key. A raw f occurrence costs n_alpha+1 private original VALUES (all K nodes and the terminal Q). A raw J_+ or J_- occurrence costs n_alpha+1 private original VALUES (the same inner points plus the baseline x+s u), with original g(x) reused from the captured record. G_H's full auxiliary output is computed and retained; it is not dropped inside its full-gradient compiler.

If N_pair and N_cov,+, N_cov,- are the actual complete raw-source occurrence counts of the admitted finite pair and covariance actions, the safe total before optional exact centering is

    M+1+(n_alpha+1)(N_pair+N_cov,++N_cov,-)
         +known/numerical work.                         (16)

This expression includes every nested source response; a changed raw input is a new complete K graph. Identical source ancestors are reused only within one literal same-record raw occurrence. No random K value is cached across independent complete banks. Each requested first/adjoint sweep has the corresponding original HVP sites; a discarded primal graph pays its replay. No producer evaluates an HVP.

Private tape dimension is D times the full pair/action occurrence counts and clock/fill banks, with each raw source having two D-roots. It is not just 2D. All finite filter orders, quadrature versions, padding, variance shares and tolerances are frozen before differentiation. At fixed requested target grade, the original occurrence counts and coefficients are public-log polynomials, subject to the same finite-clock hypotheses as the imported compilers. In particular n_alpha=O(log²(1/alpha)); exponentially tiny s outside the admitted clock retains its literal log(1/s) bill.

Raw numerical normalized g-value errors are multiplied by 1/s in f and by L_Q in G_H/s; these are absolute finite precision weights, not inverse-s source-copy counts. The outer Delta factors and the covariance-response cancellations must be propagated through the actual complete finite graph before choosing each leaf tolerance. There is no division by a possibly zero actual source energy. Original physical sites are b_M+sqrt(A)[x+s R_i W], b_M+sqrt(A)Q, and b_M+sqrt(A)[x+s u], together with the complete finite response substitutions at those sites.

Every raw zero is literal under recorded-anchor reuse. The completed source-zero carrier consists of the actual pair's independent source-zero carrier combined with its retained incoming p, and the known eta p+zeta z carriers of the two covariance actions. It is never identified with the Gaussian reference of a law coupling. Captured x and g(x) retain their actual caller derivatives. For the normalized caller a, Dx<=2; D_a K<=C As; therefore scaled (7) has caller <=Lambda Delta A/s. In J_sigma the s balancing leaves caller <=Lambda A. For the original physical caller, replace the normalized source caller A by the imported sqrt(A) profile at each leaf: the complete reserve caller is <=Lambda Delta sqrt(A)/s<=Lambda sqrt(A)s, plus the actual finite-mode/numerical profiles. These are absolute first bounds, not small derivatives of the covariance approximation.

## 8. Finite diagnostics and independent checks

The exact matrix quadratic test should give H_star=B-s²B² and H=B(I+s²B)^(-1), retaining all shared-G products in Q, f and K. The error is <=A³s⁴ sqrt(D), with no scalar-only restriction. Genuine C2 noncommuting Hessian fixtures must verify the same-source derivative/polarization identities, true conditional covariance target, and the one-energy restoration numerically. The finite compiler implementations are imported proved VALUE programs; diagnostic Hessians used in these tests are not producer leaves.

The author script check_value_curvature_reserve.py reports 908 passing assertions. It includes 64 non-diagonal quadratic matrix/heat cases in dimensions 1,2,5,17; all shared-root covariance and balanced-polarization identities; positive full-matrix targets; heat-exponent inequalities at s down to 1e-4; and a genuine C2/non-C3 two-dimensional ridge fixture with noncommuting Hessians. At two tensor-Gaussian integration resolutions, twelve nonlinear curvature-error ratios ||H_star-H||HS/(A³s⁴sqrt(2)) lie between 0.0486 and 0.0674. These finite quadrature numbers are diagnostics rather than an asymptotic proof. The exact first-chaos coefficient is integrated through Cov(u,f), avoiding an unrecorded derivative-quadrature assumption. Actual source Hessians are used only to check C1/curl/gradient certificates; they are not producer actions. The huge imported finite pair/action compilers are not end-to-end numerically instantiated by this script.
