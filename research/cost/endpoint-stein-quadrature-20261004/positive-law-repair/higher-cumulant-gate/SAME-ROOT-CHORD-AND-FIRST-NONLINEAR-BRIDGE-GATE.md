# The same-root local chord and the first nonlinear reverse-bridge gate

2026-10-04. A concrete positive VALUE construction, a full-law fallback, and a failed next-order test. Independent audit pending. This is not a new general order-three transition or an all-order impossibility theorem.

## Result first

1. The existing conditional positive Stein packet can be inserted directly into the exact affine reverse-OU bridge. This gives a complete non-Gaussian conditional law comparison, including **all** higher conditional cumulants, at

       W2 <= C (Delta/t) A^2 s^3 sqrt(D) + (Delta/(t s^2)) eta_Q.

   Its query exponent is unchanged. This is a usable second-order full-law fallback, not an order-three upgrade of the existing endpoint theorem.

2. A literal local VALUE chord, using the same two roots and shifted original query sites, removes the packet's entire matrix-quadratic second-order covariance defect at constant extra query cost. Unlike the older g(g) amplifier, it differentiates in the correct local direction in its amplitude expansion.

3. That stronger candidate still fails the next general nonlinear order. For the anchored smooth sandwich potential

       U_A(x)=A[x^2/4 + 0.1(sin x-x)],

   its third cumulant differs from the true target by

       -0.009301002921... A^2 + o(A^2),

   even after granting exact affine mean/variance matching as a favorable diagnostic. In the legal reverse bridge r=0,t=1/2 this becomes one eighth of that defect. Uniform fourth moments turn it into an Omega(A^2) full-law W2 lower bound. Thus third-order mean and curvature services cannot repair this nonlinear current by themselves.

The first failed gate is now an explicit carrier-retained higher-current identity, not positivity, a missing scalar covariance, or an unexecuted local chord. No recurrence improving the general full-law grade is proved here.

## 1. Direct full-law gluing retains every descendant

Use the exact reverse-OU notation from the bridge ledger:

    a=rz, s^2=1-r^2, Delta=t^2-r^2,
    A_rt=r(1-t^2)/(t s^2), B_rt=Delta/(t s^2),
    v=(1-t^2)Delta/(t^2 s^2).

At the exposed caller z, the exact conditional target is

    X~nu_(r,z), nu(dx) proportional to exp(-|x-a|^2/(2s^2)-U(x))dx,
    Y_t=A_rt z+B_rt X+sqrt(v)N,

where N is independent of X conditional on z. For r<t<=1, s>0. The existing complete positive packet Q has

    W2(Law(Q|z),nu_(r,z)) <= C A^2 s^5 sqrt(D)+eta_Q,

including its conditional finite-mode, original finite-mode tilt, and numerical restorations. Execute

    Y_packet=A_rt z+B_rt Q+sqrt(v)N.                    (1)

Conditional affine coupling proves the displayed result immediately. At t=1 the N coefficient is zero; no division by that vanished reserve is used. All conditional cumulants of order k>=3 are exactly B_rt^k times those of the **actual** Q. They are retained and coupled to the true cumulants through a full-law comparison; they are not declared equal to the target cumulants.

For an exact Gaussian mean/covariance replacement, the remaining defect can instead be first order in A. The concrete test in Section 4 shows this. Consequently “order-three velocity plus order-three covariance” is not “order-three bridge.”

The complete original query count of (1) is the Q count, M+1+n_alpha, plus known arithmetic and one fresh D-dimensional Gaussian block when v>0. Its private tape consists of Q's two complete roots and this new root. Its same-caller comparison retains only z,r,t and captured finite versions, not Q's private roots. Every requested first/adjoint sweep charges the actual Q VALUE sites; discarded primal records are replayed. This elementary affine composition supplies no stronger retained/protected port than Q itself.

## 2. Execute a local chord with the correct shared-root direction

First work in a standardized anchored conditional problem with force f=grad V, f(0)=0 and 0<=Df<=alpha I. Let the admitted positive quadrature have weights w_i and nodes r_i, put q_i=sqrt(1-r_i^2), and freeze its finite version. For W=(u,v), execute

    K(W)=sum_i w_i f(r_i u+q_i v),
    K1(W)=sum_i w_i f(r_i[u-K(W)]+q_i v),
    E(W)=K1(W)-K(W),
    Y_chord=u-K(W)-lambda_Q E(W),                       (2)

where

    beta_Q=sum_i w_i q_i,
    c_Q=1/4+beta_Q^2,
    lambda_Q=2(1-c_Q).

All roots in K1 are the same original u,v, and K1 reads the **actual entire** K. Every changed force site is a new original VALUE. No expectation, Hessian, covariance, tensor, density, or inverse map is queried. Equation (2) is a positive pushforward, despite a signed deterministic chord coefficient. Signed deterministic arithmetic is not a signed probability.

Since 1/2<=beta_Q<=sqrt(3)/2, one has 0<=lambda_Q<=1. The actual bounds are

    Lip K<=alpha, ||K||p<=C_p alpha sqrt(D),
    ||E||p<= (alpha/2)||K||p<=C_p alpha^2 sqrt(D),
    Lip E<=2alpha+alpha^2/2,
    Lip(Y_chord-u)<=C alpha.                            (3)

For the E energy, each shifted node moves by r_i K and sum_i w_i r_i=1/2. Its first is not falsely bounded by alpha^2: a difference of the local Hessians can be of order alpha under C2.

### Full-gradient lift of this actual chord

Write R_i=(r_i I,q_i I), P=(I,0), Pi=P*P, and

    J(W)=sum_i (w_i/r_i)R_i* f(R_i W),
    L_Q=sum_i w_i/r_i,
    E_J(W)=J(W-P*K(W))-J(W).                            (4)

J is a genuine full gradient, PJ=K, and PE_J=E exactly. The source uses the same K and nodes as (2), retaining the auxiliary output of J. Let B0=DJ(W), B1=DJ(W-P*K(W)). Then

    DE_J=B1-B0-B1 Pi B0,
    DE_J-DE_J*= -B1 Pi B0+B0 Pi B1.                     (5)

Since P DJ=DK, both P Bj and Bj P* have operator norm at most alpha. Therefore

    Lip E_J<=2alpha L_Q+alpha^2,
    ||E_J||p<=C_p alpha^2 L_Q sqrt(D),
    ||DE_J-DE_J*||op<=2alpha^2.                         (6)

All raw zeros are literal zero with recorded-anchor reuse. The admitted dyadic rule has L_Q<=4m^2, so the full-gradient lift only adds a public-log coefficient bill. This is a real near-gradient original-VALUE source, not an HVP-valued producer. To call an admitted buffered **mean law** on it, enforce that consumer's actual normalized-radius, variance-share, padding, finite-clock, and numerical guards at the full radius 2alpha L_Q+alpha^2. Alpha<=1/2 alone is insufficient. Such a mean law does not return the carrier-dependent current needed below.

## 3. What the chord really improves on matrices

For f(x)=Bx, 0<=B<=alpha I, all matrices are polynomials in the same B and

    K=B u/2+beta_Q Bv,
    E=-(B/2)K,
    Y_chord=(I-B/2+lambda_Q B^2/4)u
               +(-beta_Q B+lambda_Q beta_Q B^2/2)v.

Its exact covariance, including every shared-v cross-node product, is

    C_chord=I-B+B^2-lambda_Q c_Q B^3
                         +(lambda_Q^2 c_Q/4)B^4.         (7)

The B^2 coefficient is exactly one because c_Q+lambda_Q/2=1. The fixed-gap Gaussian-root comparison gives

    W2(Law(Y_chord),N(0,(I+B)^(-1)))<=C alpha^3 sqrt(D).

This is a genuine one-order quadratic gain with n_alpha extra original VALUES and no new Gaussian root. It does not establish the same gain for nonlinear forces.

## 4. Gaussianization already misses infinitely many first-order cumulants

Take the scalar anchored potential

    U_A(x)=A[q x^2+epsilon(sin(kx)-kx)],
    q=1/4, epsilon=1/10, k=1.                           (8)

Its gradient is A[2qx+epsilon k(cos(kx)-1)], zero at zero. Its Hessian lies between .4A and .6A globally. This smooth potential belongs to the requested general C2 class.

Let X_A have density proportional to exp(-x^2/2-U_A(x)). For r=0 and t<=1/2, the exact reverse bridge is

    T_A=tX_A+sqrt(1-t^2)N.                              (9)

At A=0 it is standard Gaussian. Differentiating its characteristic function at zero amplitude and comparing with the Gaussian having the **exact** mean and variance gives, for k=1,

    E exp(i xi T_A)-E exp(i xi G_A)
      = i A epsilon exp(-(xi^2+1)/2)
                     [t xi-sinh(t xi)]+O(A^2).          (10)

The O(A^2) is for each fixed xi,t. At xi=1, sin is 1-Lipschitz; hence

    W2(T_A,G_A)>=A epsilon e^(-1)[sinh(t)-t]-C_t A^2.

For every odd j>=3, the first-order cumulant coefficient is

    kappa_j(T_A)= -A epsilon sin(j pi/2)e^(-1/2)t^j
                    +O(A^2).                          (11)

Thus a cubic-only repair still leaves rank five at first order for fixed t; a fixed finite moment match is not an amplitude-order theorem. The existing full packet handles this whole first-order non-Gaussian current at once. This is why the next test must examine its second-order nonlinear current rather than Gaussian accuracy alone.

## 5. Exact third-cumulant failure of the quadratic-calibrated local chord

For the finite quadrature define

    m2=sum w_i r_i^2, m4=sum w_i r_i^4,
    m11=sum w_i r_i q_i, m31=sum w_i r_i^3 q_i,
    K_Q=m2+2 beta_Q m11,
    J_Q=m4+2 beta_Q m31.

Use (2) with f=A f0, f0(x)=2qx+epsilon k(cos(kx)-1). Write h=K/A and h_u=D_u h. In this smooth fixture the actual VALUE graph has the expansion

    Y_chord=u-Ah+lambda_Q A^2 h_u h+O_Lp(A^3).            (12)

This is an analytical expansion of the executed chord, not an HVP instruction. The remainder is uniform over all positive quadrature versions because f0 has bounded first and second derivatives and total quadrature mass is one.

For tau=exp(-k^2/2), the target's exact cumulant expansion is

    kappa_3(X_A)=A epsilon k^3 tau
       +A^2 q epsilon tau(k^5-6k^3)+O(A^3).             (13)

One derivation expands log E exp(hX-AU_1(X)) under the standard Gaussian. The order-A term is minus the third h-derivative of E_(N(h,1)) U_1. At order A^2, only the even/odd cross term contributes. Gaussian integration by parts gives

    Cov_(N(h,1))(X^2,sin(kX)-kX)
      =tau[2hk cos(kh)-k^2 sin(kh)]-2kh.

Its third derivative at h=0 is tau(k^5-6k^3), proving (13).

For the candidate, direct expansion of the **centered** third moment gives

    [A] kappa_3(Y_chord)=3 epsilon k^3 tau m2,
    [A^2] kappa_3(Y_chord)
      =3q epsilon tau[lambda_Q(k^5 J_Q-3k^3 m2)
                                       -2k^3 K_Q].      (14)

No shared-root term is dropped. The useful exact identities are

    E[(u^2-1)h_u h]=q epsilon tau(k^5 J_Q-3k^3 m2),
    E[u h^2]=2q epsilon k[tau(1-k^2 K_Q)-1],
    E h=epsilon k(tau-1), E[u h]=q.

For example the second identity includes the cross term between q u+2q beta_Q v and the nonlinear part of h, all on their original common u,v record.

Choose the admitted dyadic rule to have delta_Q<=A^3. Then the first-order coefficient discrepancy in (14) is O(A delta_Q)=O(A^4). Its coefficients converge to

    beta=pi/4, lambda=3/2-pi^2/8,
    K=(1+pi/2)/3, J=1/5+pi/15.

At the parameters (8), subtracting (13) from (14) gives

    kappa_3(Y_chord)-kappa_3(X_A)
       =delta_3 A^2+o(A^2),
    delta_3=q epsilon e^(-1/2)
       [3lambda(J-1)-2(1+pi/2)+5]
       =-0.009301002921... !=0.                         (15)

The displayed O(A^3) expansion at a fixed fine rule and the o(A^2) statement along growing rules are deliberately distinguished.

### Grant exact mean and variance and the defect remains

Let mu_Y,sigma_Y and mu_X,sigma_X be the actual exact means and standard deviations. For this argument only, define the favorable comparison

    Y_hat=mu_X+(sigma_X/sigma_Y)(Y_chord-mu_Y).            (16)

This is not an admitted producer: its exact moments are used only to give the candidate every possible affine mean/covariance repair. The candidate already matches the target variance to first order, so sigma_X/sigma_Y=1+O(A^2). It follows that the third cumulant changes only by O(A^3). Thus (15) holds for Y_hat as well.

Both centered fourth moments are bounded uniformly for small A. For Y_chord this follows from (3), its exact zero, and Gaussian moments. For X_A it follows directly from .2 A x^2<=U_A(x)<=.3 A x^2 and its Gaussian normalization. The affine factors in (16) stay bounded. For any coupling of two centered scalar laws,

    |E X^3-E Y^3|
      <= ||X-Y||2 ||X^2+XY+Y^2||2
      <= C ||X-Y||2

when those fourth moments have a common fixed bound. Therefore

    liminf_(A->0) W2(Law(Y_hat),Law(X_A))/A^2>0.          (17)

Independent Gaussian convolution in (9) preserves equal means/variances and multiplies their third-cumulant difference by t^3. Its fourth moments remain bounded. At t=1/2 the corresponding legal reverse-bridge law also satisfies (17), with a positive constant reduced by the factor 1/8 in the cumulant witness. The Gaussian buffer does not erase this current.

## 6. Exact missing current and consumer boundary

For a fixed shape U=A U_1 and a smooth test phi, the second-order test current of the chord is

    C_chord(phi)=lambda_Q E[(D_u h)h dot grad phi(u)]
                    +(1/2)E[h tensor h : Hess phi(u)].   (18)

The true tilt coefficient is

    C_target(phi)=(1/2)E[phi(u)
                ((U_1-EU_1)^2-Var(U_1))].              (19)

These are analytical coefficient identities. No potential value, formal tensor, or derivative-valued producer is executed. A next-order positive repair must cancel their difference for the actual carrier-dependent test, while retaining the full common-root genealogy. The cubic witness (15) proves that quadratic calibration plus complete first/second-moment repair does not satisfy that identity.

Under the explicit inherited guards after (6), the lift does produce a near-gradient **global mean** source for the local chord. A law theorem for a completed N(E E_J,I) integrates its u,v source bank. The remaining output in (2) still reads the same u and K separately. Appending those observations after that comparison is not licensed.

Freezing u as a caller does in fact preserve a useful private-v small-curl certificate. Let Ku1,Kv1 be the two D-by-D first blocks of K at (u-K(u,v),v), and Kv0 the v block at (u,v). Each individual K block is symmetric, and

    D_v E=Kv1-Kv0-Ku1 Kv0,
    D_v E-(D_v E)*=-Ku1 Kv0+Kv0 Ku1,
    ||Curl_v E||<=alpha^2 beta_Q.                       (20)

Its v first can still be O(alpha). The fixed-u origin E(u,0) is an executable caller-only two-packet record, and its centered and anchored energies carry the actual profile C alpha^2(|u|+sqrt(D)). Thus a caller-retained mean service is a legitimate next route once its origin, captured first, share, and radius bills are paid. It does not preserve the same old K(u,v) or v used separately by (18); those private labels have been integrated. A fresh independent reserve can pay the mean program's Gaussian buffer, but the source-law replacement must still account for the lost joint current, beginning with Cov_v(K,E|u). No such current is set to zero here.

This identifies a specific next source/consumer gate: a positive carrier-retained repair for the difference of the rank-one feedback and same-root rank-two currents (18)-(19), including its regenerated third-Hermite coefficient and higher descendants. It is not a request for a strong conditional-mean oracle. A different complete current compiler might solve it; the present separator rules out only (2), even with ideal affine moment matching.

## 7. Conditional caller, zero, original-query, and replay ledger

For the original conditional problem use the existing captured finite mode x=x_M and anchored force

    f(y)=s[g(x+s y)-g(x)], alpha=A s^2,
    Q_chord=x+s Y_chord.

All x and g(x) dependencies remain in the actual caller graph. For a=rz, the finite conditional mode has ||D_a x||<=2. Thus ||D_a f||<=4As, ||D_a K||<=4As. Differentiating K1 only once gives its additional feedback factor alpha/2, so

    ||D_a E||<=C As,
    ||D_a Q_chord||<=C,
    ||D_(u,v)Q_chord-sP||<=C A s^3.                     (21)

These are absolute first bounds, not derivatives of the small VALUE energy. Multiplying by the exact bridge's B_rt and adding its known A_rt caller row gives the literal first of the proposed conditional bridge. In particular B_rt s alpha=(Delta/t)A s, and B_rt<=1 because Delta<=t s^2 for 0<=r<t<=1. Every r,t version is fixed before differentiation; no derivative in a clock parameter is supplied here.

For the auxiliary full lift E_J, its captured-caller row instead retains the actual L_Q multiplier: ||D_a E_J||<=C L_Q As, and its original physical caller row is <=C L_Q s sqrt(A). The projected bound (21) is not silently assigned to this larger auxiliary source. Its finite precision likewise retains L_Q, as below.

For an original physical caller, the imported normalized source caller is O(sqrt(A)); the conditional anchor subtraction contributes O(s sqrt(A)), with the already priced finite-mode derivative. This gives ||D_y E||<=C s sqrt(A), and the corresponding physical caller row is obtained by the displayed s and B_rt readouts. No small derivative of a law error is asserted.

The complete original g-VALUE count for Q_chord is

    M+1+2 n_alpha,                                     (22)

before any separately requested terminal force, downstream consumers, or numerical/known-row work. K's already recorded n_alpha values are reused in K1-K; K1 executes n_alpha genuinely changed sites. Full lift E_J uses the same sites and only known linear output arithmetic. A changed raw input inside an imported consumer is a complete new execution of both packets. If that consumer makes N occurrences, the safe bill is M+1+2n_alpha N, not M+1+2n_alpha+N. Only the identical captured caller-only mode/anchor graph can be shared across complete independent occurrences.

Each requested first/adjoint uses original HVPs at these recorded original VALUE sites, including the feedback through K. No HVP is a VALUE leaf, and no saved HVP is differentiated. A discarded primal incurs its full replay. Within one occurrence K, K1, their common roots, every old source ancestor, and both evaluations of the same anchor remain aliased exactly as displayed. Source-zero output is u (zero at u=v=0); at original all-zero caller and private roots the physical state is the captured anchored mode. Re-evaluation at an identical zero site uses the saved original anchor.

The bridge adds one independent D-root when v>0. With the conditional packet this is 3D total; at t=1 it is 2D. A completed mean/pair consumer would instead have its entire replay/root/clock tape, not these two raw roots alone. All original physical sites are b_M+sqrt(A)[x+s(r_i u+q_i v)] and b_M+sqrt(A)[x+s(r_i(u-K)+q_i v)], plus the caller-only mode/anchor sites. Their displacement from x is at most C s||(u,v)||, uniformly in the query count.

Finite numerical g errors in (22) propagate through the actual weight sums and Lipschitz feedback; the full lift additionally carries L_Q. Assign absolute leaf tolerances through this graph, retain the already supplied conditional-mode and physical finite-mode floors, and freeze the quadrature, coefficient lambda_Q, precision version, and target grade before differentiation. No floor is divided by a possibly vanishing actual chord energy. The imported finite clock condition is still needed to call log(1/(As^2)) a public-log bill.

## 8. Diagnostics and scope

The author checker executes the actual scalar VALUE chord, computes the true target by independent one-dimensional Gaussian integration, and tests the favorable exact-moment comparison only as a diagnostic. It also verifies the analytic shared-root coefficient identities, full non-diagonal matrix quadratic covariance calibration, inherited reverse-bridge third cumulants, and the actual full-gradient lift/curl on a genuine C2/non-C3 two-dimensional noncommuting fixture. It reports 439 passing assertions. Its 337-node rule has beta_Q=0.7853981633982807; the calculated defect coefficient is -0.009301002921166684.

No numerical output is a proof of an all-order result. Equations (10)-(19) are the analytical separation. The new constructive result is the complete full-law affine fallback (1) and a correctly priced positive chord with a real matrix-quadratic grade gain. The nonlinear next-order repair remains open at the explicit current (18)-(19); no improved general c(P) recurrence is asserted.
