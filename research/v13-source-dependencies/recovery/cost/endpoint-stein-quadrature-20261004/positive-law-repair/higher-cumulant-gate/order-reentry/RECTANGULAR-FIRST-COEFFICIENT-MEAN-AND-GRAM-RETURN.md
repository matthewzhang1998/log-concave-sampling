# A rectangular first-coefficient mean and positive true-Gram return

2026-10-04. Finite original-gradient VALUE construction, independent review requested. It fills the true-orientation action for the finite j coefficient kernel, subject to the separately proved covariance-mixture lemma. It does not fill the mixed K current or the same-carrier P3 mean.

## 1. Input source and target

Use the finite j source from `../same-carrier-feedback/RESOLVENT-COVARIANCE-STABILITY-AND-FINITE-J-KERNEL.md`. All its versions, nodes and anchors are frozen. At captured X, with c_i=sqrt(1-r_i^2), it is

    H_Q(x,H)=sum_l v_l g(tau_l x+sqrt(1-tau_l^2)H),
    J_Q(x,H)=g(x-H_Q(x,H)),
    V_X(u,H)=sum_i (w_i r_i/c_i)
                   [J_Q(r_iX+c_i u,H)-J_Q(r_iX,H)].

The target is the SAME finite nonsymmetric matrix

    B_Q(X)=E_(u,H) D_u V_X(u,H).

No Jacobian is executed as a producer. The source has private first O(A), complete square-lift curl O(A^2) under P_u=(I,0), captured-X first O(A), exact V_X(0,H)=0 with coherent same-site reuse, and energy O(A sqrt(D)). Its raw complete bill is 2N_r(N_in+1) original g VALUES. The anchor J_Q(r_iX,H) depends on the PRIVATE H and cannot be globally captured across banks.

There is a literal split

    V_X=V_g,X+V_e,X,
    V_g,X(u)=sum_i(w_i r_i/c_i)[g(r_iX+c_i u)-g(r_iX)].

V_g,X is a genuine gradient in u. Put E_Q(x,H)=J_Q(x,H)-g(x). Then V_e is the identical anchored positive-weight difference of E_Q. Its useful stronger ports are

    full first O(A), square-lift curl O(A^2),
    H first O(A^2), X first O(A),
    |V_e,X(u,H)|<=C A^2(|X|+|u|+|H|),
    V_e,X(0,H)=0.

The small VALUE estimate follows from |E_Q(x,H)|<=A|H_Q(x,H)| and the bounded sums sum w r/c and sum w r^2/c. It does not make its u or X derivative order A^2.

## 2. The partial-variable first-coefficient filter

Choose the admitted K-node first-chaos rule b_l in (0,1/2), c_l=sqrt(1-b_l^2), with

    sum d_l=1, sum d_l b_l^(2h)=0 for 1<=h<K,
    sum |d_l|<=4, min b_l>=1/(2K).

At captured X and a fresh standard D-dimensional incoming p, define the actual source on two private D-roots z,H by

    F_X(p;z,H)=sum_l d_l
       [V_X(c_l z+b_l p,H)-V_X(c_l z-b_l p,H)]/(2b_l).     (1)

The SAME H is left UNCHANGED in every sign and every filter node. The SAME z is used at every filter node. This is a partial-u filter; it is not the full square-space Mehler response with public input (p,0), which would wrongly shrink H's Gaussian variance.

Condition on H first. The ordinary D-input vector Hermite calculation for u gives

    ||E_z F_X(p;z,H)-[E_u D_u V_X(u,H)]p||_(L2(p))
        <=C 4^(-K) ||V_X(.,H)-E_u V_X(.,H)||_(L2(u)).

The conditional energy on the right is at most C A sqrt(D) by the uniform u first. Jensen in H therefore proves, uniformly in captured X,

    ||E_(z,H)F_X(p;z,H)-B_Q(X)p||_(L2(p))
        <=C A sqrt(D)4^(-K).                           (2)

This is a conditional-mean calibration, not a strong random-matrix estimate. It needs no symmetry of B_Q and keeps every H_Q ancestor on its original record. Taking 4^(-K)<=A^3 divided by the chosen public-log allowance puts (2) below A^4 sqrt(D); K is public-logarithmic. Any higher absolute filter grade is equally available at a fixed requested order.

## 3. Genuine-gradient and near-gradient mean sources

Split (1) on the SAME z,H record as F_g+F_e, using V_g and V_e respectively. For fixed X,p, F_g is a genuine gradient in z. Indeed

    D_z F_g=sum_l d_l c_l/(2b_l)
       [D V_g(c_l z+b_l p)-D V_g(c_l z-b_l p)]

is symmetric. The weights need not be positive for this gradient fact. Its actual private first is at most C K A, and its captured origin F_g(X,p;0) is executable.

For F_e use the square lift P_z*F_e on (z,H). Its z-z skew comes from the difference of the two V_e skews and is at most C K A^2. Its H block is

    D_H F_e=sum_l d_l/(2b_l)
       [D_H V_e(c_l z+b_l p,H)-D_H V_e(c_l z-b_l p,H)],

also at most C K A^2. Consequently

    first(F_e)<=C K A,
    curl(P_z*F_e)<=C K A^2.                            (3)

The VALUE estimate of V_e, the finite 1/b bound, and Gaussian moments give at fixed X,p

    ||F_e||_Lp(z,H)<=Lambda A^2(|X|+|p|+sqrt(D)).        (4)

The origin F_e(X,p;0,0) is NOT automatically zero when p is retained nonzero. Compute and save it with the literal same finite program; subtract it before the near-gradient compiler and add it back afterward. It has the same deterministic profile C Lambda A^2(|X|+|p|), so the anchored source retains (4). This is a captured-X,p origin only; private H-dependent raw V anchors remain inside each complete F_e occurrence.

The actual X first of F_g,F_e is O(Lambda A), and the p first is O(Lambda A), by the finite chain rule and the corresponding raw V ports. No derivative of (4) is used. At p=0 the paired response is exactly zero by identical-query reuse. At the whole original source zero both sources vanish and only the completed mean's known Gaussian carrier remains.

Use independent COMPLETE gradient and near-gradient mean banks after X,p and every version/origin are captured. Fixed positive variance shares 1/2+1/2, gradient mean order at least four, and padding mu=A give an actual positive mean M_X(p) with

    W2(Law(M_X(p)|X,p),N(E_(z,H)F_X(p;z,H),I))
        <=Lambda A^4(|X|+|p|+sqrt(D))+absolute floors.    (5)

Every imported radius, finite-clock, active-dimension and covariance-gap guard is imposed at the ACTUAL normalized first C K A and curl C K A^2, not at A with K silently erased. The known carrier and origins are those of the finite mean programs. The claimed unit covariance belongs to the reference law, not the actual sample covariance.

Combining (2) and (5), then integrating fresh p,

    ||W2(Law(M_X(p)|X,p),N(B_Q(X)p,I))||_(L2(p))
        <=Lambda A^4(|X|+sqrt(D))+floors.               (6)

The inner p-conditional comparison retains p and X but integrates every z,H/filter/mean root. No old private root can be reattached afterward. Requested first/adjoint sweeps remain original HVPs at all recorded original VALUE sites, never producer leaves.

## 4. True orientation from a Gaussian mean law

Condition on X, then integrate standard p in the Gaussian reference of (6). Since its independent unit Gaussian is drawn AFTER p,

    B_Q(X)p+N ~ N(0,I+B_Q(X)B_Q(X)*).                  (7)

This is exactly the transposed Gram, even when B_Q is nonsymmetric. No forward square B_Q^2 and no generic adjoint oracle has been substituted. The source M_X(p) is used only through this whole conditional LAW comparison. The Gaussian reference noise is not identified with one of its actual source-zero roots.

Its actual returned private/caller first is the imported O(Lambda A) residual bound together with its known Gaussian carrier and the recorded X,p origins. This follows from the finite mean graph applied to (3)-(4), not from (6).

## 5. Positive clock and variance-budget join

Let the desired covariance target be

    C_Q(Z)=sum_j a_j E_[X_j|Z][B_Q(X_j)B_Q(X_j)*],
    X_j=q_jZ+sqrt(1-q_j^2)G_j,
    a_j=2v_j q_j>0, sum_j a_j=1.

Choose a fixed total baseline variance v0>0, put v_mean=v0/2 and v_free=v0/2, and allocate nu_j=v_mean a_j. At each node, execute the mean service for the scaled source F_X/sqrt(v_mean), on independent COMPLETE banks, with its own fresh standard p_j. Return

    R(Z)=sum_j sqrt(nu_j) M_(X_j,v_mean)(p_j)
                                      +sqrt(v_free) Z_free.       (8)

Here M_(X,v_mean) targets N(B_Q(X)p/sqrt(v_mean),I). Every scaled source radius/curl and floor is checked at the fixed 1/sqrt(v_mean) multiplier. No inverse a_j source radius appears because a_j/nu_j=1/v_mean. The final Z_free is a literal untouched Gaussian independent of every owned G_j and mean/action bank.

After the conditional mean comparisons, (8) is a centered Gaussian mixture with covariance

    v0 I+sum_j a_j B_Q(X_j)B_Q(X_j)*.                  (9)

The component errors sum with their actual sqrt(nu_j) factors. The conservative loss sum_j sqrt(a_j)<=sqrt(N_q) is a public-log factor and is explicitly included in Lambda. X_j is NOT captured globally; it is an owned coarse root inside this covariance service, and its path into every F_X source is counted.

A mixture does not become Gaussian just because (9) has the right expectation. The following separately proved covariance-mixture current is required:

    W2(mixture,N(0,v0 I+C_Q(Z)))
       <=C_v0 Lip_HS(C_random) ||C_random-E C_random||_(L2;HS),

for the actual owned Gaussian bank, with a fixed independent buffer. For the present source its hypotheses hold uniformly in Z. Indeed j_Q is globally L_j-Lipschitz, and gradient symmetry is not needed for

    D_u[P_r D j_Q](X)
       =(r/sqrt(1-r^2)) E[(D j_Q(rX+cG)u)G*].

Commuting the two spatial derivative slots and Gaussian Bessel gives the HS bound r L_j/sqrt(1-r^2). Positive clock summation therefore proves

    ||D_X B_Q(X)[u]||HS<=C A|u|,
    ||D_X(B_QB_Q*)[u]||HS<=C A^2|u|.                  (10)

For C_random=sum a_j B_Q(X_j)B_Q(X_j)*, Cauchy–Schwarz over the independent root blocks yields Lip_HS(C_random)<=C A^2. Its centered HS energy is at most C A^2 sqrt(D) by ||B_Q||op<=C A and ||B_Q||HS<=C A sqrt(D). Hence that current costs O(A^4 sqrt(D)), not an unproved covariance-matching shortcut.

Together with (6), the completed positive reserve (8) targets

    N(0,v0 I+C_Q(Z))

with conditional error Lambda A^4(|Z|+sqrt(D)) plus floors, and integrated fresh-Z error Lambda A^4 sqrt(D). The original resolvent-stability/clock theorem then restores C_Q to Cov(H_j|Z) at the same integrated grade. This last restoration retains its exact L2(standard Z) scope.

## 6. Complete original-query and private-root counts

Let Q_g=2N_r for the raw V_g source, and Q_V=2N_r(N_in+1) for V. Then a conservative complete raw occurrence count for the filter sources is

    Q_Fg<=2K Q_g,
    Q_Fe<=2K(Q_V+Q_g).                                 (11)

Same-key cancellations/anchors inside one response can reduce these safe counts. They do not allow sharing a private H or an ancestor across distinct complete mean banks.

Let N_B,N_E be the literal complete occurrence counts of the fixed-order gradient and near-gradient mean programs at their actual normalized sources. One conditional mean call costs

    Q_M<=Q_captured(X,p,origins)+N_B Q_Fg+N_E Q_Fe
                                 +known/numerical/replay work.    (12)

The origin captures in (12) are full original-VALUE response executions at z=H=0 with the SAME X,p versions. Each q bank has a new owned X_j and its own captures. The reserve costs the sum of (12) over j, plus the separate original/conditional mode and original-g anchor records. The normalized p-source scale 1/sqrt(v_mean) does not vary with a_j.

All original ancestors of every H_Q and J_Q in (11)-(12) are included. Every changed raw source argument is a complete replay. Actual Gaussian dimensions include every z,H,mean/filter/clock,p_j,G_j and fill block, not only the two schematic roots of raw V. At fixed accuracy order all counts and dimension/D are public-log polynomials, so this service adds no inverse-A query exponent. Known arithmetic/storage and all absolute inverse-width precision multipliers remain separately charged.

## Scope

This is a source-qualified positive true-Gram covariance service for the finite j kernel, contingent on the stated imported mean/filter guards and the separate covariance-mixture proof. It does not compute B_Q as a strong matrix estimate, does not make V_X a gradient, and does not fill the mixed K current or the P3 same-carrier mean. The all-order typed closure remains a separate theorem.
