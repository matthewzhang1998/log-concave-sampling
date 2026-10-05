# Radial inflation also refutes the actual centered P3 covariance correction

2026-10-05. Separate SAME-g corollary; the frozen single-history proof is unchanged. Independent exact-text review requested. This note refutes a precise covariance-derivative repair of the literal F3_Q mean. It does not refute the full resummed Gaussian-backbone VALUE consumer or the exact live shifted current.

## Result and exact candidates

Let m3=R1 psi_H, where

    psi_H(x)=E[g(x-F2_H)|X0=x]

uses the TRUE second substitution on the continuous OU history. Let I2_Q be the actual two-level common-root source inside the literal finite F3_Q, and put

    psi_Q(x)=E[g(x-I2_Q)|x],
    E[F3_Q|Z]=Q_out psi_Q(Z).

Define the actual second-substitution covariances

    C2_H(x)=Cov(F2_H|x),   C2_Q(x)=Cov(I2_Q|x),

and the actual first-displacement covariances

    C1_H(x)=Cov(F1_H|x),
    C1_Q(x)=Cov(sum_j v_j g(tau_j x+d_j H)|x).

There is a smooth anchored convex-gradient family with 0<=Dg<=A I and A=D^(-1/4) for which BOTH candidates

    m3 ?= Q_out psi_Q
          +(1/2) R1[(C2_H-C2_Q):D2g(x-mu_*(x))]
          +O(Lambda A^4 sqrt(D)),                        (C2)

    m3 ?= Q_out psi_Q
          +(1/2) R1[(C1_H-C1_Q):D2g(x-mu_*(x))]
          +O(Lambda A^4 sqrt(D))                         (C1)

have discrepancy at least c_* A for all sufficiently large D. The allowance is A^4 sqrt(D)=A^2. Here mu_* may be E[F2_H|x], E[I2_Q|x], or any pointwise convex combination of those actual means. The conclusion also covers the theta integral of such centered derivative currents.

The middle rule has multiplier tolerance delta_mid<=A^2, and the outer rule delta_out<=A^3. All three rules are positive, have mass one, and have exact first moment 1/2. The inner rule can satisfy any admitted tighter accuracy; this proof only needs its stated positivity/moment contracts. Thus the corollary applies to the literal finite raw source in `P3-MEAN-RAW-SOURCE-AND-SECOND-DECOUPLING-GATE.md`, not merely an idealized continuous cheap comparison.

## 1. Original field and the second-substitution coherent radius

Use the normalized smooth step psi from the frozen single-history radial-inflation proof. Its derivative is an even nonzero bump, strictly decreasing on (0,2), with sup psi'<1. Let

    R=sqrt(D)=A^(-2),   c=1/2,   d=1/4,
    lambda=cA,
    a2=lambda/2-lambda^2/4,
    k2=1-a2,
    R0=k2 R,
    f(x)=psi(|x|-R0) x/|x|,
    g(x)=lambda x+dA f(x),                              (1)

with f(0)=0. For D>=256 the shell is far from zero. Exactly the same Hessian calculation as in the single-history proof gives

    Dg=lambda I+dA[psi' nn*+(psi/r)(I-nn*)],
    g(0)=0,   0<=Dg<=A I,
    |g(x)-lambda x|<=epsilon:=dA.                       (2)

Every g at every level below is this same original function. The only change from the single-history family is the coherent shell radius k2 R appropriate to the actual second substitution.

## 2. Exact true OU backbone, with bounded full nonlinear remainder

On the actual OU history conditioned on X0=x define the analytical comparison variables

    H1=integral_0^infinity e^(-t)X_t dt,
    H2=integral_0^infinity t e^(-t)X_t dt.

At time t, write the genuine future first displacement as

    F1,t=lambda H1,t+R1,t,   |R1,t|<=epsilon.

The second substitution therefore has the EXACT identity

    F2_H=lambda H1-lambda^2 H2+R_H,
    R_H=-lambda integral e^(-t)R1,t dt
                   +integral e^(-t)(g-lambda Id)(X_t-F1,t)dt,
    |R_H|<=epsilon(1+lambda).                           (3)

The H2 kernel follows by integrating over t,u>=0 with t+u=s; its convolution weight is s e^(-s). This uses the complete genuine future ancestry, not independent copies.

The conditional Gaussian moments are

    E H1=x/2, E H2=x/4,
    Cov(H1)=I/4,
    Cov(H1,H2)=I/4,
    Cov(H2)=5I/16.

For example, integrating the stationary covariance gives

    K(a,b)=integral integral e^(-at-bu)e^(-|t-u|)dtdu
          =(1/(a+b))(1/(a+1)+1/(b+1));

differentiate at a=b=1 and subtract the endpoint regression products. Hence the actual Gaussian part of (3) is

    lambda H1-lambda^2 H2=a2 x+b_H N_H,
    b_H^2=lambda^2/4-lambda^3/2+5lambda^4/16,             (4)

where N_H is standard and independent of x, but generally correlated with R_H. No independence from the remainder is asserted.

## 3. Exact literal finite two-level backbone

For the finite source use its original middle and inner rules:

    y_j=tau_j x+d_j H,
    z_jk=sigma_k y_j+e_k J,
    L_j=sum_k u_k g(z_jk),
    I2_Q=sum_j v_j g(y_j-L_j),
    beta_mid=sum_j v_j d_j,
    beta_in=sum_k u_k e_k.

Here H,J are independent standard vectors, shared across their full levels. Positivity, first-moment exactness, and (2) give the EXACT decomposition

    L_j=lambda(y_j/2+beta_in J)+r_j,   |r_j|<=epsilon,

    I2_Q=lambda(1-lambda/2)(x/2+beta_mid H)
                          -lambda^2 beta_in J+R_Q,
    |R_Q|<=epsilon(1+lambda).                           (5)

Thus

    I2_Q=a2 x+b_Q N_Q+R_Q,
    b_Q^2=lambda^2(1-lambda/2)^2 beta_mid^2
                                   +lambda^4 beta_in^2, (6)

with N_Q the normalized displayed Gaussian combination. It is standard and independent of x, but it is not independent of R_Q.

Uniformly over these rules,

    b_H^2=c^2 A^2/4+O(A^3),
    b_Q^2=c^2 A^2 beta_mid^2+O(A^3).                  (7)

The middle multiplier accuracy supplies the same nonzero gap as before:

    1>=beta_mid>=2/3-delta_mid>=7/12,
    beta_mid^2-1/4>=13/144.                             (8)

No nested-mean theorem, moment-to-law inference, or path approximation is used in (3)-(8).

## 4. The actual P3 force means keep the full radial inflation

From (3)-(6), original-g Lipschitz continuity gives pointwise in x

    |psi_a(x)-E_N g(k2 x-b_a N)|
                      <=A epsilon(1+lambda)=O(A^2),    (9)

for a=H,Q. This comparison couples to the actual N_a and keeps R_a correlated with it; only the deterministic bound on R_a is used.

Let X~gamma_D, S_D=|X|-R. The radial estimate of the single-history proof applies with k2=1+O(A) and R0=k2R. More explicitly, for w=-b_aN, its radial component is O_Lp(A), its Euclidean norm is O_Lp(A R), its quadratic radial inflation is b_a^2 R/2+O_Lp(A), and the cubic geometric remainder is O_Lp(A^3R)=O_Lp(A). Equations (7) give

    E_N f(k2 X-b_H N)=psi(S_D+h_H)n(X)+O_L2(A),
    E_N f(k2 X-b_Q N)=psi(S_D+h_Q)n(X)+O_L2(A),
    h_H=c^2/8,   h_Q=c^2 beta_mid^2/2.                  (10)

The outer linear terms lambda k2 x of the Gaussian responses cancel exactly. Therefore the ACTUAL true-versus-finite P3 integrands obey

    psi_H-psi_Q
       =dA[psi(S_D+h_H)-psi(S_D+h_Q)]n(X)
                                      +O_L2(A^2).       (11)

Their actual conditional means before the outer g are

    mu_H=a2 x+E R_H,   mu_Q=a2 x+E R_Q.

Each coherent center in their pointwise convex hull has the form

    mu_*=a2 x+M_*,   |M_*|<=epsilon(1+lambda).            (12)

In particular the true and finite second-substitution means may differ; this causes only an O(A^2) outer-force term and does not create the Omega(A) separator.

## 5. Both actual covariance choices keep only the tangent

A useful elementary fact is that, for a standard Gaussian vector N and any vector R on the same probability space with |R|<=B,

    ||Cov(N,R)||HS<=B,
    ||Cov(R)||HS<=B^2.

The first inequality is the genuine conditional first-chaos Bessel bound applied separately to the components of R; it does not require N and R independent. This is not a generic vector-trace estimate.

Apply this to the actual decompositions (3)-(6), where b_a=O(A) and B=epsilon(1+lambda)=O(A):

    C2_a=b_a^2 I+B2_a,   ||B2_a||HS<=C A^2.            (13)

Similarly the literal first-displacement decompositions give

    C1_H=c^2 A^2 I/4+B1_H,
    C1_Q=c^2 A^2 beta_mid^2 I+B1_Q,
    ||B1_a||HS<=C A^2.                                 (14)

For z_*=x-mu_*(x), formula (12) places its shell coordinate at S_D+O_L2(A). The explicit radial derivative calculation gives

    sup_z ||D2g(z)||_(HS->vector)<=C A,
    Delta g(z_*)=d A R psi'(S_D)n(X)+O_L2(1).            (15)

Thus the B remainders in (13)-(14) cost O(A^3) after contraction. The difference between their isotropic coefficients is only O(A^3) by (7); its explicitly calculated trace costs O(A^3 A R)=O(A^2). It follows that, for ell=1 and ell=2,

    (1/2)(Cell_H-Cell_Q):D2g(x-mu_*)
       =dA(h_H-h_Q)psi'(S_D)n(X)+O_L2(A^2).             (16)

This is uniform over all centers (12), so it also holds for their theta-averaged current. In particular,

    || R1[((C2_H-C2_Q)-(C1_H-C1_Q))
                                 :D2g(x-mu_*)] ||2
        <=C A^2.                                       (17)

Replacing first-displacement covariance by the exact second-substitution covariance cannot fix the leading missing nonlinear radial response.

## 6. Same nonzero linear witness; explicit finite outer floors

For beta=beta_mid put

    q_beta(s)=psi(s+h_H)-psi(s+h_beta)
                            -(h_H-h_beta)psi'(s),
    h_beta=c^2 beta^2/2.

Equations (11) and (16) give the discrepancy integrand

    dA q_beta(S_D)n(X)+O_L2(A^2).                        (18)

Let S~N(0,1/2) and m(t)=E psi'(S+t). Since the bump is symmetric and strictly unimodal, m(t) is strictly decreasing for t>0. Therefore

    E q_beta(S)
      =integral_(h_H)^(h_beta)[m(0)-m(t)]dt
      >=kappa:=integral_(h_H)^(c^2(7/12)^2/2)
                                  [m(0)-m(t)]dt >0.     (19)

This is exactly the uniform positive witness from the frozen single-history proof, numerically kappa approximately 1.8145100151e-6 for c=1/2. Its positivity is analytical; numerical integration is only a diagnostic.

The Gaussian shell CLT is uniform over the compact beta range, and T_D=X/R has norm one with R1 T_D=T_D/2. Consequently (18)-(19) imply, for ell=1 or 2,

    ||m3-R1 psi_Q
       -(1/2)R1[(Cell_H-Cell_Q):D2g(x-mu_*)]||2
        >=d kappa A/4-C A^2.                           (20)

The actual literal three-root source has OWN mean exactly Q_out psi_Q: shared roots across outer nodes do not change this expectation identity. From (5)-(6), or simply the Lipschitz/energy bound,

    ||psi_Q||2<=C A sqrt(D)=C/A.

Hence the complete finite outer mean floor is

    ||Q_out psi_Q-R1 psi_Q||2
       <=C delta_out A sqrt(D)=C delta_out/A.           (21)

Choosing delta_out<=A^3 makes (21) at most C A^2. Combining it with (20) proves, for ell=1,2 and all sufficiently large D,

    ||m3-E[F3_Q|Z]
       -(1/2)R1[(Cell_H-Cell_Q):D2g(x-mu_*)]||2
        >=c_* A,
    c_*=d kappa/8>0.                                    (22)

This is a direct literal-F3_Q target comparison, including its actual finite outer clock floor.

For completeness, if one also replaces the counterterm's OUTER R1 by any operator Q_ctr with ||Q_ctr-R1||<=delta_ctr, (16) bounds its additional floor by C delta_ctr A. Thus delta_ctr<=A suffices for O(A^2), and using the same delta_out<=A^3 is more than sufficient. This grants an ideal exact current only to strengthen the negative result; no D2g/covariance producer is licensed or constructed.

The admitted raw OWN-mean compiler has error Lambda A^4 sqrt(D)=Lambda A^2, plus its restored absolute floors. Taking those floors polynomially small and choosing the stated polynomial clock accuracies cannot hide the Omega(A) target error in (22). The usual public-log qualification fixes all additional public versions/floors at bounded or polynomial scales along this family.

Any exact coherent-mean correction bounded by A||mu_H-mu_Q||2 also costs at most 2A epsilon(1+lambda)=O(A^2), so retaining such a term does not rescue (C1) or (C2).

## Scope

This refutes the centered counterpart of the proposed actual P3 covariance repair, both with its original leading first-I covariance and with the actual full F2/I2 covariance. It uses the literal original g in the true history, in every finite ancestor, and in the outer force. It preserves complete source correlations and computes the contracted radial trace explicitly.

It does not refute the finite F3_Q raw source as a source with its own well-defined mean. It does not refute the exact centered/resummed Stein identities or the separate bounded-remainder Gaussian-backbone VALUE consumer: (9) already shows why retaining the full Gaussian response succeeds on this family at O(A^2).

The unrestricted theorem allows A sqrt(D)=A^(-1); no claim is made about a separately restricted regime A sqrt(D)=O(1). No ordinary Markov path-grid algorithm or generic Gaussianization theorem is introduced.

Status: separate analytical P3 centered-covariance corollary complete; independent exact-text review pending.
