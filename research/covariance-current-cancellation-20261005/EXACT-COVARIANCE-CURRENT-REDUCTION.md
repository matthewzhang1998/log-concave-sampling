# Exact covariance-current reduction and the genuine skew-square boundary

2026-10-05. This packet corrects an analytical ceiling in the existing positive buffered consumer. It does not claim an arbitrary-order original-VALUE construction, a native rank-five admission, or closure of the two-cycle whole-bank return.

## 1. Result

In the existing finite-rank path, the moving-Hermite/root current and the covariance sampler current cancel to first order in the reference correction. Their exact sum is the single rank-two current

    t J_t Sigma J_t^T.

There is no residual term linear in Sigma and the skew map. In particular the separately estimated L^2 h3 term is unnecessary. This is an exact same-endpoint Gaussian identity, not an asymptotic cancellation or a modification of the sampler.

Consequences for the existing RAW-cumulant analytical positive reference, under its normalized small-radius guard:

    m = 3: W2 <= C_m L^3 e / v^(3/2),
    m = 4: W2 <= C_m L^4 e / v^2,
    m >= 5: W2 <= C_m L^5 e / v^(5/2).

The old quartic statement remains valid: its fifth-current Stein remainder still has L^4 e/v^2. What changes is the assertion that a moving-covariance debt prevents higher ranks from helping. The actual first obstruction to raw-rank improvement is the skew-square feedback. A scalar bounded-Hessian gradient example below proves its L^5 e scale is real.

## 2. Setup and the weighted-Hermite heat identity

Keep exactly the imported bank and positive path:

    C_t = (v/2) I + (1-t^2) Sigma,
    x_t = S_t Z1, S_t = C_t^(1/2), eta = sqrt(v/2),
    q_r,t(x) = K_r : H_(r-1)^(C_t)(x),
    p_t(x) = sum_(r=3)^m a_r(t) q_r,t(x), a_r(t)=1-t^r,
    W_t = t B + x_t + p_t(x_t) + eta Z2,
    J_t = D_x p_t(x_t).

B, Z1 and Z2 are independent as in the imported theorem. All test derivatives in Sections 2–5 are at THIS actual W_t. Let

    D_t = sum_r a_r(t) partial_t[q_r,t(S_t Z1)],

where Z1 is fixed and derivatives of a_r are excluded, exactly as in the imported current generator.

For any constant coefficient vector-valued Hermite polynomial q_C, its Gaussian density weighted form is

    q_C(x) gamma_C(x) = K : (-D_x)^(r-1) gamma_C(x).

Since covariance differentiation of a Gaussian density is the heat operator and commutes with fixed spatial derivatives,

    partial_t(q_C gamma_C) = (1/2) dot C : D_x^2(q_C gamma_C).

Subtract q_C times the same heat equation for gamma_C. Using D log gamma_C = -C^(-1)x gives

    partial_t q_C = (1/2) dot C : D^2 q_C
                   - Dq_C dot C C^(-1)x.

In this path C and Sigma commute. Thus dot S S^(-1)=(1/2)dot C C^(-1)=-t Sigma C^(-1), and the derivative with Z1 fixed is

    dot q_C(S Z1) = -t Sigma : D^2 q_C(x)
                   + t Dq_C(x) Sigma C^(-1)x.                 (1)

This argument applies component by component and does not require Dq_C to be symmetric. Summing with a_r held fixed gives (1) with q replaced by p and its left side replaced by D_t.

## 3. Exact current cancellation

Condition on B and Z2. For indices i,j,k, Gaussian integration in x yields

    E[(C^(-1)x)_k (partial_j p_i) (partial_i phi)(W)]
      = E[(partial_k partial_j p_i) (partial_i phi)(W)]
        + E[(partial_j p_i) (partial_i partial_l phi)(W)
              (delta_lk + partial_k p_l)].

Multiply by t Sigma_jk and sum. The first term cancels the -t Sigma:D^2p term in (1), leaving

    E[D_t . grad phi(W_t)]
      = t E[(J_t Sigma (I+J_t)^T) : D^2 phi(W_t)].             (2)

The imported covariance sampler current is -t Sigma J_t^T. Because D^2 phi is symmetric, J_t Sigma and Sigma J_t^T have the same contraction. Therefore

    E[D_t . grad phi(W_t)]
      - t E[(Sigma J_t^T):D^2 phi(W_t)]
      = t E[(J_t Sigma J_t^T):D^2 phi(W_t)].                  (3)

Both sides are evaluated at the same positive law. No source-bank law is Gaussianized. No derivative is transferred to B. The untouched Z2 buffer remains independent of every coefficient. Polynomial truncation and the same Sobolev approximation used by the imported consumer justify the identity for its stated C1/Lipschitz source assumptions.

### Root-gauge invariant version

The same rule extends to any differentiable positive covariance C=S S^T with any differentiable invertible actual root S. Set Omega=dot S S^(-1), A=Omega C, so dot C=A+A^T. Keep D as the fixed-root derivative of the Hermite maps with their scalar amplitudes held fixed. The weighted heat equation and Gaussian integration give

    E[D . grad phi(W)]
      = -E[(J A^T (I+J)^T):D^2 phi(W)].

The actual moving-root sampler current beyond the baseline (dot C/2):D^2 phi is +A J^T:D^2 phi. Their linear difference A J^T-J A^T is antisymmetric. Hence their exact sum is

    -(1/2) E[(J dot C J^T):D^2 phi(W)].                  (3a)

No commutation of C and dot C, or symmetry of J, is needed. One must use the actual root current A J^T; replacing it by a guessed symmetric-root formula can spoil the calculation. For the imported path dot C=-2t Sigma, (3a) is exactly (3). This provides a repeatable exact covariance-current rule for moving positive carriers, rather than a new source oracle.

## 4. Corrected all-rank current generator

Use the set-partition tensors F_(r,pi) from the imported all-rank generator. Replace its first two instructions by:

    A1 = 0,
    A2 = t J_t Sigma J_t^T.

Then, unchanged:

    for r=3,...,m and every partition pi of [r-1]:
       add -r t^(r-1) F_(r,pi) to A_(1+|pi|),
    for r=3,...,m:
       add +r t^(r-1) K_r to A_r,
    A_(m+1) includes t^m T_m.

The covariance term ADDS to A2. For m=4 the exact ranks become

    A1 = 0,
    A2 = t J Sigma J^T -3t^2 K3[E] -4t^3 K4[F],
    A3 = -3t^2(K3[H,H]-K3) -12t^3 K4[E,H],
    A4 = -4t^3(K4[H,H,H]-K4),
    A5 = t^4 T4.

No nonlinear map-derivative offspring have been dropped. In particular the skew-square terms in A2 and A3 remain.

## 5. One-energy norm, width powers and numerical constants

At v=1 let h=sum h_r and s=sum k_r, using all proper-cut norms of the actual K_r. The new covariance term satisfies

    ||J Sigma J^T||_(L2;HS) <= C_m ||Sigma||op h s
                              <= C_m L^2 h s.               (4)

Each term is a two-vertex Wick graph. Its two coefficient vertices are joined by the deterministic Sigma edge; both retain their final physical output. All Gaussian contractions join distinct initial Wick factors. Start one vertex in HS norm and attach the other through the proper cut consisting of all their connecting slots. The second physical output remains open, so the usual HS-times-operator estimate is legitimate. No self-trace or extra dimension factor is introduced.

More explicitly, let lambda=lambda_min(C), and take one q_r Jacobian and one q_s Jacobian. A valid fixed-rank coefficient is

    beta_rs = (r-1)(s-1) [sum_(p=0)^min(r-2,s-2)
       (p! binom(r-2,p) binom(s-2,p))^2
        (r+s-4-2p)!]^(1/2).

Then their covariance term is bounded by

    beta_rs ||Sigma||op min(h_r k_s,h_s k_r)
            lambda^(-(r+s)/2).

After the rank-two reserve transfer its bill is at most

    2^((r+s+1)/2) beta_rs L^2
      min(h_r k_s,h_s k_r) / v^((r+s+1)/2).                 (5)

The leading (3,3) contribution is L^2 h3 k3/v^(7/2), hence O(L^7 e/v^(7/2)). It is TWO additional powers smaller than the surviving skew-square feedback h3 k3/v^(5/2).

A conservative explicit all-rank constant is twice the imported

    C_m = 2^(10m^2) ((2m^2)!)^4.

The resulting normalized bound is

    W2 <= 2 C_m [L^m e + L^2 h s
                     + sum_(r=3)^m h_r sum_(j=1)^(r-1) s^j]. (6)

Restore v by L'=L/sqrt(v), e'=e/sqrt(v),
 h'_r=h_r/v^(r/2), k'_r=k_r/v^(r/2), and multiply by sqrt(v).
With M_m=sum_(r=3)^m c_r/r!, the guard

    L' <= min(1/2, (2 M_m)^(-1/3))

ensures s'<=1/2. The imported estimates h_r<=L^(r-1)e/r and k_r<=c_r L^r/r! then imply the rates in Section 1. For a numerical rate constant one may use 8 C_m (1+M_m)^2, which deliberately dominates all finite rank sums under that guard.

## 6. Why adding raw cumulants really stops at the next grade

At t=0 the different Gaussian Hermite degrees are orthogonal. In particular E[x p_0(x)^T]=0, and

    Cov(W0) = v I + Sigma + E[p_0(x)p_0(x)^T].               (7)

The last term is positive semidefinite, with a separate nonnegative contribution from every included raw rank. In scalar dimension its exact contribution is

    sum_(r=3)^m (r-1)! K_r^2 / C_0^(r-1).                   (8)

For an explicit admissible original gradient, fix a=1/4 and set

    g_epsilon(x)=epsilon [x+a(cos x-1)],
    U_epsilon(x)=epsilon [x^2/2+a(sin x-x)].

Then g(0)=0 and 3epsilon/4 <= U'' <= 5epsilon/4 globally. Put B=g_epsilon(G)-Eg_epsilon(G). Its third cumulant equals epsilon^3 times a nonzero constant: with c=exp(-1/2),

    κ3(B)/epsilon^3 = -3ac+a^3 E[(cos G-c)^3] < 0,

because |E[(cos G-c)^3]|<=8 and 3ac>8a^3 for a=1/4. Both L and e scale linearly with epsilon.

For each fixed m>=3 and fixed positive v, (8) has leading term

    2 K3^2 / C_0^2 = Theta(epsilon^6/v^2).

Both laws are centered, so every coupling obeys

    W2(W0,W1) >= |sqrt(Var W0)-sqrt(Var W1)|
                = Theta(epsilon^6/v^(5/2)) as epsilon->0.  (9)

This is exactly the scale L^5e/v^(5/2). Higher raw ranks cannot cancel a nonnegative variance surplus. Thus the improved raw-map ceiling is a genuine obstruction, while the previously stated L^4e moving-covariance ceiling was only a loose upper-bound artifact.

## 7. Explicit first connected repair target

Here is the exact coefficient target for the leading skew-square repair, not an executable tensor port. Set Q=C^(-1), K=K3, and define

    N_ij = sum_abcd K_iab Q_ac Q_bd K_jcd,
    M_ijbd = sum_ac K_iab Q_ac K_jcd.

For x~N(0,C), an independent additive Gaussian reserve, U=x+reserve, and p_i=K_iab H_ab^C(x), the finite Hermite product formula gives

    (1/2) E[p_i p_j phi_ij(U)]
      = (1/2) (K tensor K):E D^6 phi(U)
        + 2 M:E D^4 phi(U) + N:E D^2 phi(U).                    (10)

Consequently the positive polynomial map can include

    s_i(x) = -N_ij H_j^C(x) -2 M_ijbd H_jbd^C(x)            (11)

to remove these unwanted connected rank-two and rank-four terms at that weak expansion order around the Gaussian base U, while retaining the desired disconnected rank-six term. This is coefficient extraction in a finite Taylor jet, not an assertion about real exponential moments. In scalar standard coordinates, p=aH2 and s=-a^2 H1-2a^2 H3.

The coefficients N and M have one-energy size, bounded by h3k3 ||Q||^2 and h3k3 ||Q|| respectively. Formula (11) is therefore of the same h3k3/v^(5/2) size as the removed skew-square terms. It is a genuine positive polynomial pushforward with the reserve untouched, but positivity alone does not make its coefficients original-VALUE executable. This weak-order calculation is not an exact same-nonlinear-endpoint current cancellation; higher Taylor terms and mixed cumulant maps still require complete analysis.

This identifies the next native targets more precisely. Substitution of the literal three-force K3 trees gives six-force coefficients. N joins two physical legs of each tree, leaving two physical marks and one graph cycle; M joins one leg, leaving four marks and a tree. Coefficient histories still retain all old structural clocks, private shields, observer dependencies and actual caller sites. Neither target can be admitted merely by naming an abstract cumulant tensor.

## 8. Native accounting and stopping boundary

The exact regrouping (3) changes NO executed program. Its incremental bills are:

- original-g VALUE queries: zero;
- owned public, private, auxiliary or reserve roots: zero;
- captured callers, outside observers and replay sites: unchanged;
- source-zero carrier, positivity gaps and actual first/curl: unchanged;
- finite means, source/filter/calibration and arithmetic floors: unchanged.

Thus it is safe to use this improved analytical comparison inside an already admitted native packet without inventing a derivative source or differentiating a law estimate. It does not create the native fifth-cumulant packet needed to exploit the m>=5 rate.

The new skew-square repair (11), in contrast, still needs a positive finite original-VALUE realization of its exact six-force rank-two/rank-four coefficients AND every conditional, auxiliary and complete-old-bank return, direct first/curl, root ownership and replay bill. No such interface is proved here. The supplied two-cycle construction also keeps its whole-old-bank return and doubled-four-cycle marked-chunk re-entry gate open. This packet does not consume or close either gate.

The task's full requested repeatable all-feedback native compiler therefore remains open. The finished result is a strict and independently checkable correction to the proposed obstruction, an improved all-rank positive comparison, a sharp next raw-map obstruction, and the explicit first connected correction target.
