# Endpoint Stein packet and a polylogarithmic law-quadrature gate

2026-10-04. Constructive bounded result for the known-center lane. This is NOT an admitted all-order posterior source, a replacement of the direct-mean baseline, or a sublinear complete-cost theorem.

## Result first

There is a concrete finite original-gradient VALUE packet whose Gaussian-conditional mean approximates the first-order Gaussian Stein transport field to arbitrary heat order using only O_J(log^2(1/A)) original queries. Its positive deterministic time quadrature has no inverse-heat exponent, even under the original C2 Hessian sandwich. Endpoint singularities are resolved by dyadic panels and a literal small-panel allowance, not by assuming high derivatives of the original force.

The packet has actual bounded private/caller firsts, an exact zero after finite-mode anchoring, complete finite query counts, and no algebraic dimension/heat guard. It also exactly removes first-order quadratic random-scale fluctuations. Its remaining O(A^2) covariance and nonlinear law currents are real. An all-order positive endpoint compiler consuming those currents would be a genuinely different route past the quarter-grid bottleneck; that compiler is not proved here.

The construction deliberately uses only a Gaussian-conditional mean in analysis. The mean is never queried or estimated strongly. Approximating that mean by Monte Carlo to arbitrary accuracy would return to the independent strong-mean cost barrier.

## 1. Target, finite mode, and original queries

Fix 0<A<=1/2; the baseline may retain its stronger A<=1/8 guard. For the exact posterior mode b, define

    U(z)=V(b+sqrt(A)z)-V(b)-sqrt(A) grad V(b) dot z,
    g(z)=grad U(z)=sqrt(A)[grad V(b+sqrt(A)z)-grad V(b)].

Then g(0)=0 and 0<=Dg<=A I. The exact standardized posterior is

    mu_U(dz) proportional to exp(-|z|^2/2-U(z)) dz.

No value of V or U is queried. Execute the baseline's fixed-count mode iteration b_0=y, b_(m+1)=y-A grad V(b_m), and replace b by b_M everywhere in the displayed g formula. Record one anchor gradient grad V(b_M). The actual standardized posterior at b_M additionally has the linear term

    ell dot z, ell=[b_M-y+A grad V(b_M)]/sqrt(A).

It cannot be dropped as an exact identity. Its physical transport price is at most sqrt(A)|ell|, by strong monotonicity, and is assigned the original finite-mode/caller numerical budget. In particular |ell|<=C A^(M+1/2)|grad V(y)|. All M and quadrature versions are frozen before caller differentiation.

The nonlinear finite packet has g_M(0)=0 exactly, regardless of ell. Freeze one deterministic original-gradient numerical version and reuse the identical recorded anchor value whenever the site equals b_M, including the all-zero graph. This makes the numerical zero cancellation literal. If a primitive instead returns separately perturbed evaluations at the same site, their difference is only a separately paid numerical zero residual and is not called exactly zero. The physical output zero is b_M, with its already paid distance from b.

## 2. The analytical Gaussian Stein field

Let gamma be standard Gaussian on R^D and define the Mehler operator

    P_r f(z)=E_G f(r z+sqrt(1-r^2)G), 0<=r<=1.

The field

    v(z)=integral_0^1 P_r g(z) dr

satisfies, for smooth test phi,

    E_gamma[v dot grad phi]=Cov_gamma(U,phi).                 (S)

Proof: with L=Delta-z dot grad, the centered inverse OU potential is integral_0^infinity P_(exp(-t))(U-EU)dt. Its gradient is integral_0^infinity exp(-t)P_(exp(-t))g dt=v. Gaussian integration by parts gives (S). All these are analytical objects, not oracle calls.

The actual regularity available here is enough:

    ||Dv||_op <= A integral_0^1 r dr=A/2,
    ||v||_(L2(gamma)) <= ||g||_2 <= A sqrt(D).

For any fixed p the same moment bound holds with C_p. This does not claim pointwise analyticity at r=1.

## 3. Explicit positive dyadic Gaussian quadrature

Write t=1-r. Partition [2^-K,1] into dyadic intervals [2^(-k-1),2^-k], k=0,...,K-1. Use m-point Gauss-Legendre on every interval. Use one midpoint on [0,2^-K]. Convert each t node to r=1-t, retaining its positive integration weight. Thus

    n=Km+1, w_i>0, sum_i w_i=1, sum_i w_i r_i=1/2.            (Q0)

The constant and linear identities are exact analytical quadrature identities. Actual encodings receive numerical allowances.

An explicit sufficient bound is 8*4^(-m)+2^(1-K); in particular there are universal c,C>0 such that

    sup_(n>=0) |sum_i w_i r_i^n - 1/(n+1)|
       <= C exp(-c m)+2^(1-K).                              (Q1)

Here n in (Q1) is a Hermite degree, not the number of nodes.

Proof of uniformity, including r near 1: on one interval t in [a,2a], choose the Bernstein ellipse of parameter rho=2. Its center is 3a/2, its real semiaxis is 5a/8, and its imaginary semiaxis is 3a/8. The squared distance to 1 is a convex quadratic in the ellipse cosine parameter, so its maximum is at a real endpoint. Uniformly for 0<a<=1/2 this gives |1-t|<=1-7a/8<=1. Hence the holomorphic polynomial (1-t)^n is bounded by one on that ellipse for EVERY n. Chebyshev truncation through degree 2m-1 has uniform error at most 4*4^(-m). Exactness through degree 2m-1, positivity, and total mass a bound each quadrature error by 8 a*4^(-m). Summing all nonterminal intervals gives at most 8*4^(-m). The final tiny interval costs at most twice its length because both its integral and positive midpoint rule are bounded by its length. No derivative of g appears in this argument.

Set K,m=O(log(1/delta)). Then n=O(log^2(1/delta)) and (Q1)<=delta. For delta=A^(J+2) times the prescribed numerical/log allowance and fixed maximum J, the work is O_J(log^2(1/A)). Rapid J constants are allowed; there is no inverse-A count.

## 4. Hermite diagonalization gives a dimension-safe mean certificate

For a vector field g in L2(gamma), write g=sum_(k>=0)g_k in orthogonal Gaussian chaoses. Mehler acts as P_r g_k=r^k g_k. Thus, for

    v_Q(z)=sum_i w_i P_(r_i)g(z),

(Q1) gives exactly

    ||v_Q-v||_2^2
      =sum_(k>=0)|sum_i w_i r_i^k-1/(k+1)|^2 ||g_k||_2^2
      <=delta^2 ||g||_2^2
      <=delta^2 A^2 D.                                     (Q2)

This is an L2(gamma) certificate for the conditional mean. It is not a pathwise error bound for a fixed Gaussian bridge, and not a uniform-caller bound in z. No illegal pointwise analytic-force conclusion is used. All physical callers are handled by mode centering and the numerical residual in Section 1.

## 5. Literal finite VALUE packet and actual ports

Draw two independent full Gaussian blocks Z,G after y is exposed. Use the SAME G in every node:

    q_i=r_i Z+sqrt(1-r_i^2)G,
    H_Q(Z,G)=sum_i w_i g_M(q_i),
    X_Q=b_M+sqrt(A)[Z-H_Q(Z,G)].                            (P)

Every q_i produces one original gradient at b_M+sqrt(A)q_i. Every query and alias is retained. The field H_Q is neither an exact mean nor a call to a smoothed posterior-force oracle.

Actual finite bounds, from the literal transcript and only the original Hessian sandwich, are

    Lip_(Z,G) H_Q <= A,
    ||H_Q||_p <= C_p A sqrt(D),
    D_(Z,G) X_Q =sqrt(A)[(I,0)+O(A)],
    D_y X_Q=I+O(A),
    H_Q(0,0)=0, X_Q(y;0,0)=b_M.

For the caller bound, D_y g_M(q_i)=sqrt(A)[Hess V(b_M+sqrt(A)q_i)-Hess V(b_M)]D_y b_M, which is O(sqrt(A)); multiplying by the output sqrt(A) gives O(A). No modulus of continuity of Hess V is used. Requested original first/adjoint sweeps are HVPs only at recorded VALUE points; no saved HVP is differentiated.

The canonical force sqrt(A) grad V(X_Q) has private first O(A) and caller first O(sqrt(A)). Its complete VALUE count is M+n+2, including the anchor and terminal canonical force, up to the convention for counting the initial mode anchor. Any original first sweep has the same number of original HVP sites. Gaussian dimension is exactly 2D; deterministic nodes do not add Gaussian blocks.

For a retained finite graph, weight node i by d_i=sqrt(w_i). Its known Gaussian row is d_i(r_i I,sqrt(1-r_i^2)I); the stack has operator norm at most one, and the terminal row with coefficients d_i also has norm one. Each node remains its own genuine original-gradient primitive at recorded inverse width 1/d_i. Those widths can be inverse fixed-order heat powers and are included in the numerical precision budget. H_Q itself is a rectangular output field; this does not declare it a genuine gradient on the whole (Z,G) record or grant a native retained/proxy theorem.

All weights have total mass one, so original numerical bias is propagated by actual absolute weights rather than divided by sqrt(n). Original caller numerical profiles stay those of the baseline finite mode plus |grad V(y)|+sqrt(D A). There is no D-power domain guard.

For a fixed maximum target J, use the SAME Q, chosen for delta=A^(J+2), in every compared level. Unlike a common finest quarter grid, this changes only a public-log query factor. It is legitimate to share this identical packet in all later fine/coarse sources. This statement does not supply an adjacent pair for the unconstructed higher corrections.

## 6. Exact first-order law identity and uncancelled next current

For psi with bounded derivatives sufficient for the following identity, put Y=Z-H_Q. Taylor's integral formula and (S) give

    E psi(Y)=E psi(Z)-Cov_gamma(U,psi)
             -E[(v_Q-v)(Z) dot grad psi(Z)]
             +integral_0^1(1-s) E[H_Q tensor H_Q:
                                      Hess psi(Z-sH_Q)] ds.  (W)

The last term is an actual same-root rank-two current. It includes all cross-node correlations generated by the common G. It is NOT a covariance oracle and is NOT zero. The target mu_U has its own second- and higher-order tilt currents. Matching the first-order term in (W), even with an exponentially accurate quadrature, does not prove a W2 source grade, a dimension-safe bound on the remaining currents, or a positive all-order law compiler.

In particular a bound obtained just by |H_Q|^2 costs A^2 D. The desired one-energy source return would require using its bounded first together with its actual Gaussian/root structure and preserving the target's matching currents. That is a separate theorem.

## 7. Exact quadratic ledger: shared roots matter

Let U(z)=z^T B z/2 with 0<=B<=A I. Write

    beta_Q=sum_i w_i sqrt(1-r_i^2).

By (Q0), the actual output is Gaussian:

    Y=(I-B/2)Z-beta_Q B G,
    Cov(Y)=I-B+(1/4+beta_Q^2)B^2.                         (G)

The target covariance is (I+B)^(-1). With accurate quadrature beta_Q tends to pi/4, so the B^2 coefficient tends to 1/4+pi^2/16, approximately .866850275. It is not the target coefficient one. Thus even perfect time integration leaves a real second-order covariance defect. Treating the G at different nodes as independent would incorrectly replace beta_Q^2 by sum_i w_i^2(1-r_i^2).

Using one random uniform r instead of deterministic Q is substantially worse in high dimension: conditional covariance is I-2rB+B^2. For B=a I, its random O(a) scale remains visible in the radius as D grows. Exact first-order moment matching after averaging r does not make that a dimension-safe O(a^2 sqrt(D)) source. The deterministic Q removes that particular random-scale failure by its exact first moment, but (G) records what remains.

## 8. Outstanding positive-source/outer contract

A successful continuation must supply, rather than assume:

1. A finite positive output correcting (W) and the matching target currents to arbitrary fixed order using VALUES and declared firsts only.
2. One-energy sqrt(D) estimates including common Z/G roots, every later observer, and actual covariance/cumulant descendants.
3. Actual named fine/coarse marginals with strong adjacent gaps if the direct-mean outer is retained; alternatively an explicitly weaker endpoint-law outer contract. The conditional-mean approximation (Q2) alone is not either port.
4. Exact zeros/callers/finite-mode restoration and all additional original queries and Gaussian blocks.
5. An order-to-complete-cost recurrence with heat exponent o(J), or eventually o(J), rather than an unpriced order-dependent compiler.

The deterministic unsmoothed quarter-grid is absent from (P). The new unresolved cost is the positive higher-law correction, not the analytical quadrature of the first-order Gaussian Stein field.
