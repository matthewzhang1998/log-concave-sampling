# Fixed-dimensional existence of the separated unshifted covariance current

2026-10-05. Supplement to `SINGLE-HISTORY-EXACT-CENTERED-RESUMMED-CURRENT.md`, SHA256 `276d5873cc7d3f5ca8650955c08abe128cec7eda96680319f5c2a00b45118e9a`. This narrows its existence caveat at each fixed D. It does NOT prove the uniform one-energy A^4 sqrt(D) remainder.

## Claim

For the original C2 gradient and the exact true/finite-cheap source pair of the companion note, put J=Dg and C=C_H-C_Q. The separated current

    S=R1[C:D2g]

has a canonical L2(gamma_D) distributional definition and is the L2 limit of the corresponding smooth-regularization currents at each fixed dimension. A deliberately crude bound is

    ||S||2 <= C0 A^3(D+1)sqrt(D).                    (1)

The dimension loss in (1) is explicit. It is not an acceptable uniform one-energy estimate. The exact residual in companion equation (13) therefore exists as an L2 difference at fixed D, while its asserted A^4 sqrt(D) bound remains OPEN.

## 1. Covariance is C1 using only original first derivatives

For either source a in {H,Q}, let xi_a=I_a-mu_a and C_a=E[xi_a xi_a*]. The true source and cheap source have x Jacobians between zero and A I/2, because their positive clocks have first moment 1/2. Their mean Jacobians have the same bound. Thus

    ||D_x xi_a||op<=A/2

(the looser A bound would also suffice). Conditional Gaussian Poincare and the actual private-root Lipschitz constants give

    E[|xi_a|^2 | x] <= A^2 b_a^2 D,
    b_H=1/2, b_Q=beta_Q<=1.                         (2)

For H, the first constant follows from the OU kernel calculation int int e^(-t-u)c(t,u)dt du=1/4. These bounds hold uniformly in x.

Differentiate the original C1 VALUE functions under the conditional expectation. Their bounded x firsts and finite centered Gaussian moments justify this directly, or via finite analytical history approximations. Then C_a is C1 and

    partial_j C_a=E[(partial_j xi_a)xi_a*+xi_a(partial_j xi_a)*],
    div C_a=E[(D_x xi_a)xi_a+xi_a tr(D_x xi_a)].      (3)

Here (div C_a)_i=sum_j partial_j(C_a)_ij. Consequently

    ||C_a(x)||op<=A^2 b_a^2,
    |div C_a(x)|<= (A^2 b_a/2)(D+1)sqrt(D).         (4)

The latter bound intentionally does not claim a cancellation of the trace. Summing the two source bounds gives the same form for C with absolute constants. This is a fixed-D majorant, uniform over admissible smooth regularizations.

## 2. Distributional factorization and the Gaussian divergence bound

Because C is symmetric and C1, the row-divergence identity is

    C:D2g = div(J C)-J div C.                       (5)

This defines the left side as an order-one distribution under C2; it never asks for a pointwise D2g. Let M=J C. The bounded operator norms imply

    ||M||_(L2;HS)<=C0 A^3 sqrt(D),
    ||M x||2<=C0 A^3 sqrt(D),
    ||J div C||2<=C0 A^3(D+1)sqrt(D).               (6)

Write delta for Gaussian divergence, row by row:

    (delta M)_i=sum_j [x_j M_ij-partial_j M_ij],
    div M=Mx-delta M.                              (7)

R1 has Hermite multiplier 1/(n+1). Gaussian divergence maps a row/vector coefficient of chaos degree n into degree n+1 with norm at most sqrt(n+1) times the coefficient's L2 Hilbert norm. Therefore

    ||R1 delta M||2
       <= sup_(n>=0) sqrt(n+1)/(n+2) ||M||_(L2;HS)
       =(1/2)||M||_(L2;HS).                        (8)

The assertion extends from finite Hermite sums to every L2 matrix field by closure. It is a legitimate divergence estimate with an open physical output, not the invalid trace estimate for a coefficient times D2E.

Combining (5)-(8) gives the explicit definition

    S=R1(Mx)-R1 delta M-R1(J div C),                 (9)

and proves (1). The term J div C is where this elementary argument retains a large dimension allowance; no one-energy cancellation has been manufactured.

## 3. Consistency under smooth regularization

Mollify U and subtract its resulting constant gradient at zero, obtaining g_epsilon with the same Hessian bounds and anchor. Build BOTH source packets, their conditional means, and their covariance difference C_epsilon from that same g_epsilon. Original gradients and Hessians converge locally to g and J, while their global first/linear-growth bounds remain uniform.

For each fixed x and D, source values, centered values, and x derivatives converge in the required conditional L2 spaces. Equations (2)-(4) give uniform finite-D majorants. Hence C_epsilon and div C_epsilon converge pointwise, and (6) plus dominated convergence gives

    J_epsilon C_epsilon -> JC in L2(gamma;HS),
    J_epsilon C_epsilon x -> JCx in L2(gamma),
    J_epsilon div C_epsilon ->J div C in L2(gamma).

Use the bounded operators in (9) to conclude

    R1[C_epsilon:D2g_epsilon] -> S in L2(gamma_D).   (10)

The same definition is obtained by keeping C fixed and regularizing only g. Thus the separated current is not an ambiguous artifact of the coupled regularization.

The companion's complete resummed current already converges under C2 through its derivative-free homotopy. Its mean term also converges. Subtracting (10) therefore defines the unshifted trace residual in L2 at fixed D. Nothing in this argument bounds that residual by Lambda A^4 sqrt(D), nor constructs the original-VALUE native response consumer.

## Status

FIXED-D EXISTENCE: proved, with an explicit polynomial dimension majorant. UNIFORM A^4 ONE-ENERGY REDUCTION AND FINITE CONSUMER: still OPEN. The generic port-only trace separator remains fully compatible with this supplement.
