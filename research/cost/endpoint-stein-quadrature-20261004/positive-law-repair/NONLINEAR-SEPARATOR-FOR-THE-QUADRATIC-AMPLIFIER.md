# A genuine C2 noncommuting separator for the quadratic amplifier

2026-10-04. Scoped continuation: this disproves general nonlinear order amplification for the specific polynomial VALUE graph in Section 5 of `POSITIVE-SECOND-ORDER-LAW-AND-QUADRATIC-AMPLIFIER.md`. It does not contradict that note's general positive order-two law theorem, its arbitrary-order quadratic restriction, or any different future nonlinear current compiler.

## Result

Every fixed m>=1 in the proposed quadratic polynomial graph still has a nonzero order-A^2 variance defect for an explicit fixed two-dimensional potential with `0.4 I<=Hess u<=0.6001 I`, with noncommuting Hessians, and whose Hessian is continuous but not Lipschitz. Thus even all its polynomial levels together do not improve the general-C2 law power beyond order two. All the same-G rank-one products are retained in this calculation.

The obstruction is explicit: nested `g(g(z))` reads the Hessian at the mode, while the nonlinear second-order target includes the covariance of the potential with its local Hessian. The latter does not disappear merely because the executed graph has small actual firsts.

## 1. Exact scalar calculation with shared roots

Start with

    u_0(z)=b z^2/2 + epsilon (1-cos z),
    f_0(z)=b z+epsilon sin z,
    b=1/2, epsilon=1/10,
    g_A=A f_0.

Then `0.4<=u_0''<=0.6` and f_0(0)=0. Use ANY positive endpoint rule satisfying the admitted exact moments. Put

    q_i=r_i Z+s_i G, s_i=sqrt(1-r_i^2),
    H_Q=sum_i w_i f_0(q_i),
    beta=sum_i w_i s_i,
    c=1/4+beta^2 in [1/2,1],
    rho_ij=r_i r_j+s_i s_j,
    d_Q=exp(-1) sum_ij w_i w_j sinh(rho_ij) >= 0,
    tau=exp(-1/2).

The actual unamplified endpoint is `Y=Z-A H_Q`. Its covariance is exact through degree two in A because the endpoint is affine in A:

    E Y^2 = 1-A(b+epsilon tau)
                +A^2[c b^2+2c b epsilon tau+epsilon^2 d_Q].    (1)

The shared-G cross-node covariance is precisely rho_ij. Replacing it by independent roots would change d_Q and is not performed.

The target has even density proportional to `exp(-z^2/2-Au_0(z))`. Expanding its normalized second moment, or using Gaussian integration by parts twice, gives

    E_target z^2 = 1-A E f_0'
          +A^2 [E f_0^2+Cov(u_0,f_0')] + O(A^3)
      =1-A(b+epsilon tau)
          +A^2[b^2+(3/2)b epsilon tau
                    +epsilon^2(exp(-1)-exp(-2))]+O(A^3).      (2)

This expansion is valid uniformly for sufficiently small A because u_0 has at most quadratic growth and is nonnegative. Only the displayed fixed analytic fixture is being expanded; no higher original derivative is imported into the program.

## 2. What every finite polynomial level actually does

For fixed m>=1 execute exactly the earlier VALUE recurrence

    T_0=Y,
    T_k=(c-1)g_A(g_A(T_(k-1)))
                +c g_A(g_A(g_A(T_(k-1)))),
    Y_m=sum_(k=0)^m binom(-1/2,k) T_k.

Since g_A(0)=0 and Lip(g_A)<=0.6A, `||T_k||_Lp=O_m(A^(2k))`. Also f_0'(0)=b+epsilon, so, in every fixed Gaussian Lp,

    Y_m=Z-AH_Q + [(1-c)/2] A^2(b+epsilon) f_0(Z)+O_m(A^3).

All k>=2 terms begin at A^4 and cannot change the second-order coefficient. Thus the source-minus-target variance defect, divided by A^2, is

    Delta_Q = b epsilon[(1-c)+(c-1/2)tau]
        +epsilon^2[d_Q+(1-c)tau-(exp(-1)-exp(-2))].            (3)

For c in [1/2,1],

    (1-c)+(c-1/2)tau >= tau/2,
    d_Q+(1-c)tau >= 0.

Consequently, for EVERY admitted positive rule,

    Delta_Q >= delta_0
       := (b epsilon tau/2)
               -epsilon^2(exp(-1)-exp(-2))
        > 0.0128.                                          (4)

This is an analytic lower bound, not a numerical test. The maps and target are odd/even as appropriate, so all means vanish. The elementary inequality

    W2(Law(X),Law(Y)) >= |sqrt(E X^2)-sqrt(E Y^2)|

shows an actual order-A^2 law lower bound. The conclusion also holds if Q varies with A: all leading identities and the positive lower bound are uniform, and the VALUE Taylor remainders have uniform Gaussian moment bounds from the unit total quadrature mass and the original Lipschitz bound.

## 3. Make the fixture genuinely C2 and noncommuting

Work in dimension two. Let `v=(1,1)/sqrt(2)`, eta=10^-4, and define

    h'(t)=sqrt(|t|)/(1+sqrt(|t|)), h(0)=0,
    psi'(t)=h(t), psi(0)=0,
    u_eta(x)=b |x|^2/2 + epsilon(1-cos x_1)
                                  +eta psi(v dot x).

Here h is odd, psi is even and nonnegative, and

    |h(t)|<=|t|, 0<=psi(t)<=t^2/2, 0<=h'(t)<=1.

The potential u_eta is C2, but its Hessian is not Lipschitz along v dot x=0. Its Hessian sandwich is

    0.4 I <= Hess u_eta <= (0.6+eta) I < I.                 (5)

At the origin `B_0=Hess u_eta(0)=diag(0.6,0.5)`. At x=(1,0) its off-diagonal Hessian entry is `eta h'(1/sqrt(2))/2>0`. Hence its Hessian there does not commute with B_0. This is a literal noncommuting C2 fixture using the same original potential at every node.

For any even anchored C2 potential u with gradient f, Gaussian bounded Hessian and at most quadratic growth, the second-order target second-moment matrix is

    C_target,2 = E[f f^T]+Cov(u,Df),                       (6)

where the matrix covariance is entrywise. For the amplified packet its second-order matrix is

    C_source,2 = E[H_Q H_Q^T]
             +[(1-c)/2][B_0 E Df+(E Df)B_0].              (7)

These are analytical identities; no full Hessian matrix or covariance is executed by the producer. Equation (7) retains every cross-root product in H_Q H_Q^T. The general C1 Taylor expansion of f at zero and at the Gaussian argument gives an o(A^2) remainder, which is enough for the lower bound. No third derivative or quantitative Hessian modulus is needed.

At eta=0, the (1,1) difference between (7) and (6) is exactly Delta_Q from (3). The following explicit perturbation bound makes the noncommuting continuation rigorous:

    |Delta_Q,11(eta)-Delta_Q,11(0)| <= 5 eta+3 eta^2.        (8)

To verify it, write `f_eta=f_base+eta f_1`, with `f_1=h(v dot x)v`, and similarly `u_eta=u_base+eta u_1`. The (1,1) endpoint force satisfies

    ||H_base,1||_2<=0.6, ||H_1,1||_2<=1.

So its covariance term changes by at most `1.2 eta+eta^2`. Since Df_1(0)=0, B_0 is unchanged; the amplifier cross term changes by at most `0.6 eta`. For the target's E f_1^2 term the same `1.2 eta+eta^2` bound holds. Finally

    ||u_1||_2<=sqrt(3)/2,
    ||u_base||_2<=0.3 sqrt(8),
    ||Df_base,11||_2<=0.6, ||Df_1,11||_2<=1.

Cauchy-Schwarz bounds the two linear covariance terms by

    eta[(sqrt(3)/2)*0.6+0.3 sqrt(8)]

and the quadratic one by `eta^2 sqrt(3)/2`. Their sum with the preceding terms is less than `5 eta+3 eta^2`, proving (8). All bounds are independent of the node count and covariance between its Gaussian roots.

For eta=10^-4 the perturbation is below 0.000501, whereas (4) exceeds 0.0128. Therefore

    Delta_Q,11(eta)>0.0122.                                (9)

The one-dimensional projection onto x_1 contracts W2. Its source and target second moments tend to one, so (9) gives a strictly positive order-A^2 lower bound on the full two-dimensional W2 error for every fixed m>=1. This remains true for increasing-accuracy quadrature rules Q_A.

## 4. Scope and next genuine current

The quadratic program remains useful and correct on arbitrary symmetric matrix quadratics. Its original VALUE count, root retention, actual firsts and zero are unaffected. The present calculation rules out only its extrapolation as a general nonlinear order amplifier.

Any successful nonlinear continuation must account for both `E[H_Q H_Q^T]` and `Cov(u,Df)`, including their orientation and shared roots, or produce an equivalent same-endpoint law current that cancels them. Neither a trace-only correction nor the nested mode-Hessian action supplied by g(g(z)) does so. This does not assert that the displayed analytical covariance/Hessian expressions are legal producer oracles.
