# Independent audit of the nonlinear separator for the quadratic amplifier

2026-10-04.

## Verdict and exact scope

**PASS, scoped to the displayed VALUE graph and standardized family.** For every fixed polynomial level m>=1, the explicit fixed standardized potential u_eta in the source, with `U_A=A u_eta`, leaves a positive order-A^2 first-coordinate variance defect. Its uniform coefficient is at least

    delta_* = 0.012337794913467539... > 0.0122.

Consequently the standardized marginal W2 error has

    liminf_(A->0) W2/A^2 >= delta_*/2
                         = 0.006168897456733769...,

including when the admitted positive quadrature rule varies with A. The potential is genuinely C2, has a continuous non-Lipschitz Hessian, and has noncommuting Hessians at explicit points. No common-G product is discarded.

This refutes extrapolating the **particular quadratic polynomial VALUE recurrence** into a uniform general-C2 higher-order law amplifier. It does not contradict the earlier positive order-two standardized / order-5/2 physical theorem, the all-order matrix-quadratic specialization, or any different nonlinear algorithm.

Normalization matters: the fixed potential is u_eta in the standardized family U_A=A u_eta. An admissible original physical realization is

    V_A(x)=A u_eta(x/sqrt(A)),  y=b=0,

whose Hessian remains between .4 I and .6001 I and whose original VALUES realize `g_A=A grad u_eta`. This is a valid uniform counterexample family over the allowed original Hessian class. The note does not establish the same nonzero A^2 coefficient for one fixed **unscaled physical V** held constant while A tends to zero. No such stronger claim is admitted here.

Source:

- `../NONLINEAR-SEPARATOR-FOR-THE-QUADRATIC-AMPLIFIER.md`
- SHA256 `6fd7c048f1cf4980284d7ea196213aa1655b3cf8fabc7a18bfaf9b19963cee5d`.

The original main-law source remains separately pinned at `a085dfdef3f66f69208612034c17e4e18045fc00a21be3aed64ccfcb5e9b7085`; this audit makes no change to it or its earlier verdict.

Independent diagnostics in `check_separator_independent.py` pass **429 assertions** and write `separator_independent_checks.json`. They import no author diagnostic or earlier audit program. The separate manifest pins the source, audit, script, and output.

## 1. Scalar same-root moments and the exact coefficient

Write `b=.5`, `epsilon=.1`, `tau=exp(-1/2)`, and

    f(z)=bz+epsilon sin z,
    u(z)=bz^2/2+epsilon(1-cos z).

For every admitted positive rule, `sum w_i=1`, `sum w_i r_i=1/2`, and `q_i=r_i Z+s_i G` with one common standard G. Let `L=sum w_i q_i=Z/2+beta G`. Then `E L^2=c=1/4+beta^2`.

For correlated standard Gaussians with covariance rho,

    E[sin X sin Y]=exp(-1) sinh(rho),
    E[X sin Y]=rho tau.

The actual covariance is `rho_ij=r_i r_j+s_i s_j`. Hence

    E[H_Q^2] = c b^2+2c b epsilon tau+epsilon^2 d_Q,
    d_Q=exp(-1) sum_ij w_i w_j sinh(rho_ij),
    E[Z H_Q]=(b+epsilon tau)/2.

This derives source equation (1) exactly, not merely asymptotically. The cross term uses

    sum_i w_i Cov(L,q_i)=1/4+beta^2=c.

All off-diagonal node pairs survive. Independent node roots would alter both the linear covariance and d_Q.

For a centered Gaussian and the normalized density proportional to exp(-Au), a direct expansion gives

    E_A[Z^2]=1-A Cov(Z^2,u)
      +(A^2/2) Cov(Z^2,(u-Eu)^2)+O(A^3).

Twice integrating by parts gives

    Cov(Z^2,u)=E f',
    (1/2) Cov(Z^2,(u-Eu)^2)=E f^2+Cov(u,f').

For this fixture,

    E f^2=b^2+2b epsilon tau+(epsilon^2/2)(1-exp(-2)),
    Cov(u,f')=-(b epsilon tau)/2
       -epsilon^2[(1+exp(-2))/2-exp(-1)].

Their sum is exactly the target coefficient stated in source equation (2). This also verifies its sign and its factor 3/2.

## 2. What the actual nested VALUE graph contributes

Let `g_A=A f`, `Y=Z-AH_Q`. Since `f(0)=0`, `f'(0)=b+epsilon`, and f is globally .6-Lipschitz,

    g_A(g_A(Y))=A^2 f'(0) f(Z)+O_Lp(A^3),
    g_A(g_A(g_A(Y)))=O_Lp(A^3).

These are expansions of the actual nested VALUES, not replacements of a VALUE by an HVP. Positivity and unit total quadrature mass give a Q-uniform bound on every fixed Gaussian moment of H_Q, so replacing Y by Z in the leading term is legitimate uniformly over Q.

Each later layer contracts its Lp size by at most a constant times A^2. Thus T_k is O(A^(2k)) for k>=2, and the first binomial coefficient is -1/2. It follows that, at every fixed m>=1,

    Y_m=Z-AH_Q+(1-c)A^2(b+epsilon)f(Z)/2+O_Lp(A^3).

The sign is positive because T_1's order-two coefficient is c-1. Taking the second moment and subtracting the target coefficient gives exactly

    Delta_Q = b epsilon[(1-c)+(c-1/2)tau]
       +epsilon^2[d_Q+(1-c)tau-(exp(-1)-exp(-2))].

The expression `g(g(Z))` therefore reads the derivative at the zero/mode of f in this leading correction. It is not the local action `Df(Z) f(Z)`.

Because `rho_ij>=0`, d_Q>=0. The first bracket decreases as c increases and has minimum tau/2 on `[1/2,1]`. Therefore

    Delta_Q >= b epsilon tau/2-epsilon^2(exp(-1)-exp(-2))
            = 0.012837824913467539... .

The lower bound is uniform over the quadrature, including Q_A with increasingly accurate time integration. It is not based on the numerical value of the limiting beta.

## 3. General matrix second-moment identities require only C2

For an even anchored C2 potential u with bounded Hessian, f=grad u has linear growth and u has at most quadratic growth. The normalized target moment expansion is valid through order two with O(A^3) remainder by Gaussian integrability of powers of u. It requires no third spatial derivative.

For the matrix test `x x^T`, Gaussian integration by parts twice yields

    Cov(x x^T,u)=E Df,
    (1/2) Cov(x x^T,(u-Eu)^2)
      = E[f f^T]+Cov(u,Df).

The second covariance is entrywise. Product differentiation introduces the factor 2 in `D^2[(u-Eu)^2]`, which cancels the expansion's factor 1/2. Thus source equation (6) is correct for a genuinely C2 u.

Let `B_0=Df(0)`. Differentiability at zero gives

    g_A(g_A(Y))=A^2 B_0 f(Z)+o_L2(A^2).

A bounded derivative dominates the remainder, and uniform moments of Y allow truncation to a compact set followed by dominated convergence. The replacement of f(Y) by f(Z) uses only the global Lipschitz bound. These arguments are uniform in the admitted Q, since H_Q has Q-uniform moments of any fixed finite order.

The actual amplified map therefore has

    Y_m=Z-AH_Q+(1-c)A^2 B_0 f(Z)/2+o_L2(A^2).

Gaussian integration by parts gives `E[f(Z)Z^T]=E Df`, yielding

    C_source,2=E[H_Q H_Q^T]
       +(1-c)[B_0 E Df+(E Df)B_0]/2.

No commutation has been used: the order of both matrix products is retained. This proves source equation (7), including its orientation.

For the particular rough fixture one can make the remainder more explicit. Its rough h obeys `|h(t)|<=2|t|^(3/2)/3`, while `|sin t-t|<=|t|^3/6`. Hence

    |f_eta(x)-B_0 x|
      <= epsilon |x_1|^3/6 + (2 eta/3)|v.x|^(3/2).

The leading nested-VALUE expansion consequently has a Q-uniform remainder of order `O_m(A^3+eta A^(5/2))` in every needed fixed Gaussian norm. This reinforces, rather than replaces, the C2-only little-o argument. It is an analytical estimate; the producer does not query these derivatives.

## 4. Literal C2/noncommuting fixture and perturbation arithmetic

Let `v=(1,1)/sqrt(2)`, eta=1e-4, and take the source's odd h with

    h'(t)=sqrt(|t|)/(1+sqrt(|t|)), h(0)=0,
    psi'=h, psi(0)=0.

Then h is C1, psi is C2, psi is even/nonnegative, and

    |h(t)|<=|t|, 0<=psi(t)<=t^2/2, 0<=h'(t)<=1.

For

    u_eta(x)=b|x|^2/2+epsilon(1-cos x_1)+eta psi(v.x),

its Hessian is

    diag(b+epsilon cos x_1,b)+eta h'(v.x) vv^T.

This immediately gives `.4 I<=Hess u_eta<=.6001 I`. Along x=t v, the rank-one increment is asymptotic to `eta sqrt(|t|) vv^T`, while the cosine increment is O(t^2). The Hessian difference quotient therefore diverges, proving lack of Lipschitz continuity rather than merely a numerical failure of a chosen constant.

At zero, `B_0=diag(.6,.5)`. At x=(1,0), the off-diagonal entry is `eta h'(1/sqrt(2))/2>0`; its commutator with B_0 is nonzero. The independent numerical value of its Frobenius norm is about 3.23e-6, but noncommutation follows exactly from the displayed positive entry and the unequal diagonal entries.

For the (1,1) second-order coefficient, expand in eta. The endpoint E[H H^T] changes by at most `1.2 eta+eta^2`; the target E[f f^T] has the same bound. The origin Hessian is unchanged, so the amplifier anticommutator changes by at most `.6 eta`. These are conservative bounds; the exact v_1 factors would improve them.

For the covariance term, Cauchy–Schwarz and `Var F<=E F^2` give the valid estimates

    ||u_1||_2<=sqrt(3)/2,
    ||u_base||_2<=.3 sqrt(8),
    ||Df_base,11||_2<=.6,
    ||Df_1,11||_2<=1.

Here `u_base<=.3|x|^2` and `E|Z|^4=8` in dimension two. Thus the total linear-in-eta coefficient is bounded by

    1.2+.6+1.2+.6 sqrt(3)/2+.3 sqrt(8)
      =4.36814337969452... < 5,

and the quadratic coefficient by

    1+1+sqrt(3)/2=2.866025403784438... < 3.

This verifies source equation (8) independently. The estimates use Minkowski over the positive weights and each node's marginal Gaussian law; they make no independence assumption between nodes.

Therefore

    Delta_Q,11(eta)
      >=0.012837824913467539...-(5e-4+3e-8)
      =0.012337794913467539... > .0122,

uniformly over the admitted rule. The strict .0122 margin absorbs the Q-uniform little-o remainder for sufficiently small A at every fixed m.

## 5. Convert the moment defect into a law obstruction

The source graph is odd in its full root record because f_eta is odd. The target is even. All means therefore vanish, though the following inequality in fact works for raw second moments as well:

    W2(Law(X_1),Law(Y_1))
      >= | ||X_1||_2-||Y_1||_2 |.

It follows directly from the reverse L2 triangle inequality in every coupling. Projection onto the first coordinate contracts W2. Since both projected second moments approach one, division by their square-root sum yields

    liminf_(A->0) W2(Law(Y_m),mu_(A u_eta))/A^2
      >=delta_*/2>0.

When Q changes with A, a single coefficient limit need not exist; the uniform coefficient lower bound and uniform remainder instead give the displayed liminf. This is enough to rule out any o(A^2) general-C2 guarantee for this recurrence at a fixed m. Physical rescaling gives the corresponding obstruction at order A^(5/2) for the admissible family V_A described at the start.

## 6. Independent diagnostics and limits

The 429 independent checks cover:

- Exact analytic coefficient and perturbation constants, source hash, and six positive moment-exact rules: one midpoint, an endpoint pair, a nontrivial symmetric pair, and three dyadic rules up to 71 nodes. Endpoint/midpoint rules are additional controls; the asymptotically accurate dyadic family is included.
- The actual same-G scalar products `rho_ij` and `d_Q`, with an explicitly different independent-node-root counterfactual.
- Direct four-dimensional Gaussian integration of the actual two-dimensional rough common-root field H_Q and its full matrix second-order coefficient. Across the six rules the observed rough coefficient is between approximately .01716 and .02786, above the conservative uniform analytical lower bound.
- Two independent target-coefficient formulas, one using E[ff^T]+Cov(u,Df), the other directly expanding normalized second moments. Their Gaussian-quadrature discrepancy is about 4.92e-7; it is reported rather than treated as exact cancellation at the C2 cusp.
- Literal finite scalar VALUE programs through m=4, A down to 1/256, full three-call-per-layer counts, zero means, and convergence of the variance coefficient to its exact expression.
- Actual rough two-dimensional finite graph means and anchor zeros, explicit noncommutation, Hessian sandwich samples, and growing Hessian difference quotients. The global sandwich and failure of Lipschitz continuity are proved analytically above; samples are not their proof.

Floating-point quadrature checks have finite tolerances. They neither establish a universal theorem by sampling nor replace the analytical argument. The admitted separator is narrow: this particular finite polynomial graph cannot supply the missing nonlinear all-order compiler.
