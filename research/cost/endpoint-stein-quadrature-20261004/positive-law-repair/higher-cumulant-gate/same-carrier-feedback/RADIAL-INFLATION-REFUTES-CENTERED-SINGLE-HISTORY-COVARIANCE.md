# Radial inflation refutes the centered single-history covariance reduction

2026-10-05. Genuine SAME-g counterexample. This note concerns the proposed reduction of the exact centered/resummed current to a deterministic conditional-covariance contraction at the coherent mean. It does not refute the exact live Stein current or any finite VALUE program that retains its complete shifted target. Independent exact-text review is requested.

## Result

Replacing D2g(x) by D2g(x-mu(x)) does **not** produce a dimension-uniform O(Lambda A^4 sqrt(D)) remainder in the stated class 0 <= Dg <= A I. There is a smooth convex-gradient family with

    A = D^(-1/4),   R = sqrt(D) = A^(-2),

for which the centered covariance formula has L2(gamma_D) error at least c_* A, while its proposed allowance is A^4 sqrt(D) = A^2. The separation survives every fixed public-polylogarithmic multiplier.

The proof uses the true continuous conditional OU history and any finite positive common-root rule satisfying the source note's first-moment and mean-quadrature contracts. Every source value, covariance vertex, and outer force uses literally the same g. All source correlations are retained.

The coherent center can be mu_H, mu_Q, or any pointwise convex combination mu_theta. The same lower bound holds for the theta-averaged centered covariance current. The exact mean-quadrature term is smaller than the displayed error.

The mechanism is different from the prior unshifted separator: at this dimension the random centered displacement has isotropic norm of order A sqrt(D)=A^(-1), and produces an order-one *radial inflation* A^2 sqrt(D)=1. Centering removes the coherent inward displacement; it does not resum the nonlinear response to this fluctuation-induced radial inflation.

## 1. One admissible original convex gradient

Use exactly the smooth step from the earlier radial note:

    eta(s) = exp(-1/(1-s^2/4)) for |s|<2, and 0 otherwise,
    Z = integral_R eta(s) ds,
    psi(s) = Z^(-1) integral_(-infinity)^s eta(u) du.

Thus psi is smooth, nondecreasing, zero on (-infinity,-2], and one on [2,infinity). Its derivative is an even positive bump on (-2,2), strictly decreasing on (0,2), and has supremum less than one. Every derivative has a fixed dimension-independent bound.

Choose

    c = 1/2,   d = 1/4,   k = 1-cA/2,
    R0 = k R.

For x != 0 set r=|x|, n=x/r, P=I-nn*, and

    f_D(x) = psi(r-R0)n,
    g_D(x) = c A x + d A f_D(x).                         (1)

Set f_D(0)=0. For D>=256, R0-2 is larger than R/2, so f_D vanishes in a neighborhood of zero and is smooth everywhere. It is the gradient of the radial convex function with radial derivative psi(r-R0). Therefore g_D is the gradient of

    U_D(x) = (cA/2)|x|^2 + dA integral_0^|x| psi(u-R0) du.

Its Hessian is

    Dg_D(x) = c A I
          + d A [psi'(r-R0)nn* + psi(r-R0)P/r].          (2)

All eigenvalues are nonnegative. The radial eigenvalue is at most (c+d)A< A; the tangential one is at most (c+d/(R0-2))A< A. Consequently g_D(0)=0 and 0 <= Dg_D <= A I globally. This is a literal smooth member of the admitted C2 class.

The shell is deliberately placed at R0=(1-cA/2)R, the radius reached after the large coherent linear mean displacement. No source is changed when evaluating the terminal force.

## 2. Genuine continuous and finite sources

Take any admitted positive clock rule Q with clocks tau_j in (0,1), weights v_j>0, and

    sum_j v_j = 1,
    sum_j v_j tau_j = 1/2,
    ||Q-R1||_(L2(gamma)->L2(gamma)) <= delta <= A^2.

Write d_j=sqrt(1-tau_j^2) and beta_Q=sum_j v_j d_j. The same Hermite multiplier argument as in the prior radial note gives

    |sum_j v_j tau_j^2 - 1/3| <= delta,
    1 >= beta_Q >= 2/3-delta >= 7/12,
    beta_Q^2 - 1/4 >= 13/144.                           (3)

For the TRUE conditional OU history and the actual FINITE common-H packet, linearity of the cAx part gives exact identities

    I_a(x) = c A [x/2 + sigma_a G_a] + d A K_a,
    sigma_H=1/2,   sigma_Q=beta_Q,                       (4)

where G_a is a standard Gaussian vector independent of x. G_H is the normalized linear functional of the true conditional OU Brownian history; G_Q is the actual shared root. The nonlinear terms are

    K_H = integral_0^infinity e^(-t) f_D(X_t^x) dt,
    K_Q = sum_j v_j f_D(tau_j x+d_j G_Q),
    |K_a| <= 1 pathwise.                                (5)

G_a and K_a are generally correlated. Nothing in this proof makes them independent or replaces either source by a Gaussian law.

Let M_a(x)=E[K_a|x]. Then

    mu_a(x)=E[I_a|x]=cAx/2+dA M_a(x),   |M_a|<=1,
    ||M_H-M_Q||2 <= delta,
    ||mu_H-mu_Q||2 <= d A delta.                         (6)

The bound in (6) is exactly the mean-operator quadrature estimate on the bounded original VALUE field f_D. No path quadrature or path-grid approximation is made.

## 3. Uniform shell geometry at the larger dimension

Let X~gamma_D, S_D=|X|-R, and n(X)=X/|X|. For any sigma in [0,1], independently take G~gamma_D and set

    y_sigma = k X - c A sigma G,
    h_sigma = c^2 sigma^2/2.

Then, uniformly in sigma,

    f_D(y_sigma)
      = psi(S_D+h_sigma)n(X) + O_(L2(X,G))(A).            (7)

This is a radius/direction estimate, not a small Euclidean-displacement Taylor approximation. Here is an explicit proof.

Write rho=k|X|, n0=n(X), and w=-cA sigma G. On the event |X|>=R/2 and |w|<=rho/2,

    |rho n0+w|-rho
      = n0.w + (|w|^2-(n0.w)^2)/(2rho)
                        + O(|w|^3/rho^2),               (8)
    |n(rho n0+w)-n0| <= C |w|/rho.

For each fixed p the Gaussian shell moments ||S_D||p are uniformly bounded. Gaussian moments give

    ||n0.w||p = O(A),
    ||w||p = O(A R)=O(A^(-1)),
    || |w|^3/rho^2 ||_(Lp; |X|>=R/2) = O(A^3 R)=O(A),
    || |w|/rho ||_(Lp; |X|>=R/2) = O(A).

The centered chi-square part of |w|^2-(n0.w)^2 has Lp size O(A^2 R)=O(1); division by rho makes this O(A^2). Its conditional mean divided by 2rho equals

    c^2 A^2 sigma^2 (D-1)/(2k|X|)
       = c^2 sigma^2/2 + O_Lp(A),                        (9)

because A^2 R=1, k=1+O(A), and |X|=R+O_Lp(1). Also

    rho-R0=k S_D=S_D+O_Lp(A).

Equations (8)-(9) show that the radius relative to R0 equals S_D+h_sigma+O_Lp(A), while the direction equals n(X)+O_Lp(A). The complements of the stated events have exponentially small Gaussian probability; boundedness of f_D, psi and the usual polynomial moment bounds make their L2 contributions negligible. Since psi is Lipschitz, these estimates prove (7).

The crucial difference from A=D^(-1/2) is the constant inflation h_sigma. It is not an O(A) correction in this family.

## 4. Actual source mean difference

The true terminal argument from (4) is y_(sigma_a)-dA K_a. The global Lipschitz constant of f_D is bounded by a fixed constant, so without breaking the correlation between K_a and G_a,

    f_D(y_(sigma_a)-dA K_a)
      = f_D(y_(sigma_a)) + O_(L2)(A).

Apply (7), take conditional expectations, and multiply by the outer dA. The cAx part of the outer g has difference -cA(mu_H-mu_Q), of norm at most cd A^2 delta by (6). Thus the literal same-g means satisfy

    j_H-j_Q
      = dA [psi(S_D+h_H)-psi(S_D+h_Q)]n(X)
                                  + O_L2(A^2),           (10)

where

    h_H=c^2/8,   h_Q=c^2 beta_Q^2/2.

The O(A^2) constant is uniform in all admitted Q. It does not depend on a path approximation or on independence of the bounded perturbation K_a.

## 5. Actual covariances at any coherent center

Expanding the covariance of the actual decomposition (4) gives

    C_a(x)=Cov(I_a|x)
          =c^2 A^2 sigma_a^2 I + B_a(x),
    ||B_a(x)||HS <= C A^2.                               (11)

For clarity, this dimension-free remainder estimate is a legitimate *first-chaos coefficient* bound, not a generic contracted-trace inference. For each component K_l, the variables (G_a)_i form an orthonormal collection in conditional L2. Therefore

    sum_(l,i) |Cov((K_a)_l,(G_a)_i|x)|^2
       <= sum_l Var((K_a)_l|x)
       <= E[|K_a-M_a|^2|x] <= 1.

Also ||Cov(K_a|x)||HS <= tr Cov(K_a|x) <=1. These statements preserve all G_a/K_a correlations and prove (11) by direct covariance expansion.

Take any pointwise coherent center

    mu_*(x)=cAx/2+dA M_*(x),   |M_*(x)|<=1,               (12)

including mu_H, mu_Q, or their arbitrary pointwise convex combinations. Put z_*=x-mu_*(x). It obeys

    |z_*|-R0=S_D+O_L2(A),
    n(z_*)=n(X)+O_L2(A/R),                               (13)

on the high-probability shell; bounded global derivative estimates handle the complement.

For any matrix B the exact radial derivative formula is

    D2f_D(z):B
      = psi'' n(n*Bn)
        + a_r [n tr(PB)+P(B+B*)n],
    a_r=psi'/r-psi/r^2,                                  (14)

with scalar factors evaluated at r-R0. Since R0 is comparable to R=sqrt(D),

    sup_z ||D2f_D(z)||_(HS->vector) <= C,
    sup_z ||D2g_D(z)||_(HS->vector) <= C A.               (15)

In particular the B_H-B_Q remainder in (11), contracted with D2g_D(z_*), is O_L2(A^3). The large isotropic trace is calculated explicitly rather than bounded generically. Setting B=I in (14) gives

    Delta g_D(z_*)
      = dA [psi''+(D-1)(psi'/r-psi/r^2)]n
      = d A R psi'(S_D)n(X)+O_L2(1).                    (16)

Indeed r=kR+O_L2(1), so (D-1)/r=R(1+O(A))+O_L2(1); the O(A) change of psi' in (13), after multiplication by A R=A^(-1), is O(1). The psi'' and psi/r^2 pieces are smaller.

Combining (11), (15), (16), and A^2 R=1 now gives

    (1/2)(C_H-C_Q):D2g_D(x-mu_*)
      = dA(h_H-h_Q)psi'(S_D)n(X)+O_L2(A^2).             (17)

The error is uniform over all centers (12). Therefore (17) also holds after integrating over theta for mu_theta=(1-theta)mu_Q+theta mu_H. This covers exactly the deterministic centered covariance term obtained by freezing the shifted derivative in the exact convolution/Stein identity.

## 6. A uniform linear resolvent witness

Define

    q_beta(s)
      = psi(s+h_H)-psi(s+h_beta)
                     -(h_H-h_beta)psi'(s),
    h_beta=c^2 beta^2/2.

Subtracting (17) from (10) yields the pre-resolvent discrepancy

    d A q_(beta_Q)(S_D)n(X) + O_L2(A^2).                 (18)

Let S~N(0,1/2), and define m(t)=E[psi'(S+t)]. Since psi' is an even nonzero strictly unimodal bump, its Gaussian convolution m is strictly decreasing for t>0. An elementary proof writes psi' as a positive layer-cake mixture of centered intervals: each centered Gaussian interval mass is strictly decreasing under a positive shift.

By (3), h_beta>h_H>0, uniformly over beta in [7/12,1]. Consequently

    E q_beta(S)
      = integral_(h_H)^(h_beta) [m(0)-m(t)] dt
      >= integral_(h_H)^(h_min) [m(0)-m(t)] dt
      =: kappa > 0,                                    (19)

where h_min=c^2(7/12)^2/2>h_H.

The Gaussian shell CLT gives S_D -> S. The family q_beta, beta in [7/12,1], is uniformly bounded and uniformly Lipschitz; hence the convergence of E q_beta(S_D) is uniform in beta (for example by a finite net in beta). Since |X|/R ->1 in L2, also

    inf_(beta in [7/12,1])
        E[q_beta(S_D)|X|/R] >= kappa/2                  (20)

for all sufficiently large D.

Use the vector test T_D(X)=X/R, which has L2 norm exactly one. The Gaussian resolvent is self-adjoint and R1 T_D=T_D/2. Thus (18)-(20) imply

    < R1[ j_H-j_Q
          -(1/2)(C_H-C_Q):D2g_D(x-mu_*) ], T_D >
       >= d kappa A/4 - C A^2.

For all sufficiently large D, uniformly over the admitted finite Q and all centers (12),

    || R1(j_H-j_Q)
       -(1/2)R1[(C_H-C_Q):D2g_D(x-mu_*)] ||2
       >= c_* A,   c_*=d kappa/8>0.                     (21)

The same statement holds for the theta-integrated centered covariance term. The allowed fourth-order size is

    A^4 sqrt(D) = A^2,

so the ratio in (21) diverges at least as c_*/A. Every fixed public-log multiplier is dominated by this power separation.

## 7. Relation to the exact centered Stein identity

This is a failure of a proposed *reduction*, not a failure of the exact source-qualified identity. In the latter,

    T_theta=x-mu_theta-sqrt(theta)xi_H-sqrt(1-theta)xi_Q,

and the derivative D2g(T_theta) remains correlated with the complete source Stein matrices tau_H and tau_Q. Replacing it by D2g(x-mu_theta), then averaging tau_H-tau_Q to C_H-C_Q, produces precisely the theta-averaged current refuted above.

The exact mean term is bounded by

    A||mu_H-mu_Q||2 <= d A^2 delta <= d A^4.

It cannot absorb an Omega(A) discrepancy. Thus the fourth-order residual in the centered version of the exact identity is itself Omega(A), up to this negligible mean term.

No generic Hilbert/Bessel bound was used to justify a vector trace; that trace was evaluated explicitly. No independent Gaussian law was substituted for an original-g source. The only Gaussian approximation in the geometry is an exact linear term of the SAME g plus its pathwise bounded nonlinear remainder. No ordinary Markov path-grid producer is introduced.

## 8. Scope and next requirement

The theorem is dimension-uniform and therefore allows D=A^(-4); the source note does not impose A sqrt(D)=O(1). If a different proposed theorem explicitly restricts that product, this counterexample does not address that narrower regime.

A successful resummed target in the unrestricted dimension regime must retain more than the coherent mean and a single covariance derivative at that mean. This family already requires the nonlinear response to the covariance-induced order-one radial displacement. Keeping the full live shifted current is consistent with the example. No conclusion about its finite original-VALUE consumer follows from this counterexample alone.

The proof also explains why simply increasing a covariance Taylor truncation by a fixed order has no automatically small uniform expansion parameter here: A^2 sqrt(D)=1. It does not assert that every particular higher-order formula is refuted without identifying that formula and proving its residual.

Status: analytical SAME-g centered-covariance separator proved as above; independent exact-text review pending. General exact resummed-current consumption remains open.

## 9. Constructive non-obstruction: a full shifted Gaussian VALUE response survives

For this explicit family, define the full Gaussian response

    j_a^G(x) = E_N g_D(kx-cA sigma_a N),   N~gamma_D,
    sigma_H=1/2,   sigma_Q=beta_Q.                        (22)

This keeps the original g_D at its shifted Gaussian argument, rather than Taylor expanding its covariance response. It admits an immediate source-qualified estimate using the SAME G_a from (4):

    | g_D(kx-cA sigma_a G_a-dA K_a)
          -g_D(kx-cA sigma_a G_a) |
       <= Lip(g_D) dA |K_a| <= d A^2.

Consequently, pointwise in x,

    |j_a(x)-j_a^G(x)| <= dA^2,
    ||R1[(j_H-j_Q)-(j_H^G-j_Q^G)]||2 <= 2d A^2.           (23)

At D=A^(-4), this is exactly an O(A^4 sqrt(D)) error. Formula (7) also shows directly that (22) preserves the required nonlinear shell response psi(S_D+h_a).

The random variable g_D(kx-cA sigma_a N) needs only one evaluation of the SAME original g_D per branch, plus affine Gaussian arithmetic with the explicit family constants c,A,sigma_a. The raw comparison of the two branches needs two such VALUES; N may be shared. It uses no conditional-mean oracle, covariance oracle, Stein matrix, or Hessian-vector product. This is a family-specific analytic proxy and a finite raw shifted-VALUE source, not a general compiler for its expectation and not a general same-carrier m3 theorem.

Thus the centered derivative reduction fails while a genuinely resummed shifted Gaussian VALUE response remains viable even on the counterfamily. The counterexample does not obstruct the full shifted-current route.
