# A finite positive VALUE graph with dimension-uniform A^(11/4) target bias

## Executive result

Let g=grad U, U in C2, g(0)=0 and 0 <= Dg <= A I, with 0<A<=1/12. Let X be stationary standard OU, and define

    H_t = integral_0^infinity exp(-s) g(X_(t+s)) ds,
    B = integral_0^infinity exp(-t) g(X_t-H_t) dt,
    psi(x) = E[g(X_0-B) | X_0=x],
    m3 = R1 psi.

There is an explicit finite original-g VALUE source using J+2 original VALUES and 2D independent Gaussian coordinates, where J=O(log^2(1/A)), whose exact outer-resolvent own mean differs from m3 by at most

    5 A^(11/4) sqrt(D).                                      (E1)

Its formula retains the terminal nonlinearity and uses only positive history weights. It uses no strong mean replication, no derivatives of the original source, and no expectation oracle. A finite positive outer resolvent rule adds a separately priced operator-norm error and one D-dimensional Gaussian root. This is a target-bias theorem; completion of a Gaussian-own-mean service must still be certified separately.

The improved exponent comes from delaying the genuinely random history away from X_0, taking advantage of the Gaussian conditional variance X_0 | X_delta, and doing Gaussian integration by parts only in that conditional variable. The terminal Hessian is transferred to the resolvent test, so no bound on D2g is used.

The graph does change the old stage-two graph. It is not claimed to be exactly quadratic-exact. Its entire source tree is given below, and all descendants use the original g.

## 1. Exact weak conditional-noise lemma

Let (X,Y) be jointly standard Gaussian in R^D with covariance q I, 0<=q<1, and set sigma=sqrt(1-q^2). Given Y, let Z be an independent standard Gaussian bank, independent of X, and let b(Y,Z) be L-Lipschitz in that complete bank, uniformly in Y. Set m(Y)=E_Z b(Y,Z). Let F:R^D -> R^D be K-Lipschitz. Define

    e(x)=E[F(X-b(Y,Z))-F(X-m(Y)) | X=x].

Then

    ||R1 e||_(L2 gamma)
       <= (K L^2 sqrt(D)/2) (1/2 + 1/sigma).                  (1.1)

No derivative bound on b in Y is required. F need not be a gradient. Infinite Gaussian banks are allowed for the proof by strong finite-dimensional approximation; the executed graph below uses only a finite bank.

### Proof

It suffices first to treat smooth F and b, and smooth vector tests f; approximation preserving the displayed Lipschitz bounds gives the general result. Write h=R1 f. Hermite spectral calculus gives

    ||h||2 <= ||f||2,
    ||Dh||2 <= (1/2)||f||2,                                 (1.2)

because sup_(n>=0) n/(n+1)^2=1/4. The matrix derivative norm here is Hilbert-Schmidt, summed over output components.

Condition on Y and write e_Z=b(Y,Z)-m(Y). Interpolate lambda in [0,1]. The Gaussian covariance identity, with Z_t=exp(-t)Z+sqrt(1-exp(-2t))Z', gives

    d/dlambda E F_i(X-m-lambda e_Z)
      = lambda integral_0^infinity exp(-t)
          E sum_(j,k) C_(jk)(Y,Z,Z_t)
             partial_k partial_j F_i(X-m-lambda e_(Z_t)) dt,

where

    C(Y,Z,Z_t)=D_Z b(Y,Z) D_Z b(Y,Z_t)^T,
    ||C||op <= L^2.                                        (1.3)

The covariance identity follows by applying Cov(a(Z),phi(Z)) = integral exp(-t) E[Da(Z) dot Dphi(Z_t)]dt coordinatewise. Its use requires no derivative of Db. Importantly, C, m and e_(Z_t) are independent of X conditional on Y.

Now X=qY+sigma eta, where eta is standard Gaussian independent of Y,Z,Z'. Integrate one displayed derivative with respect to X conditional on Y. With

    M = DF(X-m-lambda e_(Z_t)) C,

one has ||M||op<=K L^2 and ||M||HS<=sqrt(D) K L^2. Thus

    |E sum_(i,j,k) h_i C_(jk) partial_k partial_j F_i|
      <= sqrt(D) K L^2 ||Dh||2
           + K L^2 sigma^(-1) E[|h(X)| |eta|]
      <= K L^2 sqrt(D) (||Dh||2 + sigma^(-1)||h||2).          (1.4)

The last step is ordinary Cauchy-Schwarz under the original joint probability law; independence between h(X) and eta is NOT assumed. Integrating lambda costs 1/2 and integrating exp(-t) costs 1. Self-adjointness of R1 and (1.2) prove (1.1).

For nonsmooth Lipschitz b, use Gaussian Sobolev approximation. For the C1 F in the application, Euclidean mollification preserves its Lipschitz bound and converges uniformly up to a deterministic O(K epsilon sqrt(D)) allowance, which tends to zero. Thus the Hessian displayed during the proof is not an executed oracle or a source hypothesis.

### Two-source corollary

If b_1,b_2 have bank Lipschitz constants L_1,L_2, the same conditional mean m(Y), and each is conditionally independent of X given Y, then

    ||R1 E[F(X-b_1)-F(X-b_2) | X]||2
      <= (K (L_1^2+L_2^2) sqrt(D)/2)(1/2+1/sigma).           (1.5)

If their means m_1,m_2 differ, add K||m_1(Y)-m_2(Y)||2. This latter term uses ordinary conditional Jensen and resolvent contraction.

## 2. The executable graph

Choose a scalar w in (0,1), put q=1-w and sigma=sqrt(1-q^2). Let Q be a deterministic finite positive OU rule

    Q f = sum_(j=1)^J a_j P_(r_j) f,
    a_j>=0, sum a_j=1, sum a_j r_j=1/2, 0<=r_j<=1,
    ||Q-R1||_(L2->L2) <= epsilon.                           (2.1)

Here P_r denotes the OU operator with correlation r. Section 5 gives a positive finite construction with J=O(log^2(1/epsilon)). The exact first moment is not needed for the target-bias proof, but is imposed for the cleaner native ports in the companion note; Section 5 supplies its positive correction with no asymptotic cost.

Given caller x and independent Z,N ~ N(0,I_D), execute

    Y = q x + sigma Z,
    U_j = r_j Y + sqrt(1-r_j^2) N,       for j=1,...,J,
    V_Q = q sum_j a_j g(U_j),
    u = x - V_Q,
    T_Q(x;Z,N) = g(u - w g(u)).                             (2.2)

Every U_j has its exact prescribed joint law with Y. The same N is deliberately used for every j. No claim is made that this finite cross-node joint law is a true OU path. Only its individual OU marginals and its whole-bank Lipschitz bound are used.

This is J+2 original VALUES: J evaluations g(U_j), then g(u), then g(u-wg(u)). In particular no terminal expectation, centroid, smoothed source, or derivative oracle is executed. The two occurrences defining the terminal composition are retained coherently.

Write mu_Q(x)=E_(Z,N) T_Q(x;Z,N) for the source's own Gaussian mean, which is a mathematical description, not an instruction to execute that expectation.

## 3. Target-error theorem with explicit parameters

For every w in (0,1) with w A<=1,

    ||R1 mu_Q - m3||2
      <= sqrt(D) [ A^3
                   + (2 sqrt(2)/3) A^2 w^(3/2)
                   + q w A^3
                   + q A^2 epsilon
                   + q^2 A^3 (1/2 + 1/sigma) ].            (3.1)

Set w=sqrt(A), epsilon=A^(3/4) for 0<A<1. Since sigma^2=2w-w^2>=w,

    ||R1 mu_Q-m3||2
      <= [2 + 2sqrt(2)/3 + (3/2)A^(1/4) + A^(3/4)]
            A^(11/4) sqrt(D)
      <= 5.45 A^(11/4) sqrt(D).                            (3.2)

For the intended range 0<A<=1/12, the bracket is below 4, so in particular (E1) holds. More generally the explicit estimate (3.1) holds whenever 0<w<1 and wA<=1; one can use 6 in place of 5 for 0<A<=1. At A=1 the limiting q=0 graph is well-defined and satisfies the limiting bound; no infinite time is executed.

### Proof

First let H_0=integral exp(-t)g(X_t)dt. Stationarity, anchoring and Lipschitzness give

    ||B-H_0||2 <= A^2 sqrt(D),

so replacing B by H_0 in the terminal original g costs at most A^3 sqrt(D), including the outer resolvent.

Choose delta=-log q and split H_0=H_short+V, with

    V=integral_delta^infinity exp(-t)g(X_t)dt
     =q integral_0^infinity exp(-s)g(X_(delta+s))ds.

Stationary OU increments give

    ||H_short-wg(X_0)||2
      <= A sqrt(2D) integral_0^delta exp(-t)sqrt(1-exp(-t))dt
      = (2sqrt(2)/3) A sqrt(D) w^(3/2).                     (3.3)

Define F(u)=g(u-wg(u)). Since 0<=Dg<=AI and wA<=1, u->u-wg(u) is a global contraction, hence F is A-Lipschitz and F(0)=0. Commuting the short deterministic prefix into F costs only

    |g(x-wg(x)-V)-F(x-V)|
       <= A w |g(x)-g(x-V)|
       <= A^2 w |V|.

Because ||V||2<=q A sqrt(D), this costs q w A^3 sqrt(D). This is a direct value comparison; no differentiation through the commutation is used.

Let Y=X_delta. Conditional on Y, V depends only on the future after delta, so V is independent of X_0. In its Gaussian future noise bank it is q A-Lipschitz: each OU time row has Cameron-Martin operator norm at most one, and the positive history weights have total q. Its conditional mean is

    E[V|Y]=q R1 g(Y).                                      (3.4)

For V_Q in (2.2), the conditional noise bank is just N. It is qA-Lipschitz, and

    E[V_Q|Y]=q Qg(Y).                                      (3.5)

The mean discrepancy has L2 norm at most q epsilon A sqrt(D). Apply the two-source conditional-noise corollary with F of Lipschitz constant A and L_1,L_2<=qA. This gives the final two terms in (3.1). Adding the preceding two coherent comparisons proves the theorem.

## 4. Source-valid bounds and finite outer rule

Under a stationary caller x, every U_j is standard. Thus

    ||V_Q||p <= q A ||G_D||p,
    ||T_Q||2 <= A(1+qA) sqrt(D),
    ||T_Q-g(x)||p <= A^2(1+wqA)||G_D||p.                    (4.1)

The latter follows by bounding the total shift V_Q+w g(x-V_Q), whose Lp norm is at most A[q+w(1+qA)]||G_D||p.

For the complete private (Z,N) Gaussian bank, each U_j has operator norm sqrt(r_j^2 sigma^2 + 1-r_j^2)<=1. Hence

    Lip_(Z,N) V_Q <=qA,
    Lip_(Z,N) T_Q <=q A^2.                                 (4.2)

Its caller first satisfies ||D_x V_Q||<=q^2 A sum a_j r_j<=q^2 A. It is symmetric positive semidefinite because all scalar row coefficients are nonnegative and Dg is symmetric. No derivative of a clock or numerical inverse is present.

At caller x=0 and private roots Z=N=0 all original VALUE sites vanish exactly. If g=0, the entire graph is identically zero. A caller origin at nonzero x must still be retained and restored in a native completion.

If every original VALUE has absolute error at most nu, then the history error is at most q nu; the u error is q nu; the first terminal evaluation has error at most (1+qA)nu; the final terminal input has error at most [q+w(1+qA)]nu=(1+wqA)nu; and the final VALUE error is at most [1+A+wqA^2]nu. A separately evaluated residual baseline g(x) adds nu. These are occurrence-local exact-arithmetic propagation bounds; clock/row/weight and replay errors remain separately priced.

For a deterministic positive outer rule Q_out with operator error epsilon_out, execute its source using one additional independent G root and common (Z,N) roots at all outer nodes. Its mean is Q_out mu_Q, and

    ||Q_out mu_Q-m3||2
       <= right side of (3.1)
           +epsilon_out A(1+qA)sqrt(D).                    (4.3)

Choosing epsilon_out=O(A^(7/4)) makes this additional error O(A^(11/4)sqrt(D)). The terminal original-VALUE bill is J_out(J+2) before legitimate exact-key sharing. A baseline/residual split uses J+3 per residual occurrence. Native completions must count every repeated occurrence, complete bank, retained origin, fill, mode, numerical source version, derivative/adjoint pass, and replay; the small source graph does not eliminate those costs.

## 5. Explicit positive polylogarithmic resolvent quadrature

For 0<epsilon<=1, set a=epsilon/8 and T=log(8/epsilon). Partition [a,T] into J0<=ceil(log2(T/a)) intervals [s,t] with t<=2s. In each interval use n-point Gauss-Legendre quadrature, where

    n >= ceil(log_4(16T/epsilon)).

For each Gauss node t_j with positive integration weight c_j, put raw OU weight c_j exp(-t_j) at correlation r_j=exp(-t_j). Add endpoint raw weights 1-exp(-a) at correlation 1 and exp(-T) at correlation 0. Normalize the entire positive weight list to mass one.

For every Hermite degree k>=0, the relevant integrand on time intervals is exp(-(k+1)t). The degree-(2n-1) Taylor polynomial at the midpoint of [s,t] has uniform error at most 4^(-n): its remainder is bounded by exp(-lambda s)[lambda(t-s)/2]^(2n)/(2n)! <=4^(-n), using lambda=k+1, t-s<=s, and m! >=(m/e)^m. Positive quadrature exactness therefore gives interval error at most 2(t-s)4^(-n). The summed interior error is at most 2T4^(-n)<=epsilon/8. Early and late replacements each contribute at most epsilon/8. Raw total mass differs from one by at most epsilon/8. After normalization, the uniform Hermite-multiplier error is at most

    (epsilon/2)/(1-epsilon/8) <= (4/7)epsilon <epsilon.

Thus (2.1) holds by Hermite orthogonality, independently of D, with

    number of nodes <= 2+n ceil(log2(T/a))
                    =O(log^2(1/epsilon)).                  (5.1)

All nodes and weights are known scalar quadrature data; no original-source oracle is used to produce them. Their finite-precision construction must be certified under the actual scalar setup rules. Endpoint P_1 is g(Y), and P_0 is g(N), so neither endpoint requires infinite-time simulation.

If exact first exponential moment is desired, start with operator error epsilon/3. Let m1=sum a_j r_j. If m1>1/2, mix Q with P_0 until the first multiplier is 1/2; if m1<1/2, mix with P_1 instead. The required extra mass is at most 2epsilon/3 and the total operator error remains at most epsilon. This preserves positivity and mass one, and adds at most one endpoint node (which can be merged if already present).

## 6. A separate exact terminal-OU commutator estimate

For completeness, suppose Kf(x)=E f(x-b(x,Z)), with ||b||2<=L sqrt(D) and Lip_x b<=L. Direct coherent coupling yields

    ||(KP_r-P_rK)g||2
       <= A L [(1-r)+sqrt(2(1-r))]sqrt(D).

Since ||R1(P_r-I)||_(2->2)<=1-r,

    ||R1 K(P_r-I)g||2
       <= (1-r)||Kg||2
           +A L [(1-r)+sqrt(2(1-r))]sqrt(D).

Thus a terminal OU smoothing of variance scale A^2 costs O(A^3sqrt(D)) weakly when L=O(A). This useful source-valid perturbation does not itself establish the finite-history accuracy theorem above; Sections 1-3 supply the actual construction.
