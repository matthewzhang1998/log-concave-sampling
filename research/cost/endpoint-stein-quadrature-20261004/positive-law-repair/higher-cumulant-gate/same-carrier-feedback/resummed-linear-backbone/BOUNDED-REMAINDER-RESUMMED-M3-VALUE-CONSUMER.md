# A finite resummed m3 VALUE consumer for a bounded-remainder linear backbone

2026-10-05. Bounded positive construction; independent audit required. This is not a generic centered-current closure.

## Result

Assume the original source satisfies

    g=grad U, g(0)=0, 0<=Dg<=A I, 0<A<=a0,
    g(x)=lambda x+r(x), 0<=lambda<=A,
    sup_x |r(x)|<=epsilon.

The scalar lambda and the certified bound epsilon are declared input parameters, frozen before differentiation. No value of r, derivative, covariance, conditional expectation, or continuous path is a producer leaf. The program uses the original g only.

Let m3 be the canonical Gaussian-history target

    m3(Z)=R1 psi2(Z), psi2(x)=E[g(x-F2)|X0=x],

where F2 is the actual second substitution on the stationary OU history, and R1=int_0^1 P_t dt. Put

    a2=lambda/2-lambda^2/4,
    b2^2=lambda^2/4-lambda^3/2+5lambda^4/16,
    psi_G(x)=E_N g((1-a2)x-b2 N).

Then, uniformly in x and dimension,

    |psi2(x)-psi_G(x)|<=A epsilon(1+lambda),              (1)
    ||m3-R1 psi_G||_(L2 gamma)<=A epsilon(1+lambda).      (2)

For epsilon<=d A and D>=A^(-4), this is at most C_d A^4 sqrt(D). In particular it applies to the same-g radial-shell family at A=D^(-1/4), despite its order-one covariance-induced radial inflation. No restriction A sqrt(D)=O(1) is used.

A finite positive original-VALUE mean-law program for this target is given below. Its completed conditional law obeys

    ||W2(Law(M_G(Z)|Z),N(m3(Z),I))||_(L2 gamma)
      <=A epsilon(1+lambda)+C delta A sqrt(D)
          +Lambda A^4 sqrt(D)+restored absolute floors.          (3)

Use delta<=A^3, the fixed-order source/compiler guards, and epsilon<=d A, D>=A^(-4) for the advertised order-four result on this source subclass. The general class with no bounded linear remainder is OPEN. The certificate must hold for the ACTUAL normalized source of a call; it is not automatically preserved at the displayed grade by every conditional rescaling.

## 1. Literal analytical comparison, with full history correlations

On the stationary standard OU process, define

    H1=integral_0^infinity e^(-t)X_t dt,
    H2=integral_0^infinity t e^(-t)X_t dt.

These are comparison variables, never executed paths. Let F1=int e^(-t)g(X_t)dt and let F2 be the second substitution at the same endpoint. Because the weights are positive and have mass one,

    |F1-lambda H1|<=epsilon                         pathwise.

At each shifted history time the identical bound holds with its genuine future ancestry. Substituting once more and using g(y)=lambda y+r(y) gives the exact decomposition

    F2=lambda H1-lambda^2 H2+R2,
    R2=-lambda integral e^(-t)(F1,t-lambda H1,t)dt
                         +integral e^(-t)r(X_t-F1,t)dt,
    |R2|<=epsilon(1+lambda)                          pathwise.   (4)

The convolution kernel for the twice-integrated linear force is t e^(-t), which proves the H2 term. No independence between R2 and the Gaussian backbone is used or asserted. Applying the A-Lipschitz original g to x-F2 proves (1) directly. Conditional expectation and R1 contraction prove (2).

## 2. Exact conditional Gaussian backbone

At the original endpoint X0=x,

    E[H1|x]=x/2, E[H2|x]=x/4,
    Cov(H1|x)=I/4,
    Cov(H1,H2|x)=I/4,
    Cov(H2|x)=5I/16.                                 (5)

For a direct check, the unconditioned Laplace covariance is

    K(a,b)=integral integral e^(-at-bu)e^(-|t-u|)dtdu
          =(1/(a+b))(1/(a+1)+1/(b+1)).

At a=b=1 it gives Var(H1)=1/2, Cov(H1,H2)=3/8, Var(H2)=3/8. Subtract the endpoint regression products 1/4, 1/8 and 1/16 to obtain (5). The covariance orientations are scalar here and exact.

Therefore lambda H1-lambda^2 H2 conditioned on x is exactly

    a2 x+b2 N, N~N(0,I),

where b2 is the NONNEGATIVE known scalar square root displayed above. The polynomial 4-8lambda+5lambda^2 has negative discriminant and is strictly positive; b2=0 at lambda=0. No source-dependent covariance root is executed. The fresh N supplies the complete conditional backbone law, not a coupling to an observed old path.

## 3. Finite shifted original-VALUE source

Freeze a positive interior Hermite-multiplier rule (w_i,t_i), mass one and exact first moment 1/2, satisfying

    ||Q-R1||_(L2(gamma)->L2(gamma))<=delta,
    N_out=O(log^2(1/delta)).

Let c_i=sqrt(1-t_i^2). At retained original Gaussian carrier Z, draw independent standard D-roots G,N, each shared across its entire level. Execute

    x_i=t_i Z+c_i G,
    F_G(Z,G,N)=sum_i w_i g((1-a2)x_i-b2 N).           (6)

This is the complete raw source: N_out original g VALUES and no hidden original ancestors. Every changed compiler argument replays all of (6). Its own conditional mean is exactly Q psi_G(Z). Since every x_i is standard unconditionally,

    ||psi_G||2<=A sqrt((1-a2)^2+b2^2) sqrt(D)
                         <=C A sqrt(D),
    ||Q psi_G-R1 psi_G||2<=C delta A sqrt(D).         (7)

All nonlinear shifts in (6) remain inside the original terminal g. There is no Taylor expansion in covariance or in the Euclidean displacement. The full Gaussian covariance inflation is retained.

## 4. Source admission, raw zero, callers, and positivity

Use the literal same-G baseline/correction split

    B(Z,G)=sum_i w_i g(x_i), E_G=F_G-B.

B is a genuine gradient in G, with first at most A. The correction has

    ||E_G||_(Lp|Z)<=C_p A(|a2|+b2)(|Z|+sqrt(D))
                        <=C_p A^2(|Z|+sqrt(D)),
    complete private first(E_G)<=C A,
    captured-Z first(E_G)<=C A,
    curl(P_G* E_G)<=C A b2<=C A^2.                 (8)

Indeed the G derivative is

    sum_i w_i c_i[(1-a2)Dg((1-a2)x_i-b2 N)-Dg(x_i)],

which is symmetric even if its difference is O(A). The only off-diagonal private block in the P_G lift is

    D_N E_G=-b2 sum_i w_i Dg((1-a2)x_i-b2 N),

whose operator norm is at most A b2. This proves the curl bound without a Hessian modulus or a differentiated HVP.

For explicit sufficient guards, let beta=sum_i w_i c_i and take A<=1/2. Then 0<=1-a2<=1, and the difference of the two displayed positive-semidefinite Hessians has operator norm at most A. Consequently

    G first(E_G)<=A beta,
    N first(E_G)<=A b2,
    full private first(E_G)<=A sqrt(beta^2+b2^2),
    captured-Z first(E_G)<=A/2,
    |E_G(Z,0,0)|<=A a2 |Z|/2.

For half variance shares, sufficient normalized radius bounds are rho_B=sqrt(2)A beta and ell_E=sqrt(2)A sqrt(beta^2+b2^2)<=1/4, together with the imported gradient-order/public-log and curl guards. This is an explicit sufficient choice, not a relaxation of those native guards.

At fixed nonzero Z the origins B(Z,0) and E_G(Z,0,0) are actual caller-only VALUE executions and are captured/restored under the same frozen version. At total Z=G=N=0 the source zero is literal from g(0)=0; b2 N and every x_i vanish. lambda=0 causes no division and no degeneracy problem.

Invoke the already admitted gradient/near-gradient OWN-mean compiler with independent COMPLETE banks for B and E_G, fixed positive shares summing to one, gradient order at least four and padding mu=A. Impose its literal actual normalized first/radius, small-curl, active-dimension, caller, quadrature and precision guards. By (8) the conditional target is N(E F_G|Z,I), with error Lambda A^4(|Z|+sqrt(D)) plus absolute floors. Conditional product coupling applies only after both complete law returns; no output is read as a strong mean, and no executed Gaussian carrier is subtracted.

All operations are finite real Gaussian pushforwards. Negative scalar arguments and the arithmetic baseline difference are ordinary deterministic operations, not signed probabilities. The final completed mean has its imported literal known source-zero unit Gaussian row.

## 5. Full root and work ledger

One raw F_G occurrence owns 2D private Gaussian coordinates; one B occurrence owns D. The completed services own every additional filter, pair, mark, mean-bank, clock and fill root required by their pinned compiler. Their COMPLETE banks are independent conditional on the same Z and actual exterior labels. The raw G,N are integrated in the output law and may not later be appended as observers.

A safe bill is

    Q_captured+N_out N_B+2 N_out N_E
                            +Q_known/numerical/replay,           (9)

where N_B and N_E are the FULLY EXPANDED original-source occurrence counts of those admitted fixed-order mean services at their actual normalized radii. No raw source leaf is counted as an uncharged function provider. If the baseline value in E cannot be cached at an identical complete key, its call is already included in the factor two. Caller-only origins each cost N_out for B and 2N_out for E. Known scalar a2,b2, quadrature coefficients, and Gaussian/arithmetic work are separately charged.

Requested first/adjoint sweeps use original HVPs only at the recorded original VALUE sites in (6) and in the baseline. A discarded primal pays complete replay. No saved HVP is differentiated. All lambda,A,epsilon bounds, scalar-root versions, shares, schedules, padding and numerical tolerances are fixed before differentiation. Original or conditional finite-mode restorations retain the imported actual caller chain rules and absolute numerical floors; they are not inferred by differentiating (3).

At the fixed grade, (9) is a public-log polynomial. The inverse-A original VALUE exponent is zero on this declared subclass; the bound epsilon and the condition D>=A^(-4) are hypotheses, not computable tests that the program infers from finite observations.

## 6. Scope

This gives a concrete positive resummed VALUE consumer for the canonical m3 target under a bounded linear-remainder certificate. It is particularly useful as a countertest companion: the high-dimensional radial family defeats a finite centered covariance truncation but is handled by this full Gaussian backbone with its original terminal nonlinearity intact.

It does not establish the general centered/resummed current theorem. In a general C2 source, g-lambda Id need not have bounded amplitude epsilon=O(A), and (4) then supplies no required absolute bound. For example, under conditional normalization f(y)=s[g(a+s y)-g(a)], the new parameters are lambda_f=lambda s^2 and a safe epsilon_f=2s epsilon. These actual parameters must be inserted into (3); one cannot assume epsilon_f=O(alpha) or D>=alpha^(-4), alpha=A s^2, from the original assumptions. No reverse-OU endpoint join is supplied by this note.

No generic one-energy trace estimate, generic conditional Gaussianization, or all-order recurrence is claimed.
