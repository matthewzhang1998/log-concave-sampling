# Gaussian-RMS residuals extend the finite resummed m3 consumer

2026-10-05. Source-qualified component theorem. The general m3 gate and fourth-order endpoint remain open.

## Result and strict enlargement

The global bounded-remainder hypothesis can be replaced by two residual norms under **explicit centered Gaussian laws**. Their use introduces no history-expectation oracle and no extra original VALUE leaf. The comparison is now integrated over the original Gaussian endpoint; it is not uniform in that endpoint.

Let

    g=grad U, g(0)=0, 0<=Dg<=A I, 0<A<=1/2,
    0<=lambda<=A, r(x)=g(x)-lambda x,
    Lr=max(lambda,A-lambda),
    sigma1^2=1-lambda+lambda^2/2,
    epsilon0=||r(Z)||_(L2),
    epsilon1=||r(sigma1 Z)||_(L2), Z~N(0,I_D).

The scalar lambda is an actual known input, fixed before fresh source roots and caller differentiation. Certified upper bounds for epsilon0,epsilon1 may replace the norms. They are analytical source qualifications, not callable expectation leaves. A computable envelope certificate and a completely explicit choice of lambda are supplied below.

For the canonical Gaussian-history target

    m3=R1 psi2,
    psi2(x)=E[g(x-F2)|X0=x],

define the same scalar backbone as in the earlier bounded theorem:

    a2=lambda/2-lambda^2/4,
    b2^2=lambda^2/4-lambda^3/2+5lambda^4/16,
    psi_G(x)=E_N g((1-a2)x-b2 N).

Then

    ||psi2-psi_G||_(L2(gamma))
       <= A[epsilon1+(lambda+Lr)epsilon0],                 (1)
    ||m3-R1 psi_G||_(L2(gamma))
       <= A[epsilon1+(lambda+Lr)epsilon0].                 (2)

The previously admitted finite positive original-VALUE program, with all of its literal fixed-order guards, therefore satisfies

    ||W2(Law(M_G(Z)|Z),N(m3(Z),I))||_(L2(gamma_Z))
      <= A[epsilon1+(lambda+Lr)epsilon0]
           +C delta A sqrt(D)+Lambda A^4 sqrt(D)
           +restored absolute floors.                    (3)

In particular the sufficient source qualification is

    epsilon1+(lambda+Lr)epsilon0 <= C_res A^3 sqrt(D),      (4)

with delta<=A^3 and all native compiler guards retained. No upper bound on A sqrt(D) is imposed. The two Gaussian norms incur no density-ratio factor exponential in D.

This is genuinely larger than the prior useful bounded-remainder class. Section 4 gives a C2 convex-gradient source for every D>=1 with an unbounded residual around its bulk scalar slope, and with epsilon0,epsilon1<=A^3/sqrt(2). Every scalar slope yielding a globally bounded residual has a bound at least of order A sqrt(D), which is too large for the previous fourth-order certificate by a factor at least of order A^(-2).

Conversely, Section 6 gives an exact two-eigenvalue quadratic obstruction: no scalar resummed backbone can approximate all admissible convex gradients at the desired order. This extension does not claim general source closure.

## 1. Genuine future ancestry and the two explicit Gaussian laws

Work on a stationary standard OU history with covariance exp(-|t-u|)I. At every actual future time t define

    F1,t=integral_0^infinity e^(-u)g(X_(t+u))du,
    H1,t=integral_0^infinity e^(-u)X_(t+u)du,
    Delta1,t=F1,t-lambda H1,t.

These remain analytical comparison variables and are not producer leaves. Since r(0)=0 and r is Lr-Lipschitz, all quantities are square integrable. The Lipschitz constant follows from

    -lambda I <= Dr=Dg-lambda I <= (A-lambda)I.

Stationarity and Minkowski, with positive mass-one exponential weights, give

    ||Delta1,t||2
      =||integral e^(-u)r(X_(t+u))du||2 <= epsilon0.       (5)

The *same-history* Gaussian comparison variable is

    Y_t=X_t-lambda H1,t.

Unconditionally, Cov(X_t,H1,t)=I/2 and Cov(H1,t)=I/2, so

    Y_t ~ N(0,(1-lambda+lambda^2/2)I).

Consequently ||r(Y_t)||2=epsilon1, exactly. This is a known Gaussian marginal, not a claim that Y_t is independent of Delta1,t or of the endpoint. The Lipschitz inequality on their actual coupling yields

    ||r(X_t-F1,t)||2
      <= ||r(Y_t)||2+Lr||Delta1,t||2
      <= epsilon1+Lr epsilon0.                            (6)

For the actual second substitution

    F2=integral e^(-t)g(X_t-F1,t)dt,

Fubini gives the exact decomposition

    F2=lambda H1-lambda^2 H2+R2,
    H1=integral e^(-t)X_t dt,
    H2=integral t e^(-t)X_t dt,
    R2=-lambda integral e^(-t)Delta1,t dt
          +integral e^(-t)r(X_t-F1,t)dt.

Using (5) and (6), again on the genuine history,

    ||R2||2 <= epsilon1+(lambda+Lr)epsilon0.               (7)

Every expectation in this proof is with respect to the original stationary law. No independently resampled F1 is inserted. In particular no conditional Gaussianity, small covariance, or independence of the residual is asserted.

## 2. Conditional Gaussian backbone and integrated target bias

Conditioning the linear comparison on X0=x gives exactly

    lambda H1-lambda^2 H2 | X0=x  =_law a2 x+b2 N.

This uses the already proved covariance block

    Cov((H1,H2)|X0)= [[1/4,1/4],[1/4,5/16]] tensor I,

and regressions E[H1|X0=x]=x/2, E[H2|X0=x]=x/4. Hence, on the original history, conditional Jensen and the A-Lipschitz property imply

    E | E[g(X0-F2)-g(X0-lambda H1+lambda^2 H2)|X0] |^2
       <= A^2 E|R2|^2.

Only after taking this conditional expectation do we replace the conditional *law* of the linear backbone by the fresh Gaussian N. This proves (1). The L2(gamma) contraction of R1 proves (2).

Unlike the bounded-remainder proof, (7) does not bound R2 at every fixed endpoint. No uniform conditional version of (1), derivative estimate obtained from (1), or retained-private-root coupling is claimed.

The older bounded class is still included at its stated grade: if |r|<=epsilon, then epsilon0,epsilon1<=epsilon and (1) is at most A epsilon(1+lambda+Lr), with lambda+Lr<=2A. The older pathwise estimate A epsilon(1+lambda) remains a sharper available bound when that stronger certificate is known. For epsilon=O(A), D>=A^(-4), both are fourth-order at the original standard carrier.

## 3. Computable certificates, not hidden expectations

The theorem does **not** assert that epsilon0 or epsilon1 can be discovered at negligible cost from arbitrary black-box g VALUES. The executing mean program never evaluates either norm. One can certify (4) from ordinary explicit source structure.

A useful general envelope certificate is:

    |r(x)|<=e_core for |x|<=R,
    r is Lr-Lipschitz globally.

Radially project x onto the ball and apply the Lipschitz inequality to obtain

    |r(x)|<=e_core+Lr(|x|-R)_+.

For either known scale s=1 or s=sigma1,

    ||r(sZ)||2 <= e_core+Lr sqrt(J_D(s,R)),
    J_D(s,R)=E[(s|Z|-R)_+^2].                             (8)

This is an explicitly specified one-dimensional chi integral:

    J_D(s,R)=1/[2^(D/2-1) Gamma(D/2)]
               integral_(R/s)^infinity
                  (s q-R)^2 q^(D-1) exp(-q^2/2)dq.

For s=0 its value is zero when R>=0. Here sigma1>0, so this convention is only for general use. Formula (8), a certified upper bound for the integral, or a directly proved radial envelope is a sufficient source certificate. A numerical upper bound, when used, must include its actual quadrature and tail error. No call to an unknown nonlinear or conditional expectation is substituted for that certification.

For example, if R=sqrt(D)+T and 0<s<=1,

    J_D(s,R)<=J_D(1,R)<=2 exp(-T^2/2).                    (9)

A short direct proof avoids any dimension-asymptotic approximation. Chernoff's bound for |Z|^2 gives, for t>=0,

    P(|Z|>=sqrt(D)+t)
      <=exp(-(x-D-D log(x/D))/2), x=(sqrt(D)+t)^2
      <=exp(-t^2/2),

because u-log(1+u)>=0 at u=t/sqrt(D). Integrating the survival function,

    E[(|Z|-R)_+^2]
       <=2 integral_0^infinity v exp(-(T+v)^2/2)dv
       <=2 exp(-T^2/2).

The constants are uniform in every integer D>=1. They are not derived by dominating one D-dimensional Gaussian density by another.

The two residual norms must generally remain distinct. The inequality sigma1<=1 does not imply ||r(sigma1 Z)||2<=||r(Z)||2 for arbitrary admissible r: the residual magnitude need not increase along rays. It does imply monotonicity for the radial hinge *envelope* in (8), and for the particular radial residual in Section 4. No generic replacement of epsilon1 by epsilon0 is used.

Lambda selection is explicit in the example below. For a generic source, choosing lambda by an unknown Gaussian regression is not a finite algorithm merely because the regression has a formula. Any empirical selection/validation adds all training VALUE calls and Gaussian roots to setup costs and requires its own error certificate. Those costs have not been concealed in the fixed-order consumer bill.

## 4. Explicit smooth tail-modified source beyond every useful global-sup certificate

Take D>=1, 0<A<=1/2, lambda=d=A/2, and set

    T=sqrt(8 log(1/A)), R=sqrt(D)+T.

Use the explicit monotone C2 step

    psi(u)=0                            if u<=0,
           10u^3-15u^4+6u^5            if 0<u<1,
           1                            if u>=1.

Its derivative is 30u^2(1-u)^2 on (0,1), and 0 elsewhere. Let

    h(s)=0                              if s<=0,
         (5/2)s^4-3s^5+s^6              if 0<s<1,
         s-1/2                          if s>=1.

Then h'=psi, h is C3, and 0<=h(s)<=s_+. Define the actual source

    g(x)=lambda x+d h(|x|-R) x/|x|     for x!=0,
    g(0)=0.

The nonlinear term vanishes on a neighborhood of the origin, so the radial expression is globally C2 (in fact C3). It is the gradient of

    U(x)=lambda |x|^2/2+d integral_0^|x| h(q-R)dq.

At radius rho>0 its radial and tangential Jacobian eigenvalues are respectively

    lambda+d psi(rho-R),
    lambda+d h(rho-R)/rho,

both in [lambda,A]. In D=1 only the radial eigenvalue is present. Thus the source is genuinely in the original convex-gradient class.

For the selected bulk slope lambda=A/2 the residual satisfies

    |r(x)|=d h(|x|-R)<=d(|x|-R)_+,

and is globally unbounded. Since sigma1<=1, (9) gives

    epsilon0,epsilon1 <=d sqrt(2)exp(-T^2/4)
                     = A^3/sqrt(2).

Here Lr=A/2, so (1) gives the entirely explicit bound

    ||m3-R1 psi_G||2 <= (1+A)A^4/sqrt(2).                 (10)

It is fourth-order against A^4 sqrt(D) for **every** D>=1, with no D>=A^(-4) condition.

This is not merely a different decomposition of a source already covered at that grade by the old global-sup criterion. For rho>=R+1,

    g(rho n)=A rho n-d(R+1/2)n.

For every scalar mu!=A, sup_x|g(x)-mu x| is infinite. For mu=A, rho-h(rho-R) is nondecreasing from zero to R+1/2, hence the exact global remainder bound is

    sup_x|g(x)-A x|=d(R+1/2).

The old bound at this sole admissible globally bounded slope is therefore

    A d(R+1/2)(1+A),

whose ratio to A^4 sqrt(D) is

    (1+A)(R+1/2)/(2A^2 sqrt(D)) >= (1+A)/(2A^2).

Thus no fixed public-polylog loss converts that old certificate into the desired order as A tends to zero. This is a strict enlargement of the order-four-certified source class. It is not a lower bound on the actual target error when using slope A; only the old global-sup certification is being compared here.

## 5. The finite positive consumer, caller graph, ancestry and costs

The executing program is unchanged from the independently audited bounded construction. Freeze a positive interior rule (w_i,t_i) of mass one, first moment 1/2, and Hermite operator error delta. Let c_i=sqrt(1-t_i^2). At the retained original standard Gaussian Z, draw independent standard D-roots G,N, each shared across its full level, and execute

    x_i=t_i Z+c_i G,
    F_G=sum_i w_i g((1-a2)x_i-b2 N),
    B=sum_i w_i g(x_i), E=F_G-B.

There are no canonical F1/F2 VALUE ancestors behind the known coefficients. One F_G occurrence uses N_out original VALUES and owns 2D private Gaussian coordinates. A B occurrence uses N_out VALUES and D coordinates; a raw E occurrence costs at most 2N_out VALUES and 2D coordinates.

Its own mean is exactly Q psi_G(Z), and

    ||Q psi_G-R1 psi_G||2<=C delta A sqrt(D).

The boundedness of r was not used for any of the following ports, so they remain unchanged:

    ||E||_(Lp|Z=z)<=C_p A^2(|z|+sqrt(D)),
    beta=sum_i w_i c_i,
    ||D_G E||<=A beta,
    ||D_N E||<=A b2,
    ||D_(G,N) E||<=A sqrt(beta^2+b2^2),
    ||D_Z E||<=A/2,
    curl(P_G^*E)<=A b2<=A^2/2.

Indeed D_G E is the sum of c_i[(1-a2)Dg(shifted_i)-Dg(x_i)], a symmetric difference of matrices in [0,AI]. The N off-diagonal block is -b2 sum_i w_i Dg(shifted_i). No Hessian modulus or differentiated HVP is used.

The fixed-Z origins are actual caller executions

    B0(Z)=sum_i w_i g(t_i Z),
    E0(Z)=sum_i w_i[g((1-a2)t_i Z)-g(t_i Z)].

They cost N_out and 2N_out VALUES before exact-key sharing. Their bounds are A|Z|/2 and A a2|Z|/2. Anchoring/restoring these values retains their actual caller first paths, and does not pretend they vanish at nonzero Z. Only the complete zero Z=G=N=0 is literal from g(0)=0.

Use the admitted fixed-order gradient and near-gradient mean consumers with independent COMPLETE banks conditional on the same retained Z and actual exterior labels, fixed positive shares summing to one, gradient order at least four, and padding mu=A. For half shares the actual sufficient normalized radii remain

    rho_B=sqrt(2)A beta,
    ell_E=sqrt(2)A sqrt(beta^2+b2^2).

Impose the actual gradient-radius threshold, ell_E<=1/4, curl, active-dimension, finite-clock/filter, numerical, and exterior caller guards of those imports. A<=1/2 alone does not discharge them. The near-gradient error is fourth-order because its energy is O(A^2)(|Z|+sqrt(D)) and its fixed-grade error bracket is O(A^2). The completed conditional targets add to N(Q psi_G(Z),I), by ordinary product coupling only after the entire law returns.

No raw G,N, training root, or executed source force is subsequently appended as an observer. The known source-zero Gaussian row is retained as required by the mean-law program; no output is read as a strong conditional mean or covariance statistic. All probabilities remain positive Gaussian pushforwards; the arithmetic subtraction in E is not a signed measure.

The complete original-VALUE bill remains

    Q_certificate_setup+Q_captured
       +N_out N_B+2N_out N_E+Q_known/numerical/replay.       (11)

Here N_B,N_E are the fully expanded native occurrence counts at the actual normalized radii, including every filter, clock, mean-bank, and source replay. The sample-free explicit envelope example has Q_certificate_setup=0 original VALUES: A,D,lambda,R and the displayed analytical bound are known arithmetic. Any alternative fitting or empirical certification pays its full additional Q_certificate_setup; it is not a free oracle. N_out=O(log^2(1/delta)), so at fixed grade and declared/certified parameters this preserves the existing inverse-A original VALUE exponent zero. This is not a dimension-independent arithmetic claim.

All additional Gaussian roots, fills, filters, pairs, marks and mean banks remain charged by the imported compilers. Original HVPs occur only at the recorded original VALUE sites in requested first/adjoint sweeps, with complete primal replay when necessary. Frozen lambda, coefficient/root versions, source certificate, shares, guards, modes, tolerances and schedules are not differentiated as if they were live estimated expectations. If a caller actually changes or computes a backbone parameter, its dependence must instead be represented and charged by a new caller graph; this note supplies no license to stop its derivative.

Absolute VALUE, clock, scalar-root, finite-mode, and replay floors are exactly those in the bounded audit and are restored additively. They are not divided by an RMS residual that may be zero. The inequalities (1)-(4) establish only target bias; caller/first guarantees come from the literal finite program above.

## 6. Exact obstruction to automatic general scalar-backbone closure

There is a simple original-source obstruction even without any nonlinear shell. Let D be even and

    K=diag((A/3)I_(D/2),(2A/3)I_(D/2)), g(x)=Kx.

This is anchored, smooth, convex, and 0<=Dg<=AI. Define a(k)=k/2-k^2/4. The genuine canonical mean is exactly

    m3(Z)=(1/2)K[I-a(K)]Z.

The fully resummed scalar-backbone mean at any chosen lambda in [0,A] is exactly

    m_G,lambda(Z)=(1/2)K[1-a(lambda)]Z.

The b2N term has zero expectation because the terminal g is linear. The positive outer rule with exact first moment 1/2 has no target quadrature error for these linear means.

For k1=A/3,k2=2A/3, minimize the exact squared norm over c=a(lambda):

    ||m_G,lambda-m3||2^2
       =(D/8)[k1^2(a(k1)-c)^2+k2^2(a(k2)-c)^2].

Its minimizer is c*=[k1^2 a(k1)+k2^2 a(k2)]/(k1^2+k2^2). Since a is increasing on [0,A] and c* lies between a(k1),a(k2), this minimizer is attained by a lambda in [k1,k2]. Therefore the exact optimum is

    inf_(0<=lambda<=A) ||m_G,lambda-m3||2
       =sqrt(D/8) [k1 k2/sqrt(k1^2+k2^2)] |a(k2)-a(k1)|
       = A^2(2-A)sqrt(D)/(36 sqrt(10)).                   (12)

This exceeds every fixed public-polylog multiple of A^4 sqrt(D) along A tending to zero. Even an oracle choosing the *best* scalar lambda cannot remove it. Adding a fixed affine intercept to the scalar comparison adds only a constant mean term, orthogonal in L2(gamma) to this centered linear error, and cannot lower (12).

The result refutes universal closure by this **scalar** Gaussian backbone; it does not refute matrix or endpoint-dependent conditional backbones, a nonlinear correction to its target, or a different finite current consumer. On this example a known matrix backbone K is exact, but declaring or recovering such a matrix for general g, implementing its actual caller dependence, and proving an adequate residual certificate are additional obligations. They are not delivered by finite scalar moment fitting or by the existence of Gaussian L2 norms.

## 7. Conditional rescaling boundary

For a genuinely normalized conditional source

    f(y)=s[g(a+s y)-g(a)], alpha=A s^2,

the corresponding scalar slope is lambda_f=lambda s^2, and its residual is

    r_f(y)=s[r(a+s y)-r(a)].

The required norms for this source are centered at the actual translated arguments a+sZ and a+s sigma1,f Z:

    epsilon0,f=s||r(a+sZ)-r(a)||2,
    epsilon1,f=s||r(a+s sigma1,f Z)-r(a)||2,
    sigma1,f^2=1-lambda_f+lambda_f^2/2.

The unshifted norms epsilon0,epsilon1 at the original source do not control these uniformly in a. The tail-modified example makes the danger concrete: a caller located in the transition/tail need not see a small residual around the original bulk slope. A conditional application needs these actual translated Gaussian certificates and its actual caller ports. No endpoint join or uniform conditional guarantee is inferred from (10).

## 8. Provenance and exact scope

The incoming continuation manifest is

    ../RESUMMED-MEAN-CONTINUATION-MANIFEST.json
    SHA256 c4cd4c154c7f8ec98a49c72953a1f99faf428232acb079f29bf880a23401a1bd.

The unchanged source/consumer theorem is

    ../resummed-linear-backbone/BOUNDED-REMAINDER-RESUMMED-M3-VALUE-CONSUMER.md
    SHA256 798d20b87124e7069e465ce4e7bcf92c5af3734d86c8a62aeb288b5fec8eaef3,

and its independent audit is

    ../resummed-linear-backbone/independent-audit/INDEPENDENT-BOUNDED-RESUMMED-M3-AUDIT.md
    SHA256 49d135ba6692c54012acd9d0d94f603c0015fca98b5d1be1164457ff45d1213a.

That audit pins the positive Hermite-multiplier rule, completed gradient/near-gradient means, coisometry adapter, complete caller origins, original VALUE replays, and source-zero Gaussian row. Those fixed-order services are imported here, not reconstructed numerically. The new proof work is the dimension-safe Gaussian-RMS comparison, computable certificate/example, and exact scalar-backbone obstruction.

Established: a strictly larger, explicitly certifiable, guarded positive m3 mean subclass using the existing finite source and costs.

Not established: a general C2 convex-gradient m3 consumer; automatic data-driven backbone selection; uniform conditional rescaling; a reverse-OU endpoint join; arbitrary-order recurrence; or c(P)/P tending to zero.
