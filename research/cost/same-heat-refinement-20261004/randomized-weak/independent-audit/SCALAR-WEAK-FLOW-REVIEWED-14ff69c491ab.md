# Scalar randomized quarter-flow: weak gain and all-layer debts

2026-10-04. New research. The scalar terminal conclusions below are not a new admitted general-dimensional source family or a change to the current c_P recurrence.

## Result and scope

An explicit finite positive scalar scheme has a real weak-law improvement over pathwise quadrature. Stratified Gaussian-generated random clocks are independent between Picard layers. Every old gradient value is stored and reused, and the incumbent is called once at the original heat. A variance-preserving split of the fresh momentum supplies a genuinely unread final Gaussian keep.

For fixed depth M, final-layer count N_M, earlier counts N_j, keep amplitude sigma, and separately paid numerical/cap errors, its scalar physical W2 ledger is

    C_M [a error_old + a^(M+3/2) + a^(3/2) sigma^2
         + sum_(j=1)^(M-1) a^(M-j+3/2) N_j^(-3/2)
         + a^(5/2)/(sigma N_M^3)].                       (A)

The sigma^2 keep term is a posterior-law statement proved below. The easier same-record coupling has only sigma a^(3/2) error. The first-layer/last-layer rank-one quadrature term cancels exactly, conditional on the complete preceding record, provided no exterior observer reads that final clock bank.

For target R>=5/2, take sigma=a^((R-3/2)/2) and M+3/2>=R. The scalar original-query exponent from (A) is

    b_R=max(0, (2R-5)/3, R/2-13/12).                    (B)

The earlier-layer term is absent for M=1. At R=5/2, one Picard layer needs N=O(a^(-1/6)), up to public logs. At R=7/2, two layers need N_1,N_2=O(a^(-2/3)). Above R=7/2, the existing bound for stochastic nonlinear feedback dominates: b_R=(2R-5)/3. Thus this does not yield o(R).

There are important additional limitations:

1. The simple conditional weak proof has a path-magnitude fourth-moment factor. In dimension d, it is generically of order d, not sqrt(d). A quadratic random-scale fixture below prevents simply changing that factor to sqrt(d).
2. Clock derivatives must be included. They grow with the path magnitude. An explicit known smooth cap gives a scalar global-first/tail remedy, but it is a new executed guard, not a free inherited primitive.
3. The late keep is genuine, with variance a sigma^2. Its shrinking width and every later inverse-width use remain charged. A complete native protected/retained-host return is not proved here.
4. The old deterministic bump bound remains confined to worst-case deterministic **first-Picard integration at fixed (z0,Z)=(0,1)**. It is not a lower bound on exact nonlinear flow, Gaussian-input Lp approximation, randomized schemes, or posterior-law compilers. Adaptive gradient/HVP queries are covered there only after the entire baseline transcript is fixed.

## 1. Normalized scalar setting

Let m be the posterior mode. For notation first use its exact value; the actual finite original-gradient mode version is restored at its separately assigned VALUE and caller budget. Put

    z=(x-m)/sqrt(a),
    U(z)=V(m+sqrt(a)z)-V(m)-sqrt(a) V'(m) z,
    b(z)=U'(z)=sqrt(a)[V'(m+sqrt(a)z)-V'(m)].

Then U(0)=b(0)=0, 0<=U''<=a, |b(z)|<=a|z|, and the normalized posterior is

    pi_a(dz)=Z_1^(-1) exp(-z^2/2-U(z)) dz/sqrt(2 pi).

Assume a is below a fixed small guard. Constants may depend on the fixed target/depth and the admitted scalar moments, with public logarithms suppressed. No third derivative of V is assumed or evaluated.

The exact normalized quarter-flow is z''=-z-b(z), at T=pi/2. It preserves pi_a after a fresh standard Gaussian momentum. Its endpoint is

    F(z0,p0)=p0-integral_0^T cos(s) b(z(s)) ds,

and it is C a-Lipschitz in z0. These facts were independently audited in `../independent-audit/INDEPENDENT-SAME-HEAT-FLOW-AND-QUERY-AUDIT.md`.

## 2. The unread keep is variance preserving

Fix 0<sigma<=1/2 and kappa=sqrt(1-sigma^2). Sample independent standards Z,K after the incumbent caller is fixed. All force paths, caps, quadrature nodes and earlier layers use only kappa Z. Emit the endpoint

    Y_sigma=F(z0,kappa Z)+sigma K.                      (K)

K is absent from every force query and all preceding computations. The true comparison uses the exact flow driven by

    Z_bar=kappa Z+sigma K,

which is standard and independent of z0. The direct momentum difference cancels exactly:

    Y_sigma-F(z0,Z_bar)
      = integral_0^T cos(s)[b(z_full(s))-b(z_omit(s))] ds.

The path difference is at most sigma |K|/(1-a). Hence the normalized strong difference is at most C a sigma |K|, or C a^(3/2) sigma physically. There is no generic sigma sqrt(a) addition error in this comparison. Replacing Z by an independent force-side record would destroy this identity.

The same O(a) starting-position contraction holds for (K), because its added K is unchanged in that comparison. Incumbent posterior error therefore still costs C a error_old.

## 3. Scalar posterior-law keep lemma: a sigma^2 gain

**Claim.** If z0~pi_a and Z,K are independent standards, then

    W2(Law Y_sigma,pi_a)<=C a sigma^2.                 (WK)

This is normalized. Its physical cost is C a^(3/2) sigma^2. The proof uses only the exact invariant phase density, first flow bounds and U in C2 with 0<=U''<=a.

Write delta=sigma^2, kappa^2=1-delta. Under the exact flow with standard initial momentum, the final pair (z,p) has law pi_a times gamma. Let p0(z,p) be its initial momentum under the inverse quarter-flow. The harmonic inverse has p0=z, and the same Duhamel estimate gives

    |p0(z,p)-z|<=C a(|z|+|p|).

Changing only the initial momentum to N(0,kappa^2) multiplies the invariant phase density by

    w_delta(p0)=kappa^(-1) exp[-delta p0^2/(2 kappa^2)].

Consequently the position law mu_delta=Law F(z0,kappa Z) has the exact density

    d mu_delta/d pi_a(z)=r_delta(z)
      =E_(p~gamma) w_delta(p0(z,p)).                   (D1)

Set g_delta(z)=kappa^(-1) exp[-delta z^2/(2 kappa^2)] and M_delta=E_pi g_delta. For delta<=1/4,

    |r_delta(z)-g_delta(z)|<=C a delta(1+z^2),
    |M_delta-1|<=C a delta.

Indeed the exponential has absolute derivative at most delta/(2 kappa^2) as a function of its nonnegative square argument, and |p0^2-z^2|<=C a(z^2+p^2). Let nu_delta have density pi_a g_delta/M_delta. Since 1/g_delta grows at most exp(z^2/6), Gaussian moments imply

    chi2(mu_delta || nu_delta)<=C a^2 delta^2.

The density nu_delta is proportional to exp[-z^2/(2 kappa^2)-U(z)] and is uniformly strongly log-concave. Its transport-entropy inequality therefore gives W2(mu_delta,nu_delta)<=C a delta. Convolution by the same sigma K cannot increase W2.

It remains to compare nu_delta convolved with N(0,delta) to pi_a. Their densities relative to gamma are respectively

    A_delta(z)/Z_kappa,  exp[-U(z)]/Z_1,
    A_delta(z)=E_G exp[-U(kappa^2 z+kappa sqrt(delta) G)],
    Z_kappa=E_(X~N(0,kappa^2)) exp[-U(X)].              (D2)

For h=exp(-U), |h'|<=a|z| and |h''|<=a+a^2 z^2. Taylor's integral formula with displacement -delta z+kappa sqrt(delta)G gives

    |A_delta(z)-exp[-U(z)]|<=C a delta(1+z^4),
    |Z_kappa-Z_1|<=C a delta.

The chi-square distance of these two position laws, relative to pi_a, is at most C a^2 delta^2, since exp(U(z))<=exp(a z^2/2). The transport-entropy inequality for pi_a proves the remaining C a delta bound and hence (WK).

This proof is scalar as stated. Its direct multidimensional polynomial-moment estimate does not have the required one-energy sqrt(d) constant. It does not bypass that issue.

## 4. Finite stratified-clock DAG

Use independent standard Gaussian clocks G_(k,j), with

    S_(k,j)=(j-1+Phi(G_(k,j))) T/N_k,  1<=j<=N_k.

Thus one node is uniform in each time stratum, and different layers have independent complete clock banks. Clock sampling and Phi precision are actual work.

Let C_B be an explicitly executed known C1 cap, equal to the identity on [-B,B], bounded by 2B, with first bounded by a fixed constant. Its zero and caller paths are retained. Start with

    q^[0](t)=cos(t) C_B(z0)+kappa sin(t) C_B(Z).

For previous layers use

    f_(k,j)=b(q^[k-1](S_(k,j))),
    q^[k](t)=q^[0](t)
       -(T/N_k) sum_j H_epsilon(t-S_(k,j)) f_(k,j).    (P)

Here H_epsilon is a known nonnegative C1 smoothing of 1_(r>=0) sin(r). It equals that function outside [-epsilon,epsilon], has uniformly bounded first, and differs by O(epsilon) only in that interval. A fixed cubic Hermite patch supplies such a finite implementation at small epsilon. Its integrated error is O(epsilon^2). It is another explicit known filter, not a derivative of V. Picard layers remain acyclic even if the smoothing reads a nearby future time from the preceding layer.

Each f_(k,j) is queried once and stored. Evaluating the preceding path at a new node reads these stored values and known kernel coefficients; it does not replay the preceding original-gradient program. The final emitted normalized value is

    Y_hat=kappa Z+sigma K
       -(T/N_M) sum_j cos(S_(M,j)) b(q^[M-1](S_(M,j))). (E)

K is sampled after, or held unread until after, every force and clock-dependent path value. The known leading carrier is kappa Z+sigma K and has variance exactly one.

Conditional on the complete earlier record R=(z0,Z,all earlier clock banks), the exact comparison center is

    B_R=kappa Z-integral_0^T cos(s)b(q^[M-1](s)) ds.

This integral is an analytical reference, never a queried mean. The final quadrature residual D=Y_hat-B_R-sigma K satisfies E[D|R]=0 exactly. Its clocks are fresh conditional on R; K is independent of both.

## 5. Exact local current and scalar buffered bound

Let h=T/N_M and F_R(s)=cos(s)b(q^[M-1](s)). Then

    D=-h sum_j [F_R(S_j)-E(F_R(S_j)|R)],
    Sigma(R)=E[D^2|R]=h^2 sum_j Var(F_R(S_j)|R).        (V)

For smooth scalar tests phi, the exact conditional law difference is

    E[phi(B_R+D+sigma K)-phi(B_R+sigma K)|R]
      =integral_0^1 (1-t) E[D^2 phi''(B_R+tD+sigma K)|R] dt. (J)

Its leading same-endpoint rank-two current is Sigma(R)/2 times E_K phi''(B_R+sigma K). The remainder is the rank-three current

    (1/2) integral_0^1 (1-t)^2
      E[D^3 phi'''(B_R+tD+sigma K)|R] dt.

Only derivatives of the test appear. No original third derivative is evaluated or assumed. Unknown Sigma is an analysis coefficient, not an executed covariance oracle.

The path and its time first obey |q|+|q'|<=C_M R_* uniformly over earlier clocks, where R_* may be taken as a constant times |z0|+|Z| after capping. Thus |F_R'|<=C_M a R_*. Within each stratum, the Gaussian clock derivative is bounded by C a R_* h^2. Summing its squared columns gives

    Lip_(G_M) D <= C_M a R_* N_M^(-3/2),
    ||D||_(Lp|R) <= C_(M,p) a R_* N_M^(-3/2).         (S)

The independent-buffer/Riesz inequality already proved for the weak-E mean applies conditional on R:

    W2(Law(B_R+D+sigma K)|R, N(B_R,sigma^2))
       <= Lip(D) ||D||_(L2|R)/(2 sigma).

Mixing these conditional couplings and using the scalar fourth moment of R_* proves normalized error C_M a^2/(sigma N_M^3), or physical error C_M a^(5/2)/(sigma N_M^3).

If an exterior observer reads the final clock bank, E[D|R,observer] need not vanish. Then the first-order term E[D phi'(B_R+sigma K)|R,observer] returns. Formula (J) applies only when every observer outside this packet ignores that bank, or when its own additional current is explicitly supplied. Earlier records and old source labels may be retained in R; the final clocks may not silently become such labels.

## 6. All preceding layers and actual source costs

For a fixed earlier path, stratified integration has centered Lp error C a R_* N_j^(-3/2). Subsequent applications of b contract path VALUE differences by C a. Therefore a noise introduced at layer j<M has the safe physical contribution

    C_M a^(M-j+3/2) N_j^(-3/2).

This bound is applied to the actual nonlinear later paths. No claim that their conditional means remain correct is made. The nearest previous layer contributes a^(5/2)N_(M-1)^(-3/2), exactly the term that eventually dominates (A).

The deterministic Picard tail is a^(M+3/2). Kernel smoothing adds C_M a^(5/2) epsilon^2, plus smaller propagated terms, so epsilon is chosen after the target budget. Finite Phi/kernel/cap/mode/gradient encodings each receive their actual amplified VALUE allowance.

The cap changes only force paths; the emitted leading kappa Z+sigma K is uncapped. Its physical error is at most

    C_M a^(3/2) [||z0-C_B(z0)||_2+||Z-C_B(Z)||_2].

If the incumbent's stated global Gaussian first and scalar concentration return are used, B can be a sufficiently large fixed-target public-log multiple, making this numerical at every requested fixed order. Alternatively a finite-moment cutoff and its explicit power cost can be used. No tail is inferred from a pointwise origin evaluation.

The clock firsts in (S) are not ignored. Without the cap they have an unbounded R_* multiplier. With the explicit cap, their global first is C a B N_j^(-3/2); firsts through earlier layers carry their additional factors a. The old-input first is O(a), the fresh Z first is kappa+O(a), and K has the exact independent coefficient sigma. A physical original-force readout therefore has the usual symmetric leading row and small remainder within these graph bounds. This is a mathematical first certificate for the stated cap/filter graph, not a completed admission of those guards into every native hidden host.

The executed original-query bill is

    Q_old(a)+sum_(j=1)^M N_j+Q_mode+Q_zero+Q_guard,

with one original HVP per queried original gradient in each directional first/adjoint sweep. Earlier gradient values are not rebuilt. Known kernel algebra can cost O(sum_j N_j N_(j-1)); it is not free in a full-arithmetic model. The actual Gaussian tape has at least sum_j N_j new clock coordinates, in addition to the old tape, Z and K. Any later prior consuming that tape pays its true dimension.

The retained graph must keep all the actual clock/filter/cap nodes unless a separate state-level deletion is proved. Old retained-state error is attenuated by O(a), but that is only one of the required host rows. The final protected variance is a sigma^2. Its compatibility with a later proxy's required protected heat, all inverse widths, and the new tape dimension remains an explicit join obligation.

## 7. Why the scalar weak estimate is not dimension safe

For a quadratic b(z)=a lambda z in d dimensions, take one Picard layer and shared scalar time nodes. Write

    A_N=h sum_j cos(S_j)^2,
    B_N=h sum_j sin(S_j) cos(S_j).

With an exact Gaussian incumbent, conditional on the clock bank the normalized output is an isotropic Gaussian with variance

    V_clock=(a lambda A_N)^2/(1+a lambda)
             +kappa^2(1-a lambda B_N)^2+sigma^2.       (G)

Here E B_N=1/2, but Var B_N is of order N^(-3). Thus V_clock fluctuates at order a N^(-3/2), before any native conditional covariance repair. This is a real positive Gaussian scale mixture.

If chi_d is a standard Gaussian radius, its output radius has variance

    E[V_clock] Var(chi_d)+(E chi_d)^2 Var(sqrt(V_clock)).

For any fixed isotropic Gaussian target, radial contraction gives the W2 lower bound

    [E chi_d sd(sqrt(V_clock))-sqrt(v_target Var(chi_d))]_+.

When d is large, this detects the order sqrt(d) a N^(-3/2) scale fluctuation. It rules out a dimension-uniform replacement of the conditional R_*^2 factor by sqrt(d) in the scalar a^2/(sigma N^3) bound. This does not rule out a covariance-calibrated clock scheme.

An exactly unbiased quadratic-calibrated alternative is

    integral_0^T cos(s)b(cos(s)z0+kappa sin(s)Z) ds
      =(pi/4)b(z0)
       +E_(U uniform(0,1))
          [b(sqrt(1-U)z0+kappa sqrt(U)Z)-sqrt(1-U)b(z0)]/(2sqrt(U)). (QC)

For every quadratic b the random bracket is the constant kappa b(Z)/2, so (QC) removes this first quadratic scale-mixture fixture. Its old-state first, however, contains

    sqrt(1-U)[Db(q_U)-Db(z0)]/(2sqrt(U)).

C2 supplies no uniform power rate for the Hessian difference. Low-clock truncation, compensation, caller derivatives and their actual prices are therefore an open gate. The finite-value cancellation in (QC) cannot be differentiated as a small error.

## 8. The reusable all-layer martingale question

The final residual is a true conditional martingale difference. An earlier residual is also centered at its own layer, but its later observer is the nonlinear force b(q+residual). Recursive conditioning alone does not make

    E[b(q+D)|past]=b(q).

Its exact mean defect is

    E integral_0^1 [b'(q+tD)-b'(q)] D dt,              (F)

after the linear b'(q) E D term vanishes. Under C2, b' is continuous and bounded by a but has no quantitative modulus. The uniform bound from (F) is O(a E|D|), which is precisely the prior-layer term in (A). A scalar smooth near-kink force b(x)=a[x/2+eta(sqrt(x^2+epsilon^2)-epsilon)] at q=0 and a symmetric D of size much larger than epsilon realizes an O(a E|D|) conditional mean defect, while its Hessian stays in a strict sandwich for eta<1/4.

This is a conditional source/caller obstruction, not a marginal posterior-law lower bound: integrating an unobserved nondegenerate Gaussian q can produce an additional cancellation/small-region factor. Any reusable all-layer current must retain the joint q,D law and the actual later observer, exploit that Gaussian structure without differentiating b' as an original oracle, and preserve one-energy dimension scaling. Replacing (F) by an unconditional covariance number or asserting that every layer stays unbiased would miss the live term.

The next concrete gate is an observer-qualified same-record current for (F) across one nonlinear subsequent force evaluation, with an actual positive implementation. Closing that gate would improve the nearest-prior-layer budget. The scalar last-layer gains alone do not change the user's general c_P objective.
