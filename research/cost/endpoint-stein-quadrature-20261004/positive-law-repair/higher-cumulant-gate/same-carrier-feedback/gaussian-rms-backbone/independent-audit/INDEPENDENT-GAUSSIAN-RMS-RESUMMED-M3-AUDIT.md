# Independent audit: Gaussian-RMS residual resummed m3 extension

Date: 2026-10-05. Independent source and mathematical audit; author files were not edited by the auditor.

## Verdict and frozen scope

**PASS for the Gaussian-RMS residual comparison, the explicit unbounded tail subclass, the strict enlargement of the useful bounded-certificate class, the scalar-backbone obstruction, and the unchanged guarded finite positive m3 mean consumer.**

The final inspected author artifact is `../GAUSSIAN-RMS-RESIDUAL-RESUMMED-M3.md`, SHA256:

    18640f9c6c567acba4fe830c7b0d64394924887ee4335a7a49f57f817435dccf

The conclusion is an integrated original-Gaussian endpoint result. It is not a uniform endpoint bound, conditional-rescaling theorem, general convex-gradient closure, strong mean statistic, or fourth-order endpoint join. The completed finite mean compilers and their guards remain imported components, not programs reimplemented by this audit.

With the declared fixed scalar lambda and residual r=g-lambda Id, the new sharp-in-form comparison is

    ||F2-(lambda H1-lambda^2 H2)||2
        <= epsilon1+(lambda+Lr)epsilon0,
    ||m3-R1 psi_G||2
        <= A[epsilon1+(lambda+Lr)epsilon0].

Both estimates are correct on the genuine stationary future OU history. No independent-ancestry substitution, Gaussian-density comparison, hidden conditional-expectation producer, or free fitting oracle is needed. The sufficient residual grade is exactly the displayed bracket O(A^3 sqrt(D)).

The independent checker passes **1,350 assertions**. It imports no author checker, verifies frozen source/import hashes, symbolically rederives covariance and obstruction formulas, evaluates tail integrals, checks exact linear residual instances, tests actual finite source first/curl/count ports on the unbounded radial fixture, and performs a genuine shared-future OU diagnostic. Maximum finite-source directional derivative discrepancy: **4.777e-11**. Maximum finite-ancestry decomposition discrepancy: **3.886e-16**.

The author accepted two precision clarifications raised during audit: the best bounded residual in the tail example is *at least* of order A sqrt(D), because its exact amplitude also contains sqrt(log(1/A)); and variance contraction does not by itself imply epsilon1<=epsilon0. The frozen artifact includes both clarifications. No remaining mathematical blocker was found in its scoped claims.

## 1. The RMS force comparison, independently derived

Assume g(0)=0, g=grad U, 0<=Dg<=AI, 0<A<=1/2 and a fixed 0<=lambda<=A. Then

    Dr=Dg-lambda I,
    -lambda I<=Dr<=(A-lambda)I,
    Lip(r)<=Lr=max(lambda,A-lambda).

This bound uses the symmetric Jacobian interval, not a triangle estimate A+lambda. In particular r(0)=0 and all required Gaussian second moments exist. Since g is A-Lipschitz and anchored, the exponential history integrals are well-defined in L2. Bochner Fubini and Minkowski are justified by their exponentially weighted finite second moments.

For the genuine future process, put

    F1,t=integral_0^infinity e^-u g(X_(t+u)) du,
    H1,t=integral_0^infinity e^-u X_(t+u) du,
    Delta_t=F1,t-lambda H1,t.

Stationarity and the positive mass-one integration kernel give, for every t,

    ||Delta_t||2
      <=integral e^-u ||r(X_(t+u))||2 du
      =epsilon0.

No independence is involved. Let Y_t=X_t-lambda H1,t. Direct covariance integration gives

    Cov(X_t,H1,t)=integral_0^infinity e^-u e^-u du I=I/2,
    Var(H1,t)=I/2.

Thus Y_t is a centered Gaussian with covariance

    sigma1^2 I=(1-lambda+lambda^2/2)I.

In particular ||r(Y_t)||2=epsilon1. Keeping the same actual Delta_t,

    ||r(X_t-F1,t)||2
      =||r(Y_t-Delta_t)||2
      <=epsilon1+Lr epsilon0.

The actual second substitution expands exactly as

    F2=lambda H1-lambda^2 H2+R2,
    R2=-lambda integral e^-t Delta_t dt
          +integral e^-t r(Y_t-Delta_t)dt,

because the double future convolution is

    integral_(t,u>=0) e^-(t+u) X_(t+u) dt du
      =integral_0^infinity v e^-v X_v dv=H2.

Consequently

    ||R2||2<=epsilon1+(lambda+Lr)epsilon0.

The factor lambda+Lr equals max(A,2lambda), hence is at most 2A. If an independently certified smaller residual Lipschitz bound is known, it can replace Lr; no such refinement is necessary.

An independent replacement of H1,t would instead give Var(X_t-lambda H1,t^ind)=(1+lambda^2/2)I, missing the -lambda term. The supplied proof correctly uses the genuine future correlation.

## 2. Exact conditional backbone and target comparison

The joint endpoint/history scalar covariance is

    Cov(X0,H1,H2)=
      [[1, 1/2, 1/4],
       [1/2, 1/2, 3/8],
       [1/4, 3/8, 3/8]].

Conditioning on X0 subtracts the endpoint regression product, giving

    E[H1|X0=x]=x/2, E[H2|X0=x]=x/4,
    Cov((H1,H2)|X0)=
      [[1/4,1/4],[1/4,5/16]] tensor I.

For B=lambda H1-lambda^2 H2, this yields the exact conditional law

    B|X0=x =_law a2 x+b2 N,
    a2=lambda/2-lambda^2/4,
    b2^2=lambda^2/4-lambda^3/2+5lambda^4/16.

The nonnegative known scalar root is legitimate, including lambda=0. There is no source-dependent covariance root.

On the original coupling, conditional Jensen and the outer Lipschitz bound imply

    ||E[g(X0-F2)-g(X0-B)|X0]||2
       <=A||F2-B||2.

Replacing only the conditional distribution of B now identifies the second conditional expectation with psi_G. This proves the author’s psi comparison. The L2(gamma) contraction of R1 proves the m3 comparison.

This argument deliberately does not try to condition the unconditional epsilon bounds at every endpoint. A small Gaussian RMS does not give a uniform endpoint error. The final artifact makes that distinction explicitly.

## 3. Computable certificates and an invalid shortcut to avoid

Suppose |r(x)|<=e_core on the Euclidean ball of radius R, and r has a certified global Lipschitz bound Lr. Projection onto that ball gives pointwise

    |r(x)|<=e_core+Lr(|x|-R)_+.

For either s=1 or s=sigma1, ordinary Minkowski therefore gives

    ||r(sZ)||2<=e_core+Lr sqrt(J_D(s,R)),
    J_D(s,R)=E[(s|Z|-R)_+^2].

The chi integral in the author artifact is exact, including its normalization. An alternative cancellation-free upper certificate is

    J_D(s,R)
       <=s^2 D Q(D/2+1,R^2/(2s^2)),

where Q denotes the regularized upper incomplete gamma function. This follows by dropping the nonnegative subtraction in the hinge and using the truncated chi second moment. Numerical evaluations of either representation need certified integration/tail/rounding errors before being used as theorem input. The checker’s ordinary floating-point quadrature is diagnostic; the explicit analytic certificate below is the proof.

For R=sqrt(D)+T, T>=0, the chi-square Chernoff bound yields

    P(|Z|>=sqrt(D)+t)<=exp(-t^2/2).

Indeed the optimizing chi-square exponent at x=(sqrt(D)+t)^2 is

    -(x-D-D log(x/D))/2
      =-t^2/2-D[u-log(1+u)], u=t/sqrt(D),

which is at most -t^2/2. Integrating the survival function proves

    J_D(1,R)<=2 integral_0^infinity v exp(-(T+v)^2/2)dv
             <=2exp(-T^2/2).

For 0<s<=1 the hinge envelope is pointwise smaller, so the same bound holds at s. This is dimension-safe and is not a ratio of Gaussian densities.

The latter monotonicity applies to the radial hinge envelope, **not** to arbitrary r. Here is an exact one-dimensional counterfixture to epsilon1<=epsilon0:

    A=1/2, lambda=eta=A/2, w=1/10,
    r(x)=eta x exp(-x^2/(2w^2)),
    g(x)=lambda x+r(x).

Writing v=x/w, the residual derivative is eta(1-v^2)exp(-v^2/2). Its maximum is eta, its minimum is -2eta exp(-3/2), and therefore 0<=g'<=A. The source is smooth and anchored. Direct Gaussian integration gives

    ||r(sZ)||2^2
      =eta^2 s^2(1+2s^2/w^2)^(-3/2).

At sigma1=sqrt(25/32)=0.883883..., the ratio epsilon1/epsilon0 is **1.062549...**, greater than one. The frozen author text correctly keeps the two general norms separate.

## 4. The unbounded radial example is valid and strictly enlarges the useful class

For lambda=d=A/2, use the author’s quintic monotone step psi, its integrated function h, and

    T=sqrt(8 log(1/A)), R=sqrt(D)+T,
    g(x)=lambda x+d h(|x|-R)x/|x|.

The step is C2 because its first two derivatives match at both endpoints; h is C3. The nonlinear radial field vanishes on a neighborhood of the origin. Thus g is globally C3, in particular C2, and is the gradient of the displayed radial potential.

For rho>0, its Jacobian eigenvalues are

    radial: lambda+d psi(rho-R),
    tangential: lambda+d h(rho-R)/rho.

Since 0<=psi<=1 and 0<=h(rho-R)<=rho, both lie in [lambda,A]. The D=1 statement needs only the radial eigenvalue. This checks all original source hypotheses.

At the selected slope lambda, r has magnitude d h(rho-R), vanishes on the core, and grows linearly at infinity. It is genuinely unbounded. Nevertheless the certified hinge comparison gives

    epsilon0,epsilon1<=d sqrt(2)exp(-T^2/4)
                     =A^3/sqrt(2).

Because Lr=lambda=A/2, the bracket is at most (1+A)A^3/sqrt(2). Hence the canonical mean bias is at most

    (1+A)A^4/sqrt(2),

which has the desired dimensional grade for every integer D>=1. No lower dimension threshold or upper bound on A sqrt(D) is needed.

This is a strict enlargement of the class certified at order four by a global scalar remainder. For rho>=R+1,

    g(rho n)=A rho n-d(R+1/2)n.

If mu!=A, g-mu Id is unbounded. At the sole bounded slope mu=A, the exact supremum is d(R+1/2). The old bias certificate divided by A^4 sqrt(D) is therefore

    (1+A)(R+1/2)/(2A^2 sqrt(D))>=(1+A)/(2A^2).

That loss is not a fixed public-polylog factor as A tends to zero. This comparison concerns the old certificate, not a lower bound on the actual error produced by the slope-A backbone. The author correctly preserves that distinction.

The diagnostic sweep covers A in {0.5,0.2,0.05,0.01} and D in {1,5,100,10000}. It checks actual one-dimensional residual integrals, the incomplete-gamma envelope, the dimension-uniform analytic envelope, both Hessian eigenvalues, and the exact alternative bounded-remainder amplitude. All tests pass. These checks are deliberately not substituted for the analytic tail proof.

## 5. Finite source ports and full costs remain unchanged

Only the analytical target comparison changes. The raw program still uses

    x_i=t_i Z+c_i G,
    F_G=sum_i w_i g((1-a2)x_i-b2 N),
    B=sum_i w_i g(x_i), E=F_G-B.

The same G and N are shared across each entire positive-rule level. Linearity of conditional expectation gives E[F_G|Z]=Q psi_G(Z), without any claim that its raw covariance equals that of the canonical force. The norm bound ||psi_G||2<=C A sqrt(D) follows solely from g(0)=0, Lip(g)<=A and the known Gaussian affine argument. The imported Hermite-multiplier rule gives the stated C delta A sqrt(D) bias.

The source-admission estimates do not involve epsilon0, epsilon1, or a global residual bound. With alpha=1-a2 and H_i=Dg(alpha x_i-b2 N), H_i^0=Dg(x_i), the exact derivatives are

    D_G E=sum_i w_i c_i(alpha H_i-H_i^0),
    D_N E=-b2 sum_i w_i H_i,
    D_Z E=sum_i w_i t_i(alpha H_i-H_i^0).

Because alpha H_i and H_i^0 belong to [0,AI], their difference is symmetric with norm at most A. This proves the stated G, N, combined-private and captured first bounds. In the P_G square lift the only skew block is D_N E, so the curl is at most A b2<=A^2/2. Hessian commutation and Hessian derivatives are unnecessary.

The original Lipschitz chord estimate gives

    |E|<=A sum_i w_i[a2|x_i|+b2|N|],

and therefore the same conditional O(A^2)(|Z|+sqrt(D)) energy. Caller-only origins are actual VALUE evaluations. Anchoring preserves the private first/curl bounds; caller derivatives of the recorded origins must still be included, as in the prior audit. One must not infer a derivative theorem from the new small L2 target-bias bound.

The independent finite-graph diagnostic uses the unbounded radial source itself, at bulk, shell and tail scales. It checks VALUE occurrence counts, same-record subtraction, deterministic origins, literal total zero, energy bounds, first/curl bounds and directional chain rules. Its Hessians genuinely fail to commute: the maximum tested commutator norm is 0.005387..., while the G-block skew remains zero to numerical precision.

The existing guarded gradient/near-gradient consumers thus apply at exactly the same ports. Half-share normalized radii are unchanged; A<=1/2 alone is not enough to discharge them. Independent COMPLETE banks are required conditional on the retained Z and exterior labels. Product coupling occurs only after the complete laws return. Original private roots cannot be reattached as observers.

The ledger is correctly enlarged only to expose possible setup work:

    Q_certificate_setup+Q_captured
      +N_out N_B+2N_out N_E+Q_known/numerical/replay.

For the explicit analytic example Q_certificate_setup=0 original g VALUES; known arithmetic is still charged. No generic Gaussian fitting or black-box residual-validation procedure is supplied for free. N_B,N_E remain fully expanded compiler occurrence counts, not callback counts. Every changed raw argument replays the original leaves. Absolute numerical floors remain additive, never divided by a possibly zero residual norm.

The source/consumer result is therefore conditional on the previously admitted compiler contracts and guards. This audit does not implement or independently reprove those complete compilers.

## 6. Exact scalar obstruction

For even D and K having equally many eigenvalues k1=A/3 and k2=2A/3, g(x)=Kx is an admissible convex gradient. The canonical and scalar-backbone means are exactly

    m3(Z)=(1/2)K[I-a(K)]Z,
    m_G,lambda(Z)=(1/2)K[1-a(lambda)]Z,
    a(k)=k/2-k^2/4.

The positive outer rule's first moment makes the linear mean quadrature exact. The scalar b2 innovation has zero conditional mean after the terminal linear g.

For c=a(lambda), the exact squared discrepancy is

    (D/8)[k1^2(c-a(k1))^2+k2^2(c-a(k2))^2].

Its minimizing c is the k_i^2-weighted mean of the two values a(k_i). Since a is increasing on [0,A], this minimizer is attainable with lambda between k1 and k2. Substitution gives

    inf_lambda ||m_G,lambda-m3||2
      =sqrt(D/8) k1 k2 |a(k2)-a(k1)|/sqrt(k1^2+k2^2)
      =A^2(2-A)sqrt(D)/(36sqrt(10)).

This is the canonical mean discrepancy, not the RMS of r itself. A fixed affine intercept adds only a constant mean component and cannot cancel its centered linear part. Thus even an oracle selecting the best scalar cannot give generic order-four closure. Known matrix or endpoint-dependent backbones remain outside this obstruction's scope.

## 7. Conditional normalization and conclusion

For f(y)=s[g(a+sy)-g(a)], the residual with lambda_f=lambda s^2 is exactly

    r_f(y)=s[r(a+sy)-r(a)].

Its required Gaussian certificates involve the translated arguments a+sZ and a+s sigma1,f Z. The original centered epsilon0,epsilon1 do not uniformly bound them. The tail example makes this visible when a moves into its nonlinear tail. The author correctly requires actual translated certificates and actual caller graphs rather than silently transporting the original grade.

No remaining blocker was found for the frozen component theorem. The general m3 consumer, data-driven backbone selection, uniform conditional rescaling, full endpoint feedback/restoration join, and all-order complexity remain separate open work.

## 8. Evidence and reproducibility

Run:

    python independent-audit/check_gaussian_rms_extension.py

Outputs:

- `gaussian_rms_independent_checks.json`: all 1,350 checks, grouped counts, numerical diagnostics and exact source/import pins.
- `MANIFEST.json` and `SHA256SUMS`: audit artifact provenance.

Directly inspected imported construction:

- `../../resummed-linear-backbone/BOUNDED-REMAINDER-RESUMMED-M3-VALUE-CONSUMER.md`, SHA256 `798d20b87124e7069e465ce4e7bcf92c5af3734d86c8a62aeb288b5fec8eaef3`.
- `../../resummed-linear-backbone/independent-audit/INDEPENDENT-BOUNDED-RESUMMED-M3-AUDIT.md`, SHA256 `49d135ba6692c54012acd9d0d94f603c0015fca98b5d1be1164457ff45d1213a`.

Their previously audited quadrature, gradient/near-gradient, coisometry, caller-origin and original-VALUE contracts remain imports. No private-root retention or endpoint closure is certified by this report.
