# Independent audit: innovation localization and cheap-source covariance separator

Date: 2026-10-04. Reviewer: independent localization audit. Scope: the two frozen extensions listed below, plus only the dependencies needed to identify their continuous and finite sources. The base same-carrier target contract and its counterexample receive a separate full review.

## Verdict

**PASS, bounded analytical claims.** The conditional innovation covariance bound, exact mixed-current identity, full covariance localization, and cheap-common-root covariance separator are correct under the stated anchored C2 potential class and the original Gaussian Markov endpoint convention.

**OPEN / NOT ADMITTED, executable order-four closure.** This audit supplies no finite original-gradient VALUE realization of j, K, the full C2 covariance reserve, or the same-carrier m3 source. It does not upgrade a conditional mean-law compiler to a strong value oracle, reuse its integrated roots, or claim a finite polylogarithmic count. The localization is an exact analytical reduction with a quantified remainder, not a completed sampler.

No failed analytical gate was found in the two pinned extensions. Several scope distinctions below are essential to this verdict.

## 1. Frozen inputs and exact meaning of the kernels

All paths in this section are relative to the parent `same-carrier-feedback` directory unless stated otherwise.

| Input | SHA256 |
| --- | --- |
| CHEAP-SHARED-ROOT-COVARIANCE-SEPARATOR.md | a336ca9e1b8e9923c70c6cbb3d88706f3e393803b5a9f4f87f0c85040d414ddd |
| CONDITIONAL-INNOVATION-LOCALIZES-THE-COVARIANCE-GATE.md | 186ccd57441273d5cdef1cb953bc0ba6f17fef3dc9c6dd2189a8dde91b965e63 |
| SAME-CARRIER-P3-MEAN-AND-FULL-COVARIANCE-GATE.md | a63c38f14206ec278f7c12a56542968a7d1675af6f550610e4539e500a653eca |
| ../THIRD-ORDER-NESTED-FORCE-MEAN-WITH-C2-DECOUPLING.md | b056e213ac28c60869727589d3ee6d1b2c53d6666908cdbe1961fa3121161ced |

Take dX_t = -X_t dt + sqrt(2) dB_t with X_0=Z and B independent of the standard Gaussian Z. This is the reversed multiplicative-clock history X_(exp(-t)); it has exactly the required conditional genealogy, not a common-G chord substituted for the process. Write

- I_t = integral_0^infinity exp(-u) g(X_(t+u)) du;
- Y_t = g(X_t-I_t), j(x)=E[Y_0 | X_0=x], e=j-g;
- H_f = integral_0^infinity exp(-t) f(X_t) dt;
- R_t=Y_t-j(X_t), S=integral exp(-t) R_t dt.

Then F2=H_j+S and E2=F2-H_g=H_e+S literally. The conditional kernel used by j is Law(x-I_0 | X_0=x), with I_0=F1. It is **not** the base contract's chi2(x)=Law(x-F2 | X_1=x), which supplies the integrand for m3. In particular m2=R1 j, whereas m3 needs the distinct second-shift kernel. None of these are posterior tilted endpoint kernels.

Cov(U,V) is E[(U-EU)(V-EV)^T]. This orientation fixes both K and the transpose in the square below.

## 2. Actual field and captured-endpoint bounds

Holding the future Brownian residual fixed,

    D_x I_0 = integral_0^infinity exp(-2u) Dg(X_u) du,
    ||D_x I_0||op <= A/2.

Because g is C1 with bounded derivative, differentiating the first-shift composition once is legitimate by dominated convergence:

    Dj(x)=E[Dg(x-I_0)(Id-D_x I_0) | X_0=x].

The leading expected Dg is symmetric. Consequently

    ||Dj||op <= A(1+A/2),
    ||Dj-Dj^T||op <= A^2,
    ||De||op <= 2A+A^2/2.

No derivative of Dg or modulus of continuity of Dg appears. Conditional Minkowski and g(0)=0 give

    ||I_0||_(L2|x) <= C A(|x|+sqrt(D)),
    |e(x)| <= C A^2(|x|+sqrt(D)).

The amplitude bound does not assert an A^2 first for e. Also j(0)=e(0) need not vanish. Thus neither a zero origin nor a small captured-x first may be inferred for a future executable source. The manuscript's explicit caution is correct.

The resulting conditional energy needed below is

    ||H_e||_(L2|Z=z) <= C A^2(|z|+sqrt(D)).

It follows directly by integrating the conditional Gaussian moments of X_t with positive exp(-t) weights. Its constant is dimension independent. Any broader caller distribution or tilt would still have to pay this actual |z| profile; the theorem here integrates the stated Z~gamma carrier only.

## 3. Conditional Clark–Ocone, orientation and the sharp constant

Let F_t=sigma(Z,B_u: u<=t). Markov conditioning gives E[R_t | F_t]=0. For s>t, X_t is fixed against the future Brownian variation, and

    D_s I_t = sqrt(2) integral_(s-t)^infinity
                    exp(-u) exp(-(t+u-s)) Dg(X_(t+u)) du,
    ||D_s I_t||op <= (A/sqrt(2)) exp(-(s-t)).

The matrix order in the terminal derivative is

    D_s Y_t = -Dg(X_t-I_t) D_s I_t,
    ||D_s Y_t||op <= (A^2/sqrt(2)) exp(-(s-t)).

At the diagonal s=t the chosen representative has no consequence for the stochastic integral. Applying Clark–Ocone on the future increments conditioned on F_t gives

    R_t = integral_t^infinity E[D_s Y_t | F_s] dB_s.

The conditional expectation is with respect to the full past F_s, and the integral starts at t. Reversing this orientation would not give the stated result. The derivative of j(X_t) with respect to future increments is zero.

Stochastic Fubini gives S=integral_0^infinity L_s dB_s with

    L_s=integral_0^s exp(-t) E[D_s Y_t | F_s] dt,
    ||L_s||op <= A^2 s exp(-s)/sqrt(2).

Conditional Ito isometry now proves, for every z and direction v,

    v^T Cov(S|z) v
      =E[integral |L_s^T v|^2 ds | z]
      <= |v|^2 (A^4/2) integral_0^infinity s^2 exp(-2s) ds
      = A^4 |v|^2/8.

Thus E[S|Z]=0, 0<=Cov(S|Z)<=A^4 Id/8, and its HS norm is at most A^4 sqrt(D)/8. This is a full operator inequality, not an entrywise estimate that accumulates an extra dimension factor.

The identities extend to the admitted C2 potential class via finite horizons and Gaussian Sobolev closure. The displayed deterministic majorants give integrable domination. They need only the original first derivative Dg. They do not supply a full-root or captured-Z Lipschitz bound of order A^2 for the same-tape S map: the small object is its predictable martingale integrand.

## 4. One-energy cross remainder and exact centered K current

Let H_e-E[H_e|Z=z]=integral Q_s dB_s be its conditional Brownian martingale representation. For every HS-unit matrix M,

    Cov(H_e,S|z):M = E integral Q_s:(M L_s) ds,
    ||M L_s||HS <= ||L_s||op.

Cauchy–Schwarz uses exactly one physical energy:

    ||Cov(H_e,S|z)||HS
      <= (E|H_e-EH_e|^2)^(1/2)
         (integral ||L_s||op^2 ds)^(1/2)
      <= C A^4(|z|+sqrt(D)).

It does not use the product of two sqrt(D)-sized energies. Integration over Z gives O(A^4 sqrt(D)). Therefore Cov(E2|Z)=Cov(H_e|Z)+O_L2HS(A^4 sqrt(D)) as claimed.

For the leading field H_g, the analogous cross is genuinely cubic. For fixed t,

    H_g = integral_0^t exp(-u)g(X_u)du + exp(-t)I_t.

The first term is F_t-measurable, so its product with R_t has conditional expectation zero. The future Markov property and the centering of R_t give exactly

    Cov(H_g,R_t|Z) = exp(-t) P_(exp(-t)) K(Z),
    K(x)=Cov(I_0,g(x-I_0) | X_0=x).

The same conditional history I_0 appears in both slots of K. Removing its centering, independently replacing the two slots, or reversing their order changes the target. Multiplying by S's exp(-t) weight and changing variables r=exp(-t) gives

    Cov(H_g,S|Z) = integral_0^1 r P_r K(Z) dr = R2 K(Z).

In particular there is no missing factor of r or 2.

The conditional Brownian first bound above even yields scalar-direction variance bounds Var(u^T I_0|x)<=A^2|u|^2/4 and Var(v^T Y_0|x)<=A^4|v|^2/4. Cauchy–Schwarz therefore gives the stronger uniform ||K(x)||op<=A^3/4, hence ||K(x)||HS<=A^3 sqrt(D)/4. The manuscript's looser bounds are valid. These are centered conditional bounds, so unbounded conditional means of I_0 do not enter. No correspondingly improved captured-x first for K follows from this amplitude statement.

## 5. General OU square and full covariance localization

For a deterministic C1 Lipschitz vector field f, define v_f=(1-L)^(-1)f=integral_0^1 P_r f dr. Ito's formula for exp(-t)v_f(X_t), followed by a finite-horizon limit, gives

    H_f = v_f(Z)+sqrt(2) integral_0^infinity
                                 exp(-t) Dv_f(X_t) dB_t.

It follows that

    Cov(H_f|Z)=2 R2[(Dv_f)(Dv_f)^T](Z).

This formula does not require f to be a gradient. The transpose is mandatory when f=j. Combining it with H_j=H_g+H_e, the exact K identity, and the remainders above yields

    Cov(F2|Z)=2 R2[B_j B_j^T](Z)
                      +2 Sym R2 K(Z)+O_L2HS(A^4 sqrt(D)),
    B_j=D(R1 j).

The retained first term includes the potentially large self-covariance of H_e. It has not been discarded or assigned an unjustified A^4 budget. The signed K correction is a covariance current, not itself a positive covariance source. No positivity or executable reserve construction is gained merely by writing this sum.

Independent scalar Brownian-filter calculations find

    k_H0(s)=exp(-s)/sqrt(2),
    k_J0(s)=exp(-s)(s+1/2)/sqrt(2),
    k_S(s)=-a^2 s exp(-s)/sqrt(2).

They reproduce Var(S|Z)=a^4/8, Cov(H_g,S|Z)=-a^3/8, K=-a^3/4, and

    Var(F2|Z)=a^2/4-a^3/2+5a^4/16,
    retained localization=a^2/4-a^3/2+a^4/16.

The residual is exactly a^4/4. A separate nonnormal linear-field check detects the erroneous B_j^2 orientation immediately; the correct covariance uses B_j B_j^T.

## 6. Cheap-source separator: source, potential and limiting covariance

The audited cheap finite source is the actual common-G,H packet from the pinned nested-force note:

    x_i=r_i Z+sqrt(1-r_i^2)G,
    H_Q(x,H)=sum_j v_j g(tau_j x+sqrt(1-tau_j^2)H),
    E_Q=sum_i w_i [g(x_i-H_Q(x_i,H))-g(x_i)].

The same G serves every outer node and the same independent H serves all inner occurrences. The limiting source is not silently given Markov roots.

For the stated counterfamily, A=D^(-1/2), v_D=(1,...,1)/sqrt(D), c=1/4, epsilon=1/2,

    g_D(x)=c A [b(x_i)]_i+c A v_D tanh(v_D^T x),
    b(t)=t+epsilon log cosh(t).

It is an anchored smooth gradient and its actual Hessian is

    c A diag(1+epsilon tanh(x_i))
        +c A sech^2(v_D^T x) v_D v_D^T.

Hence 0<=Dg_D<=c(2+epsilon)A Id=5A Id/8<A Id. The example lies strictly inside the admitted Hessian bound and uses no dimension-dependent derivative assumption.

At every inner point, the one-coordinate Gaussian marginal is standard. Positive weights of total mass one therefore fix the coordinate mean E b(N)=epsilon E log cosh(N), independently of root genealogy and node count. The coordinate averages are averages of independent triples (Z_i,G_i,H_i), before the explicitly isolated rank-one terms. Uniform fourth moments and positive weights give LLN errors O(D^(-1/2)) and the same rate for the summed Taylor remainder, uniformly over the number and positions of quadrature nodes. Thus no hidden clock-count factor spoils the joint dimension/clock limit.

The limiting projected chords are c[C+J_M] and c[C'+J_C], where

    J_M=integral_0^1 h(U_r^M)dr,
    J_C=integral_0^1 h(rZ_0+sqrt(1-r^2)G_0)dr,
    h(u)=tanh(u-d)-tanh(u), d=c epsilon E log cosh(N)>0.

Constants disappear after conditional centering. In fact the two continuum constants agree: their coordinate term uses only the pair law with correlation tau. For finite exact-first-moment rules, the relevant k also agrees because E[log cosh(Y) tanh(X)]=0 by joint sign symmetry and E[Y tanh(X)]=tau E[X tanh(X)]. Equality of the constants is not needed for the separator.

The marginal (U_r,Z_0) law is identical in the two sources. Thus their conditional means E[J|Z_0]=integral P_r h(Z_0)dr agree exactly. But for 0<r<s<1,

    rho_C(r,s)=rs+sqrt(1-r^2)sqrt(1-s^2)>r/s=rho_M(r,s).

Both correlations are nonnegative. Expanding h in Gaussian Hermite polynomials gives a sum of nonnegative terms for Var(J_C)-Var(J_M). The degree-one term is strictly positive because

    E h'(N)=E sech^2(N-d)-E sech^2(N)<0.

The strict sign follows from Gaussian convolution of the even strictly decreasing sech^2 profile, or its interval layer-cake representation. Also

    integral integral rho_M=1/2,
    integral integral rho_C=1/4+pi^2/16.

Therefore

    Delta:=Var(J_C)-Var(J_M)
      >= (E h'(N))^2 (pi^2-4)/16>0.

Since the conditional means agree, Delta equals the difference of integrated conditional variances. The additional endpoint coordinates orthogonal to v_D are independent of the scalar Gaussian history. Conditional-expectation contraction transfers the projected L2 limit to the integrated conditional variances even though the endpoint dimension varies. Consequently

    E v_D^T[Cov(E_C|Z_D)-Cov(E2|Z_D)]v_D
       =c^2 A^2 Delta+o(A^2),

and the L2(Z_D;HS) norm of that matrix difference is Omega(A^2). This is an integrated norm claim; it is not a pointwise-in-z Loewner ordering of the full covariance difference.

Positive quadrature measures with vanishing uniform moment error converge weakly to Lebesgue measure by polynomial density. The common-G scalar integrand is bounded and continuous on [0,1], so the outer sums converge in L2. Combined with the node-count-uniform LLN control above, this validates the simultaneous dimension/clock assertion. A fixed sufficiently accurate positive outer rule already preserves a strict positive gap; further accuracy cannot repair the common-root genealogy.

On A=D^(-1/2), the claimed A^4 sqrt(D) tolerance is A^3 while the source error is Omega(A^2). Fixed polynomial public-log losses cannot close that separation. This rules out only the identified cheap source as the covariance producer; it is not an impossibility theorem for other finite sources.

## 7. Independent diagnostics and remaining admission gates

`check_localization_independent.py` was written separately and does not import or run the author's diagnostic. It hashes the four inputs, derives the Brownian kernels directly, integrates their products symbolically, checks a nonnormal field's transpose, and independently evaluates Gaussian/Hermite calibrations. `localization_independent_checks.json` records its full output.

Numerical calibration gives d approximately 0.0468209009364293 and E h'(N) approximately -0.00039839913814524. The rigorous first-chaos lower-bound expression evaluates to about 5.82271628702e-8 for Delta, or 3.63919767939e-9 after multiplying by c^2. Truncated Hermite sums suggest a much larger actual gap, but they are diagnostic only and are not used in the proof or the asymptotic scope.

Before any constructive PASS, a new service must still provide:

1. Finite original-gradient VALUE targets for the actual j square and same-I centered K, or a direct full Cov(F2|Z) reserve with the stated target.
2. The distinct same-carrier m3 source conditional on the original Gaussian endpoint.
3. Full finite records, shared roots where required, independent complete consumer banks where required, actual captured-caller firsts, nonzero origins, replay accounting, HVP-only first/adjoint use, frozen guards, absolute floors, and the required fixed-public-log VALUE count.

The analytical Brownian history has zero executed cost only because it is not executed. Neither this audit nor the localization author has supplied an executable replacement for it.
