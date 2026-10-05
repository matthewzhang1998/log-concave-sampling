# Independent audit: bounded linear-backbone canonical m3 mean consumer

Date: 2026-10-05.

## Verdict

**PASS for the stated bounded source subclass, the canonical Gaussian-history mean comparison, and the guarded conditional unit-buffer mean-law consumer.** The original source must be fixed before the fresh Gaussian carrier is sampled and must satisfy the declared global certificate

    g = grad U, g(0)=0, 0 <= Dg <= A I,
    g(x)=lambda x+r(x), 0 <= lambda <= A,
    sup_x |r(x)| <= epsilon.

The new source comparison is uniform in the endpoint and dimension:

    |psi2(x)-psi_G(x)| <= A epsilon(1+lambda).

At the original standard-Gaussian carrier Z, the completed positive program consequently has the claimed integrated conditional W2 allowance

    A epsilon(1+lambda) + C delta A sqrt(D)
      + Lambda A^4 sqrt(D) + restored absolute floors,

subject to the actual fixed-order gradient/near-gradient consumer guards. If epsilon <= d A, D >= A^(-4), and delta <= A^3, the intrinsic error is O_d(Lambda A^4 sqrt(D)). This is a genuine source-qualified positive result beyond the raw three-root source's unresolved general target-bias gate.

It is not a generic C2 mean theorem, a strong mean statistic, a joint approximation retaining old private roots, a same-endpoint posterior sampler, or a reverse-OU endpoint join. The dimension condition is a **lower** bound used to compare an absolute bounded-remainder error to the requested dimensional grade. It is not a new upper bound on A sqrt(D).

### Frozen source and independent evidence

Audited source: `../BOUNDED-REMAINDER-RESUMMED-M3-VALUE-CONSUMER.md`.

Final inspected SHA256:

    798d20b87124e7069e465ce4e7bcf92c5af3734d86c8a62aeb288b5fec8eaef3

The initial request named `bb35c68c...`; the author subsequently added an explicit conditional-rescaling limitation and the sufficient first/radius/origin bounds supplied during this audit. This report certifies the final pin above. No author file was edited by this auditor.

`check_resummed_m3_independent.py` passes **1,436 assertions** without importing an author checker. It checks exact symbolic OU covariances, coefficient positivity, positive finite quadratures, source leaves and replay counts, shared-root Gaussian laws, genuine-future convolution algebra, and complete first/curl/origin bounds on a smooth convex-gradient fixture with genuinely noncommuting Hessians. Its maximum directional chain-rule discrepancy is 6.994e-11. Output: `independent_resummed_m3_checks.json`.

The tests do not instantiate the entire previously proved mean compiler. Their scope and the imported proof pins are explicit below.

## 1. The exact target and the genuine ancestry match

Use physical OU time t >= 0, with covariance exp(-|t-u|), conditioned at X_0=x. This is the logarithmic parametrization of the earlier canonical history anchored at multiplicative time 1. Define the actual shifted first force

    F1,t = integral_0^infinity e^(-u) g(X_(t+u)) du,
    H1,t = integral_0^infinity e^(-u) X_(t+u) du.

These are analytical variables on the same history, not separate independently sampled ancestors. Their difference is

    F1,t-lambda H1,t
      = integral_0^infinity e^(-u) r(X_(t+u)) du,

whose norm is at most epsilon pathwise. Bounded r and the finite Gaussian first moments justify these integrals and the necessary Fubini interchanges. The same inequality holds at every genuine future time t, without requiring independence from any other part of the history.

Now expand the actual second substitution:

    F2 = integral e^(-t) g(X_t-F1,t) dt
       = lambda H1 - lambda^2 integral e^(-t) H1,t dt
         -lambda integral e^(-t)(F1,t-lambda H1,t)dt
         + integral e^(-t) r(X_t-F1,t)dt.

The double linear integral has the exact convolution kernel

    integral_(t,u>=0) e^(-t-u) X_(t+u) dt du
      = integral_0^infinity v e^(-v) X_v dv = H2.

Thus F2=lambda H1-lambda^2 H2+R2, with

    |R2| <= lambda epsilon+epsilon.

No centeredness, small covariance, Gaussianization, or relation between R2 and the backbone is assumed. In particular the proof has not replaced F1,t by an independent force. The inner force is the one in the canonical target from `P3-MEAN-RAW-SOURCE-AND-SECOND-DECOUPLING-GATE.md`.

For every fixed x, compare the two terminal original g values on this same analytical history. Its global Lipschitz bound gives

    |g(x-F2)-g(x-lambda H1+lambda^2 H2)|
      <= A epsilon(1+lambda).

Conditioning and Jensen preserve that bound. Replacing only the **law** of the linear Gaussian backbone by its exact conditional Gaussian representation therefore proves the claimed uniform psi comparison. This replacement is not a strong coupling to a previously exposed private source record.

Applying the L2(gamma) contraction R1 proves the m3 comparison. This bypasses the general second-decoupling gate by a separate, stronger bounded-amplitude hypothesis; it does not purport to prove that earlier unproved gate.

## 2. Exact conditional Gaussian coefficients

The endpoint regressions are

    E[H1|x]=x/2, E[H2|x]=x/4,

because the conditional OU mean at t is e^(-t)x. The unconditional two-parameter covariance is

    K(a,b) = [1/(a+b)] [1/(a+1)+1/(b+1)].

Splitting the positive quadrant into t >= u and u >= t proves this identity directly. Differentiating the exponentially integrable scalar kernels at a=b=1 gives

    Var(H1)=1/2,
    Cov(H1,H2)=3/8,
    Var(H2)=3/8.

Subtracting endpoint-regression products produces

    Cov((H1,H2)|x) = [[1/4, 1/4], [1/4, 5/16]] tensor I.

The determinant of the displayed 2x2 matrix is 1/64. Applying the row (lambda,-lambda^2) gives exactly

    a2 = lambda/2-lambda^2/4,
    b2^2 = lambda^2/4-lambda^3/2+5lambda^4/16.

In particular

    b2 = (lambda/4) sqrt(4-8lambda+5lambda^2),
    4-8lambda+5lambda^2 = 5(lambda-4/5)^2+4/5 > 0.

The factored expression is a convenient numerically nondegenerate way to evaluate the same known scalar at lambda=0. It creates neither a source-dependent covariance root nor a division by lambda. With 0 <= lambda <= A <= 1/2,

    0 <= a2 <= lambda/2, 0 <= b2 <= lambda/2,
    0 <= alpha:=1-a2 <= 1.

The root N in psi_G is fresh standard Gaussian and supplies the entire conditional backbone distribution. Gaussianity, mean, and covariance are exact here; there is no central-limit approximation and no lost orientation of a nonscalar covariance matrix.

## 3. The finite source has the advertised own mean

At the retained standard Gaussian Z, draw one G and one N independently, and use them at every node of the frozen outer rule:

    x_i=t_i Z+c_i G, c_i=sqrt(1-t_i^2),
    F_G=sum_i w_i g(alpha x_i-b2 N).

One occurrence contains exactly N_out original VALUE leaves. There are no old original-g ancestor calls behind the known a2,b2 coefficients. Each summand has conditional expectation P_(t_i) psi_G(Z), so linearity of expectation yields

    E[F_G|Z]=Q psi_G(Z).

The shared roots do not invalidate this mean identity. They do change the raw covariance: for g(x)=lambda x, it is

    Cov(F_G|Z)=lambda^2[(alpha beta)^2+b2^2] I,
    beta=sum_i w_i c_i.

Using independent G_i or N_i instead would generally give a different covariance. The checker detects that distinction. The construction neither identifies this raw covariance with the canonical force covariance nor needs that identification for an own-mean law compiler.

The uniform Hermite-multiplier quadrature theorem is applied to the actual one-variable function psi_G. Since x and N are independent standard Gaussians in its L2 norm, anchoring and Lipschitz continuity give

    ||psi_G||2 <= A sqrt(alpha^2+b2^2) sqrt(D) <= C A sqrt(D).

Hence the actual multiplier certificate gives

    ||Q psi_G-R1 psi_G||2 <= C delta A sqrt(D).

This is an L2(gamma_Z) result, not a pointwise quadrature theorem for arbitrary unbounded Z. It requires Z's stated standard-Gaussian law after any original exterior labels have been fixed. No A-dependent regularity of Dg or clock derivatives of g is used.

The imported positive dyadic Gauss-Legendre rule has mass one, exact first moment 1/2, and sufficient error bound 8*4^(-m)+2^(1-K). Choosing m,K=O(log(1/delta)) supplies N_out=O(log^2(1/delta)). The exact identities refer to mathematical coefficients; their finite-precision representation remains a numerical allowance.

## 4. Complete private first, curl, energy and origins

The literal split is on one complete occurrence:

    B=sum_i w_i g(x_i), E=F_G-B.

The baseline is a genuine gradient in G. Its Jacobian is sum_i w_i c_i Dg(x_i), symmetric positive semidefinite, with first at most A beta. An analytical potential is sum_i (w_i/c_i)U(t_iZ+c_iG); this potential is never queried. Each c_i is strictly positive.

For the correction, set H_i=Dg(alpha x_i-b2 N), H_i^0=Dg(x_i). The actual chain rules are

    D_G E=sum_i w_i c_i(alpha H_i-H_i^0),
    D_N E=-b2 sum_i w_i H_i,
    D_Z E=sum_i w_i t_i(alpha H_i-H_i^0).

The two terms alpha H_i and H_i^0 are in [0,AI]. Their difference is symmetric and has operator norm at most A, even when its size does not gain a power of A. Thus

    ||D_G E|| <= A beta,
    ||D_N E|| <= A b2,
    ||D_(G,N) E|| <= A sqrt(beta^2+b2^2),
    ||D_Z E|| <= A/2.

In the square lift P_G^* E=(E,0), its skew derivative has block form

    [[0, D_N E],[-(D_N E)^*,0]].

Its operator norm is exactly ||D_N E||, so the full curl is at most A b2 <= A^2/2. No local Hessians have been commuted. The large possible leading Hessian difference has no skew contribution. No derivative of a Hessian has been taken.

The value bound is independently obtained by the original Lipschitz inequality:

    |E| <= A sum_i w_i(a2 |x_i|+b2 |N|).

For any fixed finite p this implies

    ||E||_(Lp|Z=z) <= C_p A(a2+b2)(|z|+sqrt(D))
                    <= C_p A^2(|z|+sqrt(D)).

The deterministic origins are actual values, not deductions from a Gaussian moment bound:

    B0(z)=sum_i w_i g(t_i z),
    E0(z)=sum_i w_i[g(alpha t_i z)-g(t_i z)].

They cost N_out and 2N_out original leaves, respectively, before exact-key sharing, and obey

    |B0(z)| <= A|z|/2,
    |E0(z)| <= A a2 |z|/2.

Subtracting these recorded origins makes the private source zeros literal. The anchored correction E-E0 retains the same O(A^2)(|z|+sqrt(D)) energy profile and the same private first/curl bounds. Its actual caller first is at most A; subtracting E0 need not preserve the sharper A/2 constant, but it preserves the required O(A) scale. The analogous anchored baseline has a valid O(A) caller bound. These caller paths, including any original finite-mode anchor, must remain in the actual first sweep.

At Z=G=N=0, every raw value is zero by coherent use of g(0)=0. At nonzero Z, B0 and E0 are not generally zero. When lambda=0, a2=b2=0 and E is exactly zero on the same record. This does not require the original g to vanish; bounded nonlinear monotone sources are still possible in that edge case.

## 5. The admitted positive mean programs apply at these ports

For definiteness take positive variance shares v_B=v_E=1/2. The literal source-relative construction can be written as

    M_B=B0+sqrt(v_B) N_k((B-B0)/sqrt(v_B)),
    M_E=E0+sqrt(v_E) P_G N^v(P_G^*(E-E0)/sqrt(v_E)),
    M_G=M_B+M_E.

The recorded-coisometry square-lift adapter uses the known Gaussian completion from the imported consumer. The whole B and E programs use independent COMPLETE banks conditional on the same Z and all exterior labels. Inside each E occurrence the baseline and shifted value retain their same G,N.

Sufficient normalized radii are

    rho_B=sqrt(2) A beta,
    ell_E=sqrt(2) A sqrt(beta^2+b2^2).

One must impose rho_B <= r_*(k,Delta,Lambda_k), ell_E <= 1/4, and the actual fixed-order finite-clock/filter/precision/active-dimension guards. A <= 1/2 is a convenient raw-source bound and is not, by itself, a proof of all compiler guards. The normalized square lift has private dimension 2D and curl at most sqrt(2)A b2. Its declared relative curl can be O(A); it is not divided by an unknown measured energy or an accidentally small measured first.

The original near-gradient theorem has the law-error row

    Lambda e_E {ell(a_seed+mu)+ell^3(1+mu^(-1/2))}
      + absolute numerical floors.

Use declared ell=C A, a_seed=C A, mu=A. Its bracket is O(A^2). The actual anchored energy from Section 4 is O(A^2)(|z|+sqrt(D)), giving O(Lambda A^4)(|z|+sqrt(D)). The fixed-order gradient program with k >= 4 has its separate O(Lambda A^4 sqrt(D)) allowance. Dimension enlargement inside the complete programs is included in their declared public-log/active-dimension factors.

The independent completed targets have variances v_B I and v_E I, whose sum is I, and means E B and E E. Conditional product coupling therefore gives the displayed own-mean target N(E F_G|Z,I). This uses only the finished outputs. It does not keep an old G,N, source force, or actual source-zero carrier as an additional observer of the law comparison.

The programs are finite Gaussian pushforwards. The subtraction in E and any signed deterministic arithmetic inside a compiler are vector operations, not signed probability weights. At source zero the imported mean programs return their literal known Gaussian rows; the independent positive shares combine to a row with covariance I. That actual row is not supplied by an optimal coupling and is not removed from the output.

Finally, Gaussian laws with identical unit covariance differ in W2 by the norm of their means. Combining the own-mean law comparison with the two target-bias terms and integrating over standard Z proves the claimed bound. The factor |||Z|+sqrt(D)||2 <= 2sqrt(D) is absorbed in Lambda.

## 6. Original work and numerical floors

The safe expanded VALUE count is

    Q_captured + N_out N_B + 2 N_out N_E
      + Q_known/numerical/replay.

Here N_B and N_E count every occurrence after the imported fixed-order programs are fully expanded. Every changed raw input regenerates every original leaf of that occurrence. One cannot count an E callback once while omitting all its later filter/clock/covariance/mean replays. Conversely there are no hidden canonical F1 or F2 calls in this raw source: their analytical effect has been replaced by known scalar coefficients and a fresh Gaussian row, with the target error proved in Section 1.

A raw baseline owns D private Gaussian coordinates. A raw shifted source or correction owns 2D. Those are per-occurrence counts; they are not the complete M_G tape dimension. All additional root blocks of both mean programs, their positive buffers, and numerical work retain their existing full ledger. G and N are shared across their own outer level and integrated when the corresponding completed law is used. Independent complete banks must not accidentally share a cached random source ancestor.

At fixed order the admitted clock/filter occurrence counts are public-log polynomials. With delta a fixed power of A, multiplication by N_out leaves inverse-A VALUE exponent zero. This is an oracle-call claim, not dimension-independent arithmetic cost, a growing-order cost recurrence, or a strong-mean Monte Carlo claim.

The first/adjoint implementation needs original HVPs only at recorded original VALUE sites. It follows the finite graph's actual affine maps and the captured-origin paths. No derivative of a saved HVP is a producer. Discarded primals must be replayed. All scalar versions, counts, positive shares, pads, and tolerances are frozen before caller differentiation.

For an absolute g-VALUE error nu at exact sites, the positive weights give raw errors at most nu for F_G or B and 2nu for E. Anchoring includes the separately computed origin error, hence safe 2nu and 4nu anchored allowances. Exact-key reuse is required for a literal numerical zero; separately perturbed copies would leave a priced zero defect. Errors in a2,b2 add the actual displacement term A(|Delta alpha||x_i|+|Delta b2||N|). Weight, node, Gaussian-row, finite-mode, and compiler errors retain their own original absolute allowances. None is divided by e_E, which can vanish.

The result certifies the anchored original problem. A physical finite-mode residual or a changed conditional normalization needs its separate original caller-law and numerical restoration. These follow their pinned imported contracts; they are not obtained by differentiating the small m3 law error.

## 7. Grade, radial member, and rescaling boundary

If epsilon <= d A and D >= A^(-4), then sqrt(D) >= A^(-2), so

    A epsilon(1+lambda) <= d(1+A)A^2
                         <= d(1+A)A^4 sqrt(D).

The quadrature term is of that grade when delta <= A^3. The consumer term already has that grade. Thus the dimension condition genuinely suffices for the stated bounded result.

The centered radial-inflation separator has the same original field

    g_D(x)=c A x+d A f_D(x), |f_D(x)|<=1,
    A=D^(-1/4), c=1/2, d=1/4.

It belongs to this subclass with lambda=c A and epsilon=d A. The location of its nonlinear shell does not affect the bounded-remainder argument. The current result keeps the full terminal response to alpha x-b2 N, including the nonlinear response to a dimension-sized Gaussian displacement. It does not Taylor truncate that response at its coherent mean. This explains why the counterexample to a centered covariance derivative does not contradict this positive bounded-source construction.

The final author text also correctly limits conditional rescaling. For

    f(y)=s[g(a+s y)-g(a)],

one has alpha=A s^2 as the new Hessian bound, lambda_f=lambda s^2, and only the safe certificate epsilon_f <= 2s epsilon. The original epsilon=O(A) does not uniformly imply epsilon_f=O(alpha) as s becomes small. Likewise D >= A^(-4) does not imply D >= alpha^(-4). These transformed parameters must be inserted into the general displayed error. The bounded result cannot be silently reused at every later conditional OU node.

## 8. Imports, provenance, and remaining obligations

This audit directly checked the relevant source ports against the following existing contracts. Full SHA256 entries and absolute locations are also in the JSON output and manifest.

1. Original LOW30 source, SHA `7767ac389d6b677eb9517b2585f7a24cb9f85ede47d099a653c635b7be750fb8`: `t30:lem:value-mean`, the actual source-zero carrier and VALUE-only program; `b27:compiler:mean`, its fixed-radius threshold, source-only VALUE interface, and fixed-order finite filters/serial counts.
2. `ENDPOINT-STEIN-PACKET-AND-POLYLOG-QUADRATURE.md`, SHA `4f0c182c36bf9dfbacc0f189899f689de52a2bb9bc37839e307eb76c4bda519f`: positive rule and uniform Hermite operator certificate.
3. `COMPLETE-MEAN-FIRST-ORDERING-AND-REENTRY-RECIPE.md`, SHA `a7cde258f1c86defbfdc681f473076293cfe65253479493c6ef5baa8b08a88a0`, and independent audit SHA `d9b4d45c276e72adc4b2590c1c17140bed2fe94dd2a6de5d8f7cc1b3a0a7d21a`: complete independent-branch composition, coisometry adapter, anchored-energy and caller obligations, and full replay cost.
4. Prior raw canonical m3 gate SHA `87b61a6a4ca03f81f7d910cdaa83b3783cd4792b3a428b89c8f0f629062e4467`: target identity and separation of own mean from target bias. Its unresolved generic second-decoupling claim is not used.
5. Prior nested mean SHA `b056e213ac28c60869727589d3ee6d1b2c53d6666908cdbe1961fa3121161ced`, independent audit SHA `c2fbd8b02e1d225226ef5cf822fafda31f3dfa05e505d747d0b8772cb9d2daf7`: the already admitted analogous source/consumer port. Its A^3 nonlinear decoupling theorem is not required to prove the new bounded-remainder comparison.

The new result is therefore complete as a **guarded component theorem relative to these pinned finite compilers**. This audit is not a new implementation or independent reproof of their entire recursive queue. The numerical diagnostics do not claim otherwise.

No mathematical blocker remains at the frozen source pin for the scoped bounded mean component. The general centered/resummed current, unrestricted conditional normalization, full covariance/feedback/restoration join, and any endpoint or all-order complexity claim remain separate obligations.
