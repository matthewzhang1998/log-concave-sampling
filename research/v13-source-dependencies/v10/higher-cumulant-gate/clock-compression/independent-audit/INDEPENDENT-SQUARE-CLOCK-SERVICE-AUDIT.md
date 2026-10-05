# Independent audit: positive square-clock VALUE covariance service

2026-10-04.

## Verdict and exact pin

**PASS for the stated fixed-positive-buffer, cubic covariance service, including its owned-root absorption and shared-incoming-p aggregation.** This imports the admitted finite VALUE square action and proves its specialization; it does not assert a new all-order action or a complete endpoint sampler.

Final audited file: `../POSITIVE-SQUARE-CLOCK-COVARIANCE-SERVICE.md`.

SHA256: `108af2124e3da0e2a04af92e9f0403bfce2f82f7e2912ba6b7e1975ef1c8b18b`.

The earlier pinned version `6de7f54610c60049c3323e2eaf7823852e208c96550d170fd923930ec0d19e99` was reviewed first. The final version adds the explicit basic query/tape census and the terminal-only no-gradient-reentry qualification. Those additions were checked and are correct.

The admitted final estimate is

    || W2(Law(R|Z), N(0,v0 I+Cov(H_cont|Z))) ||_L2(Z)
       <= Lambda A^3 sqrt(D) + C delta A^2 sqrt(D)
          + restored filter/numerical floors,

where Z is the standard Gaussian endpoint for the quadrature restoration, `v0=eta^2+zeta^2`, eta and zeta are fixed positive constants, and all actual small-radius/filter guards are imposed. The service comparison to `v0 I+C_Q(Z)` is uniform in Z; only the final analytical clock replacement uses L2 of the standard endpoint.

## 1. Analytical square and finite quadrature

Sections 8–9 of `INDEPENDENT-CLOCK-GRADIENT-AUDIT.md` provide the independent derivation. With `B=Dv=integral r P_r[Dg] dr`, the exact identity is

    Cov(H_cont|Z)=2 integral q P_q[B^2](Z)dq.

The positive finite rule supplies

    B_Q=sum_i w_i r_i P_(r_i)Dg,
    C_Q=2 sum_j v_j q_j P_(q_j)[B_Q^2],
    0<=B_Q<=A I/2, 0<=C_Q<=A^2 I/4,
    ||C_Q-C_cont||_L2(HS) <= (3/2) delta A^2 sqrt(D).

Noncommuting Hessians are handled by `B_Q^2-B^2=B_Q(B_Q-B)+(B_Q-B)B`; their order is not interchanged. P, Dg, B_Q and C_Q are analytical targets only.

## 2. Genuine D-dimensional source with all anchors

For each captured X, execute

    f_X(u)=sum_i (w_i r_i/c_i)[g(r_i X+c_i u)-g(r_i X)].

Its displayed analytical potential differentiates to f_X. The actual private derivative is

    D_u f_X=sum_i w_i r_i Dg(r_i X+c_i u),
    0<=D_u f_X<=A I/2.

Thus `ell=A/2` is a valid declared radius, `|f_X(u)|<=ell|u|`, the anchored L2 energy is at most `ell sqrt(D)`, and Gaussian Poincare gives the same upper bound for its centered L2 energy e. These bounds are uniform in X. Its origin is exactly zero when each anchor is reused at the identical query site.

The actual caller derivative includes both the shifted and anchor Hessians:

    D_X f_X=sum_i (w_i r_i^2/c_i)
                     [Dg(r_i X+c_i u)-Dg(r_i X)].

The PSD Hessian sandwich implies that each bracket has operator norm at most A. For every stated dyadic rule,

    sum_i w_i r_i^2/c_i <= sum_i w_i/c_i <= 1+sqrt(2).

Indeed `c_i>=sqrt(1-r_i)`, the panel `[a,2a]` has mass a and contributes at most sqrt(a), and the terminal midpoint of panel length h contributes at most sqrt(2h). Summing gives at most `1+sqrt(2)-sqrt(h)`. The final note's conservative constant four is valid, uniformly in Gauss order and panel count. The absolute VALUE coefficient sum `sum w_i r_i/c_i` is bounded the same way.

Finally, at the exact Gaussian law of u,

    E D_u f_X=B_Q(X).

This derivative is used only to identify the target of a VALUE action, not to execute it. Full joint gradient status in (X,u) is neither needed nor claimed.

## 3. Exact imported action and owned-root firsts

The relevant original LOW30 passage is `t30:lem:actions`, approximately lines 232–311, with `t30:eq:first-response` and `t30:eq:square-value`. Its zero-clock square variant has target `(s0 E Df)^2`, is pointwise odd in incoming p, and supplies

    calibration <=Lambda ell e mu,
    action energy <=Lambda ell e,
    complete p/private first <=Lambda ell^2/sqrt(mu),
    original-parameter first <=Lambda ell L_X/sqrt(mu).

For the present source, `ell<=A/2`, `e<=A sqrt(D)/2`, `L_X<=C A`. After dividing the action by fixed s0 squared, the bounds are respectively `Lambda A^2 mu sqrt(D)`, `Lambda A^2 sqrt(D)`, and `Lambda A^2/sqrt(mu)`.

For each q_j, `X_j=q_j Z+sqrt(1-q_j^2)G_j` is created before its complete action-private bank. Its G_j is now an owned private root of the final action. Differentiating through X_j costs the imported original-parameter first multiplied by a row of norm at most one. This supplies the claimed complete first in G_j and captured Z; it does not assume that the conditional law error can be differentiated.

Every same-X anchor is part of that owned caller graph. If G_j changes, those anchors change and must be differentiated/replayed. They are not globally deterministic cached values.

## 4. Same incoming p and positive aggregation

Let `a_j=2v_j q_j`, so `a_j>0` and `sum a_j=1`. Use one common p, independent owned G_j and independent complete response banks, and form

    D(p)=sum_j a_j C_sq(f_(X_j);p)/s0^2.

For each fixed Z, average each imported conditional calibration over its owned G_j, use Jensen, then sum with positive weights. This yields

    ||E_private D(p)-C_Q(Z)p||_L2(p)
         <=Lambda A^2 mu sqrt(D)+absolute floors.

The triangle inequality gives the same weighted bound for the integrated energy. The actual chain rule, again using total weight one, gives complete first `Lambda A^2/sqrt(mu)` in the shared p and all owned/private roots. The shared-p derivative is the sum of its literal paths. It is not an independent-bank variance calculation and receives no unproved sample-count improvement.

Using separate p_j and summing their actions would not give the displayed target matrix times one retained incoming Gaussian. The final construction correctly uses the same p. Its independence claims concern only the private banks conditional on p and Z.

The aggregate D need not be a gradient and is never passed back to a gradient compiler. The terminal buffer argument uses only its actual Lipschitz radius, energy, calibration and oddness.

## 5. Single positive reserve: complete proof of the error ledger

Fix eta,zeta>0 independently of A and D, `eta^2+zeta^2=v0`, and use an independent standard z. Let `c=1/(2eta)` and execute

    R=eta p+cD(p)+zeta z.

### 5.1 Absorb every owned root before the private comparison

At fixed p,Z let `X_p=D(p)-E_private D(p)`. It is a centered function of the complete standard Gaussian private tape, including every G_j. Its private first is at most `L=Lambda A^2/sqrt(mu)`. Centering reduces energy, so

    (E_p ||X_p||_L2(private)^2)^(1/2)
         <=||D||_L2(p,private)<=Lambda A^2 sqrt(D).

For a centered Gaussian-Lipschitz map X=h(W), a Gaussian Riesz/Stein matrix tau satisfies `||tau||_L2(HS)<=Lip(h)||X||_2`. Along `zeta z+t c X`, use its Stein identity and then Gaussian integration by parts in the independent z. The resulting continuity velocity has L2 norm at most

    t c^2 Lip(h)||X||_2/zeta.

Integrating from zero to one gives

    W2(Law(zeta z+cX),Law(zeta z))
        <=c^2 Lip(h)||X||_2/(2zeta).

Apply this at each p,Z, retain the translation `eta p+cE_private D(p)`, and square-integrate over p. It gives precisely the service's private-action allowance

    Lambda A^4 mu^(-1/2) sqrt(D).

The owned G_j cannot be left frozen in this comparison and later read as an extra output. The final note integrates them into the complete action, as required.

### 5.2 Calibration and the positive quadratic remainder

The conditional mean calibration costs `Lambda A^2 mu sqrt(D)`, by a same-p coupling and the fixed multiplier c. Its reference is

    (eta I+C_Q(Z)/(2eta))p+zeta z.

Since C_Q is symmetric, its exact covariance is

    v0 I+C_Q(Z)+C_Q(Z)^2/(4eta^2).

Both this covariance and `v0 I+C_Q(Z)` have the fixed gap v0. The final positive term has HS norm at most `C_eta A^4 sqrt(D)` because `||C_Q||op<=A^2/4` and `||C_Q||HS<=A^2 sqrt(D)/4`. Gapped Gaussian square-root stability therefore costs `C_(eta,v0) A^4 sqrt(D)`.

The resulting conditional bound is

    Lambda [A^2 mu+A^4(1+mu^(-1/2))]sqrt(D)
        +restored absolute floors.

For `mu=A`, `0<A<=1`, this is `Lambda A^3 sqrt(D)`. The constants retain their fixed eta/zeta dependence; a shrinking reserve cannot be substituted without repricing those factors. The source already requires a fixed positive variance share.

### 5.3 Continuous covariance restoration

At every z, both analytical covariance matrices have the same positive base v0. Apply the pointwise Gaussian-root inequality and then the separately proved L2(standard Z;HS) clock bound. This costs `C_v0 delta A^2 sqrt(D)`. The final note correctly distinguishes this L2 endpoint result from its uniform-Z service calibration.

## 6. Exact mean, carrier, and zero

The finite first filter is pointwise p-odd. Its value changes sign under p negation while all owned G_j, caller anchors and response roots remain fixed. The nested outer signed difference therefore also changes sign. Positive summation preserves this oddness.

Since p is independent standard Gaussian and the other roots are independent of p, `E[R|Z]=0` in the exact finite real-valued program. Coherent numerical symmetry must be preserved or its error charged, as the source states.

Setting the nonlinear source to zero makes D identically zero. The known source-zero carrier is exactly

    eta p+zeta z.

It is the actual carrier, not a Gaussian chosen in a comparison coupling. With p=z=0 the reserve is exactly zero even for nonzero endpoint Z and arbitrary other private roots. The residual first after removing the carrier is `Lambda A^2/sqrt(mu)`; the total program still has its fixed carrier first.

## 7. Literal original queries and Gaussian tape

For the unexpanded printed zero-clock action with K_f first-filter nodes:

- The inner first filter uses `2K_f` raw f VALUES, all on its shared z1.
- The two outer responses use four raw f VALUES in total; both signs retain the same I,z0,z2.
- Thus `N_sq=2K_f+4` before any separately added numerical restoration.

One complete bank has N_r original same-X anchors and `N_r N_sq` shifted original g VALUES when anchors are captured once under the exact caller/source/version key. Across N_q banks the basic cached bill is

    N_r N_q (2K_f+5).

Without relying on this cache, `2 N_r N_q N_sq` is a safe upper bound, exactly as stated in the final source. External global mode/anchor work, numerical restoration, discarded-primal replay and requested HVP/adjoint sweeps are added separately. Every recorded original first site is a legal derivative sweep; no HVP is a producer.

The basic tape contains one owned G_j plus z1,z0,z2 per q bank, one shared p, and one final z. Its dimension is

    D(4N_q+2),

before adding the external host and any numerical-restoration tapes. Inner deterministic clock nodes do not add Gaussian roots. Conditional source anchors may be cached within each bank, but not across changed X_j or independently complete banks.

## 8. Independent checks and qualifications

`check_square_clock_value_action.py` executes the literal finite zero-clock VALUE responses and passes **675 assertions**, with seed 620041005. It checks the finite Chebyshev filter identities, uniform dyadic caller/value coefficients, original-query counters, pointwise oddness, zero at p=0, source-zero carrier, and non-diagonal matrix quadratic outputs in D=1,3,7. The quadratic source produces `B^2 p/4` exactly to the numerical tolerance. Its reserve covariance is `v0 I+B^2/4+B^4/(64eta^2)`, including the paid fourth-degree term. Nonlinear fixtures are executed to check parity, records and counts, not to infer the universal calibration bound from Monte Carlo.

Together with the two earlier independent checkers, the present folder has **1,828 passing assertions**: 1,060 protected-lift/diagonal checks, 93 conditional-covariance-square checks, and 675 literal VALUE-action checks. The earlier source remains independently valid; the final square service is the current implementation result.

Remaining limits are explicit and material:

1. The high-order finite-action calibration is imported from its stated LOW30 port, with its true finite guards and floors. Finite tests do not re-prove that port.
2. eta and zeta are fixed positive shares. No small-variance covariance service is admitted by suppressing inverse reserve factors.
3. The owned action roots are integrated in a terminal complete-output comparison. No additional observer of those roots is licensed afterward.
4. Original physical-caller and finite-mode restoration remain the host's separately charged graph obligations.
5. This is a cubic covariance service. It does not supply arbitrary-order covariance correction, a corrected nonlinear rank-one drift, a full endpoint-law composition, or an all-order cost recurrence.
