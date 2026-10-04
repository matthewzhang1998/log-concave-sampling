# Independent audit: scalar randomized quarter-flow

Publication copy: nonmathematical provenance wording and/or local paths were sanitized. Original/public SHA-256 values are recorded in `INVENTORY.json`; historical source/audit pins refer to the original versions.

Date: 2026-10-04. New independent review, bounded to the displayed scalar construction.

Audited report: `../SCALAR-RANDOMIZED-QUARTER-FLOW-AND-ALL-LAYER-DEBTS.md`.

An exact reviewed snapshot is retained as `SCALAR-WEAK-FLOW-REVIEWED-14ff69c491ab.md`.

Audited report SHA256: `14ff69c491ab5cbcc3dfcfa196156cefced1c685f4905e5df663a6f8e48f8f27`.

## Verdict

**Scoped PASS.** The variance-preserving late keep, scalar posterior-law keep estimate, conditional final-clock current, scalar buffered `N_M^-3` bound, previous-layer strong-error exponents, and additive stored original-query census are valid under the stated small-heat and incumbent-moment assumptions. I found no material algebraic or analytic error in the displayed scalar ledger.

This is a scalar terminal law/query certificate, not a new admitted general-dimensional family or a change to `c_P`. The report correctly leaves dimension-uniform one-energy scaling, complete observed/retained-host ports, protected-width compatibility, and finite numerical/source restoration as separate obligations. The exact CDF graph's conditional centering must not be attributed without correction to a rounded CDF implementation. The posterior-law keep lemma also does not certify a joint law with every old record retained.

The accompanying independent script passes **757 deterministic checks**. These are algebraic and numerical fixtures; they supplement, and do not replace, the arguments and scope qualifications below.

## 1. Variance-preserving keep and the two different error notions

Write `delta=sigma²`, `kappa=sqrt(1-delta)`, and use the same independent `(z0,Z,K)` in

    Y_sigma = F(z0,kappa Z)+sigma K,
    F_full  = F(z0,kappa Z+sigma K).

The linear momentum carrier cancels exactly. Only the force integral sees the omitted momentum, giving normalized same-record error `C a sigma |K|`, hence physical `C a^(3/2) sigma`. No extra `sigma sqrt(a)` error belongs in this comparison. The combined momentum is exactly standard and independent of the initial state.

This comparison requires preserving the original `Z` alias; independently resampling a force-side copy destroys the identity. The `K` primitive must be absent from every force, cap, clock, and prior path computation. It may be sampled in advance but held unread, or drawn after the force computation.

The starting-position contraction remains `C a` because the same added `sigma K` cancels when two initial states are compared. Thus an incumbent physical `W2` error is multiplied by `C a`.

The stronger `C a sigma²` estimate is a **marginal posterior-law** statement for an exactly posterior-distributed initial state. It is not obtained by differentiating, squaring, or recentering the strong error.

## 2. Density proof of the scalar posterior-law keep estimate

The inverse phase-density argument is sound with only `U in C2`, `U(0)=U'(0)=0`, and `0<=U''<=a`. The vector field is `C1`, globally Lipschitz, divergence-free, and Hamiltonian. Its phase flow therefore preserves the posterior-times-standard-Gaussian density.

Under that invariant final phase density, changing the initial momentum variance from `1` to `kappa²` has exact likelihood ratio

    w_delta(p0) = kappa^-1 exp[-delta p0²/(2 kappa²)].

Integrating out the final standard momentum proves the report's `(D1)`. Inverse Duhamel bounds give

    |p0(z,p)-z| <= C a (|z|+|p|),
    |p0(z,p)²-z²| <= C a (z²+p²).

For `delta<=1/4`, the exponential is Lipschitz as a function of its nonnegative square argument with constant at most `C delta`. Consequently `|r_delta-g_delta|<=C a delta(1+z²)`. Since `r_delta` is normalized, integrating this bound also proves `|M_delta-1|<=C a delta`. This normalization estimate does not need an additional potential-value oracle.

For `nu_delta=pi_a g_delta/M_delta`,

    chi²(mu_delta || nu_delta)
      = integral pi_a M_delta/g_delta
          [r_delta-g_delta/M_delta]²
      <= C a² delta².

The inverse weight is at most a fixed constant times `exp(z²/6)`, and `U>=0` gives the required Gaussian integrability uniformly at small `a`. Normalizing constants are bounded away from zero because `0<=U(z)<=a z²/2`. The potential of `nu_delta` has second derivative at least `1/kappa²`, so the usual strongly log-concave transport-entropy inequality, together with `KL<=chi²`, yields `W2(mu_delta,nu_delta)<=C a delta`.

Convolution by the identical independent `sigma K` is contractive for `W2`. In `(D2)`, the displayed affine Gaussian argument is exactly the conditional law of `X` given `X+sigma K=z`, when `X~N(0,kappa²)`:

    X | (X+sigma K=z) = kappa² z+kappa sqrt(delta) G.

For `h=exp(-U)`, the derivatives used are only

    h'=-U' exp(-U),
    h''=[(U')²-U''] exp(-U).

The bounds `|h'|<=a|z|`, `|h''|<=a+a²z²` and integral Taylor expansion yield the stated `C a delta(1+z^4)` error. Division by the target density introduces only `exp(U(z))<=exp(a z²/2)`, which remains Gaussian-integrable under the fixed small-heat guard. This gives the second `chi²<=C a² delta²` and the remaining `W2<=C a delta`.

Thus `(WK)` is established, physically `C a^(3/2) sigma²`. Original `V` values appearing in the density proof are analytical quantities only; the executed graph need not query them.

The independent fixtures check `(D1)` algebraically on quadratic flows and numerically for the nonquadratic admissible potential

    U(z)=a[z²/4+0.2(1-cos z)],
    U''(z) in [0.3a,0.7a].

Both density normalizations and the two scaled chi-square estimates remain consistent with the proof.

## 3. Final conditional quadrature current and buffered gain

Freeze exactly the preceding record `R=(z0,Z,earlier complete clock banks)`. The last clock bank is fresh after that conditioning, and `K` is independent of it. With the exact Gaussian CDF, every last node is uniform within its stratum. Therefore the displayed residual has **exact conditional mean zero** relative to the integral defining `B_R`.

Integral Taylor expansion after averaging the independent keep proves `(J)`. Expanding its second derivative around the same endpoint gives the rank-two coefficient `Sigma(R)/2` and precisely the report's rank-three coefficient

    (1/2) integral_0^1 (1-t)²
      E[D³ phi'''(B_R+tD+sigma K) | R] dt.

The sign and `1/2` are correct. These are analytical test currents; neither an original third derivative nor an executed covariance oracle is present.

The full clock first is important. For a stratum of length `h`,

    dD/dG_j = -h F_R'(S_j) h phi_Gauss(G_j).

The first contributes two factors `h`, not one. Since `|F_R'|<=C_M a R_*`, summing squared columns gives `Lip(D|R)<=C_M a R_* N_M^-3/2`. Gaussian concentration/Poincare supplies the same conditional fixed-`p` moment scale. This remains true with the known `C1` filter because its **first**, rather than second, is uniformly bounded.

For completeness, interpolate `X_t=B_R+tD+sigma K`. Gaussian Riesz duality on the complete last clock bank gives the current

    d/dt E phi(X_t) = t E[(Riesz D)·(grad D) phi''(X_t)].

One integration by parts in the independent keep yields velocity norm at most

    t Lip(D|R) ||D||_(2|R) / sigma.

Integrating `t` gives the exact factor `1/(2 sigma)`. Mixing the conditional couplings in squared `W2` uses `E R_*^4`, not merely `E R_*²`, and proves normalized `C_M a²/(sigma N_M³)`. The physical factor is therefore exactly `a^(5/2)/(sigma N_M³)`.

This argument permits earlier records to remain in `R`. It does not permit an exterior observer to read the final clock bank: then its conditional rank-one current need not vanish. It also does not say the separate posterior keep lemma holds conditionally on arbitrary old records. Those two different uses of conditioning must remain separate in any later join.

## 4. Previous layers, smoothing, caps, and exponents

A clean way to avoid an unnecessary random-path supremum claim is to use

    e_k = sup_(t in [0,T]) ||q_actual^[k](t)-q_reference^[k](t)||_Lp.

Each new bank is independent of all earlier path errors. Evaluating an earlier error at a new random node is therefore bounded by its `sup_t Lp` norm. The bounded absolute kernel row masses and `Lip(b)<=a` give the recurrence

    e_k <= C a e_(k-1)
           + C a N_k^-3/2 + C a epsilon²,

with the scalar moment envelope restored in the constants. A layer-`j` quadrature residual starts with normalized size `a N_j^-3/2`; it passes through `M-j` further force applications before the final output. Its physical contribution is exactly

    C_M a^(M-j+3/2) N_j^-3/2,  j<M.

No assumption that nonlinear later layers preserve its mean is used. The nearest previous layer therefore costs `a^(5/2) N_(M-1)^-3/2`. Replacing this by the last-layer weak cost would be unjustified.

The continuous Picard tail is normalized `a^(M+1)`, hence physical `a^(M+3/2)`. A filter error has local integral magnitude `O(epsilon²)` and must pass through at least one additional force before the final unsmoothed `cos(s)` readout. Its leading physical price is `a^(5/2) epsilon²`; there is no such interior filter price when `M=1`.

The nonnegative cubic Hermite patch exists with the stated bounds. With `u=(r+epsilon)/(2epsilon)`, one valid patch is

    H_epsilon(r)
      = (-2u³+3u²) sin(epsilon)
        +(u³-u²) 2epsilon cos(epsilon), |r|<epsilon.

It joins `0` and `sin r` with their first derivatives, is nonnegative for small `epsilon`, has a uniformly bounded first, and has integrated error `O(epsilon²)`. Its second derivative may grow like `epsilon^-1`; no second filter derivative is used in the stated first/adjoint interface.

A concrete `C1` cap is identity up to `B`, equals `sign(z) B[1+t-t²/2]` for `t=(|z|-B)/B in [0,1]`, and is constant `sign(z) 3B/2` thereafter. Its first lies in `[0,1]`. This confirms an executable known cap with the stated force-path and global-first properties.

The cap tail error is incurred inside an `a`-weighted force, so the report's physical `a^(3/2)` coefficient is correct. A polylogarithmic cap radius is justified by an actual incumbent Gaussian-first/concentration return, not by posterior `W2` accuracy alone. If only finitely many moments are available, its power cost must be retained. Both limitations are stated in the report.

For `sigma=a^alpha`, choose `alpha=(R-3/2)/2`. The sufficient count exponents are

    beta_j = max(0, (2/3)[R-(M-j+3/2)]), j<M,
    beta_M = max(0, R/2-13/12).

The maximum earlier exponent, when an earlier layer exists, is `(2R-5)/3`. The report's formula `(B)` follows. At `R=5/2`, `M=1` has no previous-bank term and the count exponent is `1/6`. At `R=7/2`, `M=2` gives `2/3` for both banks. Above `7/2`, the earlier-layer strong certificate dominates. The exponent is a sufficient upper certificate for this execution, not an impossibility theorem for all randomized flow methods.

## 5. Actual C2 firsts, caller dependence, and stored-query cost

At a queried point, the only original differentiation is an action of `V''` on the actual incoming tangent or accumulated adjoint. Node motion differentiates the known CDF and kernel. It does not require `V'''`, a derivative of an HVP output, or a quantitative modulus for `V''`.

The caller derivative deserves an explicit calculation, since it is not the same as a derivative in `z0`. Fix heat and all schedules/caps before differentiating the entering physical caller `y`. For the exact mode,

    m_y = [1+a V''(m)]^-1 = 1+O(a).

For a finite mode iteration, bound its actual derivative directly from its own recurrence; do not differentiate its small mode VALUE error. If the incumbent has `X0_y=1+O(sqrt(a))`, then

    (z0)_y = (X0_y-m_y)/sqrt(a) = O(1).

The explicit mode dependence of normalized `b`, at fixed `q`, is

    partial_y b(q)
      = sqrt(a)[V''(m+sqrt(a)q)-V''(m)] m_y
      = O(sqrt(a)).

That is a bounded-Hessian estimate, **not** a small-Hessian-difference estimate. The finite paths have `q_y=O(1)`, so differentiating the physical emitted endpoint `m+sqrt(a)Y_hat` gives `1+O(a)`. This closes the actual caller calculation for the declared finite cap/filter graph under the same admitted incumbent first assumptions.

Its normalized old-input first is `O(a)`, fresh `Z` first is `kappa+O(a)`, and `K` first is exactly `sigma`. The leading row on the full tape is `(0_old,kappa,sigma,0_clocks)` and is coisometric. Each bank's complete clock first is `C a B N_j^-3/2`, with additional factors `a` when propagated through later layers. Under scalar moment concentration, `B` is an explicitly charged public-log factor. The standard symmetric leading original-force row conclusion follows in these first units; a complete native source/host return still does not follow from these firsts alone.

There is no hidden multiplication of the incumbent cost or the previous gradient count. At layer `k`, computing each `q^[k-1](S_(k,j))` reads cached earlier `f` values and known coefficients. It does not re-evaluate their original gradients. The gradient count is therefore one incumbent call plus `sum_j N_j` new force queries, shared anchor/mode/zero/guard work, exactly as printed. One directional forward or reverse sweep performs one original Hessian action per stored gradient after accumulating its outgoing direction. This does not mean that a dense Jacobian with many independent requested directions costs one sweep.

Known dense cross-layer algebra can cost `sum_j N_j N_(j-1)` and must be included in a full-arithmetic model. The actual Gaussian tape includes every clock coordinate. A future source that consumes that tape pays its true dimension and cannot remove clocks merely because their terminal law comparison integrated them out.

At owned Gaussian zero, `Phi(0)=1/2`, so zero-clock nodes are midpoints, not absent nodes. The complete graph's true zero trajectory and its caller path must be executed or explicitly specialized. No old origin or zero-cost claim follows from an `Lp` estimate.

## 6. Finite encoding and observer qualification

The conditional mean identity and `(J)` are exact for the real-valued graph using the exact known CDF. If a finite CDF approximation is substituted, uniform stratification is generally lost. That coherent conditional bias is a numerical restoration term. It does not receive the `N^-3` gain. The same warning applies to independently rounded cached anchors or inconsistent reused source versions.

A sound finite implementation freezes one common precision/version, bounds its VALUE effect through the actual graph, and supplies the actual bounded firsts of its known approximants. The report explicitly allocates numerical restoration, but does not construct every native finite encoding. Small VALUE error alone is not a first certificate. This is an implementation/return obligation, rather than a defect in the exact scalar estimate.

There are three separate boundaries:

1. The final-clock weak comparison preserves the preceding record but requires no exterior observer of the new last bank.
2. The posterior keep comparison uses the invariant **marginal** initial law and does not certify all observed old-source labels jointly.
3. The protected primitive has physical variance `a sigma²`, which shrinks with the target. Later inverse-width prices and required protected heat remain real join costs.

Retaining all the actual clock/cap/filter states while using an observer-free terminal marginal is not automatically a contradiction. It is a reminder that an externally accessible retained-state theorem requires additional proof, which this note does not claim.

## 7. Dimension fixture and nonlinear feedback obstruction

The quadratic shared-clock fixture `(G)` is correct. Conditional on the clock bank, all coordinates have the same scalar random variance. Its radial variance identity is an application of total variance:

    Var(radius)
      = E[V] Var(chi_d) + (E chi_d)² Var(sqrt(V)).

The reverse triangle inequality for standard deviations and radial contraction give the printed lower bound. For sufficiently large dimension, it detects `sqrt(d) a N^-3/2`. Thus one cannot uniformly replace the scalar conditional fourth-moment factor by `sqrt(d)` while keeping the `a²/(sigma N³)` coefficient. The independent fixtures make this comparison using even the coarser bounds `Var(chi_d)<=1` and `E chi_d>=sqrt(d-1)`.

The calibrated identity `(QC)` is also correct: set `u=sin²(s)`, then add and subtract `sqrt(1-u)b(z0)`. The added integral is `(pi/4)b(z0)`, and for linear `b` the bracket is constant. Its old-state derivative nevertheless contains the displayed Hessian difference divided by `sqrt(u)`. `C2` offers no uniform power modulus to remove that low-clock cost. The algebraic cancellation does not authorize differentiation as a high-order small error.

For an earlier centered residual `D`, the nonlinear defect `(F)` is the exact first-order integral remainder after `b'(q) E[D|past]=0`. The near-kink force satisfies the Hessian sandwich and, at `q=0` with `D=±d`, has mean defect

    a eta [sqrt(d²+epsilon²)-epsilon].

For `epsilon/d -> 0`, this approaches `eta a E|D|`, verifying that the safe prior-layer first-power estimate can be attained conditionally. This is not a marginal posterior-law lower bound: integrating a genuinely nondegenerate unobserved base `q` can create further cancellation. Any stronger reusable current must keep the actual joint law and later observer.

## 8. Adoption boundary and reproduction

No new general `c_P` family is certified by this audit. The successful conclusions are the scoped scalar law/error/query ledger and its explicit legal first-order graph. General-dimensional energy, all retained/observed source contracts, later protected-width compatibility, and full numerical/tape returns remain required before recursion.

Reproduce the fixtures with:

    OPENBLAS_NUM_THREADS=1 python check_scalar_weak_flow.py

The manifest pins the reviewed report, this audit, the script, its results, and the relevant preexisting audit/source files.
