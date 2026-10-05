# Independent retained-carrier native-current audit

2026-10-05. This audit supplements `INDEPENDENT-MARKED-CURRENT-AUDIT.md`. It independently checks the revised p=1 finite source, conditional native substitution, target semantics, source anchoring, and execution bill. No author file was edited by the auditor.

## Verdict and snapshot

**CONDITIONAL PASS for the single-block, retained-carrier, linear-current port, under the explicitly imported raw-split theorem and its actual numerical guards.** This is more than a change-of-measure identity: the proposed input is a bounded, finite original-VALUE source with a separately proved full private/caller first and square-lift curl. It is not a retained-original-mark-law sampler, a small-energy output mark, an unrestricted-dimensional history theorem, or a native RAW-history source.

The final checked snapshot is `REVIEWED-NATIVE-SUPPLEMENT-FINAL.md`, SHA-256 `e5f24b4f1df648767856586cb1cdbbdbda7280088e2a594d169a50fb9b56c6ae`. Its count qualification is now correct: the original-gradient count **per source evaluation** is independent of selected Hermite order p, but the total native count can change through p-dependent source firsts, guards, and precision. This audit certifies only the explicitly compiled p=1 case. Higher-p native recompilation is not inferred. The earlier reviewed snapshot is preserved separately for traceability.

The main note's revised snapshot is `REVIEWED-CURRENT-REVISED.md`, SHA-256 `51ab911bbefef7554041e6cf885479dcbdb9e181320ce4a0b6170d9e293f263b`. Its new assumptions correctly require full-bank Gaussian independence, finite retained-label second moments, scalar-clock/buffer costs, and the single-block likelihood caveat.

## 1. Retuned finite producer and source first

The original pair-source proof gives a scaled secant error `C a_k delta_k A^3 n/e_k` and complete Gaussian first `C a_k(A/delta_k+A^2)` with clocks frozen. The proposed `delta_k=e_k sqrt(tau_k)`, with `a_k sqrt(tau_k)<=C`, therefore gives:

- Mean-field error `C A^3 n`
- Covariance error `C[A^5 n^(3/2)+A^6 n^2]`
- Complete source first `C A^(1/2)d_K^(-5/4)`

At `d_K~A^(10/9)`, the last quantity is `C A^(-8/9)`. The earlier dyadic/history errors are unchanged and dominate the new fixed-block covariance secant error. The source still executes fourteen original VALUES, with the retuned version propagated to all replays and floors. Its conditional-inner first remains `C A^2` because that argument never used the old small secant scale.

The Gaussian root and clock distinction is essential. The accepted revision retains S, containing selectors and clock values, through a native call, and treats only its Gaussian source tape Omega as private differentiable input. No derivative across a discrete inverse-CDF selector boundary is claimed. The uniform frozen-S first bound suffices for the conditional construction.

## 2. Bounded source and complete first/curl

For `E=H/v+H^2/(4v^2)`, direct multiplication verifies the stated trace formula

    tr(H^2) = [|F|^2 |F'|^2 + (F dot F')^2]/2.

Thus `N=Mbar (Gbar* E Gbar-tr E)/2` needs no executed tensor, expectation, or derivative oracle. Every matrix action and trace takes O(n) arithmetic.

With fixed v and frozen S, a unit directional derivative satisfies

    ||DH|| <= 2R L_F,
    ||DE|| <= (1+h)||DH||/v,
    image(DE) subset span(F,F',DF,DF').

The last inclusion proves rank(DE)<=4, including source-span rotation. The trace derivative is therefore bounded by `4||DE||`, not by n times that norm. The proposed first bound follows by the three product-rule paths: mark first, matrix first, and retained-carrier first. It is a valid conservative bound for both the full Gaussian private source and continuous retained callers.

The deterministic energy majorant `E_N=R_M e(B_G^2+2)/2` is valid uniformly in Y,G,S. A fixed known coisometric lift to the full private Gaussian input has first at most L_N and skew-Jacobian norm at most 2L_N. Padding when needed is an explicit extra Gaussian row and must remain in the root count. There is no assumption that a small VALUE energy produces a small first.

Hard radial caps provide a Lipschitz source with an almost-everywhere first. For a native theorem requiring C1, use the stated smooth radial caps with explicit inner identity/outer cap thresholds and fixed Lipschitz constants; increase deterministic bounds and choose Gaussian tail thresholds accordingly. This is a finite deterministic modification, not an appeal to extra smoothness of the original gradient.

## 3. Clipping and truncation

Conditioning on the complete old bank gives `E_G Q0^2=tr(E^2)/2<=e^2`. Since the G cap is radial, the trace terms cancel in `Q0-Q`; its remaining quadratic difference has magnitude at most `e|G|^2 1_(|G|>B_G)/2`. These facts prove (3.1), with no independence between M and H.

The likelihood's p=1 Hermite tail is `2sqrt(2)e^2` in L2 under the stated rank-two guard. Multiplication by M, then Cauchy-Schwarz, gives (3.2) for bounded callbacks. No global law comparison or sign claim about `1+Q0` is used.

The Gaussian mark-cap estimate is an integrated-standard-Y statement with the actual affine mean/covariance/offset profile. The executed source remains bounded for arbitrary retained Y, but the restoration of the unclipped source is not uniform over arbitrary unbounded Y. The supplement states that distinction correctly.

## 4. Native padding and exact origin restoration

Set `l=L_N/sqrt(u)`, `r=sqrt(l)`, `a_curl=2r`, and `mu=r`. For l<=1, the true normalized first l is no greater than r, and its curl is no greater than `r a_curl`. These are legitimate padded deterministic majorants.

Inserting `delta=E_N/(sqrt(u) r sqrt(d_private))` in the five imported physical raw-split terms gives exactly

    E_N[r^2+r mu+r a_curl+r^3/sqrt(mu)+r^3].

At mu=r and a_curl=2r, this is bounded by `C E_N r^2=C E_N L_N/sqrt(u)`, provided r<=1. The genuine-gradient baseline is identically zero under a valid coisometry, so its Gaussian component can be emitted exactly. The omitted baseline approximation term must be zero because of that actual construction, not merely ignored as inconvenient.

If source anchoring is required, execute `n0=N(Y,G,S;0)` and compile `N-n0`. The accepted revision correctly pays `2E_N` for its energy, leaves the private first/curl unchanged, and doubles the continuous caller-first majorant. Restore the same n0 after physical readout. The target is then the mean of the original N. Capture costs one full source replay per distinct retained key unless exact-key reuse is demonstrated. No inferred zero or source-energy bound replaces the actual origin.

These substitutions do not themselves prove every imported guard. The literal radius, relative-curl, source/caller, smoothing/filter, self-reserve, dimension, and numerical precision windows remain conditions of the theorem. They cannot be certified by differentiating the resulting Wasserstein error.

## 5. Correct current target and what the native bank consumes

The accepted native target is conditional on the same `(Y,G,S)`:

    R_N approximately N(mu_N(Y,G,S),u I),
    mu_N = E_Omega[N | Y,G,S].

After S is averaged, the reference is a mixture. There is no unjustified replacement by `N(E_S mu_N,uI)`. This correction is essential: treating sampled selectors as part of one globally differentiable Gaussian native source would be invalid.

The retained G must likewise stay outside the private source. Before G clipping, `E_G[M Q0 | Y,S,Omega]=0`, so an own-mean compiler integrating G would discard the desired current. The revised construction does not do that.

Let `M0` have the original conditional finite-mark law at `(Y,S)`, with a fresh mark bank independent of the entire native bank given `(Y,G,S)`. The executing output is `(sqrt(v)G, M0+R_N)`. For a bounded callback `a(Y,G)` reading no old mark/native bank or selector S, conditional W2 implies

    |E[R_N | Y,G,S]-mu_N(Y,G,S)| <= epsilon_native.

The reference Gaussian marker noise has zero conditional mean. Pairing with a and averaging S therefore gives precisely the extra epsilon_native term in (5.2). Adding the independently proved clipping and Hermite truncation errors gives the stated current contract.

This only uses linearity in the mark. The actual marker noise need not cancel pathwise and need not be small. Nonlinear reuse of J, reading consumed roots, subtracting a completed-law carrier, or claiming the original joint mark law would require a different theorem. The supplement excludes those promotions.

## 6. Grade arithmetic and block scope

At fixed physical block size, the three first paths have A exponents 5, 37/9, and 7, all divided by v. Thus the dominant bound is `L_N<=Lambda_b A^(37/9)/v`, while `E_N<=Lambda_b A^7/v`. The resulting native error is `Lambda_b A^(100/9)/(v^2 sqrt(u))`. The Hermite truncation is `C_b A^11/v^2`, with the displayed physical-size profile.

For u comparable to v and `v>=c A^2 K^2`,

    r <= Lambda_b A^(37/18) v^(-3/4)
      <= Lambda_b A^(5/9) K^(-3/2).

At v of order A^2, the native error exponent is 55/9 and the truncation exponent is 7. These calculations are correct. They give a nonempty fixed-block/public-log small-A window, not a dimension-uniform theorem for b=D or a proof for arbitrary growing logarithmic demands.

Only block-local current/readout contracts can be combined by Euclidean root-sum-square error. Arbitrary globally coupled callbacks require a genuine global joint construction. Multiplying rank-two likelihoods across many blocks does not preserve the rank-two energy constant.

## 7. Costs and final boundary

One conditional N source call has at most nineteen original VALUES and the full old pair/mark Gaussian bank. Because G is retained, the pair's Gaussian coordinate count is at most 16n, before genuinely new mark/padding roots; selectors and clock randomness remain separately charged. M0 costs five more VALUES when unavailable. Every native occurrence owns a complete fresh Gaussian source tape at the same retained S and its own native auxiliaries.

The honest cost expression is the stated imported recurrence `5+19 N_raw+Q_capture+Q_restore+Q_numeric`, with the actual complete Gaussian dimension and scalar clock/root bill. Source origin replays, caller/source-version changes, and precision amplification must be included. This is a finite conditional mathematical construction based on the imported compiler, not an executed numerical implementation for an unspecified original gradient.

The revision repairs the three initially material port issues: unbounded exact likelihood is not the native source; G is retained instead of averaged away; and stochastic selector/clock modes are frozen through the conditional native call. Exact origin capture/restoration is now charged. The remaining true-history base/singleton and genuine Delta3 source obligations are unchanged.

The author reported a new executable module and 47 author-side tests after this mathematical review was complete. Those tests and implementation are not independently certified by this audit; its independently executed diagnostics are the ones listed in the companion audit and `current_resummation_checks.json`.
