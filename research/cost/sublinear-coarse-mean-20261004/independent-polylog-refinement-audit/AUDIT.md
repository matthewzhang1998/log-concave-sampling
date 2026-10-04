# Independent audit of the polylogarithmic scalar law refinement

2026-10-04 UTC.

## Verdict and scope

PASS for the stated bounded scalar mean-law construction in its exact Gaussian/uniform and known real-arithmetic model. There is no fatal gap in the executed pilot anchor, conditionally independent clipped batches, growing-rank Hermite bounds, failure-event coupling, additive provider census, or scalar retained-Z compatibility. The fresh private batches must be integrated out. The full pilot and the original exposed caller may be retained in a same-caller coupling.

This establishes a scalar buffered marginal-law transformer. It does not establish a smooth first/adjoint map, paired-caller port, multidimensional analogue, posterior-law identification, or complete recursive source theorem. A finite-precision port still needs the expressly reserved numerical work.

Primary source: `../POLYLOG-SCALAR-MEAN-LAW-REFINEMENT.md`.
Audited source SHA-256:

    f135f8583def97eafbee5f0f4a9b52c9465de61329577e8145366450a3dfa239

The relevant definition and identities in `../../endpoint-stein-quadrature-20261004/ENDPOINT-STEIN-PACKET-AND-POLYLOG-QUADRATURE.md` were also read to check Section 5. That companion's inspected SHA-256 was:

    4f0c182c36bf9dfbacc0f189899f689de52a2bb9bc37839e307eb76c4bda519f

The present audit checks use of the companion's literal packet, Lipschitz bound, and first quadrature moment; it does not substitute for a separate complete audit of its quadrature theorem. No external peer manuscript was accessed, used, or contacted. Neither source was modified by this audit.

## Executed anchor and independence

The pilot b is a sample average that is actually computed. It depends only on its M records and the already exposed caller. All subsequent M-record blocks are independent of the complete pilot and one another at that caller. Therefore the clipped X_i are independent and identically distributed conditional on the full pilot, not merely conditional on its average. Their common conditional mean is m = (mu_clip(b)-b)/sigma and satisfies |m| <= epsilon.

The comparison mean mu is absent from execution. The good event involving mu is used only for analysis. Clipping to the known interval [b-a,b+a] is executable. Its output lies between A_i and b, so it also stays inside the original source interval, even when an endpoint b +/- a falls outside that interval. The bound of 2S on the clipping displacement is consequently valid.

The factory guard follows directly from epsilon R = 1/16 and epsilon^2/2 <= 1/512 for R >= 1. In fact L > 3 and R = 8 sqrt(L), so the parameters have considerable additional margin.

If S = 0, W = b0 almost surely. The stated exact Gaussian-translate shortcut is correct; the large pilot is unnecessary in that special case.

## Hoeffding constants and clipping bias

A source record lies in an interval of length 2S. Hoeffding gives

    Pr(|average-mu| > t) <= 2 exp[-M t^2/(2S^2)].

Using t = a/4 yields exactly the source's denominator 32. Since a^2 = sigma^2/(256 R^2), the pilot exponent is

    M a^2/(32 S^2) = M/(8192 R^2 rho^2)
        >= [C0/8192] (1+rho^2)/rho^2 L
        >= 128 L

for C0 = 2^20 and rho > 0. Thus the weaker bound 2 exp(-8L) is safe.

On a good pilot, |A_i-b| > a implies |A_i-mu| > 3a/4. Multiplying this probability by the maximum displacement 2S gives

    |mu_clip(b)-mu| <= 4S exp[-9M a^2/(32 S^2)].

This proves the displayed bias bound with substantial slack. Both Hoeffding applications require full joint independence of the M original records, not a weaker pairwise-independence condition. The stipulated independent complete replays provide it.

## Rank uniformity is genuine

The absolute Hermite majorant is independent of K:

    |S_m(z)| <= exp(epsilon |z| + epsilon^2/2) - 1.

Dropping the minus one after squaring and completing the square gives the explicit bound

    ||1_(|Z|>R) S_m(Z)||_2^2
        <= 2 exp(3 epsilon^2) Phi_bar(R-2 epsilon).

For R >= 1 and epsilon <= 1/(16R), the Gaussian Chernoff bound then gives

    ||1_(|Z|>R) S_m(Z)||_2
        <= sqrt(2) exp(epsilon^2/2 + R epsilon - R^2/4)
        <= 2 exp(-R^2/8).

There is no rank-dependent polynomial constant. This is the crucial step that makes growing K legitimate.

Hermite orthogonality gives the exact remainder estimate quoted in the source. Since epsilon <= 1/16,

    remainder <= 2 * 16^(-K).

The cutoff correction is bounded by the same tail norm. The complete profile error is therefore at most

    2 * 16^(-K) + 4 exp(-R^2/8)
        <= 6 exp(-8L),

because K >= 8L and R^2 = 64L. This is a concrete universal bound, with no hidden C_K. The earlier maximal-common-density argument then gives an uncapped standardized W2 bound at most sqrt(12 sqrt(3)) exp(-4L).

The cap contributes at most sqrt(5) exp[-B log(3)/2]. Since B >= 8L and log(3) > 1, this is at most sqrt(5) exp(-4L). Thus the complete factory error relative to N(mu_clip(b),sigma^2) is below 7 sigma exp(-4L) for every pilot.

## Failure events and retaining the pilot

For every pilot, the standardized factory output V, including fallback, has E[V^2 | pilot] <= 3/2. Let N be an independent standard Gaussian. Using |b-mu| <= 2S,

    E[|b+sigma V-mu-sigma N|^2 | pilot]
        <= 2(b-mu)^2 + 2 sigma^2 E[V^2 | pilot] + sigma^2
        <= 8S^2 + 4 sigma^2.

Consequently the bad-pilot contribution to the squared coupling cost is at most

    sigma^2 (16 rho^2 + 8) exp(-8L).

On a good pilot, use the conditional law coupling and the clipping shift. A single coupling preserving the entire pilot exactly gives

    W2^2 / sigma^2
        <= [7 exp(-4L) + 4 rho exp(-8L)]^2
           + (16 rho^2 + 8) exp(-8L).

This inequality proves the requested C sigma delta result directly. Indeed exp(-L) = delta/(e+rho), and the displayed parameters even imply the conservative stronger bound 13 sigma delta^4. The extra slack is not needed for the advertised statement.

Retaining the pilot means that this integrated same-pilot coupling pays zero movement in the retained pilot coordinates. It does not mean that the output is within C sigma delta of N(mu,sigma^2) conditional on every individual bad pilot. That stronger uniform conditional statement is neither needed nor proved. Conditional comparison to N(mu_clip(b),sigma^2), by contrast, is uniform over all pilots.

The fresh K batches cannot be appended as retained observers after Hermite averaging. The previous scalar audit's retained-record counterexample still applies. The pilot is permitted because it was exposed before those fresh factors and the analysis explicitly pays its rare failure event.

## Literal query census and arithmetic

One execution draws M pilot records and K separate M-record factor batches. Each polynomial coefficient is a prefix product of the already generated X_i; computing degree k does not recursively invoke the provider k additional times. Hence the provider count is literally

    M + KM = (K+1)M.

With K = ceil(8L), R^2 = 64L, and M = ceil(C0(1+rho^2)R^2L), this is O((1+rho^2)L^3). A nonfree initial anchor producer remains additive. Rejection proposals do not replay the provider. There is no M^K or recursively branching sampling tree within this construction.

This arithmetic does not give a free recurrence for arbitrarily many later compositions. Each future source use would still have to pay this complete call count and establish its own ports. The source appropriately declines a complete hidden-source induction.

The correct additional known scalar operation count is O((1+rho^2)L^3+L^2), including record sums, clipping, prefix products, and at most B evaluations of a K-term Hermite recurrence. The author incorporated this clarification into the final statement. It is polynomial-logarithmic when rho is fixed or bounded by the public logarithms; it is not polynomial-logarithmic in an arbitrarily large rho itself.

The chosen C0 is extremely conservative. The polylogarithmic claim is an asymptotic order statement, not a demonstrated practical speedup. In the deterministic diagnostic cases, the literal-count coefficient relative to (1+rho^2)L^3 was as large as approximately 5.80e8. No enormous batch was actually simulated.

## Numerical model qualification

The exact theorem is in an explicitly ideal Gaussian/uniform and real-arithmetic model. Its normalization and exact rejection failure probability belong to that model. A finite-precision implementation needs separately certified proposal accuracy, source-record and anchor accuracy, and acceptance-threshold errors; those errors may also change the exact rejection mixture identity.

The author revised the numerical paragraph during this audit, and the final version was re-read in full. It now explicitly distinguishes the known factorial/polynomial intermediate scales from arbitrary exact source reals, permits rounding small products only at a selected absolute accuracy, and reserves provider/anchor precision, Gaussian generation, and propagated rounding. This resolves the earlier overstatement about finite exact bit lengths. Large absolute anchors may still require additional dynamic-range or cancellation precision under the provider's numerical model. The final note expressly leaves a full bit-complexity and numerical caller port open; no such port is inferred here.

A useful bound already established in the preceding independent audit is that a uniform acceptance-threshold error eta <= 1/6, with exact Gaussian proposals retained, contributes at most sigma sqrt(6B eta) in W2 for the B-capped program. Since B = O(L), eta can be selected exponentially small in L using O(L) threshold bits. This addresses that one numerical component but not the source, anchor, or Gaussian variate implementation, and not smooth first-action behavior.

## Retained Z and the endpoint packet

In the inspected scalar packet, H_Q(Z,G) = sum_i w_i g_M(r_i Z + sqrt(1-r_i^2)G), with positive weights summing to one and ||Dg_M|| <= A. Therefore its G-Lipschitz constant is at most A, uniformly in the retained Z. The actual finite node list can be executed at G=0 to obtain b0(Z) = -H_Q(Z,0); no conditional expectation is queried.

For W = -H_Q(Z,clip_(L0) G),

    |W-b0(Z)| <= A L0,
    |E[W | Z] + v_Q(Z)|
        <= A E|G-clip_(L0)G|
        <= 2A phi(L0).

The final bound is uniform in Z. This verifies the source's looser O(A exp(-L0^2/2)) clipping estimate. To obtain an arbitrary fixed heat order, the constant multiplying sqrt(log(1/A)) in L0 must depend on that requested order. Choosing L0 without that order-dependent coefficient would not automatically supply every order.

Independent fresh G blocks generate the literal conditional provider records. Retaining the original Z is allowed because it is exposed before the pilot and factor banks. The error coupling keeps Z fixed. Adding Z to the factory output thus yields the claimed approximation to Z-v_Q(Z)+sigma N. The fresh factor banks remain unretained.

If sigma is comparable to A and L0 = O(sqrt(log(1/A))), then rho = AL0/sigma is O(sqrt(log(1/A))). For fixed target order the transformer uses only public logarithmic factors beyond the complete packet cost. Its actual scale factor is still charged. Each W record pays the packet's whole declared node list and its other complete-provider ancestors; the b0(Z) anchor is charged separately. This verifies the scalar conditional compatibility only with the buffered v_Q target, not with the exact Stein field or posterior without their additional errors.

## Covariance debt remains explicit

For g(z) = Bz, the packet's exact first quadrature moment gives v_Q(Z) = BZ/2. The buffered target is (1-B/2)Z + sigma N, with variance

    (1-B/2)^2 + sigma^2.

The exact difference from the posterior variance is

    sigma^2 - B^2(3-B)/(4(1+B)).

This recovers the claimed second-order debt; with sigma comparable to A and B of order A, both contributions have second-order size and may only cancel under an additional tuned argument. The mean-law transformer itself does not correct them. In particular it cannot establish arbitrary posterior-law accuracy by silently identifying these distinct targets.

## Reproducible diagnostics

`check_polylog_refinement.py` checks deterministic analytic upper bounds and exact census formulas across 30 combinations, with rho from 1e-6 to 1e6 and delta from approximately 0.134 down to 1e-100. The output is in `diagnostics.json` and `diagnostics.out`.

- Minimum checked pilot exponent divided by L: approximately 128.
- All rank-uniform tail, remainder, guard, rejection-cap, and combined-error bounds passed.
- The largest checked explicit W2 upper bound divided by sigma delta was approximately 0.00012591.
- Twelve scalar quadratic covariance identities passed.

These are floating-point checks of the independent analytic bounds, not simulations or a replacement for the proof above. They expose the large conservative absolute constant while confirming the claimed parameter dependence.
