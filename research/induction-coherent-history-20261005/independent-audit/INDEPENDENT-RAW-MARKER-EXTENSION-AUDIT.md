# Independent audit: direct raw marker extension

2026-10-05. This checks the author's explicit final proposal to emit the directly executed `Mnew=M+N`, rather than the auxiliary-noise native output. It uses the finite p=1 source N audited in `INDEPENDENT-RETAINED-CARRIER-NATIVE-AUDIT.md`. No prior file was changed.

## Verdict

**PASS for preservation of the stated local finite raw-marker source bounds and its linear-current comparison**, assuming the original M already has the qualified complete-row lift, first, curl, and actual conditional/integrated energy profile. The proposal does not construct such an old mark if it is missing. It does not give equality of the old mark law, a true-history RAW approximation, a nonlinear marker-reuse estimate, or an arbitrary-depth recurrence.

The output uses the original unclipped M in its baseline and the clipped correction only:

    Mnew = M + Mbar Q,
    Q = (Gbar* E Gbar-tr E)/2,
    X = sqrt(v)G.

The entire graph is executed directly from original VALUES and Gaussian roots. There is no native mean LAW output, auxiliary marker noise, or carrier subtraction in this construction.

## 1. Old curl survives; clipping affects only the correction

Let the old physical mark have a qualified coisometry P_old from its complete Gaussian input bank, and let its square lift be `P_old* M`. Enlarge the Gaussian bank to include every added pair-source root and the fresh G. Extend the coisometry by zero columns:

    P_ext = [P_old 0],   P_ext P_ext* = I.

The old square-lift Jacobian is its previous Jacobian with zero extra rows/columns, so its first and full skew-Jacobian bounds are unchanged. This statement includes every previously owned root; it does not omit auxiliary blocks from the old curl.

For the correction, the already proved full first bound L_N gives

    Lip(P_ext* N) <= L_N,
    ||D(P_ext* N)-D(P_ext* N)*|| <= 2L_N.

This includes all mixed old-bank/new-G blocks. Hence

    Lip(P_ext* Mnew) <= L_M + L_N,
    curl(P_ext* Mnew) <= kappa_M + 2L_N.

There is no need to claim that radial clipping preserves the old mark's sharper curl. It generally need not. Here the old baseline M is **not clipped**, and all derivatives of its clipped copy inside N are already charged through the conservative L_N bound.

## 2. A direct clipping-safe energy bound

Assume both radial caps are contractions: `|Mbar|<=|M|` and `Gbar=s(G)G` with `0<=s(G)<=1`. Ordinary radial projections and suitably chosen smooth radial caps have these properties. Conditional on the complete old source bank, E is fixed, rank(E)<=2, and `||E||op<=e`. Gaussian quadratic isometry gives

    E_G[(G* E G)^2] = 2 tr(E^2)+(tr E)^2 <= 8e^2,
    |tr E| <= 2e.

Since `|s(G)^2 G* E G|<=|G* E G|` pointwise, Minkowski yields

    ||Q||_(2|old bank) <= (sqrt(2)+1)e.

Multiplying by the fixed conditional Mbar and integrating gives the useful uniform bound

    ||N||_2 <= (sqrt(2)+1)e ||M||_2,
    ||Mnew||_2 <= [1+(sqrt(2)+1)e] ||M||_2.

No clipping tail, mark fourth moment, independence of M and H, or mark supremum is needed for this energy conclusion. The Gaussian must be fresh relative to the entire old bank. The same bound holds conditionally on any retained labels for which that freshness and the original mark energy are valid.

Thus an old energy `O(A^3 sqrt(n))` stays at that order under the rank-two numerical guard. Preserve its actual retained-Y/S profile: an integrated-standard-Y bound does not become a uniform bound at arbitrary retained Y. If constants of a chosen smooth cap exceed one rather than making it a contraction, charge those constants explicitly.

## 3. First and curl grades

For the retuned finite source at fixed block size,

    L_N <= Lambda_b A^(37/9)/v.

If `v>=c A^2 K^2`, then

    L_N <= C Lambda_b A^(19/9)/K^2.

Consequently the needed condition for `L_N=O(A^2)` is the explicit public-log window

    Lambda_b A^(1/9)/K^2 <= C_0.

At fixed b and fixed logarithmic exponents this holds for sufficiently small A, but the logarithmic factor cannot be discarded outside that window. Under it, the old assumptions `L_M=O(A)` and `kappa_M=O(A^2)` give the same first and curl orders for Mnew. The energy bound above gives the same old `O(A^3 sqrt(n))` energy order.

All continuous caller derivatives use their actual source firsts. The selectors/clock modes S remain frozen exactly as in the accepted supplement. Their inverse-CDF sampling is not silently included in a claimed globally differentiable Gaussian source. Varying v, cap thresholds, or clock values requires their separately declared paths/replay rules.

## 4. Current error and cost

The previous likelihood identity and p=1 Hermite truncation apply directly to `(X,Mnew)`:

    |E[Mnew dot a(Y,G)]-E[M dot a(Y,Z/sqrt(v))]|
      <= 2sqrt(2)e^2 ||M||_2 + clipping_error,

for the specified norm-one callbacks, unread old private bank, and actual retained-label law. The clipping error is exactly the supplement's (3.1). There is no epsilon_native term because R_N is not used.

For the local finite mark, `e=O_b(A^4/v)` gives the current truncation order `O_b(A^11/v^2 sqrt(n))`, with its explicit cap/row profiles. At v of order A^2 this has grade 7. It does not supply a Wasserstein comparison retaining the old mark law and does not control arbitrary nonlinear functions of Mnew.

The mark M and its clipped copy share their executed ancestors. One draw therefore needs no extra original VALUES beyond the same finite at-most-nineteen-call pair/mark source, before exact-key reuse. The new G is part of its complete raw Gaussian bank and retains its existing coordinate-root charge. Extra mark roots, scalar clock randomness, basis transforms where relevant, and all numerical floors remain charged. No native occurrence recurrence is needed for this direct output.

This proves a valid **local source-class preservation statement plus a linear-current approximation**. Further iterations would need explicit fresh-carrier ownership, accumulated errors/energy/constants, and the requisite finite mark suppliers. None of those missing true-history or arbitrary-depth obligations is removed by this one-step audit.
