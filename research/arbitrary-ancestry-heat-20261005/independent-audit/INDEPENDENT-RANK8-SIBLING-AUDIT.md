# Independent audit: rank-eight sibling heat obstruction

2026-10-05. Local mathematical audit and runnable checks. No external messages or uploads.

## Verdict

**The proposed rank-eight history, Gaussian extraction obstruction, original-gradient source family, and growing one-node covariance are correct.** The claim must remain a failure of cutoff-uniform nodewise or absolute/squared-weight history budgets and of the blind old native normalization. It is not a lower bound for an arbitrary complete shared-bank clock/cumulant sum, and is not an all-order impossibility theorem.

There is a useful strengthening: a positive part of the divergence occupies a fixed, non-endpoint portion of the exact symmetric OU bridge. It comes from the difference of the siblings' independent coarse noises, so smoothing the broad common structural parent cannot remove it.

The implementation and numerical results are in:

- `check_rank8_sibling_obstruction.py`
- `rank8-selected-history.json`
- `rank8-sibling-audit-checks.json`

Run `python check_rank8_sibling_obstruction.py` from this directory. The script reads the supplied analytical AST generator, builds only the selected branch, and verifies finite-generator membership through rank six. Legality at ranks seven and eight is checked by the same exact recurrence construction, without enumerating all histories. Numerical tests are sanity checks; the arguments below provide the proof.

## 1. Exact legal history and clocks

Use the ordered split `(n-1,1)` whenever attaching a fresh leaf. Start from the rank-two force pair, hit vertex 1 for the first three attachments, and hit vertex 2 for the last three. At each step choose that force in the product-rule derivative and give the fresh leaf the other end of the new edge.

This is literally a branch of the supplied recurrence. Its coefficient is

    2 · 3 · 4 · 5 · 6 · 7 · 8 = 8! = 40320.

Its graph has two adjacent degree-four vertices, with three degree-one leaves on each side. Thus the degree multiset is `(4,4,1,1,1,1,1,1)`, with seven edges and six excess derivatives.

The seven structural resolvents all finish at index 8. The two high primitives have index 5, while the six leaf primitives have index 2. Total clock count is fifteen. Both high vertices have exactly the same seven structural ancestors and hence the same deepest structural caller X. The occurrence labels and derivative-edge labels remain distinct.

For a concrete legal scalar bank take all seven structural correlations equal to 1/2 and the retained caller equal to zero. Then

    X ~ N(0,v),    v = 1 - 2^(-14) > 0.

Take all six ordinary primitive correlations equal to 1/2. Their private variances are 3/8. These choices are interior clock choices. One can instead use fixed positive interior quadrature nodes; their finitely many masses change only a positive multiplicative constant.

No claim is made that rank eight is the earliest possible obstruction. Its virtue here is a symmetric pair of identical C3 centers and a particularly transparent covariance formula.

## 2. Exact Gaussian geometry

Let `x = k^(-2)` and `r = sqrt(1-2x)`, with k sufficiently large. The high queries, before the original private shields are averaged, have coarse centers

    q_i = r X + sqrt(x) W_i,     i=1,2,

where W1,W2 are independent standard Gaussians independent of X. Each query also has its separate original `sqrt(x)` private shield. Consequently their structural-plus-coarse covariance is exactly

    Gamma = v r² 11^T + x I_2.

At OU bridge parameter rho, write `h²=1-rho²`. If an exact independent diagonal shield extraction uses nonnegative increments delta1,delta2, the residual must satisfy

    h² Gamma - diag(delta1,delta2) >= 0.

Testing the difference vector `(1,-1)` gives

    delta1 + delta2 <= 2 h² x.                              (A1)

This condition remains necessary when these queries are part of the full eight-force covariance: every principal submatrix of the full residual must be PSD. The original private variances therefore become at most

    t1² + t2² <= 2(1+h²)x <= 4x,
    t_i² <= (1+2h²)x <= 3x.                               (A2)

The equal extraction `delta1=delta2=h²x` is feasible and saturates (A1). A slightly larger equal extraction is not feasible. If only one shield is enlarged, its exact maximal increment is

    delta1_max = h² x (2 v r²+x)/(v r²+x) < 2h²x.

Broad common structural variance is a common-mode resource. It does not provide an order-one *independent* shield to either sibling while maintaining the original query covariance. Arbitrary correlated heat or a different regrouping is not excluded by (A1).

## 3. Original-gradient source and exact old coefficient

Take

    g_k(t) = a t + (b/k) sin(kt),
    0 < b < a,    a+b <= 1.

Then `g_k(0)=0` and

    0 < a-b <= g_k'(t) = a+b cos(kt) <= a+b <= 1.

Thus this is a smooth admissible scalar original gradient with a uniform bounded positive Hessian; no arbitrary tensor field has been substituted. Its fourth derivative is

    g_k''''(t) = b k³ sin(kt).

Private heat of variance x at either high vertex gives exactly

    E[g_k''''(q_i+sqrt(x) Z_i)]
      = b k³ exp(-1/2) sin(k r X+W_i).                    (A3)

Its normalized C3 target is `x^(3/2)` times this derivative, and is uniformly bounded. For each ordinary leaf with private variance 3/8,

    E[g_k'(q+sqrt(3/8) Z)]
      = a + b exp(-3k²/16) cos(kq).

If G_k is the product of all six leaf factors, irrespective of their correlated structural/coarse callers,

    |G_k-a^6| <= 6b exp(-3k²/16).                         (A4)

This is uniform, not an independence approximation.

Let B be the product of the other thirteen clock masses, including any fixed nonzero bounded factors, and retain the positive history coefficient C=40320. If each high primitive has mass w_k, the literal one-node coefficient is

    L_k = C B w_k² b² k^6 exp(-1)
          sin(T+W1) sin(T+W2) G_k,
    T = k r X,    Var(T)=V=k²r²v.                       (A5)

For `w_k/x -> c>0`, this becomes

    L_k = [C B c² b² a^6 exp(-1)+o(1)] k² F_k,
    F_k=sin(T+W1) sin(T+W2),

in bounded-error and L2 senses after division by k². The exact variance is

    Var(F_k)
      = 1/4 + exp(-4)(1+exp(-8V))/8
            - exp(-2)(1+exp(-4V))/4.                    (A6)

Hence

    Var(F_k) -> 0.2184556340519386 > 0,

and `Var(L_k)` is a positive constant times k^4. Restoring eight service amplitudes gives the corresponding `alpha^16 k^4` covariance scale, in the stated dimensionless A=1 normalization.

An elementary derivation of (A6) uses

    F_k = [cos(D)-cos(S)]/2,
    D=W1-W2,       S=2T+W1+W2.

D and S are independent centered Gaussians of variances 2 and `4V+2`. In particular `E F_k = exp(-1)(1-exp(-2V))/2`, and the independent cosine variances give (A6).

The program evaluates the entire actual fifteen-coordinate AST bank, including all six broad leaf factors. With 220,000 fixed-seed samples and k=8,16,32,64, the normalized measured variances are approximately 0.21944, 0.21822, 0.21844, and 0.21880, consistent with (A6).

## 4. The symmetric OU bridge preserves a bulk obstruction

The difference mode `F_low=cos(W1-W2)/2` has exact variance

    Var(F_low) = (1-exp(-2))²/8
               = 0.0934556340519386.                   (A7)

It is independent of the structural frequency k. Its gradient is taken *before* the Mehler smoothing, as required by the supplied symmetric bridge identity. Thus

    E |P_rho D F_low|²
       = exp(-2(1-rho²))(1-exp(-4rho²))/4.              (A8)

The two cosine modes occupy orthogonal Gaussian directions. Their gradient contractions have no cross term, so (A8) is a nonnegative summand in the exact bridge integrand for F_k. Integrating `2rho` times (A8) over the fixed interval `[1/4,3/4]` gives

    0.023501377273783458.

Multiplication by the square of the coefficient prefactor in (A5) yields a k^4 contribution already on this fixed bridge interval. Broad-leaf errors and their first derivatives are exponentially small in k; including them preserves the growing contribution.

This is stronger than merely saying that an exact identity preserves the total covariance. The divergent amount cannot all be confined to a vanishing near-rho=1 strip and cheaply removed by that cutoff. A faithful finite positive bridge quadrature must approximate this positive, smooth bulk term as part of its coefficient accuracy. Deliberately omitting it changes the target covariance.

## 5. What positive dyadic quadrature does and does not imply

The inequality `w <= C Delta` alone is not a lower bound for any chosen node. The lower-bound assertion requires actual positive panel mass and a bound on node count.

For an N-point positive rule on a dyadic primitive panel of length Delta near r=1, with total base mass Delta, some node has base weight at least `Delta/N`. The R5 multiplier is `r^4`, bounded below by a positive constant there. Select the node of largest effective R5 weight and use that *same legal node* for both high primitives. Its half-variance x is comparable to Delta. Choosing `k=x^(-1/2)` gives

    w >= c x/N,
    coefficient magnitude >= c B k²/N²,
    one-node covariance >= c B² k^4/N^4.               (A9)

For a fixed-order mapped Gauss-Legendre panel rule, the ratio w/x tends to a fixed positive constant. For example, the central node of the three-point rule has `w/x -> 8/27` in the panel convention used by the checker.

If N grows only polylogarithmically in k, it cannot absorb k^4. If the other thirteen clocks are also refined, their selected masses should not silently be called fixed constants. Selecting a largest positive node in a fixed interior panel gives a lower mass of order the reciprocal of that panel's node count. Those thirteen masses add a further polylogarithmic factor to (A9); the polynomial growth still dominates. With one common node count N, the coarse lower bound is of order `k^4/N^30`.

A single rule with a fixed strictly positive endpoint cutoff has no k-to-infinity sequence of included nodes. The obstruction concerns uniformity over increasingly fine legal rules or vanishing cutoffs, not divergence inside one fixed finite program. Arbitrary positivity alone, with arbitrarily huge numbers of tiny nodes and unconstrained joint refinements, is insufficient to assert (A9).

## 6. Native radius and first-path interpretation

For the literal pair of normalized C3 targets and six normalized C0 targets, the original coefficient-matching normalization is

    beta = B w_k²/(sigma1³ sigma2³).

Since `sigma1=sigma2=1/k` and `w_k` is of order k^-2,

    beta is of order B k².                             (A10)

Independent exact reshielding cannot change this power: from (A2),

    t1³ t2³ <= (1+h²)^3 x³ <= 8x³.

Thus the old prescription “all nonroot amplitudes alpha; all normalization on one physical-leaf root” has

    |rho_root| = alpha beta / R,

where R is the positive readout product. With readout factors at most one this is at least `alpha beta`. It violates a cutoff-uniform `Lambda alpha` radius budget. For fixed alpha and fixed positive other masses, it eventually exceeds even a fixed absolute amplitude cap. This is a direct arithmetic failure of the prescribed amplitude, not an inference from a high ideal coefficient.

There is also the familiar source-zero first-path diagnostic. For the admissible linear source g(t)=t, a normalized C0 leaf equals one and a C3 center target vanishes. In the stationary ideal graph a root attached to that center still has

    Y_root = rho_root Z_center + sqrt(1-rho_root²) Z_root.

Relative to its source-zero carrier Z_root, the derivative in Z_center is exactly rho_root. A zero target high cumulant does not erase this root leakage. This ideal formula identifies the missing first budget; an executed bounded-filter native program still requires its calibrated absolute errors and cannot be declared exactly equal to the ideal graph without those certificates. The explicit overlarge amplitude already suffices to reject the blind implementation when outside its native window.

Moving weights to other vertices or using more programs changes the amplitude and feedback ledger. It is not refuted in general. No lower bound here is automatically a W2 law lower bound.

## 7. Heat powers and exact limits of the verdict

The displayed constants use zero *added* heat; every original private shield is nevertheless strictly positive. Added heat comparable to k^-1 preserves all polynomial growth with modified constants. Fixed strictly positive added heat suppresses arbitrarily high k and is outside a cutoff-uniform claim.

For `k` of order `tau^-1`, this history has the diagnostic powers

    beta ~ tau^-2,
    blind root radius ~ alpha tau^-2,
    old-coefficient covariance ~ alpha^16 tau^-4.

With `tau=alpha^gamma`, the last power is `alpha^(16-4gamma)`. This disproves dropping the width losses or blindly reusing the old O(alpha) radius ledger. It does **not** disprove every power-law heat choice. A sufficiently slow heat choice can make some of these quantities small at some fixed target grades; its approximation, radius, first, readout, and descendant terms need fresh analysis.

Finally, a large node variance does not lower-bound the variance of an arbitrary sum of correlated nodes. Positive scalar quadrature weights do not imply nonnegative cross-node covariances. The result does directly invalidate a uniform sum of nonnegative nodewise covariance allowances, any corresponding absolute-history ledger, and an independently banked assembly whose variances add. A complete shared-bank cumulant sum might exploit correlations or regrouping, and its exact nonlinear target can have cancellations absent from an isolated history. No contrary claim is established here.
