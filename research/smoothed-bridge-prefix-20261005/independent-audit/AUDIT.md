# Independent audit: smoothed bridge prefix and narrow endpoint

2026-10-05. Independent mathematical and finite-diagnostic review. This audit does not execute the native mean, Gram, mixed-K, or final raw-split compiler. All conclusions involving those services remain conditional on the imported actual native guards and path certificates.

## Verdict and reviewed sources

The construction supports the claimed guarded grade **17/5**. I found no mathematical obstruction in the new smoothing, bridge replacement, positive quadrature, narrowed endpoint branch, terminal-carrier alignment, or root accounting.

I read the full draft `SMOOTHED-BRIDGE-PREFIX.md` at SHA256 `7eb8f12b828a1a04a3920d06543e8f90f3b91b9b5fd4c1070ed24319520c1003`, the imported LAW-only shrinking-buffer theorem, and the prior delayed conditional-noise theorem. The imported LAW-only manifest hash in the draft was independently checked and is correct.

The two wording clarifications below were subsequently incorporated and checked at draft SHA256 `c26e44e932f0626cd7101d092d679f7b416d4d1e71adaf09e3cf4a48941e4258`: the Gaussian Stein operator is now explicitly defined, and the precision ledger now refers to positive interior bridge gaps. These comments are resolved.

## Prefix and smoothing checks

1. The terminal commutator is applicable to the genuine nested force. Its caller first is bounded by A/2+A²/4 and its stationary L2 energy by A(1+A)√D. Coherent coupling and the Hermite multiplier bound therefore give C(A²h+A h²)√D. The executing terminal is an original VALUE with an extra Gaussian root.
2. Omitting H_s only from the short prefix costs A³w√D after the terminal Lipschitz factor. The resulting ordinary bridge is independent of the genuine future tail conditional on (X,Y). The nested prefix itself need not have that independence.
3. Both true and common-N prefixes have conditional complete-bank first at most A w√(w/(2−w)). With the separate hG root, mean interpolation and one Gaussian integration by parts consequently give C A³w³/h√D. The conditional-mean quadrature discrepancy adds A²wε_bridge√D.
4. Positive weights of mass exactly w imply 0≤D_uP_Q≤Aw I. Thus I−D_uP_Q is a contraction when Aw≤1. No commutation between distinct Hessians, or between the terminal Hessian and prefix derivative, is needed. Commuting the tail into every X-site costs another C A³w√D.

For complete precision, the Stein field in the proof should be written as `D(−L_G)^−1(p−μ)(Dp)^T`, with L_G the Gaussian OU generator on the complete private bank. This avoids confusing that operator with the outer resolvent. Its operator norm is at most L_p², and its Hilbert–Schmidt norm at most L_p²√D. The independent smoothing root is what permits the stated vector estimate without another dimension factor.

## Positive bridge quadrature

The row l_s and the complex norm identity are exact. The inequality sin²(y)≤y² and the displayed lower bound on the real bridge variance establish contraction for y²≤x(δ−x). The parameter-2 Bernstein ellipse around a panel [a,b] with b−a≤a stays in the stated wedge; reflection treats the right half. Operator-valued analyticity gives a dimension-independent positive Gauss error, and endpoint removal/replacement costs only its scalar mass. Normalization to mass w preserves the error scale and is essential for the stated contraction bound.

The endpoint nodes s=0 or δ have σ=0 and require no inverse zero gap. Statements about positive bridge gaps in the precision ledger should refer to the interior nodes.

## Narrowed endpoint and actual native scales

For 0<η≤w, the exact reciprocal variance gives v≥η/2 on t≤1−η. For p≥1, split v^(−p) using q²/(1−q²)+1/(1−t²). The endpoint contributions are logarithmic for p=1 and O(η^(1−p)) for p>1. For the RAW branch, the square-root reciprocal is bounded by C[w^(−1/2)+(1−t)^(−1/2)], giving C[η/√w+√η].

These are also valid for the actual finite outer rule. After inserting the branch boundary, a dyadic panel of gap d has weight O(d) and node gaps comparable to d. Its contribution is O(d^(1−p)); the early positive midpoint contributes O(√d) on the RAW side. Geometric summation proves the finite bounds. Splitting at q and 1−η changes neither the order nor positivity. The early cutoff A³ is below η=A.

At w=A^(3/5), η=A, h=A^(7/5):

- A²h, A³w³/h, A⁴/w: exponent 17/5.
- Near RAW: exponents 7/2 and 37/10.
- Nested-prefix omission and coherent commutation: exponent 18/5.
- Terminal smoothing drift A h²: exponent 19/5.
- Integrated mean/Gram bill: exponents 41/10 and 9/2.
- Integrated mixed-K native bill: exponents 53/10 and 35/6.

The logarithmic A⁴ contribution and optional final ΛA⁴ bill fit the declared remainder. Absorbing ΛA⁴ into the leading term requires ΛA^(3/5)≤1.

At the minimum actual buffer u≥cA, normalized mean/Gram radius and raw self-reserve first both have exponent 1/2 (with their actual public-log factors); mixed-K radius has exponent 2/3. The direct residual first Λ[A+A^(3/2)/√u] is O(ΛA). These checks validate the declared parameter choice, not an unrestricted claim that the same O(ΛA) port holds for η≪A.

## Readsets, stronger carrier alignment, and exact roots

Each tail LAW comparison can retain (z,N,G_h) conditional on Y because those labels are unread by fresh complete service banks. The A-Lipschitz consumer is applied only after the comparison. Aligning different nodes' known carrier rows is then a direct graph operation: one-node marginal laws are unchanged, but no joint cross-node LAW comparison is asserted or needed.

The stronger baseline in the draft is correct. The known terminal Gaussian row has squared norm r²(1−t²)+h²=1−r²t² and is strictly positive. Rotating this whole row into the common G leaves every off-G terminal path passing through the O(ΛA) displacement. The remaining main Hessian difference multiplies a scalar carrier coefficient and is symmetric. Consequently the curl is O(ΛA²), the remainder energy is O(ΛA²)(|z|+√D), and no h≤A restriction is needed. This is stronger than keeping the old unsmoothed baseline.

For the sealed unaligned bulk dimension d_old=D+d_M+d_H+d_K and near dimension d_old=4D, the new per-node dimensions before alignment are d_old+2D. Sharing only the aligned D-carrier gives exactly

    d_T = D + Σ_bulk(d_old+D) + 5D N_near.

This is the draft's formula. It adds two D-roots per node relative to the old aligned source, preserves all perpendicular native coordinates, and does not discard a service bank. Every final compiler reentry must still instantiate and replay the entire d_T tape and full VALUE graph.

## Independent diagnostics

`check_bridge.py` and `bridge-diagnostics.json` accompany this audit. All checks passed:

- Bridge mean, variance, complex norm, and carrier-row identities agreed to at most 6×10^−16.
- Sampled contraction wedges and parameter-2 ellipses stayed inside the certified domain.
- 300 dimension-seven noncommuting PSD matrix trials satisfied the contraction and A-Lipschitz terminal bounds.
- Literal positive dyadic Gauss sums for A from 0.1 down to 10^−8 obeyed all stated bulk and near bounds. The largest observed bulk/bound ratio was below 0.537; the largest near/bound ratio was below 1.380.

These finite checks support the algebra and scalar quadrature discussion. They are not a numerical instantiation or certification of the full imported native compiler.
