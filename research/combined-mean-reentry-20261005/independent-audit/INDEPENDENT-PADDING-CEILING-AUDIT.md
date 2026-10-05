# Independent audit: padding, reentry ceiling and next bridge port

2026-10-05. **PASS as conditional source-contract and native-scaling analysis.** This does not independently admit the combined prefix/skew producer or numerically execute LOW30. No file in the sealed reentry packet was changed.

Reviewed companion: `../POST-SYNTHESIS-REENTRY-AND-NEXT-PORT.md`, SHA256 `53748adb366902bc875e5d74871bcd1a9cda74ec22524bba699f38fc4c918f49`.

Verified source pins:

- LOW30 `30_low_acc.tex`: `7767ac389d6b677eb9517b2585f7a24cb9f85ede47d099a653c635b7be750fb8`, particularly `b27:raw:self-reserve` and `b27:raw:split` at lines 6980–7090.
- The sealed finite-depth reentry manifest: `106350c3f20d06f6eaae1fb9b46677e0ccc7bfe8d21ec5f3e62f72ffdbaf515a`.
- For matching the proposed cost/root bookkeeping only, I inspected the then-current combined graph `../../combined-bridge-skew-20261005/COMBINED-BRIDGE-SKEW-JOIN.md`, SHA256 `11b9241e22cdac27ceab6c6629b63420358ca0f877cf666f6b61a8de8acf7b1b`. That inspection is not an independent admission of its complete producer theorem.

## 1. Variable-buffer padding substitution

Suppress only declared public-log factors. With `r=rho=A/sqrt(u)`, `delta=a=A`, baseline order 4, physical readout sqrt(u), and terminal Lipschitz A, multiply each native error by A sqrt(u). The six terms are exactly

    A^5 u^(-3/2),  A^5 u^(-1),
    A^4 mu u^(-1/2),  A^5 u^(-1/2),
    A^6 mu^(-1/2) u^(-3/2),  A^6 u^(-3/2).

Thus mu=A^(1/2) gives the claimed A^(9/2) third term and A^(23/4) fifth term. The literal native first theorem gives

    physical residual first <= Lambda[A+A^(7/4)/sqrt(u)],
    normalized reserve first <= Lambda A^(7/4)/u.

At eta=A^(3/2), every executed u is at least a fixed positive fraction of eta. Mean and reserve normalized radii consequently both have power A^(1/4), with their own actual logarithmic factors and numerical admission thresholds. Mixed-K has power A^(1/2). In fact the physical residual first is O(Lambda A) **pointwise**, since A^(7/4)/sqrt(eta)=A. Weighted closure is also valid and gives `Lambda[A+A^(4/3)]` at w=A^(5/6).

For positive dyadic integration, `sum omega u^(-p)` is bounded by C times `w^(-p)+eta^(1-p)` when p>1, by `w^(-1)+log(1/eta)` when p=1, and by `w^(-1/2)` for p=1/2. At the stated scales the six rows give, respectively:

- 15/4 and 17/4;
- 25/6 and a grade-5 logarithmic term;
- 49/12;
- 55/12;
- 9/2 and 5;
- 19/4 and 21/4.

Every one is at least 15/4. The near-RAW square-root term has grade 15/4, and its eta/sqrt(w) companion has grade 49/12. The known remaining combined scalar rows, including A^6 integrated against v^(-2) and A^7 against v^(-5/2), are higher grade at these scales. The native guards remain execution-time obligations at the changed actual normalization; the source theorem allows mu to change but does not waive them.

Three important qualifications are retained in the corrected companion:

1. This cannot upgrade the old unmodified join by padding alone. Its original prefix and Gaussianization terms would have grades 13/4 and 19/6 here. Their actual replacement requires the separately admitted new producers.
2. The leading service term is Lambda A^(15/4). A numerical-constant A^(15/4) guarantee does not follow by silently absorbing that logarithm. Any strictly lower grade needs its explicit absorption window.
3. The optional final unit-variance own-mean compiler has no subsequent terminal A factor. It must retain mu=A for its A^4 allowance. Using mu=A^(1/2) there would instead introduce A^(7/2).

The changed internal padding also does not authorize relaxing independent Gram coefficient calibration, covariance target restoration or fresh-bank/caller requirements.

## 2. Exact three-term ceiling

Let

    x=A^2h, y=A^3w^3/h, z=A^5w^(-3/2).

Then x+y>=2A^(5/2)w^(3/2), and a second AM–GM inequality yields

    x+y+z >= 2sqrt(2) A^(15/4).

Equality holds at `w=2^(-1/3)A^(5/6)` and `h=2^(-1/2)A^(7/4)`. Equivalently, the weighted average with weights 1/4, 1/4, 1/2 of the three exponent rows is identically 15/4, so their minimum cannot exceed it. This proves the ceiling even if the scalar scales are not exact powers. With a fixed displayed coefficient Lambda multiplying z, the minimum is `2sqrt(2Lambda) A^(15/4)`; optimizing does not eliminate its logarithmic dependence.

This is a lower bound on a **sum of positive certified allowances**, not on actual approximation error. It remains the ceiling of this ledger when the same three fresh allowances recur at each depth. More accurate means from earlier depths do not cancel or shrink those independently charged terms. The assertion does not exclude a tighter certificate, cancellation, or a changed producer.

## 3. Conditional finite-depth recurrence

Given the independently admitted combined source ports, the stated recurrence `e_(k+1)<=A e_k+B_(k+1)+floors` is correct. The actual existing bridge consumer is A-Lipschitz, the mean service alone is reentered, and conditional means propagate with coefficient at most A. The true-force covariance and third-cumulant comparisons remain those of coherent analytical histories, not the raw S_k graph. Uniform O(A) true-force first/energy bounds also preserve the short-prefix and terminal-smoothing comparisons at every fixed depth.

The strict 37/10 and logarithmic 19/5 ledger displayed in the companion is arithmetically consistent with the specified strict scales, conditional on the full joined producer being independently admitted. Its finite geometric recurrence and limiting-mean comparison are correct. Nothing in this audit supplies an endpoint-distribution theorem or an order-raising iteration.

## 4. Next covariance-qualified bridge port

The original proposed RAW port needed an additional explicit requirement, which is now present: its **complete private first** must be at most C A w^(3/2). Mean/covariance matching by itself would not control the replacement's own third Stein defect.

Use a common declared first majorant L*=C A w^(3/2) for both true and replacement prefixes. Conditional on the original endpoints (X,Y) and the independent future bank, Gaussianize both RAW sources using the independent h-Gaussian. Paying both third-field defects and the buffered Gaussian mean/covariance restoration gives

    A delta_m + C A delta_Sigma/h
      + C A (L*)^3 sqrt(D)/h^2,

hence the claimed A^4 w^(9/2)/h^2 term. The declared endpoint first O(Aw) justifies doing the moment comparison at the original Gaussian endpoint pair and then commuting the future force at cost O(A^3w). No certificate is silently transferred to a shifted caller law.

For a whole-LAW implementation, its separate certified law error must also be multiplied by terminal A, and its complete carrier must be allocated from the available h^2 variance exactly once. The corrected companion states both obligations. This remains a producer contract, not a construction. In a future implementation, an endpoint-first bound alone need not make u-minus-prefix nonexpansive; preserving an exactly A-Lipschitz terminal consumer requires its additional monotonicity/contractivity port. A bound A(1+C Aw) follows from firsts alone. This does not affect the current recurrence, whose existing P_Q has the required symmetric positive-semidefinite X derivative.

## 5. Literal costs and roots

The proposed combined count formula matches the inspected current graph. A bulk node uses four whole service graphs, J original bridge VALUES, and one terminal VALUE. A near node uses F_Q, the same J bridge VALUES and one terminal VALUE. Capture/restoration costs are separately included. Each source reentry expands S_k and its baseline at every occurrence, with fresh complete tapes and no changed-argument cache.

A bulk node has Y, untouched keep, bridge and terminal-smoothing roots outside its four complete service banks, totaling `4D+d_M+d_H+d_K+d_3`. A near node has 6D. Aligning one D-dimensional carrier across all nodes removes exactly `(N_out-1)D`, yielding

    D+sum_bulk[3D+d_M+d_H+d_K+d_3]+5D N_near.

The formulas are therefore correct provided the named service dimensions retain their full declared banks and the external keep root is not counted again inside them. This is consistent with the inspected graph. At fixed depth, the polynomial-log cost conclusion is conditional on the actual native guards and literal occurrence counts, not a uniform-in-depth or sublinear-in-order theorem.

## Diagnostic scope

`check_padding_ceiling.py` passes **1,827 exact rational assertions** covering the six scaling substitutions, integrated grades, radius exponents, final-completion qualification, and exponent-ceiling identity. `padding-ceiling-checks.json` pins the reviewed companion. It does not execute a native compiler or establish the combined producer. The one substantive missing RAW-port bound found in review, and the padding/log/caller qualifications, were incorporated before the pinned version above.


## Application audit: the combined-source premise is now sealed

2026-10-05. **APPLICATION PASS under the actual native guards at every fixed depth.** Reviewed application note: `../APPLICATION-TO-SEALED-COMBINED-SOURCE.md`, SHA256 `8275754fbb614d41540061fd244c29474973db06440855b4ddfa19872a80e5d3`. The conditional companion's reviewed hash remains `53748adb366902bc875e5d74871bcd1a9cda74ec22524bba699f38fc4c918f49`.

I verified the admitted producer's manifest SHA256 `04b93b69a562f49c0ab1ed2fc49da2208db3204649f84f8a95964f333a71e974`, its exact theorem hash `11b9241e22cdac27ceab6c6629b63420358ca0f877cf666f6b61a8de8acf7b1b`, all eight listed member-file hashes, and all four imported manifest pins. Its independent audit explicitly admits that final theorem hash and distinguishes its 2,597 diagnostics from native execution; the author record contains 8,598 diagnostics with the same limitation.

The theorem's Sections 1 and 8 export the required actual finite gradient-baseline/near-gradient source, own-mean certificate, complete private and retained-caller first, square-lift curl, energy and origin records. Its actual bridge consumer has Lipschitz constant at most A because its X derivative is symmetric positive semidefinite and bounded by Aw. Sections 3–6 preserve the exact bridge/future ancestry and complete independent whole-LAW comparisons; Sections 9–10 provide the full replay/root/floor bill and actual normalization guards. These are precisely the companion's antecedents. The fresh F_Q near branch and fixed covariance/third-tensor services need no promotion of S_k to a RAW force law.

Therefore the previously conditional implication now applies to this admitted source: for each fixed depth, with every newly expanded native graph admitted at its actual constants and dimensions,

    e_(k+1) <= A e_k
      +C A^(37/10)sqrt(D)
      +Lambda_(k+1) A^(19/5)sqrt(D)
      +absolute floors.

The coherent higher-depth covariance and third-cumulant restorations were already proved independently; their bills are dominated by this same displayed forcing ledger. Full private-bank integration and fresh-tape replay obligations are preserved. The resulting finite-depth construction is a genuine same-grade own-mean recurrence. It does not improve the grade through iteration. The limiting-mean comparison remains valid, and no endpoint-distribution conclusion follows.

This application does **not** admit the retuned-mu 15/4 construction or a new covariance-qualified bridge producer. Those discussions retain the algebraic/contract scope audited above. Only the previously missing premise for the strict 37/10 finite-depth application has changed.
