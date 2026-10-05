# Independent audit: typed repair algebra

2026-10-05. Scope: the formal framework and exact algebra in `REPAIR-GRAPH-ALGEBRA.md` and `check_repair_algebra.py`. No native sampler was run. The audit does not certify an unsmoothed all-order positive finite-VALUE program.

## Verdict

**The revised framework is sound as a conditional, fixed-boundary formal inverse architecture.** Its central exact identities pass an independently implemented finite test. The conditions on old ports, realizable leading currents, ownership, and local finiteness are substantive inputs, not consequences of graph notation. No all-order heat or cost theorem follows.

The audit prompted the following necessary clarifications, incorporated in the final revision:

1. For heterogeneous genuine root copies, `d_child = sum_i d_i` and `Psi_child = sum_i Psi_i + 2 Gamma (h-1)`. The expressions `h d` and `h Psi` apply only when root types have the same surplus. Retained side occurrences keep their own allocation and ownership records.
2. A right inverse satisfying `L h = Id` must also be filtration-compatible on the actual realized typed representatives. Every normalization, width, injection, and caller loss must enter their grades before the contraction argument. Positive source grades alone do not prove this compatibility.
3. The displayed primitive-order bound is proved for the conditional/bridge grammar. An additional positive-increment operation requires its own finite-branching and order-growth certificate.
4. A newly recorded derivative loss does not automatically preserve reserve superadditivity; the fully charged bridge inequality must be tested, or the loss must have been prepaid. The heterogeneous bridge bound requires `gamma_v <= Gamma` for every hit occurrence. A heterogeneous occurrence-count bound uses a fixed positive `beta_min`.
5. Rank-zero normalization currents remain in the current/test pairing until an exact cancellation is shown. Quotient identities do not confer uncharged resource or filtration equivalence.
6. Heat insertion is grade-neutral in the ordinary ledger and in the clock-mass ledger after the local charge saturates at `k >= 2`. The low-order clock credits are finite.

## What the independent script establishes

`check_independently.py` passes **159 assertions**. It does not call or import the author's checker.

- It obtains the scalar spine cumulants from the full labeled set-partition formula and Gaussian monomial moments through amplitude order six, and compares them with the closed Gaussian-integral expression. In particular the second-order branches have ranks `4, 6, 6`, with separately retained `a^2` and `b^2` coefficients.
- It checks the whole-bank tree identity for 11 polynomial examples with up to four arguments. Every labeled tree is enumerated. The Gaussian moment uses its actual min-path covariance, and every tree-parameter ordering sector is integrated exactly over its simplex. This tests the tree weights and derivative multiplicities independently of the reserve arithmetic.
- It checks mixed, heterogeneous conditional surpluses with independently omitted sides and ordinary/owned-clock bank hits. It also checks that the naive identical-root formula fails for unlike root surpluses.
- It tests the new trace-to-bridge heat-jet identity through heat order four by direct independent-versus-common Gaussian-shift polynomial moments for `(N,D)=(1,2),(2,2),(3,1)`. A nonconstant multiplier counterexample verifies why the frozen-coefficient condition matters.
- It solves a two-current formal system with old/new and new/new interactions through total degree seven; the successive difference valuations are `1,2,3,4,5,6,7,8`.
- It checks every even observer cumulant from rank 2 through 32. For `g(x)=c x+epsilon sin(x)`, with `c=1/2, epsilon=1/4`, the derivative lies in `[1/4,3/4]`. The exact first-order tilted term is `c theta^2 + epsilon exp(-1/2) theta sin(theta)`, so arbitrarily high ranks occur at the same alpha grade.
- It checks the keep integration-by-parts identity with a nonunit keep scale and a shifted endpoint, plus exact counterexamples when coefficients depend on the keep or unconditional equality is mistaken for retained-record equality.

Finite tests complement the explicit algebraic proofs; they do not establish the imported native/current tables or analytical remainder bounds.

## Formal inverse and local finiteness

Work in the positive-reserve, separated completion of an explicitly specified compiler whose residual is well-defined and reserve-locally finite. A fixed finite primitive/current-jet system with finitely branching source realizations is sufficient. A finite old/new dictionary by itself does not prove this for the unrestricted all-order grammar. All reachable types must have a fixed admissible leading-current section. Every hit old port needs reserve at least `delta`; an old `O(alpha)` first bound does not imply this. Telescoping a mixed residual monomial singles out one difference and leaves at least one additional positive-reserve root/argument. The conditional and bridge rules therefore raise difference reserve by at least `delta`.

For the core grammar, a retained node with `G<P` has bounded reserve and occurrence count. Reserve bounds both ancestry depth and operation arity. A tracked occurrence acquires at most the operation arity minus one new derivative hits per bridge stage. This gives the conservative `k_max <= k_seed + ceil(P/delta)^2` estimate. Finite primitive/current expansions and deterministic source realizations then give finite branching. Heat traces and a source-zero observer violate the needed structure and must stay explicit escapes.

The universal graph algebra need not be locally finite in pure `Psi`: one may vary the number of side occurrences while setting `a = beta N + gamma K + d` with fixed `d`, leaving `Psi` unchanged. The finite `G<P` bounds above therefore do not supply this missing universal reserve-local-finiteness claim.

The equation is solved in the reserve filtration. The operational `G<P` queue instead bounds high-`G` residuals without compiling their source programs. It is not automatically the `G` truncation of the full reserve-completed inverse. A high-`G` term may be discarded only as an analytically certified remainder of the actual final simultaneous program. Later changes must update every affected mixed, first, and replay ledger. An unexecuted high-`G` correction has no executed descendants. `G` can decrease when conditional offspring omit side occurrences, so `G` should not be silently substituted for the reserve ordering or treated as a universally rewrite-stable quotient.

## Equality and realization boundaries

The full-bank tree formula is an equality after the complete bank expectation. It does not imply pointwise equality in a retained old bank. The observable-current quotient is linear at one endpoint and retained-record contract; it is not automatically an algebra quotient respected by arbitrary graph products. The keep reduction requires coefficient/keep independence or its complete product-rule corrections, and its inverse keep powers remain in the resource ledger.

The admissible-gradient observer counterexample is exact and appropriately scoped: it disproves inference of rank local finiteness from small perturbation size. It does not prove that every enlarged functional-generator or block-inverse representation fails.

The owned-clock contract supports the charge `q(k)=(k-2)_+` only for an actual retained clock weight whose amplitude has been relocated. A derivative hit cannot mint a second weight. Its two finite heat reentry steps are compatible with this algebra; indefinite heat closure is not established.

The new trace-to-bridge construction is exact for the full uniform local-trace sum with frozen scalar coefficients. It leaves a common-translation generator and distinct-vertex same-graph edge insertion. Those are not automatically supplied by the inter-packet bank-tree port. Signed inverse cross-heat is not a positive probability kernel, and the finite formal identity neither improves saturated reserve nor proves a heat remainder. This section correctly narrows the missing-generator question without closing it.

The unresolved obligations remain an exact source/current right inverse for every reachable escape, complete simultaneous native residual coefficients rather than a norm-only floor, positive joint realization with ownership and first/curl guards, a target-qualified original-heat remainder, and a complete replay/cubature cost bound. The manuscript correctly separates these obligations from finite algebra and explicitly declines the global `c(P)=o(P)` conclusion.

## Reproducibility

Run `python independent-audit/check_independently.py` from the parent directory. Results are in `independent-checks.json` and `check-run.txt`. The result JSON records SHA-256 hashes of the exact two audited main files. The final manifest separately pins this audit and the independent script.
