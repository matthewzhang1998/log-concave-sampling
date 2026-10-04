# Log-concave sampling: mathematical research snapshot

Snapshot version: 2026-10-04-public-v5-relative-known-center-linear, updated 4 October 2026. This is a curated research archive of proofs, scoped audits, source certificates, and diagnostic code. Inclusion is not an assertion that every proposed construction is proved.

## Progress in brief

The source-relative exact **1/32** baseline remains the best end-to-end reference in this collection. A new independently audited **relative standalone known-center posterior source family has a linear query exponent**. Its complete source assumptions are stated below. The old uniform hidden-family cost and generic no-copy matrix endpoint remain separate open problems; no subpolynomial-dimension or improved end-to-end sampler bound is claimed.

### Baseline

[Exact denominator 32, version 2](research/exact-slack/exact32-proof-v2.tex) is preserved byte-for-byte, SHA-256 `c7487ca24ffa44ea73b52e42e4b6f715ad964171fa5bd0a61cf9323eb646cbff`.

Its theorem is relative to the pinned LOW30/LOW31 foundational constructors, assumes a supplied initialization, and counts original gradient and directional Hessian-vector queries. It does not claim dimension-free arithmetic/storage or a differentiable return for the new packet. Those foundations are referenced rather than silently treated as included, independently reproved sources. The historical componentwise R7 and fourth-VALUE source certificates are under `research/exact-slack/`.

### New milestone: a relative linear known-center source family

The [protected short-refresh construction](research/cost/same-heat-refinement-20261004/protected-short-refresh/DETERMINISTIC-KNOWN-CENTER-SOURCE-WITH-PROTECTED-REFRESH.md) has [independent relative admission](research/cost/same-heat-refinement-20261004/protected-short-refresh/extension-audit/INDEPENDENT-WEIGHTED-GRAPH-AND-KNOWN-CENTER-ADMISSION.md) for a standalone known-center law/source/query induction:

    k_R ≤ max(k_P,R−3/2),  R=P+1/2.

Starting from the completed CW7 known-center seed k_7=0, the displayed rebuilt family has k_R=R−3/2 for R≥7.5. This is an actual relative source-family result, including first, retained-origin/graph and fixed weak protected-proxy ports, rather than an arithmetic-only recurrence proposal.

Admission depends on the checked seed's genuine-gradient-reference and arbitrarily accurate finite-restoration ports at the actual consuming widths, consistent caller/anchor treatment, and matching finite row/coisometry encodings. The extension audit resolves the weighted-global-normalization premise left open in the [original short-refresh audit](research/cost/same-heat-refinement-20261004/protected-short-refresh/independent-audit/INDEPENDENT-PROTECTED-SHORT-REFRESH-AUDIT.md). It does not independently reprove the seed.

The old uniform known/hidden-family c_P ~ 2P² is **not replaced**: hidden-service calls retain their separate explicit bill. The simple known-center-only outer application has exponent (R−1/2)/(2(R−1)) on Dim/ε², tending to 1/2. It does not beat the denominator-32 baseline. Sublinear-order cost remains open.

The [consumer-depth clarification](research/cost/same-heat-refinement-20261004/protected-short-refresh/CONSUMER-DEPTH-AND-ACTUAL-REPLAY-AUDIT.md) distinguishes the new proximal chain's external hidden depth zero from the selected fixed kinetic consumer's seven levels. The same growing-order kinetic reference needs J≥ceil(R−1); its hidden cost can consequently become quadratic again. Fixed inner RAW7/A8 histories and actual replay remain charged. The optional depth-adapter formula is conditional, not a new hidden-family theorem.

### Closed bounded modules

Each claim below retains the assumptions and numerical/source costs in its proof and audit.

- Known affine common-input heat and one-pair skew adapters: `research/p-affine-heat/`, `research/c-common-input/`, and their independent audits.
- Explicit nonlinear shear sources, constant-reference covariance corrections, and local retained-graph qualifications: `research/p-shear-heat/`, `research/r-shear-consumer/`, `research/c-shear-audit/`.
- Exact finite two-node third-iterate source algebra and the scalar direct true-Gram reserve: [scalar construction](research/r-native-twin/SCALAR-K3-DIRECT-TRUE-GRAM-RESERVE.md), with separate P/C audits. The scalar mean corollary is local; it is not a complete posterior host.
- Auxiliary variance-preserving rotation, genuine-gradient constraint offspring, and rank-aware restoration under the stated Hessian bounds: the `auxiliary-heat/` audit set and `ROTATION-CONSTRAINT-DEFECT-COMMON-GRADIENT-LIFT.md`.
- The small-full-gradient constant mixed-covariance service now has a [scoped independent PASS](research/c-native-twin/auxiliary-heat/INDEPENDENT-SMALL-GRADIENT-MIXED-DESCENDANT-AUDIT.md), binding final source `f805e8d3…`. Absolute finite-prior and numerical floors remain explicit.

### New scoped progress

- [Retained mixed-field extraction](research/no-copy-rank-20261004/RETAINED-MIXED-KERNEL-AND-TERMINAL-FLAT-GATE.md) has an [independent scoped PASS](research/no-copy-rank-20261004/INDEPENDENT-RETAINED-MIXED-EXTRACTION-AUDIT.md): a positive conditional covariance reserve retaining the same roots, with no additional heat-power replication cost at fixed target under its stated clock/compiler assumptions. The companion retained-root observer obstruction and terminal-flat K3 witness show why the full ambient endpoint/current identification remains necessary.
- [Same-heat quarter-flow refinement](research/cost/same-heat-refinement-20261004/SAME-HEAT-POSTERIOR-REFINEMENT-AND-FLOW-GATE.md) has an [independent scoped PASS](research/cost/same-heat-refinement-20261004/independent-audit/INDEPENDENT-SAME-HEAT-FLOW-AND-QUERY-AUDIT.md) for exact posterior invariance, contraction by a/(1−a), and a legal finite VALUE/HVP DAG. Its original unprotected law-level/query certificate is linear in order and does not alone admit a recursive source family. The protected short-refresh extension above now supplies a different, explicitly source-qualified standalone known-center induction. The deterministic information lower bound is limited to fixed-input first-Picard integration and adaptive baseline transcripts; it is not an exact-flow, randomized, Gaussian-input Lp, or posterior-law lower bound. Prefix-sum arithmetic does not change the query exponent.

These results do not improve the currently admitted general c_P family or prove an all-order/subpolynomial sampler.

### Latest interface and source results

- **The selected CW7 external host needs a smaller interface.** The [host-separator audit](research/no-copy-rank-20261004/host-observer-audit/ACTUAL-CW7-WEAK-E-HOST-SEPARATOR-AND-COMPARISON-ORDER.md) and [independently audited reentry recipe](research/no-copy-rank-20261004/host-observer-audit/INDEPENDENT-COMPLETE-MEAN-REENTRY-AUDIT.md) show that its physical law consumer reads captured caller θ and one complete statistic T. A complete conditional mean law can be proved first and the host attached afterward. An arbitrary joint-fine-observer theorem is not an external prerequisite for this selected host. Internal packet aliases and any baseline still reading integrated fine variables remain genuine obligations; the actual retained-source/graph proof is separate.
- A [positive projected-gradient/V-column construction](research/no-copy-rank-20261004/endpoint-next/v-column-escape/PROJECTED-GRADIENT-V-COLUMN-ESCAPE-AND-JOINT-RETENTION.md) succeeds on the terminal-flat fixture, with actual-energy control and an exact physical orientation selector. Its zero-auxiliary-root Gram identity is not a full retained-root endpoint theorem. Earlier midpoint, fixed-frame and baseline obstructions remain limited to their stated operations. The [cache partition](research/no-copy-rank-20261004/endpoint-next/cache-partition/WEAK-29-CAPTURED-CALLER-AND-COMPLETE-PRIVATE-CACHE-PARTITION.md) allows exact shared captured work but keeps complete private bank work charged.
- [Scalar randomized quarter-flow](research/cost/same-heat-refinement-20261004/randomized-weak/SCALAR-RANDOMIZED-QUARTER-FLOW-AND-ALL-LAYER-DEBTS.md) has an [independent scoped PASS](research/cost/same-heat-refinement-20261004/randomized-weak/independent-audit/INDEPENDENT-SCALAR-WEAK-FLOW-AUDIT.md) for a real terminal weak-law gain and legal stored queries. Earlier-layer feedback still gives a linear-in-order certificate, and general-dimensional/recursive source ports remain open. The [all-layer audit](research/cost/same-heat-refinement-20261004/randomized-weak/independent-audit/INDEPENDENT-ALL-LAYER-CURRENT-AND-RESERVE-AUDIT.md) retains the paired currents and conditional/path/observer qualifications.
- The [paired scalar VALUE test](research/cost/same-heat-refinement-20261004/paired-value-test/TWO-LAYER-PAIRED-VALUE-AND-COVARIANCE.md) supplies an endpoint-scoped positive scalar covariance companion but leaves a nonlinear mean current. Its [matrix screen](research/cost/same-heat-refinement-20261004/paired-value-test/MATRIX-RANK-ONE-GAP-AND-HEAT-WINDOW.md) exposes dimension/variance-gap costs of that empirical rank-one construction, without ruling out native Jacobian/Gram reserves or other algorithms.

These are scoped source/interface improvements. They do not close a generic all-rank compiler or improve the admitted general c_P family.

### Internal packet comparison remains a distinct gate

The [concrete K3 packet audit](research/no-copy-rank-20261004/internal-packet-audit/INTERNAL-K3-ANCESTOR-TERMINAL-COMPARISON-CYCLE.md), its [paired same-potential witness](research/no-copy-rank-20261004/internal-packet-audit/PAIRED-K3-COHERENT-DRIFT-REGENERATION-GAP.md), and [carrier-rescaling test](research/no-copy-rank-20261004/internal-packet-audit/LITERAL-K3-CARRIER-RESCALING-CUT-FAILURE.md) show why particular internal primitive comparisons and complete-bank substitutions fail. The limiting paired word gap is about 5.36×10⁻⁵ at clock t=3/4 and survives positive clock averaging; no W2 lower bound is inferred from it.

The narrower missing producer retains only captured caller θ, coarse root R and incoming passive p, and may integrate its fine records internally. It must still establish the joint selected-word/skew law. This is not a blanket obstruction to a different complete mean or covariance repair.

### Main remaining generic no-copy gate

The [matrix same-query word contract](research/p-native-twin/MATRIX-K3-SAME-QUERY-TWO-SIDED-WORD-CONTRACT.md) and [ambient gradient lift](research/p-native-twin/JOINT-GRADIENT-LIFTS-AND-CONSTRAINT-DESCENDANTS.md) isolate the missing interface: an actual positive VALUE construction/current on the original nonlinear constrained-input graph, with the actual within-packet shared-root/fine-record identities, two-sided matrix words, retained/caller returns, and complete endpoint composition. The selected CW7 external physical host does not add an arbitrary-private-observer requirement. A constant covariance correction alone does not establish this interface.

Counterexamples and obstruction audits are retained because they rule out shortcuts such as replacing common nonlinear queries by independently heated factors, ignoring companion fields, or inferring marked energy from an ordinary first bound. The cost reconstructions also retain the unresolved growing-order efficiency boundary.

### Separate hidden-family cost and next target

The old uniform known/hidden-family cost certificate remains **c_P asymptotic to 2P²** under the source-relative suffix-splice recurrence, with the work convention Q_P(A) ≤ Λ_P A^(−c_P). See [cost optimization](research/cost/COST-OPTIMIZATION-PROOFS-RECONSTRUCTED.md) and [growth analysis](research/cost/GROWTH-AND-SUBCRITICAL-TARGETS-RECONSTRUCTED.md). This is quadratic growth in the cost exponent; it does not prove a subpolynomial-dimension sampler. The next target is to remove the linear-in-P replication cost per upgrade while preserving the complete source, retained/caller, and endpoint contracts. Restricted repairs and the scalar K3 closure do not yet achieve that general target.

The [fixed-seed rank and sharing audit](research/cost/recurrence-audit-20261004/FIXED-SEED-RANK-CORRECTIONS-AND-SHARING-AUDIT.md) sharpens the priority: rebuild a complete all-rank no-copy compiler from the cheap c_7=0 seed. Conditional on that still-open source construction, finite fixed-target duplication and the existing finer-heat terminal preserve c_P=0. Same-heat additive posterior refinement is a separate alternative. Perfect sharing of the old empirical core alone still gives c_P=P−5.4 in the optimistic certificate, which is linear and misses the sublinear target. These are cost implications and research directions, not a completed all-order theorem.

## Provenance and integrity

- `INVENTORY.json` records each mathematical file's public SHA-256, original SHA-256, original/public byte counts, version type, and whether it is a sanitized publication copy.
- `unchanged_source_copy` means the public bytes equal the source snapshot bytes. `sanitized_publication_copy` means local paths, nonmathematical context, and/or explicitly identified independent-audit scope clarifications were edited. Mathematical formulas and historical source/audit pins were retained; those pins do not silently become hashes of the edited files.
- `historical_bytes_verified`, `reconstructed_new_version`, and `new_research_or_diagnostic` describe the original version. A reconstructed document is new work rather than a byte-identical historical source.
- `PROVENANCE-CHECKS.json` records 308 source-index bindings, with original and public hashes. One older binding names the pre-correction mixed-descendant source. Its final source and scoped audit are explicitly identified. Source-index identifiers are checksums of external provenance records, not additional bundled sources.
- `SHA256SUMS` covers every bundled file except itself. Run `python verify_integrity.py` before use.

Private correspondence, external confidential manuscripts and derivative commentary, personal data, credentials, and unrelated material are excluded. No publication license is assigned.

## Diagnostics

This is mathematical research, not a production sampler. Install `requirements.txt` in an appropriate Python environment, then run `python run_diagnostics.py`. The runner executes each script in its own disposable copy and never overwrites research source or saved results. Its default invocation only prints results. Use `--report /your/chosen/report.json` to save an additional report.

`verification/DIAGNOSTIC-RERUN.json` records an independent fresh run of the publication copies, environment versions, exact outcomes, and any differences from saved outputs. `verification/REVIEW.md` states the public-release and verification limits. Passing arithmetic and numerical checks supports the stated identities and test cases; it does not replace analytic proof or discharge any open interface.
