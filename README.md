# Log-concave sampling: mathematical research snapshot

Snapshot version: 2026-10-04-public-v4-host-separator-randomized-weak, updated 4 October 2026. This is a curated research archive of proofs, scoped audits, source certificates, and diagnostic code. Inclusion is not an assertion that every proposed construction is proved.

## Progress in brief

The source-relative exact **1/32** baseline remains the established reference. Several finite, bounded source and covariance modules are now independently checked. The general matrix constrained-input comparison and its complete nonlinear endpoint join are still open. No new all-order, subpolynomial-dimension, or improved end-to-end sampler complexity theorem is claimed here.

### Baseline

[Exact denominator 32, version 2](research/exact-slack/exact32-proof-v2.tex) is preserved byte-for-byte, SHA-256 `c7487ca24ffa44ea73b52e42e4b6f715ad964171fa5bd0a61cf9323eb646cbff`.

Its theorem is relative to the pinned LOW30/LOW31 foundational constructors, assumes a supplied initialization, and counts original gradient and directional Hessian-vector queries. It does not claim dimension-free arithmetic/storage or a differentiable return for the new packet. Those foundations are referenced rather than silently treated as included, independently reproved sources. The historical componentwise R7 and fourth-VALUE source certificates are under `research/exact-slack/`.

### Closed bounded modules

Each claim below retains the assumptions and numerical/source costs in its proof and audit.

- Known affine common-input heat and one-pair skew adapters: `research/p-affine-heat/`, `research/c-common-input/`, and their independent audits.
- Explicit nonlinear shear sources, constant-reference covariance corrections, and local retained-graph qualifications: `research/p-shear-heat/`, `research/r-shear-consumer/`, `research/c-shear-audit/`.
- Exact finite two-node third-iterate source algebra and the scalar direct true-Gram reserve: [scalar construction](research/r-native-twin/SCALAR-K3-DIRECT-TRUE-GRAM-RESERVE.md), with separate P/C audits. The scalar mean corollary is local; it is not a complete posterior host.
- Auxiliary variance-preserving rotation, genuine-gradient constraint offspring, and rank-aware restoration under the stated Hessian bounds: the `auxiliary-heat/` audit set and `ROTATION-CONSTRAINT-DEFECT-COMMON-GRADIENT-LIFT.md`.
- The small-full-gradient constant mixed-covariance service now has a [scoped independent PASS](research/c-native-twin/auxiliary-heat/INDEPENDENT-SMALL-GRADIENT-MIXED-DESCENDANT-AUDIT.md), binding final source `f805e8d3…`. Absolute finite-prior and numerical floors remain explicit.

### New scoped progress

- [Retained mixed-field extraction](research/no-copy-rank-20261004/RETAINED-MIXED-KERNEL-AND-TERMINAL-FLAT-GATE.md) has an [independent scoped PASS](research/no-copy-rank-20261004/INDEPENDENT-RETAINED-MIXED-EXTRACTION-AUDIT.md): a positive conditional covariance reserve retaining the same roots, with no additional heat-power replication cost at fixed target under its stated clock/compiler assumptions. The companion retained-root observer obstruction and terminal-flat K3 witness show why the full ambient endpoint/current identification remains necessary.
- [Same-heat quarter-flow refinement](research/cost/same-heat-refinement-20261004/SAME-HEAT-POSTERIOR-REFINEMENT-AND-FLOW-GATE.md) has an [independent scoped PASS](research/cost/same-heat-refinement-20261004/independent-audit/INDEPENDENT-SAME-HEAT-FLOW-AND-QUERY-AUDIT.md) for exact posterior invariance, contraction by a/(1−a), and a legal finite VALUE/HVP DAG. Its law-level/query certificate c_R ≤ max(c_P,P−1) is still linear in order and is not an admitted family recurrence while the protected/retained-host join remains open. The deterministic information lower bound is limited to fixed-input first-Picard integration and adaptive baseline transcripts; it is not an exact-flow, randomized, Gaussian-input Lp, or posterior-law lower bound. Prefix-sum arithmetic does not change the query exponent.

These results do not improve the currently admitted general c_P family or prove an all-order/subpolynomial sampler.

### Latest interface and source results

- **The selected CW7 external host needs a smaller interface.** The [host-separator audit](research/no-copy-rank-20261004/host-observer-audit/ACTUAL-CW7-WEAK-E-HOST-SEPARATOR-AND-COMPARISON-ORDER.md) and [independently audited reentry recipe](research/no-copy-rank-20261004/host-observer-audit/INDEPENDENT-COMPLETE-MEAN-REENTRY-AUDIT.md) show that its physical law consumer reads captured caller θ and one complete statistic T. A complete conditional mean law can be proved first and the host attached afterward. An arbitrary joint-fine-observer theorem is not an external prerequisite for this selected host. Internal packet aliases and any baseline still reading integrated fine variables remain genuine obligations; the actual retained-source/graph proof is separate.
- A [positive projected-gradient/V-column construction](research/no-copy-rank-20261004/endpoint-next/v-column-escape/PROJECTED-GRADIENT-V-COLUMN-ESCAPE-AND-JOINT-RETENTION.md) succeeds on the terminal-flat fixture, with actual-energy control and an exact physical orientation selector. Its zero-auxiliary-root Gram identity is not a full retained-root endpoint theorem. Earlier midpoint, fixed-frame and baseline obstructions remain limited to their stated operations. The [cache partition](research/no-copy-rank-20261004/endpoint-next/cache-partition/WEAK-29-CAPTURED-CALLER-AND-COMPLETE-PRIVATE-CACHE-PARTITION.md) allows exact shared captured work but keeps complete private bank work charged.
- [Scalar randomized quarter-flow](research/cost/same-heat-refinement-20261004/randomized-weak/SCALAR-RANDOMIZED-QUARTER-FLOW-AND-ALL-LAYER-DEBTS.md) has an [independent scoped PASS](research/cost/same-heat-refinement-20261004/randomized-weak/independent-audit/INDEPENDENT-SCALAR-WEAK-FLOW-AUDIT.md) for a real terminal weak-law gain and legal stored queries. Earlier-layer feedback still gives a linear-in-order certificate, and general-dimensional/recursive source ports remain open. The [all-layer audit](research/cost/same-heat-refinement-20261004/randomized-weak/independent-audit/INDEPENDENT-ALL-LAYER-CURRENT-AND-RESERVE-AUDIT.md) retains the paired currents and conditional/path/observer qualifications.
- The [paired scalar VALUE test](research/cost/same-heat-refinement-20261004/paired-value-test/TWO-LAYER-PAIRED-VALUE-AND-COVARIANCE.md) supplies an endpoint-scoped positive scalar covariance companion but leaves a nonlinear mean current. Its [matrix screen](research/cost/same-heat-refinement-20261004/paired-value-test/MATRIX-RANK-ONE-GAP-AND-HEAT-WINDOW.md) exposes dimension/variance-gap costs of that empirical rank-one construction, without ruling out native Jacobian/Gram reserves or other algorithms.

These are scoped source/interface improvements. They do not close a generic all-rank compiler or improve the admitted general c_P family.

### Main remaining gate

The [matrix same-query word contract](research/p-native-twin/MATRIX-K3-SAME-QUERY-TWO-SIDED-WORD-CONTRACT.md) and [ambient gradient lift](research/p-native-twin/JOINT-GRADIENT-LIFTS-AND-CONSTRAINT-DESCENDANTS.md) isolate the missing interface: an actual positive VALUE construction/current on the original nonlinear constrained-input graph, with the actual within-packet shared-root/fine-record identities, two-sided matrix words, retained/caller returns, and complete endpoint composition. The selected CW7 external physical host does not add an arbitrary-private-observer requirement. A constant covariance correction alone does not establish this interface.

Counterexamples and obstruction audits are retained because they rule out shortcuts such as replacing common nonlinear queries by independently heated factors, ignoring companion fields, or inferring marked energy from an ordinary first bound. The cost reconstructions also retain the unresolved growing-order efficiency boundary.

### Best general cost and next target

The strongest general cost certificate in the current record has **c_P asymptotic to 2P²** under the source-relative suffix-splice recurrence, with the work convention Q_P(A) ≤ Λ_P A^(−c_P). See [cost optimization](research/cost/COST-OPTIMIZATION-PROOFS-RECONSTRUCTED.md) and [growth analysis](research/cost/GROWTH-AND-SUBCRITICAL-TARGETS-RECONSTRUCTED.md). This is quadratic growth in the cost exponent; it does not prove a subpolynomial-dimension sampler. The next target is to remove the linear-in-P replication cost per upgrade while preserving the complete source, retained/caller, and endpoint contracts. Restricted repairs and the scalar K3 closure do not yet achieve that general target.

The [fixed-seed rank and sharing audit](research/cost/recurrence-audit-20261004/FIXED-SEED-RANK-CORRECTIONS-AND-SHARING-AUDIT.md) sharpens the priority: rebuild a complete all-rank no-copy compiler from the cheap c_7=0 seed. Conditional on that still-open source construction, finite fixed-target duplication and the existing finer-heat terminal preserve c_P=0. Same-heat additive posterior refinement is a separate alternative. Perfect sharing of the old empirical core alone still gives c_P=P−5.4 in the optimistic certificate, which is linear and misses the sublinear target. These are cost implications and research directions, not a completed all-order theorem.

## Provenance and integrity

- `INVENTORY.json` records each mathematical file's public SHA-256, original SHA-256, original/public byte counts, version type, and whether it is a sanitized publication copy.
- `unchanged_source_copy` means the public bytes equal the source snapshot bytes. `sanitized_publication_copy` means local paths, nonmathematical context, and/or explicitly identified independent-audit scope clarifications were edited. Mathematical formulas and historical source/audit pins were retained; those pins do not silently become hashes of the edited files.
- `historical_bytes_verified`, `reconstructed_new_version`, and `new_research_or_diagnostic` describe the original version. A reconstructed document is new work rather than a byte-identical historical source.
- `PROVENANCE-CHECKS.json` records 285 source-index bindings, with original and public hashes. One older binding names the pre-correction mixed-descendant source. Its final source and scoped audit are explicitly identified. Source-index identifiers are checksums of external provenance records, not additional bundled sources.
- `SHA256SUMS` covers every bundled file except itself. Run `python verify_integrity.py` before use.

Private correspondence, external confidential manuscripts and derivative commentary, personal data, credentials, and unrelated material are excluded. No publication license is assigned.

## Diagnostics

This is mathematical research, not a production sampler. Install `requirements.txt` in an appropriate Python environment, then run `python run_diagnostics.py`. The runner executes each script in its own disposable copy and never overwrites research source or saved results. Its default invocation only prints results. Use `--report /your/chosen/report.json` to save an additional report.

`verification/DIAGNOSTIC-RERUN.json` records an independent fresh run of the publication copies, environment versions, exact outcomes, and any differences from saved outputs. `verification/REVIEW.md` states the public-release and verification limits. Passing arithmetic and numerical checks supports the stated identities and test cases; it does not replace analytic proof or discharge any open interface.
