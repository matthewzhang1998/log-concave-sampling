# Independent review of the collective-heat Task B packet

Date: 5 October 2026.

## Verdict

**PASS for the revised, explicitly bounded claims. No remaining blocking mathematical error was found.** The packet establishes an exact positive-covariance heat regrouping, a finite positive quadrature certificate, and a guarded finite original-VALUE representative of the heated **leading current**. It identifies the complete Gaussian-adjoint boundary and gives a fixed convex C2 counterexample to a universal extra power for the saturated nodewise heat generator.

This verdict does **not** supply an exact all-orders native right inverse, strict reserve gain, a terminating arbitrary-order original-target repair, or nonmembership in the full executable current image. The original Task B success alternative in the algebra packet demanded more than the limited outcome obtained here. Its strict-gain/full-image question remains open. The present report correctly isolates that boundary rather than claiming to settle it.

The review concerns the author files at the hashes recorded in `AUDITED-FILES.json`. All 15 sealed input pins were independently verified. Author inputs were not modified by this audit. The author made the corrections described below during review. No native sampler, cold-start program, or paused numerical target program was run. No external sharing was performed.

## 1. Operator identity, finite jets and covariance interpolation

For the independent formal centers,

- `H_com = H_ind + B_cross` has the correct factor and sign.
- All three operators have constant coefficients and commute. The finite jet factorization is an exact identity in the truncated polynomial algebra; it does not need convergence of an infinite heat series.
- Each choice in `B_cross` connects distinct original occurrences. Repeating a pair adds parallel edges, not a self-loop. Repeated hits and multinomial multiplicities are retained.
- `Q_theta = (1-theta)I + theta 11^T` is positive semidefinite, including the singular endpoint theta=1. The original independent private shields leave a nondegenerate smoothing layer throughout this interpolation.
- The covariance generator derivative is `s B_cross`, because both symmetric off-diagonal entries contribute to the factor 1/2. Consequently `P_ind = P_com - s integral P_Qtheta B_cross` has the correct sign.

The proof works at a fixed old coefficient record with positive private shields and sufficient bounded derivative envelopes from the smoothed C2 sources. It does not differentiate the square-root representation at theta=0 or theta=1. The covariance differentiation formula is the appropriate argument there.

The revised unequal nonnegative-rate extension is also correct: weighting common directions by square roots of the rates gives the claimed diagonal and off-diagonal covariance entries. Its quantitative replacement p_lambda=sum_(u<v)sqrt(lambda_u lambda_v) is the correct total edge weight in every repeated-bridge bound, and each actual added caller row has norm sqrt(s lambda_v). If p_lambda=0, there is no bridge contribution or quadrature. A signed rate list is not a positive covariance path. Section 9 now correctly applies `2H_com-2B_cross` to the complete equally weighted sum of local Laplacians, not to an individual local Laplacian.

## 2. Quadrature certificate and actual node cost

Equation (8) follows termwise from the local all-proper-cut source constants, the retained marked spanning tree and the loopless cyclic contraction theorem. There are `p^a` ordered pair choices, and the local derivative orders sum to `K+2a`. Bounding the product of their factorials by `(K+2a)!` is valid. A one-versus-rest proper cut also supplies a local Hilbert bound with one factor sqrt(D); the graph attachment proof propagates only that one Hilbert factor.

For completeness, the factorial step can be proved without an asymptotic estimate. Writing `A=K+2`,

`(A+2j)!/(A!(j!)^2) = binom(A+2j,A) binom(2j,j) <= 2^(A+4j)`.

Together with the source factor `2^((K+2j+2)/2)`, this gives exactly the displayed `M` and `R0=tau^2/(8ps)`.

On each panel of length at most R0, Taylor expansion of degree `2m-1` around the midpoint has norm error at most `M (1/2)^(2m)=M 4^(-m)`. Gauss positivity and exactness bound the difference between the integral and quadrature of the remainder by twice this value times the panel mass. Summing panel masses gives the advertised `2sM4^(-m)` coefficient error. This is a real-variable Banach-norm argument; no complex extension is required.

The rounding estimate is valid provided weighted node error is defined using the weights in the chosen decomposition of the perturbed rule. Rounded nodes remain in their panels, and positive normalized weights preserve each exact panel mass. `sup ||K|| <= M` and `sup ||K'|| <= M/R0` give the stated bound directly.

At fixed graph and bounded `s/tau^2`, the panel count is bounded independently of alpha. The number of points per panel is logarithmic when Q, inverse shields and the requested tolerance are fixed powers of alpha. This certifies zero **additional theta-node** inverse-alpha exponent. It does not certify zero full query exponent, uniform-in-order constants, an uncharged large heat jump, or a native-law error of the same size as the coefficient error.

The revised text now fixes Q as the **original restored target envelope before node selection**. New native adapter/readout inverses that cancel in that target belong in the guard, error and execution bill. They cannot be put back into Q to create circular node selection. This matches the pinned execution-bill source.

## 3. Native same-graph edge admission and ownership

The same-graph assertion is supported by more than a proper-cut estimate. The relevant imported chain is:

1. `MARKED-SPANNING-TREE-EXTENSION.md`, especially Sections 2–5, 7 and 9: an original marked spanning tree survives a nontree chord or parallel edge; each extra edge has its own Gaussian cut color; selected-slot symmetry produces the literal promoted C_k source; the finite root cluster gives the graph.
2. `REVIEWED-MARKED-SPANNING-TREE-EXTENSION.md`, Sections 4, 5 and 7: original-VALUE realization, source-qualified conditional offspring, selected-spine/frame and finite-current obligations.
3. `FROZEN-COEFFICIENT-WIDTH-ZERO-PORT.md`, Sections 1–4: all three dynamic derivative frames, repeated-color product rules and Hölder control, the genuine same-endpoint root-quadratic return, precision qualifications, and exclusion of captured-caller derivatives from the zero-extra-width claim.
4. The root-seed theorem and ownership interface: addition with a shared new source-zero carrier is covered only when the chosen root has no external dynamic probe and its selected mean is measurable in the frozen old record alone.
5. The weighted/heat-shell sources and execution bills: independently owned-group addition, positive physical shares after the complete census, original source versions, one untouched keep, full old-bank observer firsts and full replay remain available under their literal guards.

The packet observes these boundaries. Correlated heat increments are coefficient-center rows, not native cut colors. Their displayed row/stacked operator bounds are correct. An edge incident to the chosen root can defeat the root-seed hypotheses; the packet properly uses the independent-share alternative instead of asserting a stronger carrier theorem. Old retained observers are not erased by introducing a fresh increment bank.

The active source count for one pair and one theta node is correct: precisely the two hit occurrences double their `2^(k_v+1)` counts. The shell count `(N+1) sum_v 2^(k_v+1)` is the correct uncombined count when each of N telescope terms evaluates one doubled marked local source. Neither formula includes native/filter repetition or caller/ancestor replay, and the packet correctly charges those separately.

This is a contract audit of admitted fixed-graph ports, not an independent construction of omitted all-order native constants. Actual radius, three-frame, covariance-gap, positive readout, source numerical-first and precision inequalities must pass for the instantiated program.

### Intrinsic return versus accuracy floor

The original wording about a right inverse to arbitrary finite floors was too easy to read as removal of the native intrinsic error. The revised Section 6 is correct and materially safer: coefficient quadrature and local preparation/filter floors can be tightened, but the intrinsic same-endpoint root-quadratic current remains. Above that floor, its current coefficients must be supplied and repaired or explicitly retained as unresolved. Thus a finite representative of a leading coefficient is not an exact all-orders finite native-law right inverse.

## 4. Full Gaussian adjoint and conditional boundary

For constant directions `D_i=a_i dot D_B`, the standard Gaussian adjoint is `D_i*=beta_i-D_i`. Squaring it gives

`(D_i*)^2 = beta_i^2-|a_i|^2 - 2 beta_i D_i + D_i^2`.

Applying this to the observable at the same endpoint yields exactly the score, `W_i/2-beta_i V_i`, and `V_i tensor V_i/2` rows of (12). Both signs and coefficients are correct. No assumption that the directions are orthogonal or unit is necessary.

Equation (13) is the full scalar product rule for a fixed geometry. Widths, heat rows, normalizations, scalar clocks, contractions and caller geometry require their own chain/product terms if differentiated. The packet explicitly excludes them from the frozen operator and requires their restoration in a live realization. A live-width diagnostic in this audit has a nonzero missing term even at a simple caller value, illustrating that this warning is substantive.

Applying the same adjoint to `chi(E) D^r phi(X)` gives exactly (14). Its mixed coefficient has factor one, not one-half. Conditional equality cannot be inferred from the unconditional identity by merely retaining the differentiated Gaussian bank. The sufficient condition `D_i E=0`, or a proved equivalent annihilation/current identity, is correctly stated. The first-integration formula avoids W but leaves a real score-times-first-coefficient row.

A common-translation lift must actually solve the injection-row equations. The example `z1=B`, `z2=2B` has no allowed scalar direction satisfying both lifts. A fresh common score representation for positive b is valid and has the stated score norm `sqrt(2D)/(2b)`; it does not create a free source or width gain.

The untouched-keep formula (16) correctly lowers physical test rank. It requires coefficient independence of the keep and a fixed affine endpoint derivative `D_Z X=eta I`. It neither erases coefficient heat nor certifies score/endpoint-tensor source bounds. Rank-zero rows are outside this operation.

## 5. Fixed C2 counterexample and scope

The series for g and g' converge uniformly. With the stated epsilon, g' is continuous and strictly between 1/4 and 3/4, so U is one fixed convex C2 potential with bounded Hessian. The construction does not change g with alpha, and all higher derivatives used in the proof belong to positively heated coefficients.

At zero every term in f_n is negative and every term in f_n'' is positive. The n-th frequency alone therefore gives the trace lower bound. The cross coefficient at that point is zero by oddness. In the full finite independently heated product, both `|f_n|-|f_n,warm|` and `|f_n|+|f_n,warm|` are positive sums, so their product gives the lower bound for the **complete finite heat difference**. This is stronger than an isolated trace-norm calculation.

The rank-six current is nonzero in the fixed-endpoint observable quotient: one compactly supported cutoff of x^6 with sufficiently large support has nonzero expected sixth derivative. The six-vertex decorated variant has correct valences, K=4, five internal edges, six physical marks and marked tree leaves. Its affine C0 leaves contribute fixed nonzero constants and no heat derivatives.

The sequence `n^-4 2^(delta n^2)` diverges for every delta>0. Hence this fixed source defeats a universal extra width power at saturated clock charge. The obstruction remains after fixed public-log attenuation.

The audit requested a distinction between weights merely allowed by an upper node-mass inequality and weights of the actual sealed Gauss grid. The revised packet makes it explicit. `omega=kappa t^2` is first a contract-admissible sequence. The additional actual-node construction is valid under the stated hypotheses: a positive m-node panel of mass d supplies at least one weight at least c d/m, while `h^2<=d` implies `t^2` is comparable to d. Frequencies fixed as inverses of those selected widths give a lower bound weakened by `m_n^-2`, still larger than every width power for public-log m_n along a sufficiently fast width sequence. This argument does not assert a lower bound for every Gauss weight or silently assume that the necessary panels were included.

The conclusion is nodewise/fixed-record and universal-contract negative. It does not rule out an additional complete-history or structural-clock cancellation for a specified original target. Nor does nonzero current pairing separate the example from the full source image. Section 11 and the JSON correctly disclaim that stronger inference.

## 6. Independent diagnostics

`independent_checks.py` is independently written and does not import the author test implementation. It checks:

- All 15 sealed input hashes and sizes.
- Exact two-dimensional, three-occurrence heat identities, commuting operators, finite jets and covariance homotopies, including unequal and zero rates.
- Exact two-dimensional nonlinear-endpoint Gaussian adjoints with nonorthogonal directions, scalar product rules and retained-test defects.
- A live-width second-chain-rule defect.
- Exact integer factorial bounds beyond the author's parameter range.
- All proper cut norms and one-Hilbert bounds for deterministic multigraph tensor fixtures with cycles and parallel edges. These are tensor diagnostics, not native-program tests.
- High-precision Fourier homotopies, including vanishing total common frequency and singular endpoint covariance.
- A three-occurrence nonpolynomial positive composite-Gauss certificate, including varying private width and panel count.
- Positive-term lacunary lower-bound diagnostics and log-domain witnesses against several positive powers.

The independently implemented suite passes **157 checks**. In the nonpolynomial Gauss examples, observed errors are at most approximately `1.45e-17`, with prescribed certificates at most `7.91e-9`. The retained-observer test deliberately produces a nonzero omitted pairing, `-1866887/64`, before the exact defect is restored. One lacunary trace comparison uses a `1e-65` tolerance at 70-digit precision because the omitted positive terms are smaller than working precision; the mathematical lower bound follows exactly from the same-sign proof, not from that floating comparison.

The author's script was separately copied into `author-reproduction/` and rerun there, so its output writes did not touch author artifacts. Its **879 checks** reproduced successfully. Diagnostics supplement the proofs and scoped imported contracts; they do not certify a native sampler, full target cancellation or universal current-image exclusion.

## 7. Changes verified during review

The reviewed author revision incorporates four material clarifications:

1. The local-Laplacian replacement is a complete-sum identity; the unequal nonnegative-rate extension is stated separately.
2. The quadrature Q is fixed from the original restored target, with cancelling new adapter/readout factors kept outside it.
3. Finite heated leading-current realization is separated from the intrinsic root-quadratic return and from an exact all-orders native right inverse.
4. Contract-admissible saturated weights are distinguished from actual Gauss nodes, and the latter receive explicit panel, heat and public-log count hypotheses.

No further correction is required for the scoped mathematical result. The remaining common-boundary source/control theorem and any strict-gain original-target closure are genuinely unresolved work, not consequences of this audit.
