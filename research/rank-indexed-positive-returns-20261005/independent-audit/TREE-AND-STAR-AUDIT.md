# Independent audit: labelled trees and the rank-five star return

2026-10-05. Companion to `AUDIT.md`.

## Verdict

**PASS with the boundaries stated in the source note.** The labelled analytical tree generator is exact; its heat and positive-clock target-error constants are conservative; the proposed blind native extension loses its displayed uniform small-radius/first certificate; and the formal five-star return really emits a six-force/four-mark two-cycle current.

The star calculation is a finite formal moment-cumulant calculation. It is not an existence claim for the unbounded polynomial's real moment-generating function. The endpoint budget T9 is an upper allowance, not a lower bound on native Wasserstein error. Neither the ideal stationary guard test nor the formal return rules out a differently designed finite compiler.

## 1. Exact connected recurrence and AST audit

The imported discounted-OU generating equation gives, after taking the finite logarithmic jet,

    (n−L_OU) κ_n = sum_(i=1)^(n−1) binom(n,i) Dκ_i · Dκ_(n−i),

with `κ_1=R_1 h`. The OU commutation identity is `D R_k = R_(k+1) D`. Therefore a derivative chooses each force occurrence in turn and increments every resolvent on the path to it. This is exactly the author's AST operation.

Relabelling is correct: the left rank-i history uses force labels 1,…,i and old edge labels 1,…,i−1. The shifted right history uses force labels i+1,…,n and old edge labels i,…,n−2. The new edge is n−1. Every edge consequently has two distinct endpoints, and joining two disjoint trees makes a tree.

A useful stronger clock invariant is

    clock index = number of forces in its subtree
                  + number of spatial edges crossing that subtree boundary.

It holds at initialization and is preserved by each product-rule derivative. The independent checker verifies it for every one of the 2,448 resolvent occurrences in the 272 serialized rank-five histories. This is a check of the actual genealogy, not merely the total number of clocks.

The clock count is n leaf resolvents plus n−1 structural resolvents, namely 2n−1. Each force's derivative count equals its spatial tree degree. Thus the total extra heat derivative count is `sum_v(d_v−1)=n−2`.

Every individual unmerged history has coefficient exactly n!: the product of the binomial factors at internal nodes telescopes to n!. Hence the stronger identity `S_n=n! T_n` holds. Independently computed counts through rank eight are:

    T_n = 1, 1, 4, 28, 272, 3312, 47872, 794880,
    S_n = 1, 2, 24, 672, 32640, 2384640,
          241274880, 32049561600.

The serialized rank-five shapes are exactly 112 paths, 144 forks, and 16 four-leaf stars. Evaluating all serialized histories using an independent Hermite-eigenbasis OU resolvent agrees with the moment-cumulant recursion for two different scalar polynomial fixtures. The fixtures test finite algebra only.

## 2. Dimension safety of the tree and heat bounds

For a degree-d primitive, full symmetry of `D^d g` allows one tensor slot on each side of any desired proper cut to serve as the two slots of `Dg`. Split the private heat into two independent Gaussian halves and place the remaining derivative slots on the two sides. Gaussian Hermite isometry, the pointwise `||Dg||op<=A` bound, and Cauchy–Schwarz give a proper-cut bound no larger than

    2^(d−1)(d−1)! A tau^(-(d−1)).

For the full Hilbert–Schmidt bound, retain `||Dg||HS<=A sqrt(D)` as the one energy factor; the Gaussian derivative tensor is controlled through its covariance operator. The same displayed constant is a safe majorant.

An all-marked tree has a dimension-free proper-cut bound equal to the product of its local proper-cut bounds. One way to see the operator statement is to color its physical marks by the two sides of the desired cut. Contract the monochromatic components abstractly and orient every cross-color edge from input to output. Within each input component, orient toward one outgoing boundary edge; within each output component, orient away from one incoming boundary edge. Every local tensor then has a nonempty input and output cut. Because the underlying graph is a tree, the orientation is acyclic, and topological composition of these tensor operators costs only the product of their operator norms.

For a Hilbert–Schmidt bound, instead start with one coefficient in HS norm and attach vertices one at a time through their proper cuts. This is the same one-energy mechanism as the positive-consumer proof, and every physical mark survives.

Multiplying the primitive constants gives

    product_v 2^(d_v−1)(d_v−1)!
       <= 2^(n−2)(n−2)! = B_n.

Positive resolvent averaging preserves the required pointwise cuts and contracts Gaussian Hilbert L² with its positive mass. The theorem uses these facts only at the actual common caller and with the original clock genealogy.

## 3. T1 and its scope

A scalar Hermite-multiplier error δ₀ gives a Hilbert-valued Gaussian L² operator error δ₀ at any one clock. Telescope the actual 2n−1 clocks. At the clock carrying the error, retain the Hilbert–Schmidt bound. On other factors use pointwise proper cuts, and on other clocks use the positive rule's mass bound 2. This yields a bound no greater than

    (2n−1) 2^(2n−1) S_n B_n δ₀ A^n tau^(-(n−2)) sqrt(D).

This is T1. It is an integrated standard-caller Gaussian L² statement; it does not infer a pointwise quadrature estimate or an Lp multiplier bound from an L² spectral supremum.

At rank five the explicit constants are

    B5 = 48,
    D5 = 9 × 2^9 × 32640 × 48 = 7,219,445,760.

The supplied positive dyadic-rule assumption gives the stated fixed-rank logarithmic node count after choosing `δ₀=epsilon tau^(n−2)/D_n`. It is a finite analytical coefficient approximation, not a proof that the surrounding native packet satisfies its feedback, first, curl, or observer requirements.

T2 follows from the audited Appell inequality and same-bank heat coupling. Hermite coupling of the coefficient target has the stated `alpha^n` normalized scaling. This does not eliminate the separate native inverse-shield costs.

## 4. Blind rank-five amplitude test

For a five-star, the center supplies three inverse private-shield powers. With root-weight normalization T3, a panel of mass Δ has `sigma²` comparable to `Δ+tau²` and root scale `alpha Δ/(Δ+tau²)^(3/2)`. Its maximum over endpoint scales is order `alpha/tau`. The dyadic sums for root first and center caller first give the stated `alpha/tau` and `alpha²/tau²` powers.

Thus the blind choice `tau=alpha` cannot certify an O(alpha) root/first radius by these bounds. The note correctly distinguishes the panel power certificate from a lower bound at an arbitrary individual quadrature node.

The scalar linear-source ideal test is exact: C3 vanishes, a C0 target is constant one, and the center becomes a fresh Gaussian. The residual derivative of

    rho Z_center + sqrt(1−rho²) Z_root

relative to its prescribed source-zero carrier `Z_root` is exactly rho in the center-root direction. A positive readout does not remove it. This tests the stated stationary/native-graph extrapolation. It does not prove that every finite realization must choose this carrier or amplitude allocation.

## 5. Exact star return and original-gradient witness

With independent standard P,S,T,

    V(theta)=lambda theta²(P+b theta)(S+u theta)(T+v theta),

independence gives

    E V = lambda b u v theta⁵,

    (E V² − (E V)²)/2
      = lambda²/2 [theta⁴ +(b²+u²+v²)theta⁶
          +(b²u²+b²v²+u²v²)theta⁸].

The degree-six term inside the product bracket is the all-mean contribution; after multiplication by the external theta⁴ it equals `(EV)²` and cancels. It does not cancel the theta⁴ term. Setting u=v=0 leaves exactly

    lambda²(theta⁴+b²theta⁶)/2.

Hence both the six-force fourth current and the six-force sixth current survive. They are even in lambda. Side Gaussian probes continue to contribute when their deterministic reverse shifts vanish.

For arbitrary n, the same independent-factor second-moment calculation gives T6. For every fixed k, the shifted Gaussian moment formula and ordinary finite moment-cumulant recurrence give T8. The independent checker compares that recurrence with the full partition formula through k=6.

The explicit original-gradient witness is valid:

    g(x)=Ax/2 + A sin(kx)/(4k),
    A/4 <= g'(x) <= 3A/4,
    g''''(x)=A k³ sin(kx)/4.

Its private Gaussian smoothing gives exactly the stated C3 target

    (sigma k)³ exp(−(sigma k)²/2) sin(kq)/4.

Choosing k=1/sigma gives a nonzero order-one normalized target at suitable ordinary captured centers while retaining the global Hessian sandwich. This justifies a nontrivial coefficient in the admitted original-gradient class. It remains separate from a calibration proof for the finite clipped native field.

No ordinary MGF is needed: all these statements use finite polynomial moments and their formal logarithmic jet. The warning about possible failure of real exponential integrability is important and correct.

## 6. Two-cycle graph, grouping, and the endpoint allowance

Two three-force spines supply six force vertices and four internal spine edges. Pairing their own-center public and their two side probes adds three parallel center-center edges. Thus

    V=6, E=7, beta1=E−V+1=2, M=4.

The two center physical marks are consumed; the four leaf physical marks survive. The primitive derivative orders are 4 at both centers and 1 at each of four leaves, giving

    sum_v(j_v−1)=6=M−2+2 beta1.

This is genuinely a two-cycle/four-mark return, not an all-marked rank-six tree. Grouping each intact spine preserves two physical outputs and excludes an initial same-group Wick self-contraction: each independent Gaussian color occurs only once in that group. The resulting connected conditional diagrams admit the restricted one-Hilbert merging bound. This does not establish arbitrary graph self-traces, whole-old-bank covariance controls, or native source admission.

Under the stated root-weight prescription, the squared endpoint sum has the bound

    alpha^6 sum omega²/(Delta+tau²)^3 <= C alpha^6/tau².

At each dyadic panel, positive weights satisfy `sum omega² <= (sum omega)² <= C Delta²`; summing `Delta²/(Delta+tau²)^3` gives the displayed inverse-square heat allowance. Additional original structural-clock weights stay squared and are not replaced by independent resampling.

T9 is correctly used as a conservative allowance. Balancing this allowance alone with `alpha^5 tau` gives `tau=alpha^(1/3)` and grade 16/3. No lower bound on W2 and no completed rank-five packet follows. The note explicitly retains the other first, source, caller, auxiliary, and old-bank obligations.

## 7. Check artifacts

`check_tree_star.py` imports no author generator. It passes 6,014 exact assertions and writes `tree_star_checks.json`. It verifies serialized labels and clocks, coefficients and shapes, independent polynomial OU evaluations, finite star cumulants, graph topology, and the original-gradient derivative witness.

Together with the Appell/current audit, the two independent checkers pass 8,674 exact assertions. General dimension-free results rest on the proofs above and in `AUDIT.md`, not on these finite checks. Final input hashes are recorded separately.
