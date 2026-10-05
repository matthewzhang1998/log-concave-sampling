# Shared-auxiliary coefficient frames and the same-endpoint gate

2026-10-05. Addendum to the amplitude-surplus construction. This supplies the specific finite polynomial comparison step that cannot be inferred from a scalar cumulant check or a final marked-graph HS bound. It concerns the new conditional Gaussian bank at fixed complete old Q. It does not admit whole-old-bank derivative currents.

## 1. The extra induction invariant

Fix the finite center tree, covariance-root Taylor orders, finite filters, and all original coefficient roots Q. Expand the ideal forward map only to a fixed amplitude degree. Work at the collapsed-spine level. For each matrix-valued coefficient A, retain:

1. Its operator Lp moments, for every member of a finite, preselected Hölder list.
2. BOTH rectangular Gaussian-series frames of its first derivative across the complete relevant auxiliary and reset bank, and the transposed row frame needed when a matrix factor is transposed.
3. A record of an independent terminal or completion Gaussian X when a vector term has the form A X.

Here is an explicit frame convention. If A is D by D and the differentiated Gaussian bank has N D coordinates, define B_j, of shape D by N D, by

    B_j(i,k) = partial_k A(i,j).

The two derivative frames are the operator norms of

    sum_j B_j B_j^T,       sum_j B_j^T B_j.

The first equals sum_k (partial_k A)(partial_k A)^T. The square root of the second is exactly the operator norm of D A from the differentiated bank into matrix HS space; it is NOT the transposed row frame. Retain that third flattening separately:

    sum_k (partial_k A)^T(partial_k A).

Thus all three proper flattenings of the output/input/derivative slots of D A are recorded. This makes transposition a legal closure operation.

The finite number N and all needed Hölder exponents may depend on the fixed graph and cutoff, but not polynomially on D. Public dimension logarithms are charged. The desired bounds are fixed powers of public logarithms for operator moments and the square roots of these frame moments.

At a native C0/C3 coefficient, all proper tensor cuts and the required first-derivative frames are supplied by the literal privately shielded original-gradient source. A direct auxiliary port is a Gaussian slot, not a new physical derivative. Each direct auxiliary color occurs at most once at that local center. Thus the ordinary rectangular Gaussian-series bound can be iterated across its distinct local Gaussian slots. Its two variance matrices are precisely proper flattenings of the local tensor. When a slot is differentiated, retain it as an input index in those flattenings. There is still a selected/output index on the other side, so the required flattening is proper. All local constants are uniform in fixed Q.

## 2. Three closure operations

### Matrix multiplication and numerical covariance-root powers

If A and B have the invariant, AB has it. Operator moments follow from fixed Hölder, without any independence assumption. For the derivative term (D A)B, the resulting series matrices are linear combinations of the old ones with coefficients B. Both frames obey

    sum_j Bnew_j Bnew_j^T <= ||B||op^2 sum_j Bold_j Bold_j^T,
    sum_j Bnew_j^T Bnew_j <= ||B||op^2 sum_j Bold_j^T Bold_j.

For A(D B), left multiplication gives the corresponding bound with ||A||op^2. Apply the triangle inequality to the square-root frame norm for the finite product-rule sum, then use Hölder. The third, transposed row frame obeys the same left/right operator-factor inequalities. Transposition exchanges the two row frames and preserves the HS-input flattening, so it is included without a missing cut.

Finite powers of M M^T and all covariance-root Taylor coefficients therefore preserve the invariant. They are treated as whole matrix coefficients. Their individual cyclic contractions are NOT declared to be executable permanently marked graphs. This distinction is essential.

### Multilinear side substitution

A fixed coefficient of a side-spine vector is a finite sum of terms A_i X_i. Here X_i is that spine's terminal innovation or a completion innovation, independent of A_i; its private tape is disjoint from the other sibling subtrees. The A_i may share all auxiliary Gaussians G with their parent and siblings.

Expand the finite multilinear substitution one term at a time. For a chosen side slot, condition on all shared auxiliaries and all other private variables before applying the Gaussian-series step in X_i. Its coefficient is independent of X_i. Inserting A_i on that Gaussian port transports both rectangular frames with a factor ||A_i||op, by exactly the two inequalities above. Repeat over the other distinct side-private Gaussian ports. The local source proper cuts supply the uncontracted tensor frames at each stage.

The resulting bound is a finite product of local frame/operator factors and side-coefficient operator factors. Remove their possible dependence on shared G by fixed Hölder, not by an independence assertion. After the side matrix factors have been bounded in this way, the remaining direct auxiliary Gaussian ports are estimated from the local tensor cuts. Shared occurrences at DIFFERENT centers only cause Hölder factors. At a single center the colors are distinct, as required for the local Gaussian-series step.

For a derivative, use the finite product and chain rules. A derivative either hits a direct local Gaussian port, one local coefficient, or one side A_i X_i. The first two cases use the retained local derivative cuts. In the third case apply the child's retained frames to (D A_i)X_i plus A_i D X_i. Equivalently, first finish the child's full vector Jacobian DY_i, then compose the parent's local side-derivative map with DY_i. Each Gram sum is bounded by ||DY_i||op^2 times the corresponding local Gram sum, and the HS-input flattening composes with ||DY_i||op. Local multilinearity removes the differentiated side slot, so the remaining local factor does not depend on that same child Y_i; there is no circular step. Transport other left/right factors exactly as above. Differentiating one shared G that occurs at many centers produces a finite SUM of such one-hit terms. It never requires falsely treating the repeated G values as independent.

The induction is on the finite side-spine tree, followed by the finite degree of the covariance-root expansion. Sibling private tapes are disjoint at every step. Products that repeat a side X_i, such as M M^T, are handled by the matrix-product closure and Hölder; they do not invoke a Gaussian step conditional on a coefficient that depends on that same X_i.

### Finishing a vector coefficient

For a vector term A X with X independent of A, Gaussian conditioning gives

    ||A X||_Lp <= C_p sqrt(D) ||A||_Lq(op).

Its complete Gaussian-bank derivative is (D A)X + A D X. The first term is a rectangular Gaussian matrix series conditional on A and its derivative frames; its operator Lp norm is bounded by the two frame square roots and a fixed power of log(2+N D). The second is bounded by ||A||op. Summing the finite terms gives

    ||B_N||_Lp <= Lambda_P sqrt(D),
    ||D_(G,reset,physical publics) B_N||_Lp(op) <= Lambda_P.

This is a whole-coefficient estimate. It does not use Poincare on the old Q roots, does not introduce a second sqrt(D), and does not equate an energy bound with a derivative bound.

## 3. Weighted amplitude degree and the finite remainder

Set t=alpha^(1/6) and zbar=tau z. The latter is uniformly bounded by the actual positive-clock endpoint estimate. Freeze alpha-dependent shields, clocks, filters, Q and zbar while taking the finite t-amplitude comparison. Then the actual coefficient-source amplitudes are

    C3 center: t zbar,
    nonroot C0 leaf: t^3,
    root C0 leaf: H t^(6d+3).

All exponents are positive integers for the present queue d=4,8,12,... . A returned graph with n centers and surplus d has t-degree 4n+6d. The fixed-target cutoff alpha^P is the ordinary t-degree cutoff 6P. Choosing covariance-root Taylor orders from this minimum positive degree is therefore finite.

The native calibration/value hybrid and bounded numerical covariance gap are the same ones used in the sealed local producer. First compare bounded fields at their actual Gaussian caller distributions, root first; use the covariance-root Taylor expansion at bounded/gated child values, not at an unbounded polynomial child. The preceding invariant proves the energy and complete new-bank Jacobian moment bounds for the resulting fixed-degree polynomial comparison, including shared auxiliaries. The Taylor and gate/calibration remainder has one sqrt(D) and is assigned its propagated absolute floor.

For the same-endpoint law step, keep Q fixed. Divide the remaining Gaussian block into P, the source-zero physical readout variables, and U, the auxiliaries and all zero-readout reset innovations. There is an untouched independent positive keep. Apply the finite same-endpoint Riesz argument to the actual t-path:

- A centered U coefficient is transferred by U-Riesz onto the endpoint derivative. There is no zero-degree base identity in U.
- The U-mean is a finite polynomial in P. Its centered part is transferred by P-Riesz. A base-Gaussian branch lowers visible polynomial degree. A nonlinear branch has strictly positive t-degree.
- The preceding complete coefficient Jacobian moments control these products through one Hilbert current norm and an operator derivative norm, using finite Hölder exponents.
- Stop at degree 6P. There can be only finitely many positive-degree steps, and finitely many degree-lowering base steps between them.
- Constant operators below the cutoff are canceled by the finite characteristic-coefficient identities at that SAME endpoint. Any deliberately uncorrected currents remain explicit; they are not called a Gaussian approximation.
- The independent keep transfers boundary currents to transport velocities with their actual rank-dependent inverse-width factors.

This is the weighted-degree, shared-auxiliary version of the admitted same-endpoint marked remainder interface. It supplies a finite remainder certificate for the conditional sector after its complete finite census is corrected. It is not a convergence assertion for the polynomial moment-generating series. A preexisting unmatched current outside this family is still present and must be handled separately. In particular, a retained nonzero polynomial log generator is a finite CURRENT target, not a claimed standalone probability law. A Gaussian/W2 endpoint conclusion is asserted only when an eligible joined program cancels every required lower coefficient against its Gaussian target and retains the stated independent gap. In that case the conditional remainder is at most C_P alpha^P sqrt(D), plus every propagated native/filter/clock/arithmetic error and every intentionally retained current transferred with its actual keep-width score factor. Without those joined cancellation hypotheses the output is the finite current ledger and its same-endpoint boundary estimate, not an invented positive comparator law.

## 4. Actual old-row guard and computable constants

The direct source path formulas in the main note must be multiplied by the ACTUAL original-bank injection rows. If a row norm J_old depends on alpha, that dependence is part of the guard. The assertion O(alpha) is valid only after verifying, for every recorded caller row and packet path,

    Lambda_P J_old alpha^[d+(m+1)(beta-gamma)] <= c alpha

and the corresponding root and leaf inequalities. An arbitrary unbounded old row is not covered by merely writing its norm as a factor.

At fixed P, enumerate all packet coefficients, tree choices, node counts, pair/filter orders, readout inverse powers and Gaussian moment constants. The resulting finite upper bounds K_P and public-log exponents m_P are computable from that finite list and the imported fixed-order constants. Check each actual source radius against its fixed pair/source guard, each first path against its retained-observer guard, and each propagated absolute error floor. This is the threshold definition; unspecified cancellations are not used.

When the admitted old rows have their stated public-log bound, a scalar sufficient threshold is found by increasing an integer x until all inequalities of the form

    K_P (L0+6x)^m_P exp(-x) <= required fixed radius gap

and the analogous stronger first/error inequalities hold, and setting alpha<=exp(-6x). Eventual success follows from the exponential dominating the fixed polynomial. If the actual old rows fail their guard, report that blocker instead of asserting smallness.

## 5. Limits that remain

This addendum is conditional on the complete old Q. It does not construct native sources for new whole-old-bank cumulants, physical-tilt derivative histories, or changed C_k orders. Nor does it prove a terminal genuine-gradient/curl identity for the whole joined program. The rank-five coefficient smoothing mismatch at tau=alpha^(1/3) remains alpha^(16/3); closing this conditional cycle sector does not erase that debt or establish an arbitrary-order algorithm for the original unsmoothed target.
