# Independent audit: positive sector-Gauss clock compression

2026-10-05. This audit did not modify the author's main files or the sealed inputs. It reviews the mathematical theorem and the rational scalar quadrature builder. It does not execute a native sampler or certify an original VALUE oracle's Sobolev accuracy.

## Verdict and scope

**PASS for the fixed-family positive coefficient quadrature theorem**, with its explicit original-private-heat floor, complete finite history census, external-caller derivative scope, and separate native-admission/replay guards. The theorem replaces an inverse-tolerance power in the two-edge coefficient node count by logarithms. At an inherited floor eta >= c alpha^q with fixed q, its new node-count exponent is at most 2q, independently of the requested coefficient accuracy grade. It does not prove that the total original-VALUE complexity has that exponent.

The exact-rational builder is independently checked below. Root proposals are not treated as proofs: successful calls are certified by exact signs, degree counting, rational interval bounds, positivity, mass, and first moment. The exact parameter-selection and adaptive precision-refinement portions also pass, as detailed in the implementation addendum below.

The imported manifest SHA256 was verified exactly:

09ef144321423e22831c3cc720c3b8e3c48cc2bd090c891bc675d10ca7841574

The reviewed main-file hashes are recorded in the accompanying check output and review pin. Hashes identify the precise reviewed versions; later edits require review of the changes.

## 1. Primitive constants and arbitrary Price order

The normalized privately heated derivative tensor has proper-cut bound 2^(k/2) sqrt(k!) sigma^(-k), and HS bound (k+1)^(k/2) sigma^(-k) sqrt(D). The enlarged constant C_k=2^(k+1) sqrt(k!) dominates both. For example, the ratio (k+1)^k/k! grows by a factor ((k+2)/(k+1))^(k+1) < 4 at each step, proving the slightly stronger inequality (k+1)^k <= 4^k k!.

At Price order l, precisely 2l primitive hits are added. With K the old total order and N the force count,

    product_v C_(k_v+d_v)
      <= 2^(N+K+2l) sqrt((K+2l)!).

The two binomial inequalities

    (K+2l)! <= 2^(K+2l) K! (2l)!,
    (2l)! <= 4^l (l!)^2

give sqrt((K+2l)!) <= 2^(K/2+2l) sqrt(K!) l!. Charging one eta^(-1) and at most P actual row-weighted hit choices for each Price edge yields exactly

    B_K l! (16P/eta)^l,
    B_K = W 2^(N+3K/2) sqrt(K!) eta^(-K/2).

There is no derivative of a selected covariance square root or a history-dependent heat choice in this argument. Price differentiation applies to the analytical Gaussian expectation before optional exact heat extraction. Positive original shields give uniform local derivative bounds. Adding an auxiliary covariance regularization for the proof and then taking its limit does not add target heat in execution.

P >= 1 and eta <= 1/2 are explicitly imposed. This matters when differentiating the Jacobian r: the bound r+eta/(16P) <= 2 would otherwise require a qualification for arbitrarily small row bounds.

An external affine caller derivative adds its actual force-hit and caller-row sum to W and increases K by one. The final W is explicitly a full-family sum over replica trees, new force-hit assignments, and original histories/cross-histories. B_* >= 1+sum_T B_T therefore pays the total family error once. A per-tree or per-history bound could not simply be reused at the full-family tolerance without this census or a tolerance allocation.

An initially omitted P^l in an intermediate display was corrected; the final high-derivative bound always had the correct factor.

## 2. All proper cuts survive the additional cycles

Higher Price edges can make the force graph cyclic. The earlier active-source note's tree theorem and HS-only cyclic extension would not, alone, justify all proper cuts. The separately sealed MARKED-SPANNING-TREE-EXTENSION.md, sections 2–3, supplies the needed stronger theorem and its proof is valid.

For any proper external partition, orient the marked spanning tree using a balanced terminal flow with no zero tree-edge flows. Every local vertex then has at least one incoming and one outgoing index. Choose a topological order of this directed tree and orient every additional edge forward in that order. The full graph is still acyclic as a directed circuit, and every local cut remains proper. Composition and tensor products of these local operators prove the product all-cut estimate. The argument permits parallel edges and excludes local self-loops. Price edges join distinct replicas, so they satisfy that restriction.

For HS, start at a marked spanning-tree leaf with one local HS factor and attach vertices parent first, contracting all their edges to the built set. Every remaining local cut is proper. This spends precisely one sqrt(D).

The independent tests explicitly construct these orientations for all 510 proper partitions of a nine-mark triple-cubic tree and all 4,094 proper partitions of a twelve-mark mixed tree, while allowing every possible additional inter-replica edge. These finite checks supplement the general proof rather than replace it.

## 3. Order sectors, Jacobian, and the tail

For each of the two orderings, t_large=r and t_small=rs maps the unit square to the corresponding triangular sector with Jacobian r. A three-replica covariance has one off-diagonal entry r and two entries rs. It is separately affine in each transformed coordinate. Consequently the pure-coordinate Price bounds apply at every order even though the original min-path integrand is not globally smooth across its ordering wall.

For G=rF, pure r differentiation gives r times the l-th derivative of F plus l times its (l-1)-st derivative. Pure s differentiation has no Jacobian derivative. Both are bounded by 2 B_T l! (16P/eta)^l. No mixed-derivative theorem is needed by the tensor-product error argument.

The omitted strips r>1-zeta or s>1-zeta have transformed area at most 2zeta in each sector. Since ||rF|| is at most B_T, summing sectors and trees gives at most 4 B_* zeta. The s strip includes a neighborhood of the ordering diagonal, so it must not be described as solely an original endpoint tail. The final text correctly calls it a sector-strip quadrature tail and preserves the full original finite coefficient as the approximated target.

## 4. Gauss remainder and node count

On an interval of length ell, the m-point Gauss-Legendre remainder coefficient is

    ell^(2m+1) (m!)^4 / [(2m+1)((2m)!)^3].

Combining this with the 2m-th derivative bound gives the displayed estimate

    2 B_T ell [ell(16P/eta)]^(2m)
      (m!)^4 / [(2m+1)((2m)!)^2].

Since ell <= eta/(32P), this is at most 2 B_T ell 4^(-m). Scalar duality proves the same norm estimate for the finite-dimensional tensor-valued integral. Summing intervals and using the positive tensor telescope gives 8 B_T 4^(-m) for both sectors. Summing B_T gives the stated full-family bound.

The choices 2^J >= 32 B_*/epsilon and 4^m >= 32 B_*/epsilon give tail <= epsilon/8 and exact-Gauss error <= epsilon/4, leaving more than half the tolerance for certified arithmetic and any separately allocated coefficient floors.

For dyadic gap panels, the sum of subdivision counts is at most ceil(1/h)+J. Two sectors and m nodes on each interval in two dimensions therefore yield

    Q_tree <= 2 m^2 [ceil(32P/eta)+J]^2.

This is a genuine executed node count. It is not the number of tensor entries, and no tensor expectation is evaluated to construct it.

## 5. Positivity and preservation of heat majorants

In transformed gap coordinates x=1-r and y=1-s, the two original edge gaps are x and x+y-xy. A subinterval's length is at most every gap in its dyadic panel. Every positive interval weight is at most its interval length. Therefore

    w = r lambda_r lambda_s <= xy <= x(x+y-xy).

This holds also for the symmetric positive rational replacement, because its weights have exact interval mass and nodes stay in the same interval.

Within a dyadic product cell, both original edge gaps vary by at most two; so does each replica innovation gap, and adding eta preserves this comparison. A positive monomial majorant in these regularized gaps thus has a fixed-rank sup/inf comparison. The quadrature integrates r exactly on every cell. Its sum of a majorant is consequently at most that fixed comparison factor times the original sector integral, independently of how many quadrature nodes the cell contains.

This is sufficient to transfer the imported integrated S0/S1 estimates and their literal node maxima. It does not assume the transformed bridge rule factors in the original edge coordinates. Literal squares then follow from maximum times sum. In particular, the mixed caller square allowance is eta^(-2chi), while the nonsquared caller allowance is eta^(-chi). No independent-clock variance identity is substituted for the repeated-history weights.

## 6. Exact target and native boundary

The construction retains the complete original coefficient bank, its actual source-zero physical boundary, all repeated clock labels, and old six-packet residual roots. At each node the replica covariance is exact, original derivative contractions keep their actual injection rows, and one simultaneous m_star extraction is frozen for all main/caller estimates. In the mixed family only the new independent innovation is used for the new extraction. The original private heat is preserved as a summand; eta is a lower bound, not replacement heat.

The proved approximation concerns the bank-mean coefficient and fixed external affine callers. It is not a pointwise-in-old-bank identity or an environment-stable native LAW. Internal tapes are not exposed as new unattenuated observers. Actual captured-caller first paths, source/filter Sobolev floors, readout normalizations, same-endpoint mixed currents, and independent keeps remain separate obligations.

The shrinking sector strip can create small coefficient-bank innovation gaps. The displayed disintegration only uses their nonnegative square roots with bounded center rows; it does not invert those gaps. Nevertheless any chosen native implementation's inverse gap, shared-readout dilution, source-radius restriction, or ancestor/replay cost remains charged. The final text states this explicitly. Thus no global c(P)=o(P) conclusion follows from the coefficient node count.

## 7. Rational root and weight implementation

The author's Legendre recurrence is exact over Fraction. For each nonzero root, a rational sign-changing interval supplies at least one root. Disjointness and degree counting supply exactly one root in each interval and exhaust all roots; the odd-degree zero root is checked exactly. The rational Horner interval for P'_m, squared denominator interval, and positive weight formula 2/[(1-x^2)P'_m(x)^2] enclose the exact Gauss weight whenever the successful call asserts a positive denominator lower bound.

Paired root intervals and equal paired positive weights preserve reflection symmetry exactly. Normalizing the rational weights to total mass two preserves positivity and gives first moment zero. Affine interval mapping therefore preserves its exact mass and first moment; the product rule integrates r exactly. The returned node and total-weight error bounds compare the rational rule with the exact Gauss rule, even when normalization moves a proposed weight outside its original weight enclosure: the code measures the maximum distance to both enclosing endpoints after normalization.

Independent successful calls at degrees 1, 2, 3, 5, 7, 16, and 24 verify exact mass, first moment, positive weights, and reported error envelopes. Finite precision does not imply native Sobolev or curl accuracy, and the author correctly excludes that inference. Every changed width, center, source version, requested original-VALUE floor, or affected ancestor requires its actual replay.

## Reproducibility and limits of the diagnostics

- check_positive_gauss_audit.py and positive-gauss-audit-checks.json: independent integer factorial checks, all-cut orientation checks, 160 positive-rule geometric cases, Gauss budget algebra, and an analytically integrable sector fixture.
- author-builder-independent-checks.json: independent calls to the exact-rational builder at additional degrees.
- REVIEW-PIN.json: precise reviewed source hashes and final implementation review status.

The numerical fixture exhibits decreasing error down to floating-point roundoff on the retained sectors, while recording the genuinely nonzero removed-strip integral separately. It does not equate the truncated quadrature with the full integral. These tests validate scalar formulas and finite diagnostics; the dimension-sharp tensor and private-bank assertions rest on the proofs and explicitly imported contracts above.

## Final implementation addendum

The final builder selects positive root proposals by sorted index, avoiding an absolute floating threshold. exact_B uses an integer upper bound on sqrt(K!), rounds only upward in powers of two, and rounds the eta exponent upward when K is odd. Its result is a rational upper bound for the displayed B_K. The extra conservative eta half-power occurs only inside the logarithmic degree/panel budget; it does not change the leading eta^(-2) node-count power.

exact_parameters constructs B_*=1+B_K+B_(K+1), compares rational numbers with integer powers of two, and chooses m=ceil(J/2). It therefore cannot under-round a logarithm at an integer wall. The default caller_weight=N W is explicitly qualified as using caller rows at most one; other callers require their actual scalar hit-weight sum. Unsigned inputs are checked for nonnegativity.

certified_legendre_to_tolerance raises precision until both certified error envelopes meet the requested rational floors. prepare_certified_rule combines this with the exact parameter chooser and panel construction. The Cartesian node iterator is lazy; the one-dimensional panels and nodes are still genuinely constructed and billed. The floating theorem_census remains explicitly only a displayed estimate, rather than the production certificate.

For a normalized reference-rule node perturbation d and reference total weight perturbation e_w, the full-family difference is at most

    64 B_* (P/eta) d + 2 B_* e_w.

To see the absence of a count loss, interval mapping multiplies both perturbations by half the interval length. Summing positive interval masses gives at most one, and a two-axis telescope gives the corresponding tensor bound. The code's floors d=epsilon eta/(512 B_* P), e_w=epsilon/(512 B_*) cost at most epsilon/8+epsilon/256, below epsilon/4. The author's final text now states these constants explicitly.

Independent exact checks cover 108 combinations of N, K, eta and epsilon, including non-dyadic eta; they verify the squared rational B upper bound and both integer-exponent budgets. Integer-wall checks cover powers 2^j and nearby rational values for j=0,...,64. Adaptive rules at degrees 2, 7 and 16 meet 10^(-40) node and weight floors. A full preparation call returns a degree-nine rule and 75 intervals, with both floors certified. All tests pass. This completes the scalar implementation review; native execution and native Sobolev certificates remain outside its scope.
