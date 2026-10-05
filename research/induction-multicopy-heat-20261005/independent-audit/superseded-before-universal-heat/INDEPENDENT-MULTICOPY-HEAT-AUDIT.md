# Independent multicopy heat audit

2026-10-05. Independent review of the three-cubic, mixed 3/3/6, and conservative finite-cubature constructions. This review did not modify the sealed inputs or the three primary proof files.

## Verdict and limits

**PASS for the stated fixed-family coefficient identities, exact Gaussian heat redistribution, weighted main/one-hit bounds, literal node-square bounds, and conservative finite positive cubature certificate**, relative to the imported marked-spanning-tree, normalized source, native conditional-return, and same-endpoint contracts.

The numerical constants attached to actual native sources, readout products, node counts, finite errors, and complete replay remain admission guards. This is not a claim that those guards hold at every requested tolerance or inherited cutoff. It does not close the full induction or its feedback queue. The absence of a public-log higher-cumulant quadrature bound is not treated as a defect: that bound is expressly not claimed.

Two concrete errors in the initially reviewed text were reported to the author and repaired before this verdict:

1. A repeated diagonal leaf can have aggregate inverse-width power **five**, not four, when both old covariance hits and three new hits land there. Its squared clock mass leaves at most eta^(-1/2). This agrees with the stated conservative diagonal result, but required a different proof sentence.
2. The affine **bank** row is bounded by Gamma_vv; the concatenated retained-Y-plus-bank row includes d_v^2 as well. For the literal cubic, d_v^2+Gamma_vv=(1+r_v^2)/2<=1. The revised text distinguishes these statements.

The reviewed SHA256 values are recorded in `exact_exponent_ledger_checks.json`; review snapshots are stored beside this audit.

## 1. Replica covariance and preservation of the actual old bank

For each threshold u, remove tree edges below u. The resulting connectivity matrix is a sum of all-ones block matrices and is PSD. Its integral is K. The full connected block is present for u<=c=min t_e, while vertex i is a singleton exactly on an interval of length delta_i=1-max_{e incident to i} t_e. Removing those disjoint contributions leaves only non-singleton block matrices, proving

    H=K-c 11^T-diag(delta_i)>=0.

Consequently the proposed realization with the actual retained B, an H-root, and independent delta_i innovations has covariance K tensor I. It does not replace the old service's B with a new bank. The rows used by the already-differentiated coefficient remain the original A rows: multiplying those contractions by sqrt(c) would be incorrect, and the construction does not do so.

Within replica i, extracting a diagonal D_i from its independent innovation changes its residual query covariance to delta_i(Gamma_i-D_i) and its private shields by delta_i D_i. Their sum is unchanged. Cross-copy covariances are unchanged because the moved part was replica-independent. The coefficient mean and the deterministic retained-Y row are unchanged. The new bank-row squared norm is Gamma_i,vv-delta_i(D_i)_vv.

The identity holds after the prescribed Gaussian averages. It is not pointwise cancellation after exposing arbitrary old/private tapes, and the primary notes correctly retain the resulting old/new mixed currents.

## 2. Hit-adapted cubic floor

For target force j=1, retain only the common Y0 randomness and the independent W_v terms. Its variance is z_1=1-q^2, with common coefficients r_1,r_2 s,r_3 s. For target j=2 or 3, regress Y0 onto X. The regression coefficient is

    s(1-q^2)/(1-s^2 q^2),

which is between zero and one because

    1-s^2 q^2-s(1-q^2)=(1-s)(1+s q^2)>=0.

The retained common variance Var(X) is at least 1-s^2. If r_j>=1/2, writing the common linear form as r_j x_j plus the other two coefficients gives

    |x|^2<=13(z^2+sum_{v!=j}x_v^2).

The independent W_v for v!=j and the common variance therefore give Gamma>=m_j I/13. If r_j<1/2, sigma_j^2>=3/8 and m_j<=1/2, so the W-bank alone gives Gamma>=3m_j I/4. Thus the advertised weaker floor Gamma>=m_j I/16 is valid for all three possible hit locations.

D=m_j I/32 is safely below Gamma. The floor's omission of the badly hit force is essential: using only a floor containing that same private variance would not prove the displayed weighted bounds.

## 3. Three-cubic old-clock sums, bridge integrals, and literal squares

For total cubic inverse-width degree 1+d<=4, there can be at most one force with p_j>2. Set a=p_j/2-1. Summing that force's one positive clock mass gives

    C m_j^(-a)(delta+eta)^(-a).

This follows from the ordinary dyadic split and sigma_j^2+delta m_j/32>=const*m_j*(eta+delta). The expansion of m_j^(-a) charges exactly one remaining private or structural clock. Since the total degree is at most four, every resulting endpoint exponent is at most its available single mass. Equality costs only an endpoint-panel count. This proves J_d with a_d=max(0,(d-1)/2). The same calculation at a literal node proves its stated nodewise version.

Repeated labels remain repeated in the target. A repeated-label sum can be bounded above by the nonnegative Cartesian product of single-copy sums because the diagonal is a subset of that product. This is a triangle majorant, not an equality and not a variance/independence claim. No independent-clock bound is imported as an exact squared-clock identity.

For a three-vertex replica tree the independent variances are s_1,s_2 at the leaves and min(s_1,s_2) centrally. The main envelope has exponent 1/2 on the minimum, with integral 8/3. The extra central hit gives

    integral_[0,1]^2 (min(s_1,s_2)+eta)^(-1) <=2 log((1+eta)/eta).

For an extra leaf hit, splitting at s_leaf=s_other gives at most log((1+eta)/eta)+2. Thus the claimed main and complete one-hit allowances follow, including their logarithm.

The scalar node-square deductions are valid. If S_h=sum_v 1/t_v, the two separate nodewise bounds sup b_h<=C and sup b_h S_h<=C imply

    sum b_h^2 <=(sup b_h) sum b_h,
    sum b_h^2 S_h^2 <=(sup b_h S_h) sum b_h S_h.

There is no replacement of literal squared weights by an independent-clock theorem. The bridge weights pay the minimum-gap singularities before these maxima are taken.

The graph census is consistent: nine forces, nine physical marks, eight edges, local width degree seven, and maximum main C3. There are 3*3^4=243 tree/hit histories before old histories and physical permutations.

## 4. Full mixed 3/3/6 residual extraction

The six-packet's recorded residual block in one constituent is s(Gamma-m_old I/32). The recorded floor gives

    Gamma-m_old I/32 >= Gamma/2 >= m_j I/32.

After the next replica decomposition, the *new independent* delta_6 innovation of this recorded block can therefore contribute delta_6 s m_j I/64 to fresh private heat while leaving a PSD remainder. The two constituent blocks are separately available and independent. Old private heat is retained as a summand.

This is a valid redistribution inside the new innovation of L6(B_6). It does not alter the recorded conditional function L6(B) at its original boundary. Modifying that old function and then claiming its old mixed currents were unchanged would fail, but is not the construction being audited.

The complete original bank can contain both old common roots and old residual roots. L3 simply has zero rows on the latter. All new derivative contractions use the actual original rows on this joint bank. No derivative of frozen clock roots is needed.

## 5. Full versus diagonal six-packet powers

### Full cross-history family

A constituent has its original center degree one, one old covariance hit, and d_a additional hits, hence total degree 2+d_a. Through d_a=2 the cubic lemma gives gap power d_a/2 without an inverse-cutoff loss. For d_a=3, total degree is five. If p_j>2, extracting its clock costs m_j^(-(p_j/2-1)). Every remaining endpoint exponent is at most 3/2. Since the expansion of m_j^(-a) charges only one clock per term, the extra loss is at most eta^(-1/2), not a product of multiple such losses.

This step needs the mixed note's eta to bound **both private and structural original endpoint variances**. A private-only minimum would not justify that cutoff charge. The note states the stronger requirement explicitly.

Multiplying the two literal constituent sums produces the stated effective gap delta*s. Integrating the old bridge gives the displayed square-root, logarithmic, and 3/2-power bounds. The worst d=3 case combines the local eta^(-1/2) with the old-bridge eta^(-1/2), yielding the safe bound

    J6(3,delta)<=C L^C eta^(-1)(delta+eta)^(-1).

These are upper bounds only. No assertion of sharpness or a lower-bound obstruction follows.

### Repeated diagonal family

Let P_v be the aggregate inverse-width powers on the repeated three clocks. Then P_center>=2 and sum P_v=4+d.

If P_center>4, summing its squared mass costs

    a=P_center/2-2<=d/2.

Its canonical old floor omits the center. Each other squared-clock exponent after charging that floor is at most

    P_leaf/2+a <=(4+d)/2-2=d/2<=3/2,

and the structural exponent a is likewise at most 3/2. All are below the available squared mass of two. Only the old bridge remains; its integrated power a leaves at most eta^(-max(0,a-1)), hence eta^(-1/2) for d=3.

If P_center<=4, every clock except possibly one leaf is paid directly by its squared mass. A leaf can reach five only at d=3. That leaves eta^(-1/2); then the center is exactly its baseline two and there is no simultaneous center deficit. This proves all nine old-hit cases, including the corrected leaf-five case.

The independent exact enumeration checks 9,54,324,1944 ordered patterns at d=0,1,2,3, respectively, a total of 2331. Maximum center powers are 4,5,6,7; maximum leaf powers are 2,3,4,5. The corresponding maximum diagonal inverse-cutoff powers are 0,0,0,1/2. The script uses exact integer and rational arithmetic.

## 6. New mixed bridge budget, source amplitudes, and return scales

With L6 central, its main gap power is one on min(s_1,s_2), so the integrated main is logarithmic. Its caller carries the safe eta^(-1) allowance. With an L3 central, the product is the cubic minimum-gap square root times the six-leaf gap square root, again logarithmic. Any one new caller adds at most eta^(-1/2) relative to that envelope, within the declared eta^(-1) full allowance. The diagonal bound similarly leaves only eta^(-1/2).

The nodewise main bounds are consistent: u_old/(delta*s+eta)<=C/(delta+eta), and the two new bridge weights pay the resulting central or leaf singularities. For a six-packet caller, the rough but sufficient estimate

    u_old eta^(-1/2)(delta*s+eta)^(-3/2)
        <=C eta^(-1)(delta+eta)^(-1)

follows by u_old<=C s and (delta*s+eta)>=s(delta+eta). Therefore the displayed sup, sum, and squared-sum caller bounds hold. The squared caller allowance is eta^(-2 chi), as printed; it is not confused with the single derivative factor in the bank-mean self current.

The mixed census is twelve forces, twelve marks, eleven edges, width degree ten, with main C4 and analytical caller C5. The exact new tree/hit total is 324+162+162=648. The graph still carries its original marked spanning tree after every bank edge is inserted.

The root amplitude alpha^(12-11 beta) times b_h and eleven alpha^beta nonroots multiplies to alpha^12. The intrinsic root-square exponent is 24-22 beta. The Riesz self-bank scale uses one coefficient energy and one bank derivative, hence alpha^24 eta^(-chi), not the square of the derivative loss. The old-mixed scale alpha^13 requires the separately stated actual complete old-bank first O(alpha) at the same endpoint. These uses match the imported contracts.

Substituting eta=alpha^(2 gamma) gives precisely the printed sufficient inequalities. For beta=1/22 the root exponent is 23/2, and the source-first condition is gamma<21/(4 chi). The arithmetic in both examples is correct. Readout/node constants that scale with alpha must still be debited, as the text requires.

## 7. Conservative positive cubature certificate

The weighted Prüfer census H_n=(product N_i)(sum N_i)^(n-2) is exact; a Prüfer occurrence increments that replica's degree by one. Independent integer checks include the 243 and 648 cases.

For cubature, use the analytical tensor integrand before optional extraction. In an open ordering sector, each derivative of a min-path covariance entry is zero or one. Price differentiation adds one bank edge between distinct replicas and one force hit at each end. These edges may duplicate existing edges, but never create a self-loop. The original marked spanning tree remains a certificate. Uniform privately heated derivative bounds therefore give the proposed known scalar B1,e and one-HS bound B1,e sqrt(D).

Positive regularization handles singular bank covariances. Positive original private heat makes the force derivative envelopes finite. Continuity across sector walls and integration along coordinate segments give the global Lipschitz bound; a generic perturbation of a segment handles the case where it lies in a wall. No derivative of a chosen covariance square root, and no numerical tensor evaluation, is used.

On a tensor midpoint cell, coordinate displacement from its midpoint is at most h/2. Summing the Lipschitz bounds proves the displayed sufficient mesh h<=2 epsilon/max(1,sum B1,e). The stated one-edge node count safely bounds the sum of panel ceiling counts, and its (n-1)st power bounds the tensor rule.

Every node is interior and its one-dimensional weight obeys u<=2s. On the final interval [0,eta], nodes need not be comparable to eta when that panel is subdivided. Its total mass is eta, and the existing denominators are regularized by eta, so the summed endpoint bounds still hold. This minor endpoint nuance does not invalidate the rule.

The numerical count can be polynomial in inverse heat and tolerance. That cost, and its effect on native readout/source admissions, remains explicitly separate. The cubature theorem does not by itself certify a public-log or eventual-sublinear return.

## Evidence and reproducibility

- `check_exact_exponent_ledger.py`: independent exact enumeration of all local cubic and mixed hit powers and several weighted Prüfer censuses.
- `exact_exponent_ledger_checks.json`: PASS output, counts, worst powers, and exact reviewed hashes.
- `REVIEWED-THREE-CUBIC-EXACT-HEAT.md`, `REVIEWED-MIXED-THREE-THREE-SIX.md`, and `REVIEWED-FINITE-CUBATURE-AND-N-COPY-RULE.md`: reviewed snapshots.

The new analytic claims were checked directly against the three sealed source documents listed in the assignment. The native finite source programs were not numerically executed. This audit's finite exponent diagnostics supplement the proofs above and are not substitutes for the imported source or endpoint contracts.
