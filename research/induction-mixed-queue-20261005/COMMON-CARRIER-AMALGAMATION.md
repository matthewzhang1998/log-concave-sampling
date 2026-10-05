# Common-carrier raw amalgamation without per-node variance dilution

2026-10-05. The Gaussian construction, first/curl algebra, and IDEAL finite-polynomial current theorem are proved below. The independent audit found that the supplied ACTUAL native prior does not export the necessary carrier-retaining hybrid floor. Section 4 isolates that missing inequality; the native triple-cubic application remains conditional on it. The ideal construction avoids inverse **physical readout-share** growth with the number of quadrature nodes. It does not remove those nodes from the query bill, or remove existing original-clock target errors.

## 1. Actual Gaussian coupling, not a marginal subtraction argument

For h=1,...,J take a fixed-rank native raw program

    X_h(r_h; B,P_h,U_h) = L_h P_h + u_h(r_h;B,P_h,U_h),
    u_h(0;...) = 0 pointwise,
    L_h L_h^T = v I_D,   v>0.                              (1)

P_h is that program's complete source-zero physical-public block, standard Gaussian. U_h comprises its remaining Gaussian private/cut/filter/reset tapes. The readout L_h is a known, B-independent scalar-block row. All original coefficient callers B are independent of these new dynamic roots. The source-zero identity in (1) is an actual frozen-version map identity, not merely equality in law or zero marginal mean. The h-th program's internal positive covariance gaps are already checked.

Choose one H~N(0,vI_D), independent of B and old tapes. Let V_h be independent standard copies of the P_h space, and let

    Pi_h = I - v^(-1) L_h^T L_h,
    P_h = v^(-1) L_h^T H + Pi_h V_h.                       (2)

Since `v^-1 L_h^T L_h` is an orthogonal projection,

    Cov(P_h)=v^(-1)L_h^T L_h+Pi_h=I,
    L_h P_h=H exactly.                                     (3)

Every packet's full marginal `(P_h,U_h)` is unchanged. Across two packets,

    Cov(P_h,P_g)=v^(-1)L_h^T L_g.

This is a known positive Gaussian geometry because (2) constructs it directly. Pi_h may be used as its own positive semidefinite square root. It is a known scalar-block operation; no unknown force covariance is inverted.

Execute the SINGLE owned raw group

    X_group = H + sum_h u_h
            = H + sum_h (X_h-H).                           (4)

Every term of (4) is computed from its original raw finite program and its known actual source-zero carrier. We are not subtracting a carrier from an abstract marginal LAW. No individual X_h, P_h or H is subsequently exposed as a separate external observer. The group, not each summand, owns the shared H and all null/private tapes.

Reserve `sqrt(kappa)G_keep` independently once, with v+kappa and all old fixed rows fitting the actual total covariance budget. The internal v-sized readout of each X_h is the same H; it is not J independent variance shares. Source covariance roots within each X_h retain their original positive gaps. The final distribution is positive because (4) plus the keep is an explicit Gaussian pushforward, even though the known linear arithmetic in (4) contains subtraction.

## 2. Complete firsts, including null-space directions

Write H=sqrt(v)G0. The map from `(G0,V_h)` to P_h is

    Q_h=[L_h^T/sqrt(v), Pi_h],  Q_h Q_h^T=I.

It is a coisometry and has operator norm one. Suppose the imported normalized native port supplies, for the complete residual and required finite moments,

    ||u_h||_Lp <= Gamma_h |r_h| sqrt(D),
    ||D_(P_h,U_h)u_h||_(Lp;op) <= Gamma_h |r_h|.             (5)

Then after (2) these same bounds hold on the complete shared/null/private dynamic bank. In particular

    D_G0 u_h = D_P u_h L_h^T/sqrt(v),
    D_Vh u_h = D_P u_h Pi_h,
    D_Uh u_h = its original complete private first.         (6)

All null-space derivatives are retained. Minkowski, without pretending the shared H is independent, gives

    S=sum_h Gamma_h |r_h|,
    ||sum_h u_h||_Lp <= S sqrt(D),
    ||D_(G0,V_1,U_1,...,V_J,U_J) sum_h u_h||_(Lp;op) <= S.  (7)

The sharper square sum for separate null blocks is optional; the displayed l1 bound is always safe. Captured B/Y firsts likewise add their actual path bounds because L_h,Pi_h are B/Y-independent. If the row geometry depends on B or Y, its derivatives must be added and this simplified result no longer applies.

The physical-H first satisfies `||D_H sum u_h||<=S/sqrt(v)`. Source private curl is preserved at each literal native source; this alone is NOT a terminal curl theorem. At an actual retained terminal port

    g(s H + z_0 - R(H,other tapes)),

with scalar s and z_0 independent of H, the literal derivative is

    D_H g = Dg(sH+z_0-R)(sI-D_H R).

Thus its first is at most `A(|s|+||D_H R||)`, and its antisymmetric part has norm at most `A||D_H R||` (or twice this under the unhalved curl convention). For `R=R_old+sum_h u_h`, insert the FULL old first plus `S/sqrt(v)` and every numerical Sobolev floor. All other bank/caller directions use the same chain rule with the full sums (6)-(7). This preserves an eligible existing O(A alpha) terminal-curl guard when the old residual first is O(alpha) and the new summed first is smaller. It does not make the grouped output itself a genuine gradient or infer a better curl from its LAW error.

## 3. Root-linear current and all cross terms

The needed analytic input is the admitted finite bounded/gapped native CURRENT interface, not merely (5). For each packet it supplies:

- a heat-uniform finite polynomial comparison with complete Gaussian derivative frames;
- its exact root-linear constant current T_h, including the complete side-source response and original coefficient labels;
- root-zero identity (1), so every dynamic nonlinear endpoint defect of that packet has r_h;
- separately propagated native/filter/source/mean/gate/arithmetic/Sobolev floors.

These are the ideal polynomial/current inputs used by the sealed frozen-coefficient port. Its separate actual native-prior floor is not automatically compatible with sharing the carrier; Section 4 states the additional requirement. Keep B fixed in this section. Original coefficient tensors and captured original query centers are B/Y-dependent and H-independent; dynamic polynomial coefficients of u_h may depend on H and the null/private roots, as the original packet already does.

At all roots zero, the endpoint in (4) is H. Because each `(H,P_h,U_h)` marginal is exactly the original `(L_hP_h,P_h,U_h)` marginal, the ENTIRE r_h-linear base-Gaussian current is unchanged: it is T_h. This follows from the original current identity, not from a zero-mean assertion. All root-linear Gaussian coefficient-derivative branches are included in that original identity.

Use one common reference field

    E_sum(B,H) = sum_h T_h(B):H_(R_h-1)^(v)(H),

with each actual rank and inverse-covariance convention. It is linear in the root amplitudes and uses the same H. The positive interpolation is

    X_t = H + sum_h u_h(t r_h) +(1-t)E_sum(B,H)
                 + W_old(B,U_old,Y) + sqrt(kappa)G_keep.    (8)

In the finite same-endpoint current extraction, a path derivative contributes some r_h. A Gaussian integration-by-parts/Riesz branch either:

1. differentiates the source-zero H row and is one of the root-linear base terms already matched by T_h; or
2. differentiates a nonlinear endpoint residual u_g or the nonlinear reference drift, acquiring a further r_g; or
3. differentiates a finite current coefficient, which is handled by its original Gaussian polynomial/frame identity; its root-linear base portion is again part of T_h, and any endpoint nonlinear branch has case 2's second root.

Sharing H permits h!=g in case 2. It introduces no new root-free endpoint defect, because every u_g(0)=0 pointwise. Null/private directions have no source-zero output row, so they only enter case 2. Old W reads none of these new dynamic roots. If a captured original coefficient tensor itself reads H, or an outside observer retains H, this argument is invalid unless its additional root-linear branches are separately included.

The full finite native-current argument therefore bounds every unmatched intrinsic current by products of at least two root amplitudes. The complete shared-root derivative frames are controlled by (7), and the fixed-current-order product/Hölder sums are bounded by powers of S. At S below the actual numerical native radius, the independent keep gives

    conditional intrinsic error <= C_* S^2 sqrt(D).        (9)

C_* depends on the fixed maximum graph/current orders, normalized local bounds, v,kappa, actual moment list, and declared public dimension logarithms. It need not contain a factor J from dividing physical variance into J shares. To see the count dependence explicitly, at each fixed multilinear order q sum the absolute amplitudes before Hölder: the full ordered sum is bounded by S^q. The fixed q! and moment constants depend on order, not the number of labels h. Any implementation-specific extra J dependence in a supplied frame/compiler constant must nevertheless be exposed; (9) is not permission to hide one.

All finite linear source/filter errors are paid separately at their summed propagated floors. They do not become root-quadratic by (9).

## 4. Exact native-hybrid gap and two sufficient repairs

The audited native prior retains the incoming private x and external labels but discards source-private tapes/internal values. Its actual-prior replacement is a conditional W2 coupling of completed spines. Its subsequent strong polynomial comparison is for the IDEAL forward forest. Sharing H in (2) makes other packets read linear combinations of those internal source-zero carrier tapes. Therefore one cannot apply the old marginal prior comparison while holding those tapes fixed.

This issue does not affect the ordinary independent-group rank-eight theorem: first replace the independently owned actual groups by their ideal counterparts, using their original conditional native floors. Only afterward introduce the proof-only common-H reference in that ideal comparison. Here (2)-(4) instead changes the EXECUTED native dependence, which is why a stronger hybrid port is needed.

A sufficient carrier-retaining floor for each packet is

    [ E_H W2^2(Law(u_h^actual | H,B,E_ext),
                Law(u_h^ideal | H,B,E_ext)) ]^(1/2)
       <= epsilon_h(B,E_ext),                              (13)

with H coupled identically and all unchanged external labels retained. Conditional on the common H,B,E_ext, couple the packets on their independent null/private tapes. Minkowski then bounds the group hybrid by the sum of the epsilon_h. This is exactly the additional joint property required, not a consequence of a marginal LAW floor.

A weaker sufficient port is an environment-stable replacement, uniformly at every finite hybrid stage:

    W2(group_before_h, group_after_h)
       <= epsilon_h + C Gamma_h |r_h| S_-h sqrt(D),         (14)

where the unchanged remainder shares only H and has its actual residual first/energy budget S_-h. Summing (14) gives the original floors plus O(S^2 sqrt(D)), which is enough for (9).

There is an elementary sufficient condition for (14). Suppose the rest of the group R_-h(H,Z) is GLOBALLY H-Lipschitz with constant L_-h for every independent Z. Couple the marginal actual/ideal packet outputs at error epsilon_h, and disintegrate their carriers conditional on those outputs. Since `X_h=H+u_h` on either side,

    ||H_actual-H_ideal||_2
       <= epsilon_h+||u_h^actual||_2+||u_h^ideal||_2.

Use the same independent Z in the rest of the group. Then

    W2(group_before_h,group_after_h)
       <= (1+L_-h)epsilon_h
          +L_-h (||u_h^actual||_2+||u_h^ideal||_2).         (15)

With L_-h=O(S_-h) and the root-linear residual energy bound, this is (14). Operator-Lp frame bounds do NOT by themselves imply the global Lipschitz hypothesis. The audit additionally located global raw residual Lipschitz bounds for ACTUAL fixed-order LOW30 compilers (proposition `b27:nc:raw-paths`, lines 6488-6506), with external labels requiring their own uniform source-radius bounds. This is useful, but does not supply a global bound for IDEAL finite polynomial residuals already inserted at intermediate hybrid stages. Those mixed environments still need a proof before (15) closes the full native hybrid.

For an elementary nonimplication example, take independent scalar standard Gaussians H,U and `X_r=cos(r)H+sin(r)U`. It has exactly the same marginal standard Gaussian law as H at every r and satisfies X_0=H. Nevertheless the H-retaining conditional coupling error to H has square `2(1-cos r)`, hence size asymptotic to |r|. An arbitrarily small marginal floor does not imply (13). This example does not rule out (14) or a redesigned carrier-preserving native implementation.

Thus the ideal common-carrier current theorem is established, but the actual native claim below needs (13), (14), or the verified global bound permitting (15). The missing port is specific and quantitative.

## 5. Conditional fixed triple-cubic application

The independent heat construction in

`/workspace/shared/induction-multicopy-heat-20261005/THREE-CUBIC-EXACT-HEAT.md`

supplies the exact whole-bank third cumulant of THREE complete finite cubic sums. It expands every ordered old-history triple, all three labelled replica trees, and all 81 ordered force-hit assignments per tree. Each main graph has `(a,N,K,R)=(9,9,7,9)` before admissible physical contractions. Its new replica heat is an exact redistribution of the original Gaussian law, not additional smoothing.

Let b_h be its unsigned scalar node coefficient including both bridge weights and all original weights and inverse widths. It supplies

    sum_h b_h <= C L^15,
    sum_h b_h sum_v 1/t_hv <= C L^15 [1+log(1/eta)],
    max_h b_h <= C,
    max_h b_h sum_v 1/t_hv <= C,                            (10)

where eta is the GENUINE minimum private variance of the already finite old program. All old original-clock cutoff/target errors remain outside this new exact-family comparison. Known fixed-rank adapter/readout constants are multiplied back in.

Choose the same fixed v-sized source-zero readout geometry for every rank-nine packet, using (2)-(4). Its per-packet inverse readout constants now depend on this fixed rank and v, not on the number of nodes. Use

    rho_nonroot = alpha^(1/16),
    rho_root,h = signed b_h * fixed normalization
                              * alpha^(17/2).              (11)

The exact main product is alpha^9 times its literal coefficient. From (10), the summed root size is `S<=C L^15 alpha^(17/2)`, and the complete captured-caller first has the same alpha grade times the stated logarithm. All per-node guards also hold at sufficiently small alpha from the maxima in (10). Thus (9) pays the complete intrinsic own return, INCLUDING cross-group terms, at

    alpha^17 * public-log factors * sqrt(D).                (12)

After integrating all replica coefficient banks, the simultaneous Riesz lemma pays the FULL old/new bank join against a literal O(alpha) old first at alpha^10, and the self-bank term at alpha^18 times the one-hit logarithm. No new target-heat mismatch is present. Finite cubature/source/native/arithmetic/replay errors must be chosen to total alpha^10 sqrt(D).

Accordingly the ideal positive family has the alpha^9 -> alpha^10 bound. It becomes an ACTUAL native complete exact-finite-target mixed repair only after the extra hybrid condition (13) or (14), the heat theorem's finite positive quadrature/error port, and the actual native/first/Sobolev guards are supplied. If eta is bounded below by a fixed power of alpha and the original endpoint-panel counts have their declared public-log bounds, the logarithms in (10) do not spend algebraic grade. If eta has worse dependence, enter it explicitly. This is a fixed-family claim; the original sampler's other currents and prior target debts remain.

## 6. Node counts still cost queries

Common-carrier amalgamation removes physical readout dilution, not quadrature work. Every node's raw packet still executes. For the triple-cubic main, max k_v=3, sum k_v=7, and there are nine original forces. At a given history/node the raw active-source bill is `sum_v 2^(k_v+1)` before wrappers; the exact maximum in the 81 force-hit assignments for one labelled tree is 44 raw VALUE calls (20 in the twice-hit copy and 12 in each once-hit copy). There are 243 force-hit histories across all three labelled trees per ordered old-history triple, before physical permutations.

If a positive bridge rule has `N_BKAR<=C alpha^(-nu_B) public_log^m`, its work contributes nu_B to the complete max-plus query path exponent. Native repetitions, original caller captures, requested higher-precision old versions and all ancestor/anchor replays contribute their own exponents. Allocating total numerical error alpha^10 across N_BKAR groups generally tightens individual floors by that count; the cost of those tighter floors is charged by the actual native compiler.

The exact cost is the sum over ALL ordered histories, bridge nodes and force hits of their full original-VALUE programs and every replay. A finite grid is an acceptable existence rule but may have a large nu_B. Nothing in the shared-carrier construction proves public-log cubature, a bounded-in-order native exponent, or `c(P)=o(P)`.

## 7. An explicit conservative cubature query exponent

The heat worker's `FINITE-CUBATURE-AND-N-COPY-RULE.md` gives a positive midpoint rule with

    Q_T <= [2+ceil(Lipschitz/(2 epsilon))
                +ceil(log2(1/eta))]^(n-1).                 (16)

Its known Lipschitz constant is the sum of Price extra-two-hit source envelopes at the ORIGINAL private widths. For the triple-cubic main, K=7 and a Price derivative adds two, so the safe bound is `Lipschitz<=C eta^(-9/2)` after summing literal old weights and fixed-rank factors. To approximate a normalized alpha^9 coefficient at alpha^10, choose `epsilon=c alpha`. If `eta>=c_eta alpha^q_eta`, then

    Q_T <= C alpha^(-(2+9 q_eta)) public_log^m.             (17)

Thus this fallback's explicit new cubature exponent is at most `2+9q_eta`. The fixed native compiler's local exponent is zero by LOW30's serial recurrence; complete caller/ancestor versions add their actual exponents through the max-plus rule. The coefficient count is genuinely finite even when expensive. Equation (17) is a conservative cost certificate, not an optimality claim, and the actual native shared-carrier hybrid remains the separate condition of Section 4.
