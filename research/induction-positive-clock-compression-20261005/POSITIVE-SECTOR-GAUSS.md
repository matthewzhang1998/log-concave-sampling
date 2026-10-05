# Positive sector-Gauss compression at a genuine private-heat floor

2026-10-05. New, constructive accuracy-cost result. This replaces the inverse-tolerance power in the sealed midpoint BKAR rule by logarithms. Its cutoff cost remains quadratic for a two-edge forest. It preserves the original finite coefficient target up to a declared deterministic cubature error. It does not introduce smoothing, a signed quadrature, a coefficient-estimation oracle, or a carrier-retaining native LAW assumption.

## 1. Result and exact scope

Consider one of the sealed finite three-argument families:

* Three full old cubic sums: N=9 force occurrences, K=7 total primitive order, target alpha^9, main maximum primitive order 3.
* Full or diagonal 3/3/6: N=12, K=10, target alpha^12, main maximum primitive order 4.

Use the COMPLETE original coefficient bank, exact history/cross-history census, original private shields, affine injection rows, and actual finite old source versions. Let 0<eta<=1/2 be a lower bound on their ORIGINAL private variances. Where imported clock envelopes require structural floors, eta also lower-bounds those recorded floors. Let W be the literal sum of unsigned non-width scalar factors over ALL replica trees, new force-hit assignments, and finite original history/cross-history choices. It is a full-family sum, not a per-node or per-tree factor. All physical/source/readout factors that are not absorbed into separately displayed constants belong in W.

For any normalized coefficient tolerance 0<epsilon<1, there is an executable positive two-edge rule with

    Q_tree <= 2 m^2 [ceil(32 P/eta)+J]^2,
    J = max(1,ceil(log2(32 B_*/epsilon))),
    m = max(1,ceil(log4(32 B_*/epsilon))).                 (1)

Here P is a known Price hit-count/row majorant, B_* a known main-and-required-caller envelope defined below, and rank constants remain literal. Thus

    Q_tree = O( [eta^(-1)+log(B_*/epsilon)]^2
                   log^2(B_*/epsilon) ).                (2)

For the triple cubic, P=27 is safe when every relevant injection row has norm <=1; for 3/3/6, P=45. The full three-tree factor and original history/hit census are additional: 243 new hit histories for a fixed ordered cubic history triple, or 648 for a fixed 3/3/6 history. Equation (1) is per tree and is a node bound, not the whole VALUE bill.

The error is <=epsilon sqrt(D) in HS and <=epsilon in every proper cut, uniformly in fixed external affine callers. The same node rule can cover any fixed finite number of external-caller derivatives by enlarging B_*. All Gauss high derivatives are analytical certificates only; the executed main sources remain C3 or C4 at most.

Crucially, if eta>=c alpha^q_eta with FIXED q_eta and epsilon=alpha^p, the NEW bridge-node inverse-alpha exponent is at most 2 q_eta, independently of p. Hence its ratio to p tends to zero for this fixed family/floor regime. This is not an any-order theorem when the inherited q_eta or complete old/replay exponent itself grows proportionally to the desired order.

## 2. A factorial Price bound with one Hilbert factor

Write the normalized force tensor as A^(-1) E D^(k+1)g(z+sigma Z). The active-source theorem gives proper-cut bound

    2^(k/2) sqrt(k!) sigma^(-k),

and full HS bound <=(k+1)^(k/2) sigma^(-k) sqrt(D). Both are bounded by

    C_k sigma^(-k),  C_k=2^(k+1) sqrt(k!),                (3)

with the single sqrt(D) only on a chosen marked root. The elementary factorial inequality (k+1)^k<=4^k k! proves the HS assertion: the ratio of consecutive (k+1)^k/k! terms is (1+1/(k+1))^(k+1)<4. Use the actual A normalizations, rather than silently setting A=1.

Fix an open ordering sector for the two BKAR parameters. Gaussian Price differentiation adds an edge between DISTINCT replicas and one force hit at each endpoint. It preserves the original marked spanning tree; parallel edges are permitted, self-loops are not introduced. For any one scalar coordinate in which every covariance entry is affine with slope at most one, let

    P >= sum_(i<j) N_i N_j a_i a_j,

where a_i bounds the actual coefficient-bank injection rows of argument i. Increase P to at least one. At derivative order l, the triangle sum has at most the weighted multiplicity P^l. All primitive orders sum to K+2l. The marked-tree all-cut/one-Hilbert theorem and (3) imply

    ||F^(l)||_* <= W P^l 2^(N+K+2l) sqrt((K+2l)!)
                         eta^(-K/2-l).

Here ||.||_* denotes any proper cut or HS/sqrt(D). Using

    sqrt((K+2l)!) <= 2^(K/2+2l) sqrt(K!) l!,

we obtain the explicit bound

    ||F^(l)||_* <= B_K l! (16P/eta)^l,
    B_K = W 2^(N+3K/2) sqrt(K!) eta^(-K/2).              (4)

Replace the fractional powers in this displayed convenient constant by any rational upper bound in execution. This is a dimension-sharp analytic radius certificate at the ORIGINAL private floor. It differentiates the exact Gaussian expectation, never a numerically chosen square root or the hit-adapted heat allocation. Singular coefficient-bank covariance is handled by independent positive regularization and its limit, as in the sealed Price proof. The strictly positive original force shields remain unchanged.

A q-fold external affine caller derivative adds q original force hits with their ACTUAL rows. Thus use K+q in (4) and replace W by its finite q-hit scalar sum, including every caller/ancestor Jacobian that the source contract actually requires. For main plus first, B_*=1+B_K+B_(K+1) is safe with these distinct FULL-FAMILY W values. Equivalently form per-tree main/first constants B_T and take B_*>=1+sum_T B_T. This convention pays the total family error once; a per-history use instead needs an explicit tolerance allocation whose sum is at most the desired error. Higher analytic Price order does not require executing higher derivatives of g.

## 3. Remove the min-path walls before using high order

In either of the two ordering sectors write

    t_large=r,  t_small=r s,  0<=r,s<=1.

The sector Jacobian is r. A three-replica tree covariance has one off-diagonal entry r and two entries rs, after relabeling. Thus each entry is separately AFFINE in r and s. Define

    G(r,s)=r F(t_large=r,t_small=rs).

Pure r and pure s derivatives, the only ones the tensor-product error proof needs, obey

    ||partial_r^l G||_*, ||partial_s^l G||_*
                      <=2 B_* l! (16P/eta)^l.           (5)

For r derivatives use the two terms r F^(l)+l F^(l-1); eta<=1/2 and P>=1 make the factor two safe. No mixed-derivative bound or nonexistent global smoothness through the ordering wall is assumed. The diagonal itself has zero integration mass.

Choose zeta=2^(-J), with J from (1), and keep only the square [0,1-zeta]^2 in each transformed sector. Both omitted strips together have area <=2zeta. From (4), both sectors' total omitted coefficient is <=4B_T zeta for tree T, hence <=4B_*zeta for the entire family, including each specified caller derivative. This is a declared quadrature tail, not deleted target mass or a new original-clock cutoff. The quadrature approximates the full original finite coefficient to its stated epsilon.

## 4. Positive dyadic panels and exponential accuracy

Partition the gap variable 1-r into [2^(-j-1),2^(-j)] for j=0,...,J-1, and likewise 1-s. Subdivide each panel into equal intervals of length at most

    h=eta/(32P).

The number of intervals on either axis is at most ceil(1/h)+J. On each interval use the m-point Gauss-Legendre rule. All nodes are interior and all weights positive. If lambda_r,lambda_s are interval weights, the actual sector weight is

    w=lambda_r lambda_s r>0.                             (6)

No negative mixture probability or cancellation is used. The total mass is the exact area of the truncated sectors: 2 times integral r dr ds over the truncated square, at most one. Its deficit is accounted for by the tail estimate.

For an interval I of length ell<=h, the classical Gauss remainder and (5) give

    ||integral_I f-Q_I f||_*
      <= 2 B_* ell [ell(16P/eta)]^(2m)
          (m!)^4 / ((2m+1)((2m)!)^2)
      <= 2 B_* ell 4^(-m).                              (7)

The formula holds for Banach-valued integrands by scalar duality. Sum the panels. For a square, telescope I_r I_s-Q_r Q_s into the two pure-axis errors. Positivity makes both integration and quadrature bounded by their masses, at most one. Both sectors of tree T cost at most 8 B_T 4^(-m); summing trees costs at most 8 B_*4^(-m). Together with the tail and the conservative choices (1), this is strictly below epsilon/2; reserve the remainder for node/weight arithmetic. No tensor entry is evaluated to choose the rule.

## 5. The old heat budgets survive this non-product bridge rule

Set x=1-r, y=1-s. In either sector the two original bridge gaps are

    s_small=x,  s_large=x+y-xy.

The names refer to gap size, not the ordering of t. Every subinterval sits in a dyadic gap panel, whose length is at most any gap value in that panel. Thus lambda_r<=x and lambda_s<=y, and

    w <= xy <= s_small s_large.                          (8)

Inside a dyadic product panel, x and y vary by at most a factor two. So does x+y-xy: monotonicity and f(2x,2y)<=2 f(x,y), with endpoint clipping if necessary, prove this directly. Every replica innovation gap delta_i, equal to one of the two gaps or their minimum, also varies by at most two. Adding eta preserves comparability.

The positive scalar heat majorants in the sealed cubic and mixed notes are finite sums of products of (delta_i+eta) powers. On each cell their sup/inf ratio is bounded by their fixed-rank power sum, independently of m, eta, and the number of cells. Furthermore the Gauss rule integrates the Jacobian r EXACTLY on every product cell. Hence its weighted sum of each majorant is at most that fixed comparison constant times its exact sector integral.

Consequently all sealed integrated main/first estimates transfer to this rule, with no Q penalty:

Triple cubic:
    S0=sum |b_h| <= C L^15,
    S1=sum |b_h| sum_v a_v/t_v <= C L^15(1+log(1/eta)),
    max |b_h| <= C,  max |b_h| sum_v a_v/t_v <= C.

Mixed 3/3/6:
    S0<=C L^C,  max |b_h|<=C,
    S1<=C L^C eta^(-chi),
    max |b_h| sum_v a_v/t_v<=C eta^(-chi),
    chi=1 for the full cross-history family and 1/2 for its certified diagonal.

Here b_h contains the actual weight (6), original weights at their literal powers, inverse widths, injection scalars, and declared scalar normalizations. Fixed source/readout factors and complete caller rows are restored as in the sealed input. Distinct repeated-label histories are never replaced by independent clocks.

The literal root-square sums follow WITHOUT an independent-clock variance argument:

    sum |b_h|^2 <= (max |b_h|) S0,
    sum |b_h|^2 (sum_v a_v/t_v)^2
       <= [max |b_h| sum_v a_v/t_v] S1.                  (9)

Thus the triple-cubic square-first allowance is logarithmic, while mixed has eta^(-2chi). The positive quadrature preserves the existing scalar admission tests, rather than multiplying them by its node count.

## 6. Actual retained bank, private heat, and native status

At EVERY retained node execute the sealed exact disintegration

    B_i=sqrt(c) B+(R V)_i+sqrt(delta_i) E_i

on the same COMPLETE old bank B and the required new complete banks. Keep the original injection rows in every Price/BKAR contraction. Use the same deterministic, simultaneous m_star heat allocation for all callers of that main history. For 3/3/6 retain the recorded old six-packet disintegration and extract only from the new independent innovation as specified in the input. Original private heat is not discarded, increased, or replaced by eta; eta is only a lower bound.

The cubature certificate concerns the bank-mean coefficient and uniform external callers. It does NOT claim a pointwise-in-B coefficient identity or equality of the bank-conditioned laws. All old/new mixed currents therefore remain owed and are paid by the actual sealed S1/caller and same-endpoint ledger. No tape is exposed as an unattenuated observer. All complete captured-caller first paths are the actual source paths with the original derivative rows and strict ancestors. Frozen clock labels and covariance-root choices are not differentiated as dynamic inputs.

For independently owned native groups, use their actual admitted readout shares and inverse normalization products; if those shares depend on Q, charge their resulting powers. The proposed common-carrier geometry would avoid this readout dilution ONLY after the separate actual carrier-compatible native replacement theorem is proved. This quadrature theorem does not infer that theorem from marginal LAW and does not alter its status.

## 7. Finite precision and complete work

A positive rational implementation can approximate Gauss nodes and weights without losing (8)-(9). Isolate every Legendre root by a disjoint rational sign-changing interval; degree counting certifies one root per interval. Evaluate its positive weight formula with rational interval arithmetic. Choose paired symmetric rational nodes and equal paired positive weights, and normalize each interval's weights to its EXACT length. Symmetry preserves its first moment. Thus constants and the Jacobian r stay exact, positivity remains literal, and all nodes remain in their original cells.

If the reference [-1,1] node perturbation is <=d and its total weight perturbation <=e_w, positivity and the pure-axis derivative bound give an additional full-family coefficient floor <=64 B_*(P/eta)d+2 B_*e_w. For example d<=epsilon eta/(512 B_*P) and e_w<=epsilon/(512 B_*) use less than epsilon/4. Mapping and summing panels does not multiply this budget by Q, because their positive masses sum to at most one. The executable uses these exact rational budgets. The required bit precision is polynomial in log(1/epsilon), log(1/eta), log Q, and the finite rank/scalar input bit sizes. The included executable uses exact rational sign/interval checks and an adaptive precision retry until both requested error bounds hold; floating/high-precision root estimates merely propose intervals. Its production parameter chooser uses rational upper bounds for B_* and integer power comparisons for J and m. The separate floating-point census is only a displayed estimate.

Small known replica/query covariance roots and all square roots are computed to separately certified scalar-center/width floors. Degenerate roots do not justify adding new target heat: use explicit positive partition factors or certified nonnegative scalar roots, and charge their approximation through actual source/caller bounds. With source width t, original VALUE precision delta_g costs at most 2 delta_g/(A t); every changed center, private root, width, requested floor, source version, or ancestor requires full affected call-and-anchor replay. Absolute VALUE errors do not imply Sobolev-first or curl errors. The native source/filter/arithmetic Sobolev floors remain independent certified obligations.

The complete original-VALUE work is

    sum_(tree,old history,new hit,node) sum_v
         2^(k_v+1) M_native(k_v,b_v,epsilon_v)
    + all old complete requested versions, caller captures,
      known scalar roots, physical rows, keeps and affected replays.

For triple cubics the raw main force maximum per hit history remains 44 VALUE calls before wrappers. The Gauss degree m does not raise that executed source order. All Q nodes are really executed; no expectation oracle or N^(-1/2) sampling estimate is hidden.

At fixed graph order the imported LOW30 wrapper has public-log dependence on its requested numerical floors under its stated guards. Consequently the NEW bridge NODE-COUNT exponent in the fixed-power floor regime is 2q_eta. Turning that into a complete query exponent still requires the actual native gap/readout/radius and replay ledger. The sector strip creates innovation gaps as small as the accuracy-dependent zeta. This design never inverts those coefficient-bank gaps: they appear only in explicit Gaussian roots with bounded center rows. If a chosen native implementation nevertheless introduces an inverse gap, that additional cost must be charged; the node theorem does not erase it. Every complete old source, retained-bank ancestor, and replay exponent still enters the full max-plus cost, and actual readout dilution may add more. No global c(P)=o(P), heat-restoration theorem, or finished any-order sampler follows from this result alone.
