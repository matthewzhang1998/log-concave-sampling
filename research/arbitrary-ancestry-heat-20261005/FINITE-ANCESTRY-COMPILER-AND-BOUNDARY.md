# Arbitrary finite ancestry: exact heat, finite VALUE allocation, and a sharp historywise boundary

2026-10-05. This extends the canonical symmetric-OU and second-Riesz machinery where it genuinely extends, and gives a source-qualified counterexample to a universal historywise clock bound. It does not claim a new all-current sampler or eventual-sublinear complexity theorem.

## 1. Outcome and scope

There is a computable exact heat compiler for every fixed finite positive interior clock record. There is also a computable original-VALUE main-coefficient allocation for that fixed record, under the already imported fixed-order native pair/filter contracts. Both may have constants polynomially large in inverse original widths. The second-Riesz positive-law remainder needs only a first-bank derivative at any fixed physical polynomial degree.

What fails is the stronger desired conclusion that arbitrary admitted fixed-rank ancestry has clock-summed, cutoff-uniform/public-log constants of the kind proved for the old three-force packet. Two high-degree original force occurrences can share the same last structural parent. Their private coarse roots are then the only noise that distinguishes them. At rank eight a literal double-star history gives an actual bounded-Hessian example with a divergent single-node old coefficient covariance after its original positive weights. Exact Gaussian redistribution preserves, rather than removes, this covariance.

This separates a constructive finite compiler from an unproved all-order uniform theorem. It also separates source-genealogy obstructions from arbitrary ill-conditioned Gaussian rows, and stored replay roots from actual retained law observers.

## 2. Literal ancestry matrix, not an arbitrary row input

Start with an exact labelled cumulant AST produced by the connected generator. A resolvent R_k has a positive clock r with density r^(k-1) dr. Keep its original order k, multiplicity and scalar weight. Traverse from the retained deterministic caller y.

At a structural resolvent node, replace the current row X by

    r X + sqrt(1-r^2) G_e.

At a primitive resolvent ending at force vertex i, put

    q_i = r_i X_parent + sigma_i W_i,
    sigma_i^2=(1-r_i^2)/2,

and keep a second, independent private Gaussian of variance sigma_i^2 in that vertex's original force average. Thus the total old OU variance is unchanged. Each structural G_e and each W_i has its own recorded column. Shared ancestors are shared columns; never duplicate them merely to make products independent.

The traversal produces q=d y+A Q, Gamma=A A^T. For this literal genealogy,

    Gamma_ii+sigma_i^2+d_i^2=1,
    ||A_i||<=1, |d_i|<=1,
    Gamma >= diag(sigma_i^2)>0.

These are exact row identities. They also retain the original derivative injection A_i. The injection in P_r(D_Q L) is A_i, not r A_i.

For a finite collection sharing roots, construct one common column dictionary and group the old coefficient before differentiation/covariance. Cross-node and cross-history terms stay present. A collection declared independently banked has disjoint dictionaries, as in the canonical independently banked node convention.

## 3. A finite certified floor that can exploit ancestry

The trivial certificate is ell_0=min_i sigma_i^2. An entirely scalar finite algorithm can improve it while displaying exactly which primitive shields it avoids.

Choose a protected row set S, let R be the other rows, and choose |S| STRUCTURAL columns J such that T=A[S,J] is invertible. Primitive W columns are not used here, so their R diagonal heat is not spent twice. Put U=A[R,J], eta=||T^(-1)||op, and x_R=min_(i in R) sigma_i^2. For nonempty S,R, a certified floor is

    ell_(S,J)=min( x_R/(1+2 eta^2 ||U||op^2), 1/(2 eta^2) ).

Indeed z=A[:,J]^T v gives v_S=T^(-T)(z-U^T v_R), hence

    ||v||^2 <=(1+2 eta^2 ||U||^2)||v_R||^2+2 eta^2||z||^2,
    v^T Gamma v >= x_R||v_R||^2+||z||^2.

For S=all rows use ell=1/eta^2; for S empty use ell_0. Enumerate all finite (S,J), discard singular minors, and take the maximum certificate. Norms and inverses are those of known rank-sized scalar matrices, never source tensors. A certificate needing no spectral computation replaces eta^2 by ||adj(T)||_F^2/det(T)^2 and ||U||op^2 by ||U||_F^2; both replacements only decrease the certified floor. Interval upper bounds for the norms and certified nonzero determinants suffice for fully rigorous numerics; the trivial floor remains available if a minor is numerically ambiguous.

This is a concrete rank-indexed algorithm. Its constants are computable from the actual clocks. It does not assert that the inverse minors have a harmless endpoint power. In the counterexample below, every structural minor trying to protect both sibling centers is singular for an exact genealogical reason.

## 4. Exact symmetric heat redistribution for any number of sites

At bridge correlation r=sqrt(c), h^2=1-c, choose any diagonal delta>=0 with Gamma-diag(delta)>=0. The floor algorithm gives the safe choice delta_i=ell/2. Set

    t_i^2=sigma_i^2+h^2 delta_i,
    B=(Gamma-diag(delta))^(1/2),
    z_i^+ =d_i y+r(AQ)_i+h(BV_+)_i,
    z_i^- =d_i y+r(AQ)_i+h(BV_-)_i.

Use independent V_+,V_- and independent original-VALUE private shields t_i Z_i^+,t_i Z_i^-. Then the within-copy covariance is Gamma+diag(sigma_i^2), and the cross-copy covariance is c Gamma. All means are unchanged. Every residual affine row has squared norm Gamma_ii-h^2 delta_i<=1. No 1/h or 1/ell factor enters the actual affine query rows.

This represents precisely

    Cov_Q L = integral_0^1 E[(P_sqrt(c) D_Q L) contracted_Q
                              (P_sqrt(c) D_Q L)] dc.

It never reheats the target. Known scalar roots B are permitted Gaussian geometry, not expectation/covariance oracles for g. Unequal shields can be optimized by a finite semidefinite calculation, but that optional optimization is unnecessary for the stated construction and cannot defeat the principal-submatrix obstruction in Section 9.

Reuse the positive dyadic Gauss bridge rule with sum u/h<=1/(sqrt(2)-1), and its certified scalar multiplier tolerance. For a frozen finite old derivative envelope M_old, delta_op<=e_bridge/max(1,M_old^2 D) gives conservative numerical tensor error <=e_bridge. This numerical estimate may use D, separately from one-HS analytical estimates. It is not a claim of cutoff-uniform M_old.

## 5. Complete observer test before calling this a joint-law return

The zero-readout hypothesis concerns ALL actual retained observers, not merely the displayed physical output row or the original Q subblock.

Form B_total from every coefficient root being integrated by Riesz: original structural/coarse Q, new residual roots, and any reset/auxiliary root whose coefficient dependence is now included. For every source-zero physical output, retained affine observer and untouched Gaussian keep, verify its actual row on B_total is zero. Also verify that the nonlinear base and nonlinear retained observers are independent of B_total, or put their bank derivatives explicitly into the same-endpoint current census.

A replay tape may remain stored internally without being a law observer. Appending Q itself to the comparison observable invalidates the zero-bank observer test: its row is the identity. Similarly, a reset innovation with a nonzero old physical readout cannot be promoted silently into a zero-readout coefficient bank.

For affine observations OQ there is a computable conditional alternative: disintegrate Q at fixed OQ and replace A by A Pi, where Pi projects onto ker O. The applicable covariance is A Pi A^T. Recompute the floor and the conditional mean target. A direction observed exactly has no residual conditional heat. This produces a different, explicit conditional target; it does not prove the original unconditional cancellation jointly with the observer. A general nonlinear observer gives a non-Gaussian conditional bank and is not admitted by this Gaussian compiler.

## 6. Literal VALUE sources and a finite root-heavy allocation

Suppose the old tree has n forces of spatial orders d_i, and a covariance derivative copy hits vertex j. At that copy's vertex i the required fixed-probe adapter is C_(d_i-1+1_(i=j)). All forces keep their physical marks. The two tree copies are joined at the differentiated ports, so the resulting marked graph is still a tree, with 2n original force occurrences. Root it at a genuine physical leaf and execute its selected path and disjoint side programs in the admitted chronological order.

At every vertex execute the original source

    f_i(x)=[g(z_i+t_i x)-g(z_i)]/(A_original t_i),

using its same-center anchor, finite signed affine C_k formula, separate finite filters, bounded-field preparation, calibration and complete selected-pair program. The source has private first <=1 and center first <=2/t_i. C_k uses 2^k original affine VALUE terms before its filter/pair replication. Original HVPs occur only in requested sweeps at recorded VALUE sites.

Write the exact required scalar product of amplitudes as

    product_v rho_v = alpha^(2n) K_j,

where K_j includes the positive bridge mass u, both literal old weights, old integer/physical permutation constants and signs, the original derivative-row contraction, inverse readout shares, inverse selected-matrix constants, and

    product_(copies,i) t_i^[-(d_i-1+hit_i)].

All factors in K_j are known scalar geometry except the separately executed normalized VALUE sources. No derivative tensor is supplied to the program. Do not reinsert old sigma denominators already canceled by the change to exact t-normalization.

For a FIXED finite clock list and n>=2 choose m=2n, root exponent 2 and nonroot exponent

    beta=(m-2)/(m-1)>0.

Set rho_root=sign(K_j)|K_j| alpha^2 and rho_v=alpha^beta at the m-1 other vertices. Their product is exactly K_j alpha^m. Zero K_j nodes are omitted. The root-heavy factor is not an amplitude approximation.

Let C_j be a computed upper envelope for all imported fixed-order source/pair constants, inverse t_i, known query/observer row norms, chosen positive readout inverse powers and every finite path sum in this graph. The root residual/private first is at most C_j |K_j| alpha^2. A coefficient caller at a nonroot has the root and that source on its actual ancestor path, hence at least alpha^(2+beta), apart from its explicitly included t_i^-1. A carrier derivative uses strict ancestors only; its first root ancestor still gives alpha^2. Actual branching contributes the enumerated finite sum, not an independence assumption.

Thus all complete radius and first guards of the finite graph hold after decreasing alpha until

    |K_j| alpha^2 <= native_root_gap/C_j,
    alpha^beta <= native_nonroot_gap/C_j,
    C_j |K_j| alpha <= requested_first_gap.

Take the minimum threshold over the fixed list. This is finite and computable. The threshold can depend badly on inverse original cutoffs; nothing in this step proves it uniform in a family whose clock list changes with alpha. For such a family substitute the actual clock powers into these inequalities; do not freeze them while taking alpha small.

A more economical allocation is a finite logarithmic linear feasibility problem: exact sum of log amplitudes, one upper bound per source radius, and one inequality per actual caller/carrier path. The explicit allocation above proves feasibility for a frozen finite list and avoids using an optimization oracle in the construction.

Every actual source prior has its actual source radius and ancestor attenuation. At fixed frozen widths, all radii have positive alpha powers, so a finite pair order suffices for any requested absolute native floor. Allocate finite means, filters, clipping, quadrature, arithmetic, anchors and replay separately after graph/readout/path constants are frozen. Old calls must be rebuilt if their exact cached version does not meet the new floor.

This is an executable main-coefficient rule conditional on the imported native contracts. Its finite reverse jet still has to be generated, including all covariance roots, side substitutions, physical tilts, old/new products and all observers. It does not assert that the emitted queue has a uniform surplus invariant. That missing assertion cannot be obtained from Gaussian factorization alone.

### 6.1 Explicit rank/width derivative and cost envelopes

For an n-force original tree with spatial degrees d_i, put e_i=d_i-1, so sum e_i=n-2. Split-heat Hermite bounds give a conservative local fixed-order cut constant 2^(e_i) e_i!, and hence a complete coefficient constant at most

    B_n=2^(n-2)(n-2)!.

A single additional old-bank derivative gives the fully explicit crude allowance

    M_(h,1)=n 2^(n-1)(n-1)! |w_h|
                 product_i sigma_i^(-e_i) max_i sigma_i^(-1).

The original row norms are at most one for the admitted AST; for another already-authorized row representation multiply by their actual recorded norms. Sum these quantities over histories before estimating a shared-bank coefficient. This is a valid finite envelope, not a good uniform endpoint estimate. At a covariance node replace sigma_i by its exact t_i and use the exact derivative-hit order rather than the maximum whenever possible. All physical/proper derivative cuts use the marked-tree contraction, retaining the one HS factor sqrt(D) at the unique energy slot. The corresponding finite physical-chaos matrix estimate has degree at most n-1 and its explicit fixed-degree/Hölder constants plus the recorded dimension logarithms.

For minimum frozen width s_min, the displayed bound is at most n 2^(n-1)(n-1)! |w_h| s_min^(-(n-1)). A same-rank covariance node requires total inverse-width degree 2n-2. These are directly computable fallback estimates; Sections 8–10 show why their inverse widths cannot generally all be replaced by public logarithms.

The covariance-node native census has one occurrence of every vertex in each of the two original histories, and one promotion C_k to C_(k+1) at the selected derivative hit in each copy. Enumerate every ordered old-node pair, every bridge node, both derivative-hit choices, and the chosen physical permutations. Its additional VALUE bill is exactly bounded by

    sum_(enumerated graphs) sum_(vertices v) Q_(C_k(v), b_v, precision_v)
       + Q_known_scalar_roots + Q_capture + Q_readout + Q_scalar_encoding
       + Q_required_replay.

Each Q is the COMPLETE finite original-g source, including its own origins/anchors, all 2^k signed affine terms, finite filters, private/pair/calibration replicas and absolute precision. Required old calls appear separately at their newly required full precision. An old cache contributes zero incremental work only for an exactly matching caller/source/order/precision version. Root-heavy allocation does not reduce any of this work by fiat.

## 7. All-fixed-degree second Riesz with explicit constants

Here the positive-law step genuinely generalizes without a second derivative of an original coefficient.

Let Q and P be independent finite-dimensional standard Gaussian banks, with P the physical publics. Let F(Q,P) be centered in Q at every P and polynomial of total physical degree at most k. The coefficients may only be first Sobolev functions of Q. Assume

    ||F||_L2(Q,P) <= E,
    sup_Q ||D_Q F||_L4(P;op) <= L.

The latter can be supplied by the complete extended coefficient cuts and fixed-degree matrix-chaos estimates. It must not be replaced by a scalar-test covariance bound or a bound after only averaging Q.

Put R_Q=D_Q N_Q^-1, tau=(R_Q F)(D_Q F)^T and C(P)=E_Q tau=Cov_Q F. With a Q-independent base b(P) in L2(P), define

    W_t=b(P)+tF(Q,P)+sqrt(kappa/2)G0
                  +[(kappa/2)I+(1-t^2)C(P)]^(1/2)G1.

G0,G1 are independent of Q,P, and kappa>0. This is a positive analytical reference path; the covariance root is not executed as a source tensor. At the SAME W_t,

    d E phi(W_t)/dt
       =t E[(tau-C):D^2 phi(W_t)]
       =t^2 E[(R_Q(tau-C)) D_Q F:D^3 phi(W_t)].

The second integration by parts differentiates only the endpoint, not tau. There is no D_Q^2 F term.

Hilbert-valued Gaussian hypercontractivity gives a factor 3^(k/2) from physical L2 to L4. Expand in a finite orthonormal physical Wick basis before each Q-Riesz application. Since Riesz is an L2 contraction on centered Hilbert coefficients,

    ||tau-C||_2 <= 3^(k/2) E L,
    ||(R_Q(tau-C)) D_Q F||_2 <= 3^(3k/2) E L^2.

There is exactly one Hilbert energy. Physical degree at the second Riesz stage is at most 2k. Centering is an orthogonal projection, so no factor two is required. Transfer the last two test derivatives through the constant independent G0. Gaussian Hermite isometry contributes sqrt(2)/(kappa/2), and integral_0^1 t^2 dt=1/3. Therefore

    W2(Law(W1),Law(W0))
       <= [2 sqrt(2)/(3 kappa)] 3^(3k/2) E L^2.

This finite-degree statement includes the cubic case k=2 and gives an explicit, deliberately conservative rank constant. It also supplies the usual mixed first-Riesz bound E_new L_old with the corresponding physical-degree hypercontractivity factor when centering a new coefficient at a base still containing the old field. The old bank has not disappeared just because the new coefficient was centered.

The observer test in Section 5 is essential. If a retained observer depends on Q, its derivative also appears in D_Q of the test endpoint, with no alpha-small coefficient guaranteed. The present formula then has additional currents. It is not a joint-law theorem with arbitrary old keeps appended.

## 8. A legal rank-eight genealogy that exposes the limit

Use the exact connected AST generator. Start from the rank-two sibling pair. Add six rank-one force leaves one at a time using the i=n-1 term in the rank-n recurrence. For the first three attachments let the product-rule derivative hit sibling 1; for the last three hit sibling 2. The resulting tree has degrees

    (4,4,1,1,1,1,1,1).

Its edges are the original sibling edge, three leaves at center 1, and three leaves at center 2. The history multiplicity in this ordered construction is 2*3*...*8=40320. Both center primitive resolvents are R5. The seven structural resolvents have index 8 after the literal derivative increments. The six new leaf primitives are R2. These labels can be checked recursively without generating every rank-eight history.

Crucially the two centers are siblings at the DEEPEST structural parent X. All later derivatives add edges and increment resolvents; they do not insert a new independent structural root between those two primitive children. Take all seven structural clocks and six leaf primitive clocks in fixed interior panels. Conditional on y=0, X is centered Gaussian of variance v bounded below. Put the two center primitive half-heats equal to x and r_center=sqrt(1-2x). Then their BEFORE-private-shield covariance block is exactly

    Gamma_pair=v r_center^2 [[1,1],[1,1]]+x I_2.

This is a legitimate original genealogy with all old injection rows bounded by one.

## 9. A Gaussian impossibility inequality for simultaneous independent shields

For ANY diagonal extraction delta_i from h^2 Gamma, positive semidefiniteness tested on the vector supported at these centers with values (1,-1) gives

    delta_1+delta_2 <= 2 h^2 x.

Including the original private shields,

    t_1^2+t_2^2 <=2(1+h^2)x.

Thus each t_i remains O(sqrt(x)). No choice of scalar square root, unequal allocation, unused outer structural variance or additional residual roots can give both centers an O(1) independent private heat. This is a principal-submatrix necessary condition and is independent of the other six sites. It is not an artifact of using the minimum eigenvalue floor.

## 10. An actual C2 source with a divergent weighted node covariance

Work in D=1 and A_original=1. Fix 0<b<a, a+b<=1, for example a=1/2,b=1/4, and set

    g_K(z)=a z+(b/K) sin(K z).

It is the gradient of a smooth convex potential, g_K(0)=0, and a-b<=g_K'<=a+b<=1 uniformly in K. Hence the original C2 assumptions hold uniformly, despite higher derivatives growing.

Take x=K^-2, K>sqrt(2). The privately averaged degree-four center factor is exactly

    P_(1/K) g_K''''(r_center X+W_i/K)
       =b K^3 exp(-1/2) sin(T+W_i),
    T=K r_center X.

The six privately averaged degree-one leaf factors have private widths bounded below, so each is a+O(exp(-c K^2)), uniformly in its old coarse center and all correlations. Their product is a^6 plus a uniform exponential error. This justifies replacing them by a^6 in the leading asymptotic without replacing the actual source.

For positive center clock masses w_1=w_2=K^-2 (or fixed comparable masses), suppress only fixed positive structural/leaf/permutation factors. The literal old node coefficient satisfies in L2

    L_K = b^2 a^6 K^2 exp(-1) S_K
                          +O(K^2 exp(-c K^2)),
    S_K=sin(T+W_1) sin(T+W_2),
    V=Var(T)=K^2 r_center^2 v.

With independent standard W_1,W_2, exact Gaussian characteristic functions give

    E S_K=exp(-1)(1-exp(-2V))/2,
    Var(S_K)=1/4+exp(-4)(1+exp(-8V))/8
                    -exp(-2)(1+exp(-4V))/4.

In particular the variance tends to

    v_infinity=1/4+exp(-4)/8-exp(-2)/4>0.

Therefore

    Var(L_K)/K^4 -> b^4 a^12 exp(-2) v_infinity>0.

There is an additional exact decomposition:

    S_K=cos(W_1-W_2)/2-cos(2T+W_1+W_2)/2.

The two Gaussian arguments are independent. The first term alone has variance (1-exp(-2))^2/8>0. Its symmetric-OU covariance integrand, in correlation rho, is

    2 rho * [exp(-2(1-rho^2))(1-exp(-4rho^2))/4].

This is strictly positive on every fixed interior bridge interval. The K^4 growth is therefore not concentrated in a removable near-identity bridge strip.

This is the actual covariance of an old coefficient on its complete old Q, not merely an inverse-width upper bound or one bad derivative-hit term. Summing the product-rule derivatives first cannot remove it. Exact symmetric OU covariance grouping has this same positive integral by identity. The second-Riesz lemma remains true, but its E and L constants cannot be bounded by a rank-only/public-log quantity in this history family.

For a literal positive dyadic rule, a width-Delta primitive panel near one has R5 mass comparable to Delta. If it has N positive nodes, at least one node has weight >=c Delta/N, and every node has sigma^2 comparable to Delta. Take K=Delta^-1/2 and choose that same rule node for both primitive centers. The displayed limit has harmless constants x K^2 in a compact positive interval, and gives coefficient at least c K^2/N^2 and variance at least c K^4/N^4 asymptotically. Public-log N cannot absorb these powers. The constants from other fixed interior rule nodes and readout shares must be retained; if their number/share changes by public logs, that only adds public-log losses. A mere upper endpoint-mass envelope would NOT prove this lower bound; positivity, exact panel mass and a bound on the actual node count are used.

This disproves cutoff-uniform nodewise, absolute-sum and squared-weight per-history ledgers of this independent-vertex heat-redistribution family. In an independently banked old-node assembly, covariance contributions add and this bad node survives. For nodes sharing a bank, positive scalar weights alone do not imply nonnegative cross covariances; no lower bound for their complete sum is claimed here. Complete shared-bank history cancellations, integration-by-parts regrouping into different composite VALUE sources, or inverse-cutoff replication are possible different approaches and are not ruled out.

This is not claimed to be the earliest obstructed rank. Nor does it rule out every power-law cutoff/amplitude choice: at K=alpha^(-gamma), the old coefficient is alpha^(8-2gamma) and its covariance alpha^(16-4gamma), up to constants/public logs. Some finite target grades can tolerate those explicit losses. What fails is omitting them or asserting a universal cutoff-free extension. An added heat comparable to K^-1 preserves the phenomenon with changed constants; a fixed positive added heat suppresses it but changes the target and carries its own mismatch debt.

## 11. Exact stopping point and next useful invariant

The finite construction returns:

1. The literal AST and common original root dictionary.
2. Exact means, covariance and derivative injection rows.
3. A certified heat split and all observer checks.
4. Complete original-VALUE vertex/anchor/filter/pair calls, scalar coefficients, a legal finite amplitude allocation and its actual guard threshold.
5. The finite conditional/current jet queue, numerical floors and full replay bill.
6. An explicit PASS or FAIL for the desired alpha-dependent/public-log clock ledger.

It does NOT label a failing ledger as an all-order theorem. The rank-eight test shows why an all-rank extension must include more than a spectral-floor procedure: it needs either a proved cancellation/grouping invariant BEFORE independent node assembly, or a changed family of composite sources whose constants control these colliding high-jet siblings. The useful three-site and repeated-three-site sectors can still pass their stronger clock tests. Nothing here invalidates their bounded audit.
