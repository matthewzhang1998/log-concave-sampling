# Deterministic same-heat known-center source with a protected short refresh

2026-10-04. Source-level construction relative to the completed CW7 known-center retained-graph input. This document separates a standalone known-center family from an inaccessible-center mean/hidden-family service.

## Result and precise scope

Appending one short same-posterior Hamiltonian refresh supplies a legal protected row for the deterministic quarter-flow sampler. The new increment is used in the emitted finite program. It is removed from nonterminal queries only in the auxiliary proxy graph. Thus its width can be chosen independently of the requested posterior order.

At weak proxy mark J=29/10, choose

    v=19/15,  h=a^(v/2),  gamma=sin(h)^2,
    epsilon_tw=a^(9/10),
    beta=(gamma^(-1)+epsilon_tw^(-2))^(-1/2).

Then gamma is comparable to a^(19/15), beta to a^.9, and the full gradient radius at physical force radius r=a is a^.1. There is no R-dependent inverse protected width.

The deterministic gradient nodes can be weighted so their stacked Gaussian row and absolute edge operator stay bounded independently of their inverse-heat count. These weights are part of the actual source normalization, not a deletion of queries.

Under the input/output signature below, the **known-center** query recurrence is

    k_R <= max(k_P,R-3/2),   R=P+1/2.                 (K)

Starting from CW7 with k_7=0 gives k_R=R-3/2 for the displayed rebuilt finite construction at R>=7.5. This is a linear-order known-center family with actual first, retained graph, origin and fixed weak protected-proxy ports. It is not an o(R) result.

An arbitrary hidden-center service is not needed to define this family or to run a standard outer proximal Gibbs chain using it. If a different application requires an inaccessible center, its Gaussian-mean statistic must still be constructed and charged. The current weak-E mean arm can leave a multiplicative previous-hidden-source bill; that separate recurrence is printed in Section 9.

The source proof relies on the completed input graph and its original-gradient grammar, not on the input's marginal law alone. The short-flow/projection formulas have a separate independent audit. The weighted graph normalization and the complete join are stated explicitly here for independent review before any broader family admission is claimed.

## 1. Input signature at a known center

Normalize the original Hessian upper bound to one. Fix the heat a, potential, target order, caller/reference/moment lists, all grid and finite iteration counts and all numerical versions before differentiating the caller y. The original queries are grad V and its directional first/adjoint HVP actions.

The input is a completed finite known-center program X_P(a,y;W_old), with:

1. Conditional posterior law error Lambda_P sqrt(Dim) a^P at the actual admitted caller, for Q_(a,y).
2. Full fresh first Lambda_P sqrt(a), ordinary source profiles at the actual caller, and an actual center first I+O(Lambda sqrt(a)) for the CW7 seed or I+O(Lambda a) afterward.
3. An actual retained STATE on the same complete tape, with state error Lambda_p sqrt(Dim) sqrt(a) a^K, common finite origin, and a literal original-gradient graph. CW7 provides K=24/5. A later retained-state replacement is not obtained by inverting a force error.
4. After eliminating only the known affine state skeleton, this graph has bounded known original Gaussian/caller rows, absolute edge operator under the admitted fixed guard, and a state readout whose nonlinear predecessor row has norm O(Lambda a). Equivalently the terminal original-force query argument is that state; its known Gaussian row is a coisometry. Already admitted implicit blocks retain their own finite versions and guards.
5. Actual moving anchors and zero trajectories are original computations with their caller paths. A varying physical center is adjoined as a full label before a gradient embedding. Its port contains none of the newly sampled Gaussians introduced below.
6. The complete VALUE and zero counts, and each original first/adjoint sweep, obey their explicitly expanded same-heat count. No discarded tape is replayed for free.
7. Every retained whole block used as an embedding primitive has its genuine-gradient reference and an available arbitrarily accurate finite VALUE/zero restoration at the actual newly enumerated query widths, with only fixed-target polynomial-logarithmic precision overhead. A finite nongradient approximation is not itself declared a gradient primitive. The old exact-reference/finite-column constructors supply this numerical port for their named blocks; it is not inferred from the old posterior law grade. Freeze the newly refined versions before executing the incumbent and every matched proxy/source occurrence. This refines numerical VALUE error, not an immutable intrinsic mean/prior law floor.

These are supplied by the named old retained certificates for the seed. The extension below preserves them, with the weighted graph interpretation in Section 5. The known-center construction does not call an old hidden service.

## 2. Actual finite quarter program

Let T=pi/2. Evaluate X_0=X_P(a,y;W_old) once. Then sample one new Z independently after the caller is exposed. Set t_i=iT/N_q and

    w_ij=cos(t_i-t_(j+1))-cos(t_i-t_j), j<i,
    w_ij>=0, sum_(j<i) w_ij=1-cos(t_i).

Initialize and execute fixed M_q Picard layers:

    x_i^[0]=y+cos(t_i)(X_0-y)+sqrt(a) sin(t_i) Z,
    x_i^[k+1]=x_i^[0]-a sum_(j<i) w_ij grad V(x_j^[k]). (Q)

The endpoint X_q=x_Nq^[M_q] has the exact affine Gaussian part y+sqrt(a)Z. Its source-dependent remainder is the literal a-weighted sum of the final-layer original gradients. Each gradient is executed and stored once per layer, and all later linear readouts reuse that exact value.

The exact continuous quarter-flow preserves Q_(a,y) and contracts its starting-position error by at most a/(1-a). Its finite product-integration/Picard error has the already audited bound

    Lambda sqrt(Dim)[a^(M_q+3/2)+a^(3/2)/N_q],        (Q-error)

after actual mode/profile and numerical rows. Mode centering is used for this estimate, not as an unexecuted replacement in (Q). A finite computable mode m may also be used in execution if its exact residual d=m-y+a grad V(m) is retained in the harmonic baseline. Because the kernel row sums are exact, the anchored and unanchored formulas then agree as finite VALUE identities. Dropping d without its numerical allowance is not allowed.

Choose M_q at fixed target R and N_q=ceil(Lambda_R a^[-(R-3/2)]), with the actual public-log/error allocation. All grids and versions are fixed before caller differentiation.

## 3. Append a short same-posterior refresh

Sample a new standard L after X_q and every earlier record is fixed. Set h=a^(v/2), u_i=i h/N_s and

    w^s_ij=cos(u_i-u_(j+1))-cos(u_i-u_j),
    w^s_ij>=0, sum_(j<i)w^s_ij=1-cos(u_i)<=C h^2.

Execute fixed M_s layers

    z_i^[0]=y+cos(u_i)(X_q-y)+sqrt(a) sin(u_i)L,
    z_i^[k+1]=z_i^[0]-a sum_(j<i) w^s_ij grad V(z_j^[k]),
    X_R=z_Ns^[M_s].                                  (S)

The exact short flow has initial velocity sqrt(a)L and preserves the same Q_(a,y), without a Metropolis step. The finite map (S) is only approximately invariant; its error is explicitly paid. Its Lipschitz factor in its starting position is bounded by a fixed constant (in fact close to cos h), so it does not remove the quarter-flow accuracy gain.

The short product-integration error is

    Lambda sqrt(Dim) a^(3/2) h^3/N_s.

Its finite Picard tail is bounded by

    Lambda sqrt(Dim) sqrt(a) (a h^2)^(M_s+1).

Choose M_s so the latter fits a^R, and

    N_s=max(1,ceil(Lambda_R a^[-(R-3/2-3v/2)] )).

At v=19/15, the quadrature numerator is a^(17/5), and N_s has exponent at most (R-17/5)_+. It is smaller than the quarter count's R-3/2 exponent.

No old posterior call is made at a^eta. Both (Q) and (S) use the same original potential; every new query is an original gradient at a literal computed point.

## 4. Actual first, state profiles and retained substitution

All derivative statements below come from the finite gradient transcript and bounded original Hessian, not from differentiating a law or VALUE error.

For (Q), differentiating the exact positive row formulas gives bounded intermediate firsts. At its endpoint the direct cos(T) dependence on the old state vanishes. Therefore

    D_y X_q=I+O(Lambda a),
    D_Wold X_q=O(Lambda sqrt(a) a),
    D_Z X_q=sqrt(a)I+O(Lambda sqrt(a) a).

The first assertion holds even from the CW7 actual I+O(sqrt(a)) input. For (S),

    D_y X_R=I+O(Lambda a),
    D_(Wold,Z,L) X_R=sqrt(a)P_new+O(Lambda sqrt(a)a),
    P_new=(0,cos(h)I,sin(h)I),   P_new P_new*=I.

The ordinary centered-path profiles follow by comparing the actual path to its same finite zero trajectory through the positive kernel resolvent. This gives Lambda_p sqrt(Dim) sqrt(a) for state fluctuations and Lambda_p sqrt(Dim) a^(3/2) for the Gaussian-subtracted centered state. It does not use sqrt(the number of deterministic gradient nodes). If an inherited tape is large, use its already supplied source moment profile; do not substitute a radial full-tape estimate.

Replace only X_0 by its actual same-record retained state and keep all new original gradients and both Z,L. The quarter endpoint contracts that state difference by C a; the short map has a bounded subsequent factor. Thus the new retained STATE error is

    Lambda_p sqrt(Dim) sqrt(a) a^(K+1).               (Ret)

The actual new graph is otherwise retained. There is no deletion of a large covariance correction or empirical signal in this deterministic program. All source-zero Gaussian rows and every original anchor remain.

At physical force radius r, apply the terminal original-gradient readout with factor r/sqrt(a). Its physical first is O(r), its caller is O(r/sqrt(a)), and its full recorded square curl is O(r a): the leading block r P_new* Hess(V) P_new is symmetric, while the complete remaining state derivative is O(sqrt(a) a). This includes all old columns.

## 5. Weighted original-gradient graph: no sqrt(N) radius loss

The unweighted stack of N identical-scale Gaussian rows has norm O(sqrt(N)); that stack cannot be inserted into the native signed embedding. Use an explicit deterministic diagonal normalization.

For each new group of N original-gradient nodes, give each node weight mu_i=1/N and set d_i=sqrt(mu_i). If its unscaled normalized primitive is f_i with Lip(f_i)<=1, write

    p_i=d_i f_i(z_i/d_i).

This is still a genuine original-gradient primitive at an explicitly recorded width: if f_i=grad U_i, then its potential is d_i^2 U_i(z_i/d_i). Every VALUE/HVP is still executed at its actual original physical point. A small d_i is never treated as a free source amplitude; its inverse input width is recorded and cancels in the first.

For a physical edge matrix M_phys, the normalized edge is

    A_ij=d_i (M_phys)_ij/d_j,

and the original Gaussian/caller row is multiplied by d_i. Terminal force nodes use weight one. Keep already admitted old blocks in their admitted coordinates.

For equal-size successive groups in (Q), |(M_phys)_ij|<=C a/N and each row and column sum is O(a). The normalization leaves an O(a/N) matrix, so its Euclidean operator norm is O(a). More generally, a map from a group with N_j nodes to a group with N_i nodes has entries bounded by C a/sqrt(N_i N_j), and hence operator norm O(a). The short-step internal block has the additional h^2 factor. The substitution of the quarter endpoint into short queries creates a cross block of the same O(a) form, not a unit edge.

The terminal row is a w_j/d_j. Since max w_j<=C/N and sum w_j<=C,

    sum_j (a w_j/d_j)^2 <= C a^2.

The short terminal part has O(a h^2) norm. The old retained-state nonlinear readout is already O(a); its replication into a normalized node group multiplies it by a vector of norm O(1), so that connection also stays O(a).

The stacked new Gaussian/caller rows obey

    ||C_new|| <= C sqrt(number of fixed Picard groups),

independently of their node counts. In particular the fresh-L columns in the short groups have weighted norm O(h), because every such row is sin(u_i)L and sum_i mu_i is one per group. The L column is exactly zero in every earlier group and every original external label.

The number of groups and the old graph constants depend on the fixed target, so the common small-a guard is chosen after enumeration. The full absolute edge operator is below the native contraction threshold by combining the old admitted block and these O(a) off-diagonal attachments. The signed twin's symmetric edge matrix then has its fixed contraction guard. No inverse node count is hidden in that guard.

This also explains why the many deterministic gradient nodes need not enlarge the next genuine-gradient prior's Gaussian input dimension: they introduce no new Gaussian primitives. Along this known-center induction, each step adds only Z and L (and a fixed number of proxy auxiliaries when actually requested). The seed's existing tape is retained. Row count, gradient-query count and Gaussian tape dimension are distinct ledgers.

### Full callers and actual anchors

For the full-joint graph, write the external physical center as y=y_ref+sqrt(a)U and adjoin U before the embedding. Expand new physical queries directly from (Q)--(S); their dependence on U is affine plus the old original-gradient graph. No moving force anchor is silently declared constant in U.

Actual zero computations are retained as their own original VALUE graphs. In particular, adding B*F_actual(U;0) is not automatically a full-joint gradient merely because it is constant in the private Gaussian tape. A moving anchored primitive is expanded by the native two-node identity, with its actual negative small predecessor edge. Alternatively keep the unanchored fixed-potential original query graph and apply the existing global-origin/flat-gauge construction. The new diagonal normalization scales both an anchor and its caller row by the same d_i. It does not differentiate a small restoration error or assert that a private zero is a full-joint constant.

Finite zero evaluation of the new program may cost its same order of deterministic gradient queries. That cost is charged below; the former polylogarithmic zero specialization is not assumed for this new program.

### Quantitative extension beyond a polylogarithmic node count

The cited concrete LOW30 retained-graph admissions enumerate only polynomial-logarithmically many nodes. We are not relabeling this larger graph as one of those unchanged concrete admissions. The extension used here is the following explicit compiled-source calculation.

Let m_nodes be the actual internal block count, n_W the dimension of its private Gaussian record, n_U the dimension of the explicitly adjoined caller port, and n=n_W+n_U. Suppose the centered original-gradient graph has

    ||C||<=C_R, ||M||<=q<1, ||A_t||<=C_R a.

For a symmetric signed system, its reference p solves p=f(C_tw u+S p), with Lip(f)<=1. Center at the FULL fixed reference origin, or retain the native portable-anchor expansion. Then

    ||p(u)-p(0)|| <= ||C_tw|| |u|/(1-||S||).

This is an operator bound from n inputs to m_nodes*d deterministic internal values. It does not contain sqrt(m_nodes*d). For Gaussian u and the admitted caller profile, the corresponding Lp bound is C_(R,p)(sqrt(n)+|U|). The output source is

    G_*(u)=(r/beta) C_tw* p(u),
    Lip(G_*) <= (r/beta)||C_tw||^2/(1-||S||).

The full joint gradient has output dimension n, not m_nodes*d. Captured caller coordinates remain fixed; they are not resampled merely because they are included in this dimension bound. Writing C_tw=[C_U,C_W], the private readout G_theta=(r/beta) C_W* p is a genuine gradient in W for each fixed caller, and the full joint reference additionally retains the C_U* p component when its port requires it. B is zero on the caller port, so the physical output is unchanged. A mean/pair compiler uses its declared frozen-external-label interface and complete private records; any full-joint adapter retains its actual caller component. Neither version runs a separate Gaussian prior at every internal weighted node. In the known-center induction, n is the seed's actual tape dimension plus two physical d-blocks per refinement and the explicitly requested twin auxiliaries. Fixed-order copies in a later compiler copy that entire n-dimensional source and pay its complete graph evaluation.

The output-nearest physical first is bounded from the exact identity B C_tw*=beta e_t*. The opposite physical column is bounded from C_tw B*=beta e_t. Thus both physical sides are O(r); full caller first is O(r/(beta sqrt(a))) and physical caller first O(r/sqrt(a)). For the physical d-output map, its Hilbert first is at most sqrt(d) times its operator first. For the full n-dimensional gradient source it is at most sqrt(n) times its operator first (and the corresponding private bound uses n_W). Neither is sqrt(m_nodes*d). The marked projection/twin VALUE error has only the actual fresh d-dimensional L and twin Gaussian as its Hilbert factor, together with the supplied old retained-state error.

Here is the original-query cancellation, including errors. With d_i=sqrt(mu_i), a new primitive has the physical form

    g_i(z)=(d_i/sqrt(a)) grad V(y_ref+(sqrt(a)/d_i)z).

A forwarded input row is d_i C_i u, and an edge is d_i M_ij/d_j. The products of d factors telescope along every literal forward or adjoint path. Its physical point is the original point from (Q) or (S), not a point whose Gaussian variation has been divided by d_i without the matching input row. Its directional first is the original HVP at that point, with the stated scaled direction, followed by d_i/sqrt(a). Large intermediate directions are charged in numerical precision; they do not create an extra oracle or additional first direction.

If an original gradient at node i has physical VALUE error xi_i, the scaled block error is d_i xi_i/sqrt(a). A stable full graph solve therefore propagates it to G with bound

    C_R r/(beta sqrt(a)) [sum_i mu_i |xi_i|^2]^(1/2).

For a common per-node tolerance, sum_i mu_i is the fixed number of new groups, plus the supplied old-block allowance. Hence no unpaired inverse minimum mu occurs in this VALUE error bound. Original locations, moving anchors and old blocks retain their own actual input-domain and restoration certificates. A block whose error is not controlled at that actual scaled query cannot be admitted by this cancellation alone.

Finite simultaneous Picard restoration has contraction q independent of m_nodes. After centering, its required iteration count is

    K >= log(C_R (r/beta)(sqrt(n)+caller envelope)/epsilon_G)/log(1/q),

with the separate deterministic-zero allowance. This is polynomial-logarithmic at each fixed target. The executed cost is K times the COMPLETE internal original-query count, including every old node replayed at a changed query. Thus the inverse-heat node count is charged once per iteration; it is not put inside Lambda_R.

The internal small row d_i is not a free new source-heat floor. No low-clock Gaussian integration or derivative restoration is performed in an internal p_i coordinate. The complete-source compiler's genuine Gaussian width is in u. If a proposed implementation instead inserts independent Gaussian records or separate calibrated priors at each internal node, it has changed this construction and must charge the resulting dimension and width factors separately.

Finally, the finite row encoding must use one matching set of positive d_i and their inverses, consistent positive quadrature rows, and an exactly specified final rotation/coisometry (or its explicitly budgeted covariance error). Elementwise arithmetic precision may need an inverse polynomial in m_nodes; its bit length is O(log m_nodes)=O_R(log(1/a)). Dense known matrix actions still cost their actual m_nodes-dependent arithmetic. Nothing here asserts a dimension-free arithmetic/bit bound or erases a performed original query.

## 6. Protected projection on the actual retained graph

Let Pi select only L. The final normalized Gaussian row satisfies

    P_new Pi P_new*=gamma I,   gamma=sin(h)^2.

All earlier queries ignore L. Project L out of the Gaussian rows of the short-step nonterminal original-force queries only; keep the terminal physical force's P_new row intact. Let the resulting auxiliary graph use the same anchors and finite versions.

At each short time the removed direct row has size at most sqrt(a) sin(h)|L|. The force-feedback row mass is a(1-cos h). The finite positive-kernel/resolvent comparison therefore gives

    |X_full-X_prot|
      <= sqrt(a) a(1-cos h)sin h |L|/[1-a(1-cos h)].  (Proj)

For a terminal force readout r/sqrt(a), this is O(r a h^3)|L|. It is zero whenever L=0, with all earlier records left arbitrary. Thus the zero-origin identity is preserved; no Gaussian Lp estimate is evaluated at a deterministic zero to infer it.

The weighted C norm from Section 5 supplies the global O(h) restricted row needed in the signed embedding. Projection does not change the emitted sampler (S); its VALUE price belongs only to the proxy comparison.

After projection the known row R_L=P_new Pi/gamma annihilates every nonterminal query row, is zero on the caller port, and satisfies R_L P_new*=I. The standard finite signed-twin construction applies to this actual projected graph:

    B=beta [R_L,-I/epsilon_tw],  BB*=I,
    B C_tw*=beta e_terminal*.

The squared terminal path bound is O(a^2), so the twin comparison is O(r a^2 epsilon_tw). Choose v=19/15 and epsilon_tw=a^.9. Both projection and twin errors have relative force grade J=2.9. The old retained-state error (Ret) has higher grade already from CW7.

Consequently the finite weak physical proxy has its actual physical row/column O(r), physical caller O(r/sqrt(a)), full first O(r/beta), and full caller O(r/(beta sqrt(a))). At r=a its full shrinking radius is O(a^.1), while its actual firsts come from the finite terminal/resolvent columns. The exact implicit reference is genuinely gradient; the finite Picard proxy is not declared exact-gradient by VALUE convergence.

Choose its finite proxy iteration count only after all original-query widths, caller moments, scalar readouts and Gaussian-row encodings are enumerated. Restore random and actual zero queries separately. The actual-to-retained state and protected projection preserve the original finite origin. For a finite twin, either execute the native same-version origin-preserving centering through the full caller graph, or retain the explicit finite zero offset F_actual(theta;0)-H_fin(theta;0), with its separately bounded origin restoration and caller paths. Do not infer this offset from a Gaussian Lp bound or manufacture a full-joint gradient by an unsupported caller-dependent translation. If a proxy/source version changes, rebuild every matched F/H/E occurrence. Original gradient arguments may have large declared normalized widths under the diagonal weights; those widths are part of the numerical/domain census, not a new uncharged oracle.

## 7. Query recurrence and the known-center induction

The actual VALUE list contains:

    one complete X_P(a,y) call,
    M_q N_q quarter original gradients,
    M_s N_s short original gradients,
    actual original mode/zero/anchor computations as required,
    known coefficient and Gaussian-row operations.

A zero computation or directional sweep has the analogous expanded list. Cached original queries are shared only at the identical point/tape/version. Old adjoints accumulate all outgoing directions before the saved old sweep. Streaming that discards a record pays its actual replay.

Thus, with zero and numerical/anchor calls included,

    Q_R(a) <= Lambda_R [Q_P(a)+a^[-(R-3/2)]].         (Cost)

N_s is lower-order. The finite signed proxy, when requested, has a polynomial-log number of evaluations of its complete weighted graph; each includes all old and new original calls. This changes Lambda_R and its log degree, not the heat exponent. Its genuine-gradient input dimension has no new inverse-a clock bank from this deterministic construction.

Known dense product-integration algebra can cost O(M_q N_q^2+M_s N_s^2), and storage has the actual deterministic node count. This is an original-query/first-action result, not a claim that full arithmetic has exponent R-3/2.

The output is again a known-center source with law grade R, actual g=1 firsts, same-origin retained state grade at least K+1, weighted original-gradient graph, the new fixed weak protected/twin port, and exact counted tape/zeros. This permits iteration from the completed known-center CW7 seed without invoking hidden-center services in the construction.

## 8. Decoder/phase normalization: what does not introduce a new heat power

If this known-center routine is later used as a completed decoder child at eta=a/2, construct its short refresh with its own eta and h_eta=eta^(19/30). Its fresh L contribution to the full hidden state's normalized carrier is

    gamma_hidden=(eta/a) sin(h_eta)^2,

up to the known preceding scalar decoder coefficient, which is one for the last child's fresh state. Hence gamma_hidden is comparable to a^(19/15), with only fixed/logarithmic factors. Earlier statistic/observation/decoder records ignore this new L, and its entering center was already exposed before L was sampled.

The child state projection costs sqrt(eta) eta^(29/10). The host force readout r/sqrt(a) turns this into O(r a^(29/10)) at eta comparable to a. A local child force normalized by r/sqrt(eta) similarly gives r eta^(29/10). Do not multiply by an unpaired inverse sqrt(eta).

For eta=A zeta/(1+zeta), where inverse zeta is an admitted public-log factor, the same relations preserve heat powers. Eta^(-k_R) changes only the fixed-R logarithmic degree relative to A^(-k_R). Every inherited smaller heat beyond such comparability would instead be charged explicitly; none is introduced by the new quarter/short routine.

This proves the protected/caller normalization for a completed hidden host **if its statistic and other source rows have been independently supplied**. It does not construct that inaccessible mean for free.

## 9. The separate inaccessible-center adapter and its bill

For an actual hidden center q=E[H|theta], the selected CW7 host needs only a complete conditional statistic law T approximately N(q,aI/(k+1)), followed by its fresh observation pool and completed known-center children. It does not require an arbitrary fine-record observer theorem. Its actual retained graph is nevertheless rebuilt on the original same tape.

If a statistic with law allowance delta_T, its actual first/caller/origin/retained rows, and work Q_T is supplied, the known-center program above may be used in that decoder. The output law is bounded by

    C delta_T + C Lambda sqrt(Dim) a^R
      + 2^(1-k) sqrt(ad),

with the actual child derivative I+O(a) and the host's geometric first/kept-state transport. Its work is Q_T plus a polynomial-log number of Q_R(a/2) calls and all literal anchors. Its protected last-child row is exactly the one in Section 8.

The current weak-E mean constructor is one available way to produce such a statistic from a supplied complete force/center provider of work Q_H. For target physical allowance a^R, its nominal count is a^[-(R-5.4)]_+, so its bill remains

    Q_hidden,R <= Lambda_R [Q_known,R
                    + a^[-(R-5.4)]_+ Q_H].           (Hidden)

If Q_H is the previous complete hidden family, this still adds an order-dependent exponent to that inherited hidden bill. If a particular application builds only a fixed number of such hidden levels, state that number and its cost explicitly. A depth that grows with target order cannot be put into an order constant when it changes an a exponent.

Thus (Cost) is a standalone known-center/source induction; it is not automatically the old uniform known/hidden-family recurrence c_P. A no-copy mean theorem, a separately efficient hidden reset, or an application that does not need that hidden service is required to promote it accordingly.

## 10. An application requiring only the known-center family

A standard proximal Gibbs chain for the original strongly log-concave target uses only known supplied centers:

    Y_n=X_n+sqrt(a) G_n,
    X_(n+1)=R_R(a,Y_n).

The exact restricted-posterior chain contracts W2 by rho=(1+a alpha)^(-1). Under a uniform or appropriately integrated per-call allowance delta=Lambda sqrt(Dim) a^R, the usual exact-kernel comparison gives

    error_n <= rho^n error_0 + delta/(1-rho).

Thus approximately (a alpha)^(-1) logarithmically many outer steps suffice, with stationary numerical allowance Lambda alpha^(-1) sqrt(Dim) a^(R-1). This comparison uses exact posterior contraction; it does not substitute that contraction for the finite sampler's actual caller derivative.

With beta normalized to one and alpha=1/kappa, the original-query bill is

    Lambda_R (kappa/a) Q_R(a),

plus the stated initialization/mode work. At k_R=R-3/2 and accuracy balancing a^(R-1) approximately epsilon/(kappa sqrt(Dim)), the dimension/accuracy contribution is

    (Dim/epsilon^2)^[(R-1/2)/(2(R-1))],

with the corresponding explicit kappa factor and logarithms. Its exponent tends to 1/2, not zero. This application establishes why a standalone known-center family is meaningful without inaccessible means, but it does not beat the supplied finite-order denominator-32 result or solve the desired sublinear-order cost problem.

Any higher-order outer construction that actually uses hidden conditional force histories must instead use the adapter and work bill in Section 9. Its contract cannot be silently imported into this simpler chain.

## Source pins and verification boundary

- Completed CW7 retained state/force and same-version census in `exact-slack/cw7/`.
- LOW30 `t30:eq:retained-graph`, `t30:eq:buffer-price`, `t30:eq:twin`, `t30:eq:proxy-readout`, and `b27:raw:protected`; source SHA256 7767ac389d6b677eb9517b2585f7a24cb9f85ede47d099a653c635b7be750fb8.
- The quarter-flow law/query audit at `../independent-audit/`.
- The actual host separator at `no-copy-rank-20261004/host-observer-audit/`.
- The short-refresh audit at `independent-audit/INDEPENDENT-PROTECTED-SHORT-REFRESH-AUDIT.md`.

No law proof here identifies an analytical Gaussian coupling with the actual retained source tape. No smaller-heat old call, discarded original query, HVP derivative, covariance oracle, or hidden mean is free.
