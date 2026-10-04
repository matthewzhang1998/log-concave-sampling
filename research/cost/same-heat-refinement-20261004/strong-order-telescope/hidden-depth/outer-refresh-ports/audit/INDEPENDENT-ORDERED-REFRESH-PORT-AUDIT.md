# Independent audit of the ordered short-refresh outer ports

2026-10-04. This is a new audit and does not modify any earlier pinned theorem.

## Decision and exact boundary

**PASS_SOURCE_AND_NATIVE_PROXY_PORTS; FULL_COVARIANCE_GRADE_JOIN_PENDING_INTRINSIC_DIMENSION.** The ordered-refresh source-port construction passes, subject to the explicit numerical, caller-profile and complete-source inputs below. The finite ordered list repairs the previously missing low-J pure and higher-J covariance ports while retaining the strong-allocation source query exponent. Projecting only the chosen short stage does **not** supply a protected row; every later nonterminal original-gradient row must also be projected. The signed coisometry must not replace the recorded physical carrier when computing curl.

This audit independently derives the mechanism from the pinned source contracts. It checks the concrete list

    J=(2.4,2.7,2.8,2.9,3,3.1,3.2,3.3,3.4,3.5,3.55).

It establishes the native source, protected-width, finite-derivative and strict-radius gates needed by the stated exact32 pure-fourth and covariance/mixed ports. It does **not yet establish the unchanged complete covariance grades**: the ambient factor in the marked parent's capped error branch remains an extra obligation, distinct from arbitrary-order gradient-child restoration. Section 7.1 records the gap. It does not assert that replacing exact32's old polylogarithmic source by this positive inverse-heat-cost source preserves exact32's total sampling complexity. It also does not extend the result to unspecified opaque histories, arbitrary callers without integrated profiles, heat-dependent hidden depth, or uniform growing-order constants.

The qualified candidate `ORDERED-REFRESH-LITERAL-OUTER-PORTS.md` was inspected in full and its final source/proxy version was verified at SHA256 `c3174b91086a819ba047435cdf6d828468e30f1a4e45b11a139b8133d53d02b1`. Its Section 8.1 now correctly withholds the complete covariance-grade join. The present decision applies to that exact qualified version. A future intrinsic-cap extension requires its own audit and is not assumed here.

## 1. Independent source inspection

The following SHA256 values were checked directly:

- Strong hidden continuation: `83267ae699003d637610f796150ec6d9a526551e8ff0e0bd026f92816e977d33`.
- Complete strong/law/zero audit: `b1943d4774f7cf77cbd9d11c0e3baac7ffdb83bc269992fb7fda14a100980fa6`.
- Protected quarter/short source: `44cfd4b2d5ef98b86df8b824c1b98217281547eb7cf19c3d3ff05f095d816036`.
- Exact32-v2: `c7487ca24ffa44ea73b52e42e4b6f715ad964171fa5bd0a61cf9323eb646cbff`.
- LOW30: `7767ac389d6b677eb9517b2585f7a24cb9f85ede47d099a653c635b7be750fb8`.
- LOW31: `3de62d349338304cce43af4e49c8862e8c260b1f5816234d63d6abfbc1d004ba`.

Directly checked LOW30's retained graph, signed twin and finite terminal-row formulas near lines 784–905; its fixed-order original-gradient interface near 5703; and the separate high-moment child restoration near 635. Checked LOW31's linear-growth and conditional-smoothing arguments near 2670–2780. Checked exact32's finite protected principal block near 650, matched widths near 683, covariance chain and strict gates near 1403–1467, and pure-fourth body near 1685–1695.

## 2. One actual sampler and separate auxiliary ports

Replace a known refinement map by its quarter map followed by the fixed ordered list of short maps. At heat eta, let

    h_i=eta^((J_i-1)/3).

The actual sampler executes every short map and uses fresh independent L_i, exposed after its entering state and all earlier records. This is one finite source shared by all ports. Port i changes only an auxiliary retained/proxy graph, not the source's law or executed sampler.

The final normalized affine row is exactly

    P=(product_l cos h_l) Z
        +sum_i (sin h_i product_(l>i) cos h_l) L_i.

Let s_i=sin h_i product_(l>i) cos h_l. Successive exact rotations show PP*=I, and gamma_i=s_i^2 is comparable to eta^(2(J_i-1)/3). Use one consistent encoded rotation/coisometry, or explicitly budget its covariance error. Independent unrelated rounded sine/cosine values do not prove exact coisometry.

For a hidden source select the short bank in the **last refinement of the last decoder child**. The child runs at eta=A/2. All statistics, observation pools, prior decoder states, old history banks, and earlier refinements then ignore the selected L_i. Its host-normalized coefficient is sqrt(eta/A) s_i, so the protected variance remains comparable to A^(2(J_i-1)/3). The fixed eta/A factor changes no power.

Each port has a separate projected graph and its own twin variable. All these variables can be padded into one complete record; an unused auxiliary coordinate contributes zero to the original source. The proxy differences in the covariance chain must use this one common padded record and the same immutable finite source/proxy versions. Independent complete records are reserved for different covariance packets, not substituted inside a same-tape difference.

## 3. All later force paths: the decisive projection check

Eliminate the known affine state skeleton before projecting. Each short-stage output is an affine carrier plus a sum of original-gradient values with physical coefficients. A later short stage contains the earlier carrier, so its nonterminal gradient arguments have a direct L_i row. Merely noting that L_i was fresh at its own stage does not remove those later reads.

For port i remove the L_i column from **all nonterminal gradient rows in stages s>=i**. Retain the final physical terminal gradient's affine row P. Keep every original force edge, including paths through earlier stage outputs, and keep all actual anchors and zeros. In an implementation that represents stage endpoints explicitly, do not accidentally remove their affine L_i carrier from the final physical output: only the nonterminal force query rows are modified.

For s>=i, the normalized direct L_i row is

    sin(u)                             when s=i,
    cos(u) sin h_i product_(i<l<s)cos h_l   when s>i.

It is at most C h_i. Earlier rows are exactly zero. The total final physical predecessor mass on force groups of stage s is at most C A h_s^2, times bounded downstream cosine/causal resolvents. This statement follows from the literal positive quadrature rows; it does not count one unit edge for a whole predecessor state.

The weighted absolute-edge resolvent therefore gives, pointwise in all other records,

    |X-X_prot,i|
       <= C sqrt(A) A h_i sum_(s>=i)h_s^2 |L_i|,
    |F-F_prot,i|
       <= C r A h_i sum_(s>=i)h_s^2 |L_i|.

At a child heat eta comparable to A the same estimate has only fixed/logarithmic conversion factors. In particular no unmatched A^(-1/2) appears.

Because J_s is increasing, h_s<=h_i for s>=i when A<1. The finite suffix length is fixed. Consequently

    A h_i sum_(s>=i)h_s^2 <= C A h_i^3 = C A^J_i.

This also proves every fixed Gaussian Lp estimate with a single sqrt(Dim) factor, supplied by the actual d-dimensional L_i. It is not a sqrt(total tape dimension) projection estimate. At L_i=0 the two finite graphs agree exactly with all other private records arbitrary. Hence the original all-zero tape is preserved as an algebraic identity.

Ordering is necessary for the asserted individual marks. If a J=3.55 stage precedes a J=2.4 stage, its later-path exponent becomes

    1+(3.55-1)/3+2(2.4-1)/3 < 3.55.

A stage-only projection is a different failure: even in the correct order its later rows do not vanish, so it fails the exact signed-annihilator identity. It may still have a small VALUE price, but it has not created the required protected port.

## 4. Weighted graph, twin and finite protected principal block

The original weighted normalization extends stage by stage. Give each fixed Picard/grid group of N nodes weights d_j=N^(-1/2). Then its original Gaussian rows have bounded stacked norm; an s-stage L_i column has weighted norm O(h_i). Entries between group sizes N_1,N_2 are bounded by C A/sqrt(N_1 N_2), with h_s^2 for its short-force predecessor mass. Earlier nonlinear state readouts keep their existing O(A) physical factor. The number of groups and stages is fixed at the requested order.

Thus the full normalized retained graph has bounded C, a small absolute edge matrix, and terminal predecessor row A_t=O(A). Small node weights are not free amplitudes: their inverse query widths must be kept. Products of weights cancel along actual forward/adjoint paths, and all new VALUE tolerances are chosen after this census.

After projection R_i=P Pi_i/gamma_i annihilates every nonterminal row and satisfies R_i P*=I. It is zero on the full adjoined caller port. Choose

    epsilon_i=A^(J_i-2),
    beta_i=(gamma_i^(-1)+epsilon_i^(-2))^(-1/2),
    chi_i=max((J_i-1)/3,J_i-2).

The suffix cosine factors only change constants, so beta_i is comparable to A^chi_i. LOW30's literal signed twin gives

    B_i=beta_i[R_i,-I/epsilon_i],  B_i B_i*=I,
    C_tw,i B_i*=beta_i e_t.

The terminal squared path is O(A^2), hence twin VALUE error is O(r A^2 epsilon_i)=O(r A^J_i). Together with projection, retention and separate finite numerical restoration, this supplies the required port-i weak error.

For finite simultaneous Picard iteration p^[k+1]=f(C_tw w+S p^[k]), p^[0]=0, its actual first satisfies

    Dp^[k+1]=H_(k+1)(C_tw+S Dp^[k]).

Only original Hessian actions at the actual gradient nodes occur. With ||H||<=1 and ||S||<=q<1,

    ||Dp^[k] B_i*||<=beta_i/(1-q),
    ||B_i DG_k B_i*||<=r beta_i/(1-q),
    G_k=(r/beta_i) C_tw* p^[k].

These hold at every executed finite iterate. They do not rely on convergence of derivatives. With P_i=B_i*B_i, Q_i=I-P_i and T_i=P_i+eta_i Q_i, the four blocks yield

    (alpha/eta_i)||T_i DG_k T_i||
      <=C alpha r (beta_i/eta_i+2+eta_i/beta_i).

This is precisely the structural gate used in exact32. Each consumer must use its child's own beta_i, rather than the smallest width in the full port list.

### Recorded physical row versus signed coisometry

The recorded physical carrier remains P_rec, padded by zeros in the twin coordinates. It is generally different from every B_i. The physical source is F_i=B_iG_k=r p_t^+[k], but its recorded curl is computed using P_rec.

Because the terminal column of the original edge matrix vanishes, the finite terminal first is exactly

    Dp_t^+[k]=H_t^[k]([P_rec,0]+A_t Dp_N^+[k-1]).

The leading square lift r P_rec* H_t P_rec is symmetric. The rest is O(r A), giving the full recorded curl O(r A), including old-label and twin off-diagonal blocks. This proves the physical row/curl independently of weak VALUE proximity. The signed identity supplies the opposite physical column and protected principal block. Swapping the two coisometries would invalidate the proof.

Moving caller-dependent zero offsets must be executed and differentiated through their actual original graphs. They are private-constant, but not necessarily full-joint constants. B_i is zero on caller coordinates, so a legitimate private zero restoration does not spoil the private principal block; its full caller path still needs its own stated certificate.

## 5. Source law, actual pair, zero and hidden induction are rebuilt

The new finite list is not appended only to an analytical proxy. It is executed in every known source occurrence, including all fine/coarse levels, base sources, decoder children, actual zeros and changed-query restorations.

A single quarter map contracts its starting-state dependence by C A. The finite short suffix has bounded starting-state first, at most product_i(1+C A h_i^2), and preserves that quarter contraction. The exact short maps each preserve the same Q_(A,y); only their finite discretization errors require payment. A short map's product-integration numerator is

    A^(3/2) h_i^3=A^(J_i+1/2).

For local state target A^R choose its own grid count with exponent

    (R-J_i-1/2)_+ <= (R-3/2)_+.

Choose the fixed short Picard counts to pay their sqrt(A)(A h_i^2)^(M_i+1) tails. Allocate the local allowance across the finite list before freezing versions. The quarter grid continues to dominate original-query exponent. Dense arithmetic and storage keep their actual larger bills.

At adjacent levels, share one complete old seed and the literal suffix of complete momentum banks, where a bank contains its quarter Z and all eleven L_i. Each marginal retains its own local grid counts and finite versions. The first unmatched map gives an O(sqrt(A)) state discrepancy; each shared later quarter/short map contracts it by C A. Two different finite grids contribute their separately bounded discrepancies from the same exact composed map. The unattenuated final coarse-grid error is retained. This reproduces the original known-child state pair bound sqrt(A)A^(j-1) and canonical force bound A^j. The same literal deterministic recursion gives the known-child zero bound; no Gaussian estimate is evaluated at zero.

Use the corrected convention K_0=K_1 for the integer ladder. Both now mean the same one-quarter/eleven-short finite program. It has exponent-zero cost, actual g=1 firsts and retained state mark at least 29/5. The quarter still removes the direct old seed row, and finite first propagation gives

    D_y K=I+O(A),  D_W K=sqrt(A)P_rec+O(A^(3/2)).

Nothing about the extra stages changes the complete empirical identities

    E M_j=E F_j,
    M_j(theta;0)=F_j(theta;0).

Each difference sample is a full named paired program. Independent complete banks are independent only conditional on the original theta. The source's law target q_(d,j)=h_base+sqrt(A)E F_(d-1,j) is unchanged in definition, though the named finite marginals are the newly rebuilt ones. Fine/coarse signal errors, the shared statistic/observation coupling, the actual two-caller decoder recursion, and the physical A history-bias factor therefore reproduce the strong pair and law recurrences of the pinned theorem. The exact-zero telescoping identity reproduces the deterministic-zero recurrence.

The once-refreshed external phase Gaussian stays captured, rather than being redrawn in empirical copies. The selected LOW30 Picard history still uses independent **full** lower-depth replays and physical original-gradient predecessor O(A), so it supplies the required causal-history port. An arbitrary legacy history would still have to supply that port separately.

The finite mode schedule and all source numerics remain relative to the original caller/domain moment profiles. A fixed finite prox has a caller-weighted residual; the extra refresh list does not make that residual uniformly small over every unbounded caller.

Every lower complete bank remains charged. Thus at integer level j the work/dimension exponent is still 2(j-1), while known-child grids have exponent at most j-1. Whole-copy weights N^(-1/2) preserve the O(A) physical aggregate row. Fixed hidden depth changes constants and logarithmic degrees, not this exponent.

## 6. Literal exterior gates that are now available

For pure fourth, J=12/5 gives chi=7/15. At canonical r=A^(3/4), its full radius has exponent

    3/4-7/15=17/60,

and the tuple-weighted body exponent is

    3/4-7/15-1/4=1/30>0.

This is the previously unavailable low-J port. It is a property of the actual multi-refresh source's finite proxy, not a regrading of its J=2.9 proxy.

If the pure queue's boundary is bounded conservatively using sqrt(n_amb) rather than a separately proved physical sqrt(Dim) Hilbert factor, its force-K allowance becomes sqrt(Dim) A^[1+K/30-c/2]. Choose a fixed K so 1+K/30-c/2>=9+strict_slack, then choose its native prior order after the expanded finite census. Keeping K=600 for every source order is not justified. This is a repairable finite-boundary loss because K is freely increasable; it does not repair the marked parent cap discussed below.

The remaining ten marks are exactly exact32's covariance chain. For transitions (x,y), using its unchanged u=(327/50-x)/4 and e=y-7/2+u+1/100, every child gate

    3/2-u+chi_child-e,
    3/2-u-chi_child+e,
    3/2-u

is positive for both adjacent children. The newer child's second gate equals 1/100. The final x=71/20,u=3/4,e=81/100 has the same positive minimum. The principal-block proof above supplies the finite derivative structure these inequalities require.

The source supplies the differences on one immutable padded tape, physical first/caller O(r),O(r/sqrt(A)), recorded curl O(r A), and the requested fixed-moment weak energies. Its retained mark 29/5 exceeds every J in the list. These are source ports. The outer coefficients, marked priors, matrix guards, common endpoints, owned-root restrictions, counterterm chronology and actual-amplitude scaling still come from their respective exterior theorems and must not be inferred from the source law alone.

## 7. Ambient dimension and high-moment child restoration

The old LOW31 argument used only a public-logarithmic number of physical-D Gaussian blocks. Here the complete source tape may have

    n<=C D A^(-c),  c=2(R-3/2) at the requested source level,

up to public logs. Substituting n=D into that argument is invalid. The larger dimension is harmless for a fixed-order source exponent only after it has been explicitly paid in prior and moment allowances.

First enumerate the finite moment list and all original query/caller paths. The exact genuine-gradient reference is globally Lipschitz in its genuine Gaussian input, and the admitted anchored compiler growth is linear in its input norms. Thus its actual finite children and parent readouts have bounds

    M_(2p)<=C_p sqrt(n) A^(-b_p)

for fixed finite b_p, allowing all physical widths, caller envelopes and numerical floors. This is a conditional/integrated profile assertion at the admitted callers, not a uniform theorem for an arbitrary enormous moving anchor. The completed graph/restoration input must supply it. Finite node count alone does not create independent Gaussian dimension; complete empirical copies do.

The genuine-gradient prior at conservative radius A^delta, delta>0, and fixed order B gives

    e<=C sqrt(n) A^(delta B).

It is applied to the exact reference; finite VALUE/zero occurrences are separately restored at every consuming width. No finite nongradient Picard output is declared exact-gradient because its VALUE is close.

For a later rough parent use fixed Gaussian damping and LOW31's conditional TV argument, retaining every outside label. Its joint mismatch probability is at most C e. The same maximal coupling and Holder inequality give for every requested p>=2

    ||Y-Y_ideal||_p <= M_(2p) (C e)^(1/(2p)).

On D<=A^(-dmax), this is bounded by

    C_p sqrt(D) A^[delta B/(2p)-b_p-c/2-(c+dmax)/(4p)].        (H)

For p=2, the exponent is delta B/4-b_2-5c/8-dmax/8. The more conservative penalty 3c/4 in place of 5c/8 is also valid. Choose one fixed B after c, delta, the finite list b_p, p*, numerical floor amplification and the desired outer allowances. This pays the larger tape and all moments. In particular the old number 64000 must not be preserved without checking these quantities.

If the parent is actually Lipschitz, one may instead interpolate a retained-label L2 coupling with its actual 2p moments. W2 accuracy by itself is not an Lp coupling claim at the same exponent. Formula (H) gives a direct route even for the rough kernels that motivated the LOW31 smoothing lemma.

There is no moment-order cycle: the finite outer moment list is fixed first; LOW31 Proposition `cpl:prop:growth` explicitly excludes an inverse-source-radius power and keeps the random-input degree one at every fixed B, with constants/logarithmic powers depending on B. The exponent c is that of the underlying complete source, and increasing B inserts a fixed-order compiler with public-logarithmic source calls rather than changing that source's inverse-heat exponent. If an implementation inserts an extra prior recursively inside every empirical bank, this cost argument no longer describes it.

## 7.1. The capped marked-parent dimension gate remains open

Exact32's cap-lift theorem explicitly has an ambient sqrt(n) factor and the fixed cap 21/2, with working cap 52/5 in the covariance ledger. LOW30's marked source is a square map R^m to R^m; its displayed local rows include unmarked sqrt(m)rho^16, sqrt(m)(r kappa)^8 and fixed-force queue boundaries. The new source's square lift is D=P_rec* d on the complete ambient tape. Its centered VALUE energy has only sqrt(Dim), but that alone is not a theorem replacing sqrt(n_amb) in every capped branch.

For the final mixed call, u=3/4 and B_work=52/5. The cap allowance is

    sqrt(n_amb) A^(B_work u)
        <=C sqrt(Dim) A^(7.8-c/2).

As a rank-four row its conservative interior grade becomes 16.6-c. Already c=2 drops below the needed 16. Increasing the arbitrary-order gradient child's B in Section 7 cannot change this substantive capped term. The smaller working-cap slack also cannot absorb an arbitrary c/2. This is an admission gap in the theorem-as-stated, not a numerical counterexample to the proposed source or its exact protected identities.

A plausible repair is a separate intrinsic-dimension marked-cap lemma. For D=P_rec* d, its range lies in a known fixed Dim-dimensional subspace. Its first has Hilbert norm at most sqrt(Dim)rho, and its curl has rank at most 2Dim with Hilbert norm at most sqrt(2Dim)rho^(1+k). Scalar-affine/sign/OU source descendants preserve the original output subspace. LOW30's basic rectangular D2/E2 helpers can be executed on a genuinely fixed physical output subspace, then lifted with an exact independent complementary Gaussian. These facts are useful inputs, but a proof must propagate them through every seed, frozen-curl helper, four-force word and one-curl repair, nonlinear side replacement, and capped queue branch. Merely observing that the final output is physical, or that a formula mentions one Hilbert factor, does not complete that census.

Until this lemma or a different fully priced repair is proved, the statement that all 269 covariance categories retain minima (1607/100,141/10,22603/1250) is **pending**. The passed source/pair/law/zero/work results, the eleven native proxy contracts, pure weighted body 1/30, and all strict anisotropic child gates are unaffected.

## 8. Genuine references and numerical ports are substantive inputs

For every port, enumerate the weighted graph, caller components, exact source-zero paths, signed twin, finite proxy, all differences and all consuming widths before fixing numerical versions. Obtain the genuine-gradient reference by the symmetric implicit system. Restore its finite VALUE and actual zero computations with the original contraction mechanism. The contraction guard is independent of the number of equally weighted quadrature nodes; a complete solve still evaluates every original query at every restored iteration.

The original scaled primitive d_i grad V(y_ref+sqrt(A)z/d_i)/sqrt(A) uses the original physical query. Its Jacobian is an original Hessian there, and its first/adjoint action is one original HVP with the recorded scaled direction. No HVP output is used as a new differentiable VALUE producer, and no third derivative is required. Small d_i may enlarge numerical directions and bit precision; their logarithms are fixed multiples of log(1/A) at fixed order/depth, while dense arithmetic is charged separately.

A common finite source version cannot be changed only inside one side of a marked difference. Refine and rebuild all matched actual/proxy/difference occurrences together. The numerical allowances for random inputs and actual zeros are separate. Contraction of the VALUE solve does not authorize differentiating a numerical error bound.

## 9. Executable diagnostics and limitations

`check_outer_refresh_ports.py` produces `outer_refresh_port_checks.json`. It reports PASS with 851 assertions:

- 33 finite eleven-refresh graph/projection cases at three heats, using explicit original-gradient nodes and their inherited cross-stage edges.
- Exact final coisometry, exact full-suffix annihilator, and explicit failure of a stage-only projection in every nonfinal stage.
- Pointwise projection and exact protected-coordinate-zero equality.
- 33 signed systems, with six literal Picard VALUE/first iterates each, checking the protected principal block, recorded physical row, and recorded curl separately.
- Rational projection/twin exponents, ordering necessity, the pure-fourth 1/30 body, all transition and final child gates and fixed-order ambient moment allowances.

The test potential has continuous bounded Hessian .5+.2 sqrt(|x|)/(1+sqrt(|x|)) and is C2 but not C3 at zero. The code implements its gradient and Hessian only. The largest observed projection ratio was 0.129761 and principal-block ratio 0.513699. They are diagnostic values, not theorem constants.

These diagnostics use a simple Gaussian starting-state fixture and small finite grids. They do not implement CW7 or numerically establish the complete hidden law. The analytic source reconstruction in Section 5 and the explicitly imported complete CW7/restoration interfaces are necessary parts of the conclusion.


### Complete hidden fixture cross-check

The companion author's `check_ordered_refresh_join.py` was inspected and independently rerun. It reuses the previously audited complete affine-Gaussian source definitions, replacing the entire known refinement and record allocator with the eleven-short bank. It passed 471 assertions on nine fully expanded source/pair fixtures, including hidden depth two, exact named pair marginals, deterministic level-one zeros, finite and ideal law comparisons, and the cost/gate algebra. The largest force-gap and ideal-law ratios were respectively 0.332461 and 0.206020. As with the original fixture, its seed is an exact Gaussian test fixture rather than CW7. This is an additional implementation cross-check, not a replacement for the analytic complete-source induction.

Publication scope: ancillary off-grid-ladder commentary and workflow-only narration are omitted. This copy retains the audited integer-ladder result; the original and public hashes are separately recorded.
