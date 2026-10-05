# Shared carrier and positive cubature interface

5 October 2026. Read-only comparison of the sealed carrier-replacement and positive-clock-compression packages. The peer statement is a research lead, not an imported theorem: only its supplied Slack excerpt was available for this comparison.

## Main conclusion

A carrier-preserving packet theorem could remove the remaining polynomial-in-dimension guard for the original common-carrier execution. It must retain the **actual physical carrier**, support the exact packet kernels below, and control an integrated conditional transport cost with that carrier held identical. An ordinary joint Wasserstein bound that permits the retained carrier to move does not, by itself, give this composition theorem. The existing positive cubature already has logarithmic accuracy dependence at fixed original heat; a positive OU covariance rule alone does not establish the mixed two-clock contract or remove its inverse-width bill.

## 1 The exact exposed Gaussian and private banks

Freeze a descriptor record E before drawing fresh dynamic Gaussians. E contains the complete original coefficient bank B, required original external callers Y, frozen history/clock/hit/readout labels, actual source versions and numerical floors, and any old randomness needed to specify the unchanged old output W_old. Conditional statements require the specified uniform envelopes, or an explicitly integrable envelope over E. The old output reads no new dynamic tape.

For packet h, the executed source-zero map and common-carrier construction are

    X_h^a = L_h P_h + u_h^a(E,P_h,U_h),     u_h^a(root=0)=0 pointwise,
    L_h L_h^T = v I_D,
    H ~ N(0,v I_D),                       G = H/sqrt(v),
    Pi_h = I - L_h^T L_h/v,
    P_h = L_h^T G/sqrt(v) + Pi_h V_h.

Here G, every V_h, and all U_h are mutually independent standard Gaussian banks conditional on E. Set Z_h=(V_h,U_h), including the actual private/cut/filter/reset tapes in U_h. Packet h depends only on E,G,Z_h; there is no other cross-packet tape sharing. P_h is standard Gaussian marginally and L_h P_h=H exactly. Its conditional law given G is generally singular and noncentral. The maps L_h and Pi_h are known and independent of E's dynamic caller coordinates; otherwise their extra derivative paths must be charged.

The original full-v group is

    W_old(E) + H + sum_h u_h^a(E,H,Z_h) + sqrt(kappa) G_keep,

with one independent, untouched keep. H and all Z_h are owned internally by this group. Its final external output does not expose H, P_h, individual packet outputs, query centers, or the keep. Internal exact-key replay storage does not make a tape an observer. Retaining E for this conditional proof does not authorize appending the whole bank B to the final output; the original admitted observer ledger remains binding.

**Map to the peer notation.** The least demanding candidate is peer P=G and peer P_i=G for every packet, so every projection is the identity and overlaps are complete. This is allowed if the peer result covers that case. Our symbol P_h is a different object: it also contains private Pi_h V_h and is not a deterministic projection of G. Taking peer P to include every V_h would demand a stronger, unnecessary internal-tape-retaining theorem.

These geometric hypotheses are present. The missing fact is that the actual native packet replacement remains valid with G retained. Freshness relative to E does not make an internal carrier an already admitted external observer. In particular, active-source response identities integrate their physical probes and do not automatically survive exposing them.

Sources: COMMON-CARRIER-AMALGAMATION §§1–4; FINITE-DIMENSION-CARRIER-HYBRID §§1,4,6; ACTIVE-PROBE-GENERATOR §§1–4.

## 2 Kernel and metric required from a carrier replacement

Define actual and reference residual kernels

    K_h^a(E,g,du) = Law(u_h^a | E,G=g),
    K_h^p(E,g,du) = Law(u_h^p | E,G=g).

The reference u_h^p is the original finite polynomial used by the current theorem, with its exact root-linear current and complete derivative frames. For the general projection form, a packet kernel may read P only through its declared P_i=A_i^T P, plus E. Descriptors cannot secretly depend on undeclared parts of P. Private packet kernels must factor conditionally on E,P on both sides.

A sufficient local bound is the **fiber-preserving** cost

    epsilon_h(E)^2 >= integral W2^2(K_h^a(E,g),K_h^p(E,g)) dgamma_D(g).

It permits integrated rather than pointwise-in-g error. Couple E and G identically, choose measurable conditional packet couplings, and take their conditional product. Minkowski then gives the residual-sum group error at most sum_h epsilon_h(E). After integrating E, the bound is at most the sum of the L2(E) norms of epsilon_h(E). Identity-overlapping projections cause no problem.

For a more general complete caller C(E,P,y_1,...,y_J), the corresponding statement requires its actual coordinate Lipschitz/path constants: if these are deterministic c_h, the bound is sum_h c_h epsilon_h. Arbitrary amplification by a caller is not free. Our displayed group caller has c_h=1 and keeps the same single H term.

An ordinary W2 bound between joint laws (P_i,Y_i) may move P_i. It is weaker than the fiber-preserving cost above and does not supply independently glueable conditional couplings on one unchanged P. A replacement theorem in a different metric could still work, but it must provide the precise complete-caller stability statement and its constants. “Integrated joint” alone does not identify which contract is meant.

**Reserve boundary.** The desired kernel is before the group's final untouched keep. A bound only after adding a private per-packet reserve cannot be deconvolved, and that reserve cannot be appended as a retained variable without a separate theorem. An alternative valid statement could instead compare the entire grouped output with its one shared final keep. Please specify the exact reserve/readout boundary.

**What is already proved here.** The sealed native input is an ordinary pre-keep marginal floor e_h, retaining E but not H. Proof-only polynomial taming gives

    W2(actual group, original polynomial group)
      <= (1+L S) E0 + (2+2L S) delta S sqrt(D) + 2L S^2 sqrt(D),

where S=sum_h s_h and E0=sum_h e_h. The admitted triple-cubic choice has S=B_amp alpha^(17/2) and a factor L containing D^((m-1)/2) for the actual finite physical polynomial degree m. A dimension-sharp carrier-preserving or complete-environment theorem could replace this loss; the peer excerpt is not yet evidence that its native packet hypothesis holds for these kernels.

Sources: FINITE-DIMENSION-CARRIER-HYBRID §§1–6; COMMON-CARRIER-AMALGAMATION §4, especially its carrier-retaining condition (13).

## 3 Exact tensor integrands and Gaussian node realization

The coefficient target is the full ordered bank cumulant kappa_B(L^(1),L^(2),L^(3)). For triple cubic each L is the entire finite old cubic sum. For mixed 3/3/6, the third L is the recorded finite symmetric-OU covariance six-packet, including its old bridge, original heat, all old hit pairs, and either the explicitly selected diagonal or the full cross-history census. These are distinct targets.

For each labelled replica tree T on three arguments, with two parameters t_e, define

    K_ii(t)=1,     K_ij(t)=min_{e on path i--j} t_e,
    D_ij = sum_a partial_(B_i,a) partial_(B_j,a),
    F_T(t;y) = E_K [ (product_{ij in T} D_ij)
                    (L^(1)(B_1;y) tensor L^(2)(B_2;y) tensor L^(3)(B_3;y)) ].

The Gaussian expectation uses covariance K(t) tensor I on the COMPLETE bank. Products/contractions preserve the recorded ordered physical slots, source normalizations, and any final physical permutations. The coefficient is sum_T integral_[0,1]^2 F_T(t;y) dt, with the one 1/3! factor only if the chosen log-generating convention requires it. Alpha^9 or alpha^12 may be factored out once for the normalized tolerance.

Every D_ij is expanded over the actual source occurrences and their original affine injection rows. It is not a derivative of a numerical covariance root. The new hit counts are 243 for a fixed ordered triple-cubic history triple and 648 for a fixed mixed history, across all three trees. Mixed old nine-pair and cross-history factors remain separate.

At every retained scalar node, keep the actual original B and realize

    c=min_e t_e,             delta_i=1-max_{e incident i} t_e,
    R_bank R_bank^T=K-c 11^T-diag(delta_i),
    B_i=sqrt(c)B+(R_bank V)_i+sqrt(delta_i)E_i.

New complete banks V,E_i are independent of B and the old physical publics. Within each cubic, extract only existing innovation covariance using the simultaneous deterministic m_star allocation D_i=(m_star/32)I. New private widths obey t_iv^2=old_sigma_iv^2+delta_i(D_i)_vv. The full mixed six-packet extracts delta_6 s_old m_star/64 only from its recorded NEW innovation residual blocks; old private heat and common-bank rows stay intact. This is exact redistribution, not added target heat.

The quadrature approximates the **bank-mean** coefficient uniformly in fixed external callers and its required caller derivatives. It does not assert a coefficient identity pointwise in retained B, nor a bank-conditioned native LAW identity. All complete firsts use the actual rows, strict ancestors and separately certified numerical Sobolev floors.

Sources: THREE-CUBIC-EXACT-HEAT §§1–7; MIXED-THREE-THREE-SIX §§1–5; POSITIVE-SECTOR-GAUSS §§1–6.

## 4 Positive scalar quadrature and executable callback

On each of the two ordering sectors, write t_large=r, t_small=r s. The integrand is G_T(r,s)=r F_T(r,r s). Its pure-coordinate derivatives have a certified factorial bound at the original floor eta. Use Price majorant p_Price>=sum_(i<j) N_i N_j a_i a_j, at least one; safe unit-row choices are 27 for triple cubic and 45 for mixed. B_star includes the full family's main and required caller envelopes, actual scalar factors and ancestor rows. With ||.||_* denoting any proper cut or HS/sqrt(D), the certificate is ||partial_r^ell G_T||_*, ||partial_s^ell G_T||_* <= 2 B_star ell! (16 p_Price/eta)^ell. One full-family main envelope is B_K=W 2^(N+3K/2) sqrt(K!) eta^(-K/2); add the K+1 envelope with its actual distinct caller-weight sum to obtain B_star>=1+B_K+B_(K+1). Execution uses rational upper bounds.

The sealed rule has

    j=max(1,ceil(log2(32 B_star/epsilon))),
    m=max(1,ceil(log4(32 B_star/epsilon))),
    h=eta/(32 p_Price),
    Q_tree <= 2 m^2 [ceil(32 p_Price/eta)+j]^2.

It truncates each transformed coordinate at 1-2^(-j), paying the omitted strip error. Dyadic gap panels are subdivided to length at most h; each gets positive m-point Gauss-Legendre weights. Actual sector weight is w=lambda_r lambda_s r. Rational root isolation, positive interval weights, symmetry and exact mass/first-moment normalization certify finite scalar execution.

The code interface is `prepare_certified_rule(N,K,P,eta,epsilon,W,caller_weight)`, returning `(params,base,panels)`, then `sector_nodes(base,panels)`, yielding `(t_large,t_small,w)` for **one** sector. Run both assignments of large/small to the tree's two edges. The weight already contains r; do not multiply the Jacobian again. N,K are 9,7 or 12,10. Code parameter P is the scalar Price majorant, not the peer's exposed Gaussian. The default caller_weight=N W is valid only under its unit-row assumption; supply the full actual caller envelope otherwise.

The callback is mathematical, not an implemented tensor-evaluation API: for each labelled tree, old history, hit assignment, sector node, and frozen source version, instantiate the exact Gaussian program above and compile its nine or twelve active VALUE sources. The scalar builder never evaluates F_T, estimates its expectation, queries a tensor, or executes a native sampler. Higher analytic Price derivatives certify accuracy but do not raise executed main order beyond C3 or C4.

At an order-k vertex the actual source is the normalized signed same-center VALUE difference

    f_k(x;P,z,t) = (A t)^(-1) 2^(-k) sum_eps (product eps_j)
      [g(z+a t(x+sum eps_j P_j))-g(z+a t sum eps_j P_j)],
    a=(k+1)^(-1/2),

with the recorded radial pullback if required. Its analytic response is a^(k+1)t^k A^(-1) E D^(k+1)g(z+tZ); those response factors are restored in the actual scalar coefficient. No derivative oracle or random raw-tensor estimator is substituted.

## 5 Cost and preserved admission ledgers

The complete original-VALUE bill is

    sum_(tree,old history,new hit,node) sum_v
      2^(k_v+1) M_native(k_v,b_v,epsilon_v)
    + original complete-version services + captures + scalar geometry
    + all affected ancestor and same-center anchor replays.

Triple-cubic main costs at most 44 raw VALUES per force-hit history/node before wrappers. Each changed center, width, private root, requested floor or source version triggers its affected full replay. Original VALUE error delta_g gives source-value error at most 2 delta_g/(A t); this does not certify first/curl error. At fixed admitted native order, LOW30's serial-work theorem gives public-log local floor dependence under its guards, not a growing-order theorem or a free original oracle.

Positive weights preserve the literal scalar S0, complete-first S1, nodewise maxima, and root-square sums without a node-count multiplier. For triple cubic, S0<=C L^15 and S1<=C L^15(1+log(1/eta)). Mixed has S1<=C L^C eta^(-chi), with chi=1 for full cross-history and 1/2 for the certified diagonal. The corresponding squared-first allowance costs eta^(-2chi).

At fixed eta the new node count is polylogarithmic in inverse accuracy. At eta>=c alpha^q_eta it still has inverse-alpha exponent 2q_eta. Independent variance partition is already a valid alternative, but equal shares across J packets cost J^(17/2) in the alpha^17 native row and J^8 in the alpha^18 null-bank row. With otherwise public-log census, the present conservative sufficient admission needs q_eta<7/17. A genuine common-carrier theorem could avoid those share factors; it cannot erase Q itself, old target/cutoff debts or complete replay.

## 6 Precise items requested from the peer

1. A versioned theorem/proof defining “joint” error: carrier-fixed integrated conditional W2, another explicitly caller-stable metric, or ordinary joint W2 with additional hypotheses and constants.
2. Its retained record and exact reserve boundary, and proof that it applies to the actual native residual kernels with P=H/sqrt(v), arbitrary allowed frozen E, private Pi_h V_h, and all mixed replacement stages. State how complete caller constants and finite native floors enter.
3. For the finite-amplitude matrix mixture, its exact target/current, positive sampler map and erased tape list. Keeping only P plus the external record can fit; retaining internal coefficient centers or final reserves changes the contract. Erased coefficient-generation tapes must not still be read by an old output or caller. In particular, our shared B cannot simply be erased before the separately paid original-bank join.
4. For OU quadrature, the precise operator integrand, norm, error and node bound, inverse-width dependence, required derivative/row assumptions, and executable callback. To replace this cubature it must cover F_T and its required caller derivatives, both ordering sectors, full mixed history labels, all proper cuts and HS/sqrt(D), and positive scalar budget preservation. A Hilbert-L2 covariance estimate alone supplies neither those extensions nor a native-law coupling.

A matching carrier result could close the current public-log-dimensional native gap. A matching stronger clock result could improve the eta^(-2) node cost or mixed caller losses. Neither is assumed here, and neither alone proves global arbitrary-order closure or c(P)=o(P).

## Source versions

All local package manifests and imported input pins were checked read-only: 44 references covering 43 distinct files, all matching. Full paths, byte counts and SHA256 values are in SOURCE-PINS.json; individual checks are in SOURCE-VERIFICATION.json. Principal exact versions:

- FINITE-DIMENSION-CARRIER-HYBRID.md: 071a978861d63cd80c7fbca125e76e6c721df14a2afae54e93453c722203beb2
- VARIANCE-PARTITION-NATIVE-REPLACEMENT.md: 535e4c6475b94ac9b930c1d6467d5909c68715f0dfc52a14bac8b83e7d99c5d5
- POSITIVE-SECTOR-GAUSS.md: a4fba9037d305540242dda09adc9dbe70823755652e4a94e585d392a0dfd7dbf
- positive_sector_gauss.py: 5646f982f27791a0026aece20f7299f42eb2b4e0fb106bfb0f4ee1b3c26b920c
- COMMON-CARRIER-AMALGAMATION.md: 75de84123a2a30996a937cc711ca53eb988197460507fd6aebb1765e5af5a5c8
- THREE-CUBIC-EXACT-HEAT.md: d365907bda2bf9d24cb302410a656f9d6612f637bc286e84a0780487e5db23dc
- MIXED-THREE-THREE-SIX.md: 495a6b333f71e5e22fac8c2d2d4fdc4bf3a53771ec414eefed525ae2c1f90b0a
- ACTIVE-PROBE-GENERATOR.md: 3707e549cf3550a8b24c26f54f1dce77ce1469c657309695369403809c58b058
- 30_low_acc.tex: 7767ac389d6b677eb9517b2585f7a24cb9f85ede47d099a653c635b7be750fb8. Relevant labels: b27:nc:raw-paths / b27:compiler:paths, lines 6488–6506; b27:compiler:serial, lines 6853–6873. Fixed-order source eligibility is imported, not independently re-proved here.

Peer excerpt provenance: Slack channel C0C5G45H6P5, thread 1791084240.768399, message 1791217304.784449, supplied by the parent. No peer source file, theorem version or proof hash was supplied. No external message or upload was made for this comparison.
