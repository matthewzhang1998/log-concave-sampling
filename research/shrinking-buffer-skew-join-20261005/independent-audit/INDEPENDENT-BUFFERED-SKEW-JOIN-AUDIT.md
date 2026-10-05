# Independent audit: bounded shrinking-buffer cubic join

2026-10-05. Audit of the proposed bounded higher-cumulant extension of the sealed LAW-only shrinking-buffer construction. This is a source-qualified theorem under the literal imported finite native guards, not a numerical execution of LOW30 and not an all-order recurrence.

## Verdict

**PASS under the declared imported native guards.** The reviewed complete draft is `../POSITIVE-SKEW-SHRINKING-BUFFER-JOIN.md`. Its Section 7 now includes the explicit positive-reference reconciliation below; it is a necessary part of this verdict. The scaling, true-F2 versus H cumulant replacement, positive variance budget, complete-bank law consumption, and canonical-m3 rate are valid. A native quadratic reference is initially output-slot oriented, whereas the skew consumer uses a fully symmetric third tensor. A third-cumulant match alone does not identify their full laws. Section 5 supplies the missing positive interpolation, at exactly the already budgeted cubic-feedback order. No new executing tensor input or packet symmetrization is necessary.

A particularly clean variance allocation leaves an external untouched buffer v/2 and gives the independent mean, Gram, mixed-K, and cubic banks v/8 each. Include all their Gaussian references, including the cubic packets' internal buffers, in the visible carrier S. Its covariance is exactly vI/2 plus the force covariance, after separately paying the admitted covariance-target restoration. This matches the skew consumer's half-buffer convention literally.

The final substantive canonical-m3 mean error is bounded by a fixed public-log polynomial times A^(7/2) sqrt(D), plus the declared absolute floors, at w=A. This claim keeps the logarithms: they occur at the leading exponent and cannot be silently absorbed into a numerical constant.

## 1. True-F2 versus H: a dimension-safe third-cumulant bound

Work conditionally on the retained Y. On the same complete Gaussian bank write u=F2-E[F2|Y], h=H-E[H|Y], and d=u-h. Suppose the complete-bank Lipschitz constants of u,h are at most CA, and the already established coherent comparison gives

    ||d||_(L2(Y,bank)) <= C A² sqrt(D).

The covariance matrices of u,h are bounded by CA²I by Gaussian Poincare. For an arbitrary deterministic matrix M, another Poincare inequality gives

    Var(uᵀMu|Y) <= C A⁴ ||M||HS²,
    Var(hᵀMu|Y) <= C A⁴ ||M||HS².

Indeed D(hᵀMu)=DhᵀMu+DuᵀMᵀh; the operator bounds on Du,Dh and the covariance bounds on centered u,h control both squared terms by CA⁴||M||HS². No dimension-sized vector energy is multiplied by another such energy.

For any centered random vector d and Hilbert-valued Q, the joint covariance block is PSD and consequently

    ||E[d tensor Q]||HS² <= ||Cov(Q)||op E|d|².

Apply this to Q=u tensor u, h tensor u, and h tensor h after centering Q. The exact tensor telescoping identity is

    u³-h³=d tensor u tensor u
           +h tensor d tensor u+h tensor h tensor d.

Permuting tensor slots is isometric. Therefore, conditionally and then in integrated L2(Y),

    ||kappa3(F2|Y)-kappa3(H|Y)||_(L2(Y);HS)
        <= C A² ||d||2 <= C A⁴ sqrt(D).                 (1)

The delayed negative force -qF2 multiplies the third tensor by -q³. This scalar, its sign, and the analogous mean/covariance scalings must be included at the actual source readout. Since 1/2<=q<=1, they change no exponent or native admission scale.

## 2. Exact shrinking-buffer specialization of the skew consumer

For a centered Gaussian-bank source with first L, energy e, covariance Sigma, and K=kappa3/6, apply the audited fixed-buffer theorem to X/sqrt(v). Its first and energy are L/sqrt(v), e/sqrt(v); its third coefficient is K/v^(3/2). Multiply its law error by sqrt(v). Under L/sqrt(v)<=1 and the admitted public-log smallness guards,

    W2(buffered true source, positive quadratic reference)
        <= C Lambda² L³ e / v^(3/2).                  (2)

The four unsimplified scales are

    L³e/v^(3/2), L²h/v², Lambda kh/v^(5/2),
    Lambda²k²h/v⁴,

where h=||K||HS<=CL²e and k is a proper-cut bound <=CL³. Each is dominated by (2) at small L/sqrt(v). This displays the v factors rather than importing a fixed-buffer constant unchanged.

With the exact half-buffer convention, the reference is

    m+S+b_K(S)+sqrt(v/2)Z,
    Cov(S)=C=vI/2+Sigma,
    b_K(S)=K:[(C⁻¹S) tensor (C⁻¹S)-C⁻¹].             (3)

K and C are analytical coefficients only. The existing native VALUE graphs execute the actual program.

For a common C with C>=cvI, Hermite isometry gives

    ||b_K(S)-b_Kprime(S)||2 <= C ||K-Kprime||HS/v.      (4)

Thus (1) costs CA⁴sqrt(D)/v at the input law, and CA⁵sqrt(D)/v after the A-Lipschitz terminal F. There is no extra dimension factor. All statements can be integrated in Y because their constants are uniform and the outer geometry gives the correct standard-Y marginal.

## 3. Native cubic packet scaling and ports

Let c0,a0,b0,eta0 be the old fixed normalized positive shares for a node, with its native root coefficient carrying 4w_node/(c0 a0 b0 sigma2). Give the complete cubic bank physical variance u. Use physical readouts sqrt(u) times the old readouts and put

    alpha=A/sqrt(u)

at every original normalized vertex, with the root additionally multiplied by its old weight/shield/share coefficient. Then the physical leading coefficient is exactly

    u^(3/2) c0 a0 b0
      * [4w_node alpha³/(c0 a0 b0 sigma2)]
      * (sigma2/A³) tree =4w_node tree.                (5)

The target is unchanged, including its 24/6 normalization. If the output is subsequently multiplied by -q, run the unscaled bank at variance u/q² and apply -q to its whole output; this supplies exactly -q³kappa3(H) while returning buffer u.

The native full-bank and retained-caller path estimate rescales to

    sqrt(u) Lambda alpha = Lambda A.

The descendant terms sqrt(u)alpha² and sqrt(u)alpha³ are no larger when alpha<=1. The original positive clock/shield sums and exact anchor rules remain in force. This is a graph derivative bound, not a derivative inferred from W2.

The substantive node feedback sums to

    Lambda sqrt(u) alpha⁶ sqrt(D)
       =Lambda A⁶u^(-5/2)sqrt(D).                     (6)

The pair prior is Lambda sqrt(u)alpha^b sqrt(D). Choosing any fixed b>=6 absorbs it into (6). At b=4 it is instead Lambda A⁴u^(-3/2)sqrt(D), which is still covered by the existing mean/Gram bill but cannot be called an arbitrarily tiny numerical floor. Finite clock/filter/calibration/clipping/numerical errors are separately assigned their actual readout and inverse-buffer multipliers.

At u comparable to A, alpha is comparable to sqrt(A). Every original source radius, clipping gap, selected-field bound, pair order, and serial graph guard must be checked at that actual normalization. In particular, the known root coefficient contains its explicit weight/sigma2 and public-log inverse shares. No realized source energy appears in a denominator.

## 4. Exact anisotropic joining with positive laws

Freeze (z,Y). Reserve a fresh independent external Gaussian of variance v/2. Allocate v/8 each to the mean, Gram, mixed-K, and cubic banks. Each bank is complete and fresh conditional on Y. The covariance service targets are the actual admitted targets, including the signed mixed-K correction, whose own strictly positive bank gap remains enforced.

After consuming the entire individual conditional LAW comparisons, write all service reference Gaussians in one vector bank. Their sum S has covariance C=vI/2+Sigma_tilde. Cubic corrections depend initially on their own independent visible pieces S_j with positive isotropic covariance C_j, and have the form

    B=sum_j K_j:[(C_j⁻¹ S_j) tensor (C_j⁻¹ S_j)-C_j⁻¹].

Other reference Gaussians, including cubic internal eta_j buffers, enter S linearly. Standard Gaussian regression gives exactly

    E[B|S]=(sum_j K_j):[(C⁻¹S) tensor (C⁻¹S)-C⁻¹].    (7)

This remains true for anisotropic C: Cov(C_j⁻¹S_j,S)=I, and the conditional covariance correction cancels the old C_j⁻¹ term. No commuting matrix assumption is involved.

Let U be the independent Gaussian complement of S after whitening and an orthogonal Gaussian coordinate change. Set E=B-E[B|S]. Along the positive path

    W_t=S+E[B|S]+tE+sqrt(v/2)Z,

the private-U Riesz identity yields

    d E phi(W_t)/dt
      =t E[(R_U E)(D_U E)ᵀ:D²phi(W_t)].               (8)

The base S+E[B|S] contains no U; there is no first-order unpriced branch. The fixed-degree Hilbert/matrix-chaos bounds used in the audited feedback lemma imply

    ||R_U E||_(fixed moment) <= Lambda h_sum/v,
    ||D_U E||_(operator, fixed moment) <= Lambda k_sum/v,

where h_sum=sum_j||K_j||HS and k_sum is the corresponding summed proper-cut scale. Individual inverse shares contribute only their actual public-log powers. One integration in the external untouched buffer proves

    W2(Law(S+B+sqrt(v/2)Z), Law(S+b_Ksum(S)+sqrt(v/2)Z))
       <=Lambda k_sum h_sum/v^(5/2).                 (9)

For the finite native clocks, k_sum<=Lambda A³ and h_sum<=Lambda A³sqrt(D), so this is precisely Lambda A⁶v^(-5/2)sqrt(D). This is a one-energy estimate. The covariance Sigma_tilde is bounded in operator norm by CA² after target restoration, and every visible Gaussian covariance has the required gap. Its dependence on Y is harmless because Y is retained and frozen during (7)-(9).

## 5. Output-slot orientation: the necessary positive reconciliation

A native packet's quadratic coefficient need not already be fully symmetric in its three physical slots. Denote its post-regression coefficient T, symmetric in the last two slots, with Sym3(T)=K. Put D=T-K, so Sym3(D)=0. Proper cuts and Hilbert norms of K,D are bounded by fixed multiples of those of T.

Use the SAME C and external buffer and interpolate genuine positive laws

    W_t=S+b_(K+tD)(S)+eta Z, eta²=v/2.

Let J_t=D_S b_(K+tD). Two Gaussian integrations by parts with covariance C give the exact identity

    d E phi(W_t)/dt
      =E[ D_iab ((I+J_t)_pa (I+J_t)_qb-delta_pa delta_qb)
                partial_i partial_p partial_q phi(W_t)]
       +E[ D_iab partial_a partial_b b_(K+tD),p(S)
                partial_i partial_p phi(W_t)].        (10)

The subtracted leading term vanishes because third test derivatives are symmetric and Sym3(D)=0. This is the full justification for discarding the output-slot asymmetry; a third-cumulant match alone would not have sufficed.

Writing h,k for fixed-multiple Hilbert/proper-cut bounds on T, Gaussian matrix-series moments and C>=cvI give

    ||J_t||_(operator, fixed moment) <=Lambda k/v^(3/2),
    rank-3 current <=C[Lambda kh/v^(3/2)+Lambda²k²h/v³],
    rank-2 current <=Ckh/v².

The rank-2 contraction is the product of two tensor flattenings with C⁻¹ on each contracted index, bounded by one proper cut times one Hilbert norm. All currents are independent of Z. Integrating rank-minus-one derivatives in eta Z gives

    W2(reference T, reference Sym3(T))
      <=C[Lambda kh/v^(5/2)+Lambda²k²h/v⁴].           (11)

Under the declared smallness Lambda k/v^(3/2)<=c, the second term is absorbed into the first. For k<=Lambda A³ and h<=Lambda A³sqrt(D), (11) is covered by the existing cubic-feedback bill. This proof adds no executing packet, derivative query, or retained root.

## 6. Covariance restoration also moves the quadratic reference

The reviewed draft restores the independent Gaussian component BEFORE the final regression. That ordering is valid: keep all cubic references fixed, couple only the independent Gaussian component, and then regress directly onto the restored C. No nonlinear-reference covariance correction is needed in that ordering. For completeness, if restoration were instead performed AFTER regression, one must not restore only the Gaussian covariance and silently leave the quadratic reference at its old C. Couple S_i=C_i^(1/2)G for C_i>=cvI and use b_K(S_i)=K:(C_i^(-1/2))^(tensor 2):(GGᵀ-I). The inverse-root Sylvester estimate and tensor proper cuts yield

    ||S_1-S_2||2 <= C||Delta C||HS/sqrt(v),
    ||b_K(S_1)-b_K(S_2)||2 <= Ck||Delta C||HS/v².     (12)

Thus restoring the complete reference costs

    C||Delta C||HS/sqrt(v) * [1+Ck/v^(3/2)].          (13)

At the admitted smallness this is the same order as the previously budgeted Gaussian covariance restoration. An integrated CA⁴sqrt(D) covariance discrepancy therefore costs CA⁴sqrt(D)/sqrt(v) before the terminal F. No extra leading term is hidden by the nonlinear reference.

## 7. Error ledger and own-mean source ports

Keep the previously audited prefix reduction, exact conditional geometry, bulk boundary t=q=1-w, positive outer quadrature, and fresh mean-matched near-endpoint branch. The bulk has v comparable to w. The complete mean estimate is bounded by sqrt(D) times

    C[A²w^(3/2)+A³w+A³sqrt(w)+A⁴
       +Lambda² A⁵w^(-3/2)+A⁵/w]
     +Lambda[A⁵w^(-3/2)+A⁶w^(-7/6)+A⁶w^(-2)
              +A⁷w^(-5/2)],                          (14)

plus the separately restored absolute floors and smaller covariance/mixture terms. The A⁵/w term pays the true-F2 versus H tensor discrepancy and the declared same-grade clock discrepancy. The reviewed draft uses b=5 and explicitly pays the displayed A⁶w^(-2) term. Alternatively, b>=6 absorbs the prior into the A⁷w^(-5/2) term. Mean/Gram source normalization remains the actual A/sqrt(u), and the rebalanced mixed-K source radius remains A u^(-1/3).

At w=A the exponents in (14) are respectively 7/2,4,7/2,4,7/2,4,7/2,29/6,4,9/2. Hence the result is Lambda A^(7/2)sqrt(D)+e_abs, with Lambda a fixed polynomial in actual public logs. The critical mean/Gram residual first remains Lambda[A+A^(3/2)/sqrt(u)]=Lambda A. Its normalized self-reserve guard has scale Lambda sqrt(A), as does its normalized first. The K radius is Lambda A^(2/3); the cubic vertex radius is Lambda sqrt(A).

The executing bulk graph has a known isotropic source-zero carrier of variance v, including the external reserve and all four bank carriers. Combine that scalar Gaussian row and use the already audited common-carrier rotation. Each service's actual carrier-subtracted residual has full private/caller first Lambda A. The terminal original-g graph consequently has the same gradient baseline and square-lift curl Lambda A², energy Lambda A²(|z|+sqrt(D)), and full first Lambda A as before. The final optional own-mean compiler adds its existing Lambda A⁴sqrt(D) term.

These source-port claims are direct actual-graph chain rules. They do not retain any reference root after a law comparison. The complete service laws are compared at fixed (z,Y), with every owned bank integrated; afterward one-node expectations are summed. The rotated source is separately evaluated on its actual complete tape. All source banks, including passives, filters, coarse roots, original anchors, and zero-coefficient roots, remain charged. The added five-clock cubic packets multiply fixed-order native/filter counts by a public-log clock count. The untouched external reserve adds one D-root per bulk node. Every reentered source occurrence replays its full graph; no covariance or tensor is an executing leaf.

## Scope of diagnostics

The sibling independent checker tests exact covariance/cross-covariance inequalities, tensor telescoping, anisotropic conditional regression, the Gaussian-row variance identities, physical cubic normalization, feedback/prior exponents, inverse-root quadratic-reference sensitivity, and the exact positive orientation-interpolation identity on polynomial fixtures. Those checks support the displayed algebra; they do not execute the imported native compiler or replace its guarded source qualification.

## Final reviewed source and evidence

The complete reviewed draft has SHA256 `3c93c8964ccdb1f1028746aae661b768eacb8870ec616ce95d613110ebcaa708`. Its actual choice b=5, explicit A⁶/v² terminal prior, half-buffer allocation, complete-bank readsets, and paid output-slot symmetrization all pass. No further substantive correction is required. The independent checker passes 6,418 assertions, including four exact rational, non-isotropic Gaussian polynomial fixtures for (10). Those fixtures also verify that the orientation feedback is nonzero: the correction cannot be omitted merely because its leading fully symmetric tensor vanishes.
