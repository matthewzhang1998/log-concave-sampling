# Independent audit: cutoff-free two-sector clock quadrature

5 October 2026. This review leaves the sealed inputs and the author's theorem/code unchanged. It audits the Hilbert-valued complex-OU proof, its application to the actual coefficient families, positive scalar quadrature, and the scope of the execution bill. The accompanying independent program checks exact arithmetic and finite force censuses; it does not execute a native sampler or certify source Sobolev floors.

## Verdict

**PASS for the stated fixed-family coefficient-quadrature theorem, with its explicit imported hypotheses and separate native/admission/replay obligations.** The argument removes the prior new-node factor eta^(-2). The displayed bound is Q_tree = 2(128 J m)^2 = O(log^4(M/epsilon)); the original cutoff, physical dimension, and fixed caller envelopes enter this new node count through logarithms.

The initial draft defined eta using only original private variances while importing heat bounds that also need original structural endpoint floors. This was reported immediately and has been corrected in section 1. Private heat alone suffices for the new analytic proof; the preserved old heat/caller bounds require the fuller sealed assumption. The final statement also explicitly fixes W_q before selecting the new node census and identifies contractions as complex-bilinear.

This is not a carrier-retaining native theorem, a conditional-in-retained-bank coefficient identity, a removal of the mixed caller losses, or a whole-engine arbitrary-order complexity theorem.

## 1. Hilbert-valued complex OU

Use the Hermitian convention conjugate-linear in the first slot solely for norms and the energy calculation. For one derivative direction set a = <u,partial u> = x+i y and A_j = |u| |partial u|^2. Gaussian integration by parts gives the energy contribution

    A_j + |u|^(-1) x(x+i y).

Consequently its real part is A_j+x^2/|u| >= A_j, and its imaginary part has absolute value at most

    |x y|/|u| <= |a|^2/(2|u|) <= A_j/2.

Summing proves exactly the claimed Re E >= A and |Im E| <= A/2. Replacing |u| by sqrt(|u|^2+delta) makes the zero-set calculation legitimate and gives the same limiting inequality. Finite-dimensional vector polynomials have sufficient Gaussian integrability for all steps.

For u(t)=exp[-t(1+i beta)L]f, differentiation of the cubed L3 norm yields -3 Re[(1+i beta)E]. This is nonpositive for |beta| <= 1, with room to spare. Neither the Gaussian dimension nor Hilbert-space dimension appears. For w with Re w >= |Im w|, integrating along its ray gives the contraction for P_exp(-w). The w=0 boundary is the identity.

Finite Hermite polynomials are dense in Gaussian L3 with finite-dimensional Hilbert values. For polynomial approximants f_n, contraction gives uniform L3 convergence of P_z f_n on every compact subset of the stated Omega, indeed a uniform difference bound over all its points. The limit is holomorphic. The disk |z|<exp(-pi), including zero, lies in Omega; the power-series approximants handle zero directly. Multiplication of z by a real number in [0,1] only enlarges -log|z| without changing its argument. This proves the holomorphic extension and contraction needed in later sections. No complex-L-infinity estimate or Gaussian density determinant bound is being assumed.

## 2. Actual factorization and bilinear tensor contraction

For the complete original bank dimension, the real Gaussian realization

    B_1 = sqrt(r)G + sqrt(1-r)E_1,
    B_2 = sqrt(r)G + sqrt(1-r)E_2,
    B_3 = s sqrt(r)G + sqrt(1-r s^2)E_3

has unit marginal covariances and pair covariances r, rs, rs. Conditioning on G therefore proves the displayed three-field OU identity for the original sector integrand. It is not an identity conditional on the separately retained execution bank B.

The f_i are the original packet fields after literal occurrence/hit/row expansion, with their original private shields. Each has its original marked force tree. The old bridge inside a six-packet remains internal to that packet's six-force tree, so the same imported one-Hilbert bound supplies ||f_i||_HS <= C_i sqrt(D). This applies to each actual full cross-history term, not just its diagonal. Any old normalization, scalar row contraction, and caller/ancestor row mass needed to sum terms belongs in the actual pre-quadrature W_q.

All coefficient contractions are complex-linear in each argument. They are not Hilbert inner products. Contracting two Hilbert tensors bilinearly along designated slots has HS norm at most the product of their HS norms by absolute-value Cauchy-Schwarz. Successive contraction along the three-packet tree gives the same product bound for its output tensor. Holder with exponents 3,3,3 and the OU lemma then gives C_1 C_2 C_3 D^(3/2).

Complex bilinearity is essential. For example f_1(x)=f_2(x)=x and f_3=1 yields F(r,s)=r; inserting a conjugate would replace the z^2 dependence by |z|^2 and destroy holomorphy. The revised theorem expressly excludes that substitution. The independent check program verifies the exact Hermite pairing exponents for degrees through ten and tests the bilinear HS inequality on complex integer chain tensors.

The analytic envelope spends three Hilbert factors, deliberately. Final absolute-HS error <= epsilon implies each proper-cut error <= epsilon and the requested weaker HS error <= epsilon sqrt(D). This does not claim a complex dimension-sharp proper-cut estimate.

## 3. Local disks and the square root

All four cases in the stated proof of |z-a| < (1-a)/64 lying in Omega are valid:

- a<=1/64 gives |z|<1/32<exp(-pi).
- 1/64<a<=1/16 gives a right-half-plane disk and |z|<5/64<exp(-pi/2).
- 1/16<a<=1/4 gives |Arg z|<=arcsin(1/4)<1 and |z|<17/64<exp(-1).
- a>=1/4 gives |Arg z|<=8b, while -log|z|>=63b, where b=(1-a)/64.

For pure r variation, a local square root of an Omega point remains in Omega because argument and negative logarithmic modulus are both halved. Multiplication by real s preserves it. For pure s variation, only P_(s sqrt(r)) varies and multiplication by the real sqrt(r) preserves the domain.

There is no hidden branch singularity at r=0. For Hermite-polynomial inputs the integrated three-field expression is an even polynomial in z=sqrt(r): changing all three OU parameters from z to -z is canceled by changing G to -G. Equivalently, all nonzero three-Hermite moments have even total degree. The resulting polynomial in r converges uniformly on small r disks by the L3 estimates. It supplies the removable extension at zero and glues the branch expressions. This argument avoids assuming an unrelated negative-z contraction near z=1.

Multiplying by the Jacobian r costs no envelope factor on the coordinate disks because |r|<=1. Cauchy therefore supplies the displayed pure-coordinate factorial bounds with M. Tensor Gauss needs these pure bounds only, not a mixed analytic polydisk or differentiability through the old min-path wall.

## 4. Positive dyadic Gauss and certified scalar error

On a dyadic gap panel [g,2g], subdivision into 128 cells gives length ell=g/128, while every point has gap at least g. Combining Cauchy with the ordinary m-point Gauss remainder gives

    M ell (64 ell/g)^(2m) (m!)^4 / [(2m+1)((2m)!)^2]
    <= M ell 4^(-m).

Scalar duality gives the same estimate in the finite-dimensional Hilbert output norm. The positive two-axis telescope and two sectors cost at most 4M 4^(-m). The omitted strips cost at most 4M 2^(-J). The chosen exact integer budgets make their sum at most epsilon/4. The removed s strip includes part of the ordering-wall neighborhood; it is paid cubature error against the full original finite target, not a change of target.

A normalized reference root displacement d moves a cell node by at most ell d/2. The local first derivative bound therefore costs at most M d/4 per changed coordinate, because 64 ell/(2g)=1/4. Two axes and two sectors cost at most M d. Reference weight L1 error e_w contributes at most 2M e_w. With both at most epsilon/(64M), this is at most 3epsilon/64. There is no hidden inverse-zeta or inverse-eta factor in this scalar perturbation calculation. Total exact quadrature, tails, and reference scalar errors use at most 19epsilon/64, leaving a genuine unspent budget for separately allocated execution floors.

The copied rational core is byte-identical to the sealed prior builder. It certifies root intervals by exact signs, disjointness and degree; encloses weights by rational interval arithmetic; and symmetrizes/normalizes positive weights to exact mass and first moment. These properties survive cell mapping. Root proposals alone are not treated as proofs. The new parameter chooser uses an exact rational M, integer comparisons, and exactly 128J cells. The Cartesian rule is lazy, but the stated Q remains an executed per-tree node count, not a claim that execution can skip those nodes.

## 5. Heat budgets and complete target scope

For x=1-r and y=1-s, the original edge gaps are x and x+y-xy. Positive cell weights satisfy lambda_r<=x and lambda_s<=y, hence w=r lambda_r lambda_s<=xy<=x(x+y-xy). Both edge gaps, and therefore each replica innovation gap, vary by at most two across a parent dyadic product panel. Adding eta preserves that comparison. Positive monomial heat majorants have a fixed-rank sup/inf ratio. Exact integration of r on each product cell transfers each majorant from its sector integral to the node sum with no Q multiplier.

Thus the sealed S0, S1, nodewise maxima, and maximum-times-sum literal square estimates are preserved under their original hypotheses. In particular the full mixed first loss eta^(-1), diagonal mixed first loss eta^(-1/2), and squared losses eta^(-2) and eta^(-1) remain. The structural-floor correction is necessary for quoting those inherited eta bounds. Repeated old clock labels retain their literal powers; no independent-clock variance identity is used.

Execution continues to use the original complete bank B and the sealed retained-bank disintegration with R R^T=K-c11^T-diag(delta_i). The OU coupling is proof-only. It neither replaces that execution geometry nor exposes a new observer. Simultaneous m_star extraction remains frozen across caller estimates. For the six-packet, new extraction uses only its recorded new innovation residual blocks, preserving old private heat and original common-bank rows. Full cross-history targets are never replaced by diagonal targets.

The only coefficient identity used is bank-mean, uniformly in the admitted fixed external callers and their actual derivative slots/rows. All old/new mixed currents, strict-ancestor captures, and same-endpoint joins remain owed. Marginal native LAW is not promoted to a carrier-fixed or environment-stable theorem. This agrees with the sealed peer-carrier interface.

## 6. Literal VALUE and replay bill

The execution cost is the sum over every tree, original history, new hit assignment, node, and force vertex of 2^(k_v+1) M_native(k_v,b_v,epsilon_v), plus old complete requested versions, captures, roots, physical rows, keeps, and affected replays. Analytic OU operators and Cauchy derivatives are not computational callbacks. Main maximal adapter order stays C3 for cubic and C4 for mixed.

Independent enumeration gives:

- Triple cubic: 243 new hit histories, total primitive order 7, maximal order 3, maximal raw VALUE count 44. A maximizing order list is (0,3,0,0,2,0,0,2,0).
- Mixed: 648 new hit histories times all nine old bridge hit pairs, hence 5,832 histories; total primitive order 10, maximal order 4, maximal raw VALUE count 72. A maximizing order list is (0,2,0,0,2,0,0,4,0,0,2,0).

The maxima are before native wrappers, original old-source versions, physical/readout factors, and replay. Any changed center, private root, width, requested floor, source order/version, or ancestor triggers its complete affected replay. The source value floor 2 delta_g/(A t) does not certify source-first, curl, or native Sobolev error.

W_q in M is an envelope of the original target fixed before choosing Q. Later Q-dependent native readout/adapter factors that cancel on reconstructing the target belong in the native bill; inserting them into M without a certified joint parameter choice would be circular. The final text correctly states this distinction.

Small retained innovation gaps are only square-rooted in the displayed Gaussian disintegration. Any chosen native implementation that introduces an inverse innovation gap must charge it separately. Original history growth, W growth, native admissibility/readout shares, mixed caller loss, old truncation debt, source/replay exponents, and rank growth can still obstruct a whole-engine c(p)=o(p) result. The coefficient theorem does not remove them.

## Reproducibility

The independent program `check_complex_ou_clock_audit.py` passes 3,175 checks, recorded in `complex-ou-clock-audit-checks.json` and `audit-check-run.txt`. It verifies all 13 sealed source pins and includes independent raw force enumeration, exact Hermite monomial identities, complex-bilinear tensor contraction checks, rational envelope and error budgets for private floors down to 2^(-80) and D up to 10^8, and additional certified Legendre degrees 1, 3, 7 and 11 at 10^(-30) tolerances. Source hashes in the check output identify the reviewed versions, including `COMPLETE-EXECUTION-BILL.md`. That bill's active VALUE formula matches the pinned original callback; its per-history census and changed-key replay requirements agree with the checks above. `REVIEW-PIN.json` records the final reviewed file hashes.

These checks supplement, rather than replace, the analytical energy, holomorphy, and exact-bank arguments. No source-native theorem or target-tensor oracle is numerically tested or presumed.
