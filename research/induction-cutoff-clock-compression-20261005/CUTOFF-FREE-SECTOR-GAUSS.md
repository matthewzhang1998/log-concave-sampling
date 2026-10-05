# Cutoff-free positive sector Gauss via Hilbert-valued complex OU

5 October 2026. Constructive clock-node theorem for the complete original three-cubic and 3/3/6 coefficient families. No target heat is added. The new bridge-node count is polynomial in logarithms of the original cutoff, dimension, and requested accuracy. The previous eta^-2 node factor is removed. The native carrier theorem, mixed caller losses, old-source/replay exponents, and original truncation debt are not changed.

## 1. Quantitative statement

Let D>=1 be the physical dimension and 0<eta<=1/2 a lower bound on the ORIGINAL private variances in the finite input program. Where the imported heat/caller budgets require structural endpoint variances or recorded floors, eta also lower-bounds those ORIGINAL quantities, exactly as in the sealed inputs. The new analytic quadrature argument itself needs only the private-variance condition. Fix the complete literal finite history, hit, source-version, affine-row and physical/readout census. For main triple cubic use N=9, K=7; for full or diagonal 3/3/6 use N=12, K=10. A required q-fold external AFFINE caller derivative has total primitive degree K+q and its own full actual scalar mass W_q, including injection rows and ancestors. These are the same original bank-mean coefficient and callers as in the sealed input.

Put

    B_q = W_q 2^(N+3(K+q)/2) sqrt((K+q)!) eta^(-(K+q)/2),
    M >= 1 + D^(3/2) sum_(required q) B_q.

W_q is the total unsigned scalar mass over ALL trees, original histories and new hit assignments, including physical/source normalizations unless they are separately restored. It is not a per-history constant. It is fixed from the ORIGINAL coefficient target before selecting the new node count. New native readout/adapter factors which cancel when reconstructing that target belong to the separately charged native bill. One cannot insert Q-dependent factors into M circularly without a separate certified parameter choice. One may use the prior rational `exact_B` upper bound and D ceil(sqrt(D)); no floating comparison is needed. For main plus first take q=0,1. The unit-row shorthand W_1=N W_0 is valid only when all actual caller rows are <=1 and no omitted ancestor constant is required.

For 0<epsilon<1 let

    J=max(1,ceil(log2(32M/epsilon))),
    m=max(1,ceil(log4(32M/epsilon))),
    H=128.

There is a finite positive two-sector rule with

    Q_tree = 2 (H J m)^2.                                      (1)

Its error for the COMPLETE summed bank-mean coefficient and every included caller is at most epsilon in absolute HS norm, hence at most epsilon in every proper cut and at most epsilon sqrt(D) in HS. The absolute-HS certificate is deliberately stronger than required; the temporary D^(3/2) analytic envelope affects only logarithms. Therefore

    Q_tree = O(log^4(M/epsilon))
           = O([log W + (K+q_max)log(1/eta)
                 + log D + log(1/epsilon) + rank constants]^4). (2)

Constants and the fixed H are literal. Execute all three trees, every original history and all 243 cubic / 648 mixed new hit choices; Q_tree is a per-tree node count, not the whole bill. A fixed finite set of callers uses one common rule by summing its envelopes in M.

## 2. An elementary dimension-free Hilbert-valued complex OU lemma

Let gamma_d be any standard finite-dimensional Gaussian and H any finite-dimensional complex Hilbert space. Let L=-Delta+x.dot grad and P_z be the Hermite multiplier z^ell on degree ell. For real 0<=z<=1 this is the ordinary OU conditional expectation. Then

    ||P_z f||_(L3(gamma_d;H)) <= ||f||_(L3(gamma_d;H))          (3)

whenever z=exp(-w), Re w>=|Im w|. The norm is independent of d and dim H. Also P_0 is the expectation and belongs to the domain below.

Here is a direct proof, not an assumed analytic-radius theorem. For a Hilbert-valued finite Hermite polynomial u, Gaussian integration by parts gives

    E(u)=integral |u| <u,Lu>
        =sum_j integral [|u| |partial_j u|^2
             + |u|^-1 Re<u,partial_j u> <u,partial_j u>].

The harmless zero set is handled by replacing |u| with sqrt(|u|^2+delta) and passing to the limit. Consequently Re E>=A, |Im E|<=A/2, where A=sum_j integral |u| |partial_j u|^2. Along u(t)=exp[-t(1+i beta)L]f,

    d/dt ||u(t)||_3^3 = -3 Re[(1+i beta)E(u(t))] <=0

for |beta|<=1 (indeed <=2 is available). The same calculation with either convention for the complex inner product changes only the sign of beta. Integrate along the ray to w. Finite Hermite polynomials are dense in L3(H), so contraction extends by completion. Polynomial approximants converge uniformly in L3 on every compact subset of

    Omega={z:0<|z|<1, |Arg z|<-log|z|} union {0}.

Thus z -> P_z f is holomorphic there. The disk |z|<exp(-pi) is contained in Omega, which includes an open neighborhood of zero. Multiplication by a real number in [0,1] preserves Omega and the contraction. This proves both analyticity and the operator bound without determinant factors, an L-infinity analytic assertion, or a dimension-dependent Gaussian density bound.

## 3. The exact integrand is a three-field OU product

Fix one replica tree, one original ordered history tuple and one new hit assignment. Rename replicas so the sector covariance has

    K_12=r, K_13=K_23=rs, K_ii=1.

Expand every original bank derivative over its literal force occurrences and scalar injection rows. The coefficient of a bridge contains the original row inner product, not a differentiated covariance root. After this expansion let f_i(B_i;y) be the Hilbert tensor field of the entire ith old packet with its assigned bridge/caller hits and its ORIGINAL private heat. Its free slots are all original physical slots plus those new edge/caller slots. Its original internal force graph has a marked spanning tree. The imported one-Hilbert bound therefore gives, uniformly in every B_i,y,

    ||f_i(B_i;y)||HS <= C_i sqrt(D).

The product of the three C_i, summed with the literal scalar weights, is at most the relevant B_q. This also applies to the six-packet: its old covariance bridge remains internal to its original marked force tree, and all original residual-bank rows remain in its descriptor. A full cross-history term is never replaced by a diagonal term. Repeated old clock labels remain repeated.

For real r,s in [0,1], conditionally independent replicas with common G can be realized by

    B_1=sqrt(r)G+sqrt(1-r)E_1,
    B_2=sqrt(r)G+sqrt(1-r)E_2,
    B_3=s sqrt(r)G+sqrt(1-r s^2)E_3.

All banks have the COMPLETE original bank dimension. Their covariance is exactly K. Thus the actual integrand, with the two bridge-index contractions denoted C, is

    F(r,s;y)=E_G C[(P_sqrt(r) f_1)(G;y),
                  (P_sqrt(r) f_2)(G;y),
                  (P_(s sqrt(r)) f_3)(G;y)].             (4)

This is a proof identity only. Execution uses the original sealed covariance disintegration and original VALUE sources, not evaluations of P_z or a tensor oracle.

All contractions here use the original coefficientwise bilinear pairing extended complex-linearly in both slots, not a Hermitian pairing; this is essential for holomorphy. The absolute-value Cauchy-Schwarz bound still applies. The bilinear contraction of Hilbert tensors along a tree is bounded by the product of their HS norms, by successive Cauchy-Schwarz. Hence Holder with exponents 3,3,3 and (3) bounds the HS norm of (4) by C_1 C_2 C_3 D^(3/2) whenever the three OU parameters lie in Omega. This intentionally sacrifices the one-Hilbert factor ONLY in the analytic envelope. We do not assert dimension-free complex L-infinity or complex proper-cut contractivity.

## 4. Proved local radii

For every a in [0,1), the complex disk

    |z-a| < (1-a)/64                                      (5)

lies in Omega. Here are elementary sufficient estimates. If a<=1/64, then |z|<1/32<exp(-pi), so the whole disk lies in Omega. If 1/64<a<=1/16, then the disk is in the right half plane, |Arg z|<pi/2, and |z|<5/64<exp(-pi/2). If 1/16<a<=1/4, then |Arg z|<=arcsin(1/4)<1 and |z|<17/64<exp(-1). If a>=1/4, put b=(1-a)/64. Then b/a<=3/64, so |Arg z|<=arcsin(b/a)<=2b/a<=8b; meanwhile -log|z|>=1-|z|>=63(1-a)/64=63b. In all cases the desired strict inequality holds.

For a complex r in (5), choose a local square root. If r is in Omega, so is sqrt(r), since its logarithmic radius and argument are both halved; multiplication by real s also preserves Omega. The expression (4), as a function of z=sqrt(r), is even: P_-z f(G)=P_z f(-G), and the simultaneous sign change disappears on integrating G. Therefore it extends holomorphically through r=0 and does not acquire a square-root singularity. This justifies using the whole disk even when it crosses the negative real axis. For a complex s in its disk and real r, s sqrt(r) lies in Omega, while the other two parameters remain real.

Define G(r,s)=r F(r,s), including the sector Jacobian. In the respective pure-coordinate disks, |r|<=1. Summing histories, trees and required callers gives the Cauchy certificates

    ||partial_r^ell G(r,s)||HS <= M ell! [64/(1-r)]^ell,
    ||partial_s^ell G(r,s)||HS <= M ell! [64/(1-s)]^ell.     (6)

These are pure-coordinate estimates, which are all tensor Gauss requires. No mixed derivative or smoothness through a min-path wall is assumed. Caller derivatives are fields with their real extra slots and rows before applying the same argument; no completed LAW statement is differentiated.

## 5. Positive dyadic rule and complete error

Set zeta=2^-J. Retain r,s in [0,1-zeta]. On each coordinate use the dyadic gap panels

    1-r in [2^(-j-1),2^-j], j=0,...,J-1.

Split each panel into exactly H=128 equal intervals. Put the ordinary m-point positive Gauss-Legendre rule on every interval. For the two sectors use both assignments of large and small original tree parameters:

    t_large=r, t_small=rs, w=lambda_r lambda_s r>0.        (7)

On a subinterval I, its length ell is g/128, where g is the minimum gap on its dyadic parent. From (6), the Gauss remainder is at most

    M ell [64 ell/g]^(2m)
        (m!)^4/[(2m+1)((2m)!)^2] <= M ell 4^-m.

The formula holds for Hilbert-valued integrands by duality. Positivity bounds quadrature/integration operator masses by one. Telescope the two axes, sum the two sectors and all family terms: interior error <=4M 4^-m. The omitted strips cost <=4M zeta. Our J,m make their sum <=epsilon/4. The remaining budget includes certified rational node/weight error and every separately allocated source/native floor; the quadrature certificate itself is below epsilon after scalar arithmetic.

No change of target is hidden in truncation: the removed strip is paid as a cubature error against the original full finite coefficient. Original old-clock cutoffs and errors stay unchanged.

## 6. Literal weights, clock heat, and caller admissions are preserved

Each interval remains inside its original dyadic gap panel. Its positive node mass obeys lambda_r<=1-r and lambda_s<=1-s, so w<=xy<=(1-r)(1-rs), with x=1-r,y=1-s. The gaps x and x+y-xy and each replica innovation gap vary by at most two within a dyadic product panel; adding eta preserves this. The rule integrates r exactly on each product cell. Thus every positive heat majorant in the sealed cubic/mixed proofs has its weighted sum bounded by a fixed-rank comparison factor times its original sector integral. There is no Q multiplier in S0,S1 or the literal root-square estimates.

Triple cubic retains S0<=C L^15, S1<=C L^15(1+log(1/eta)), and bounded nodewise main/first maxima. Full 3/3/6 retains its caller eta^-1 loss; diagonal retains eta^-1/2. Their squared-first allowances remain eta^-2 and eta^-1, respectively. These are real source/caller costs not removed by this clock theorem. Actual physical/readout inverse factors are still charged. Nothing here imports a carrier-retaining native LAW from marginal LAW.

## 7. Finite implementation and all original-VALUE work

`cutoff_free_sector_gauss.py` reuses a pinned local copy of the prior exact rational Legendre builder. It computes a rational M, chooses J,m by integer power comparisons, forms the HJ rational intervals, and returns a lazy positive node iterator. Legendre estimates only propose rational sign-changing root intervals; exact polynomial signs, disjointness and degree certify root isolation. Positive rational weights are symmetrized and normalized to exact interval mass; their exact first moment is preserved. Every approximate node remains strictly in its cell.

Let d bound the reference [-1,1] root error and e_w the total reference weight error. Mapping the cell and using (6), one changed coordinate costs at most M d/4 pointwise, because ell/g=1/128. The two axes and two sectors give <=M d. The weight telescope contributes <=2M e_w. Taking both tolerances <=epsilon/(64M) uses <epsilon/16. This local scaling is why scalar precision does not require an inverse-zeta or inverse-eta factor. Required bit precision is polynomial in log(M/epsilon), m and input bit sizes. Node count, scalar preparation, and actual execution are separate bills.

At every retained node execute the sealed exact retained-bank disintegration

    B_i=sqrt(c) B+(R V)_i+sqrt(delta_i) E_i,
    R R^T=K-c 11^T-diag(delta_i),

on the SAME complete original B and required fresh complete banks. The proof coupling (4) is not substituted into the native caller graph. Keep the original injection rows and exact simultaneous m_star private-heat extraction. For a six-packet extract only from its recorded NEW innovation residual blocks. Original private heat is neither enlarged nor discarded. The identity is bank-mean, not pointwise in retained B; all old/new mixed currents and captured strict-ancestor callers remain owed.

The complete bill is

    sum_(tree,original history,new hit,node) sum_v
        2^(k_v+1) M_native(k_v,b_v,epsilon_v)
    + every old complete requested version, caller capture,
      scalar root, physical row, keep, and affected anchor/ancestor replay.

Main executed order remains C3 for cubic and C4 for mixed. Analytic OU, Cauchy and high Price orders are never queried. Raw cubic sources use at most 44 VALUE calls per hit history/node before wrappers; mixed uses at most 72 (the exact census is checked in the executable). A changed center, private root, width, source order, requested floor or original source version triggers its complete affected replay. VALUE accuracy delta_g gives source-value error <=2 delta_g/(A t); first/curl and native Sobolev floors remain separately certified obligations. No ideal tensor, conditional expectation, Monte Carlo coefficient estimator, or uncharged precomputation is executed. Innovation gaps can be as small as the accuracy-dependent zeta; the displayed disintegration uses only their nonnegative square roots and does not invert them. If a selected native implementation nevertheless introduces an inverse-innovation-gap cost, it must be charged separately; this node theorem does not certify its absence.

For each desired final order p, eta=alpha^q_eta(p), epsilon=alpha^p, this NEW node factor is polynomial in [p+(K+q_max)q_eta(p)] log(1/alpha), log D and log W, at fixed family/rank data. It has ZERO new inverse-alpha cutoff exponent, even if q_eta grows with p; its p-dependent constants and all inherited exponents remain literal. A whole-engine c(p)=o(p) conclusion still requires bounding W, old source/replay costs, rank growth, actual native readout shares, and every mixed caller loss. If independently owned groups dilute variance across the number of packets, those powers now inherit only this polylogarithmic new-node count; existing polynomial history counts or other losses do not disappear.

## 8. Lowest case and scope of extension

For two replicas and one bridge, F(t)=E C[f_1(G),(P_t f_2)(G)]. Hilbert-valued L2 contraction holds for the entire unit disk by Hermite orthogonality. Cauchy radius (1-t)/2 and positive dyadic Gauss give Q=O(log^2(M/epsilon)); M can be D times the two packet envelopes. This already removes the cutoff power in the scalar D=1 one-bridge case and in every finite D, with D only inside logarithms.

The nontrivial extension proved above is the two-clock, three-replica product. It includes arbitrary finite packet force counts with an original marked-tree uniform Hilbert envelope, in particular the actual nine-force triple cubic and twelve-force full/diagonal 3/3/6, plus any fixed required caller order. It does not assume independent old histories, a diagonal mixed target, or pointwise-in-bank cancellation. A general n-replica hierarchical OU extension is plausible, but its complete positive sector Jacobian, Holder exponents, scalar weights and replay contract are not asserted here without that additional proof.
