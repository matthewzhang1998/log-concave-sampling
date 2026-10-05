# Independent audit: fractional B mark and paired-rotation covariance

Date: 2026-10-04.

## Verdict and exact sources

**PASS for the conditional fractional-B reference/scaling construction, and PASS for the source-qualified native paired-rotation all-cut/clock estimate.** In particular p=1/2 gives the correct B derivative mark and an A^7 sqrt(D) bound on the new pure-six coarse covariance, with all proper cuts O(A^7), up to the declared public logarithms.

This improves the formerly open coefficient bound. Fresh uniform recalibration also closes the imported source/filter/pair/clipping floors at conservative sqrt(D) scale under changed coarse labels. It does not close the conditional quartic six-force correction, unrelated original clock-target transport, complete public-tilt/current export, or all-order recurrence.

- `../INTERMEDIATE-SCALE-B-MARK-COROLLARY.md`: `070cac5569e9c6a1386bc26c5f98146bc8e25b5e00850cd146950c769e2235c7`.
- `../PAIRED-ROTATION-COVARIANCE-ALL-CUT-LEMMA.md`: `b8383590bdea35e5783c4fe86fdd68551d3d817a1ce81ee1c9004de8964ca668`.

The native derivative-cut and clock inputs are the pinned adapter and native independent audit listed in `INDEPENDENT-GROUPED-LAW-RESPONSE-AUDIT.md`. This audit uses all proper cuts of D_Q A, including root-versus-all-physical, rather than only proper physical cuts of A.

## 1. Fractional B marking and raw paths

With R the retained standard physical input and Cov(P,R)=M=A^(-p)JB, public tilting shifts P by bMtheta. Thus the whole-bank shift is s delta b Mtheta=sbJBtheta when delta=A^p. The first and third odd terms are proportional to epsilon B and epsilon B³: grades four and six, with no amplification of an old cubic kernel.

The orientation matters. Reading a source-dependent generated Gaussian as if it were an unchanged known input would not prove the displayed actual-first ledger. With the declared orientation and caller-independent executed old zero carrier, the source-dependent P path has scale

    A^p · A · A^(1−p) = A².

Its ordinary Gaussian part has scale A^(1+p), and direct old-kernel paths still have scale A. Every old inverse width and every new pair caller factor remains. This is conditional on the complete native pair/path certificate at its actual A^(1−p) radius; the exact Gaussian reference alone does not certify a finite sampler.

The native pair replacement must retain R and hide its discarded pair-private tape. Its P error passes through a 2|alpha|s delta L-Lipschitz downstream map, since the old zero carrier does not depend on q. A pair prior floor A^((1−p)k) sqrt(D) therefore becomes A^[1+p+(1−p)k] sqrt(D). For p=1/2, k>=ceil(2P−3), subject to the native minimum, reaches a fixed target grade P. This is fixed-order work for fixed P; all pair/embedding/filter/query/anchor costs still have to be included.

## 2. Exact Price identity

The pair of rotated banks has equal marginals N(x,s²I) and correlation tau=1−2delta². Its exchangeability implies that E_delta=F(q+)−F(q−) is centered and that its second tensor moment is twice the difference between diagonal and correlation-tau moments. Price differentiation gives

    E[E_delta⊗E_delta]
      = 2s² integral_tau^1 E[sum_a DF_a(x+sG)⊗DF_a(x+sG_t)] dt.

The factor 2, heat factor s² and interval length 2delta² are correct. Only first derivatives occur. Uniform derivative-cut bounds imply the necessary finite-dimensional Lipschitz and Gaussian integrability control; smoothing permits the C1 version. The identity uses the entire shared bank and does not create independent shields.

## 3. All-cut contraction proof

Fix any proper cut of the six external indices. For the two derivative four-tensors, let A_a and B_a be their slices across this cut, indexed by the contracted root a. The output matrix is sum_a A_a⊗B_a. Factoring it as an operator row times an operator column gives

    ||sum_a A_a⊗B_a||
      <= ||sum_a A_a A_a*||^(1/2) ||sum_a B_a* B_a||^(1/2).

The reversed orientation is available as well. The first orientation is valid whenever A has a left external slot and B a right external slot; otherwise the reversed orientation is valid because the global cut is nontrivial. Its factors are exactly two proper local flattening norms. Hence all output cuts are at most K1 K2, including cuts that separate whole old subtrees.

The central three-versus-three flattening is an ordinary matrix product. Consequently its Hilbert-Schmidt norm is at most min(K1 ||B||HS, K2 ||A||HS). Applying these two estimates inside Price's integral gives

    all proper cuts <= 4s²delta² K1²,
    HS <= 4s²delta² K1 h_D <= 4s²delta² sqrt(m) K1²,

where h_D is the actual marginal L2 Hilbert-Schmidt norm of DF. This proves the general tensor step, not merely the examples tested in the checker.

## 4. Native clock completion and its boundary

For one native node m=5D and K1<=C sum_i sigma_i^(-1). Its squared cubic amplitude supplies A^6 w²/sigma_2², with all known shares retained. At s<=1 the one-HS bound is therefore

    Lambda A^6 delta² sqrt(D) w²/sigma_2² (sum_i sigma_i^(-1))².

On each original endpoint clock panel, total positive mass Delta_i is O(sigma_i²). Summing squared node weights is at most the squared panel mass. Using (sum_i sigma_i^(-1))²<=3 sum_i sigma_i^(-2), the worst center factor is Delta_2²/sigma_2^4=O(1) per panel, hence logarithmic in the endpoint cutoff. Every other affected clock is summable with Delta_i²/sigma_i²<=C Delta_i. The remaining clock factors have bounded squared-mass sums. This proves the stated public-log total.

For independently owned node banks, the untilted signed differences are independent and centered, so their covariance tensors add without cross terms. A shared random outer caller must stay conditioned on, or its additional cross-node currents must be handled explicitly. The corollary now states this limitation.

At p=1/2, the bound is Lambda A^7 sqrt(D), with all proper cuts Lambda A^7. This is one dimension-sized Hilbert factor. It is not automatically the optional actual-H-energy refinement, and it is not a bound on every higher-rank term created by public tilting. The inherited conditional quartic remains grade six.

## 5. Uniform recalibration and strict spine

The directly inspected LOW30 import is `external:LOW30`, SHA-256 `7767ac389d6b677eb9517b2585f7a24cb9f85ede47d099a653c635b7be750fb8`. Its scalar-filter/Sobolev statements around lines 5748–5815, pair contract and frozen-label seed proof around 5840–5895, and vertex calibration around 5940–6000 have constants depending on the declared source radius, dimension and fixed degrees, not the value of a frozen coarse label.

For every deterministic q, the anchored Fi,q has Fi,q(0)=0 and private Lipschitz constant at most one. Thus ||Fi,q(G)||Lp<=C_p sqrt(D) uniformly in q. One fixed worst-case absolute calibration tolerance, multiplied by every actual share/readout/inverse-width/ancestor factor, therefore works for every shifted or rotated coarse law. The fixed external OU preparation provides the declared finite list of probe derivatives and moments. Uniform proper cuts and the declared HS upper bound C sqrt(D) likewise permit a uniform analytical clipping threshold. The retained-input pair envelope must still be used in its declared conditional form.

This is a valid new finite-provider recalibration, with all complete occurrences recompiled and tolerances frozen before execution. It is not an upgrade of an old marginal estimate, nor does it prove unrelated original clock-compression errors or a full connected-current export. The positive Price clock sum above is separately established.

The new strict-spine paragraph also checks out for these two three-vertex paths. A hit center has eccentricity one and a hit leaf eccentricity two. Joining hit vertices by one edge gives diameter ecc(v)+ecc(w)+1, namely 3, 4 or 5. Attaching one B leaf to either path instead gives max(2,ecc(v)+1), namely 2 or 3. Each of the nine hit pairs is strictly larger than both corresponding attachments. Coarse injection weights and the actual hit labels remain attached; this statement does not order arbitrary older trees by a pooled diameter.

## Diagnostics

The shared checker has 130 passing fixtures: exact Gaussian-pair tilt and p=1/2 exponents, shifted-caller Price differentiation, its integrated covariance constant, all 62 proper cuts of six-slot contractions for several four-tensor fixtures, a saturating diagonal example, and dyadic squared-mass envelopes, in addition to the main-note checks. It imports the native compiler contracts rather than instantiating that compiler. See `grouped_law_response_checks.json` for source hashes and numerical results.
