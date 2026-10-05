# Finite cumulant-tree generator and the first rank-five native return test

2026-10-05. This gives an exact all-rank analytical generator and a falsifiable native-extension test. It does not promote an unadmitted contraction graph to an original-VALUE source.

## 1. A concrete finite generator for the ORIGINAL gradient cumulants

Take g=grad U, g(0)=0, 0<=Dg<=A I, and g_tau=E g(.+tau Z). At one repeated physical theta, let h=theta.g_tau. Represent a history by:

- a labelled force leaf g(v;E_v), where E_v is its list of spatial derivative-edge labels;
- a positive resolvent node R_k(child);
- a product node with its actual common caller;
- a positive integer coefficient.

Start κ1=R1 h. At rank n use the exact connected recurrence

  κn = Rn sum_(i=1)^(n-1) binom(n,i) Dκi . Dκ_(n-i).

Relabel the two subhistories to have disjoint force and old-edge labels. Add one new edge label. Apply one derivative to each subhistory: it chooses EACH force leaf once, increments that leaf's derivative order, and increments by one EVERY resolvent index on its path from the root. This is D R_k=R_(k+1) D plus the product rule. Contract the two new derivative ports with the same new edge label. No derivative index is symmetrized with a physical output. Polarization at the end recovers the fully symmetric physical cumulant.

The implementation generate_rank_terms.py stores these exact labelled histories, without relying on cancellation or isomorphism merging. At rank n each term has n original force occurrences, n-1 spatial edges, exactly 2n-1 positive clocks, derivative order d_v at vertex v equal to its tree degree, and total extra heat derivatives sum_v(d_v-1)=n-2. Every force has one physical mark.

If T_n counts unmerged histories and S_n sums their positive integer coefficients, then

  T1=S1=1,
  Tn=sum_i i(n-i) T_i T_(n-i),
  Sn=sum_i binom(n,i) i(n-i) S_i S_(n-i).

Observed exact rows through rank six:

- Rank 3: 4 unmerged histories, S3=24, five clocks; one scalar-clock history after symmetry.
- Rank 4: 28 histories, S4=672, seven clocks; four scalar-clock histories, agreeing with the sealed 192,192,192,96 packet.
- Rank 5: 272 histories, S5=32640, nine clocks; 13 scalar-clock histories. Degree multisets are 112 paths (1,1,2,2,2), 144 fork histories (1,1,1,2,3), and 16 four-leaf stars (1,1,1,1,4).
- Rank 6: 3312 histories, S6=2384640, eleven clocks.

These are ordered product-rule histories, not counts of unlabelled trees. The rank-five JSON contains the actual AST, integer weight and all edge labels of every history. The source decoder can therefore recover literal ancestor Gaussian sites, rather than drawing unrelated centers for products.

The exact check independently compares the generated scalar expression to raw moments M_n=nR_n(h M_(n-1)) followed by the moment-cumulant recursion, using h=x+x²/7+x³/11 through rank five. This fixture tests finite algebra, not the bounded-Hessian theorem.

## 2. Coefficient smoothing, clocks and explicit rank scaling

At primitive occurrence with derivative order d_v and primitive correlation r_v, set

  sigma_v²=(1-r_v²+tau²)/2,
  f_v(x)=[g(center_v+sigma_v x)-g(center_v)]/(A sigma_v).

Its private first is <=1, its zero is exact, and its caller first is <=2/sigma_v. It is a genuine gradient; every actual query is original g VALUE. Structural ancestors and primitive coarse roots stay on the exact AST. Each occurrence gets a fresh private shield. Added tau heat belongs to that primitive coefficient only and never borrows common structural heat.

The literal native vertex required by a degree-d all-marked tree is C_(d-1). Its selected target is the known scalar normalization times (sigma_v^(d-1)/A) E D^d g. LOW30 explicitly defines these finite VALUE sources for every fixed probe number. However the sealed rank-four return only checked complete feedback/first admission for its C0/C1/C2 occurrences. Rank five requires a new C3 occurrence at its star; an all-rank source formula alone does not certify its complete returned packet.

Heat cuts and tree leaf elimination give a conservative history bound

  proper cuts <= B_n A^n tau^(-(n-2)),
  HS <= B_n A^n tau^(-(n-2)) sqrt(D),
  B_n=2^(n-2)(n-2)!.

For the full positive sum multiply by S_n. Factors 2^(d-1)(d-1)! dominate the fixed-order split-heat Hermite constants at each leaf. Positive R_k are contractions on Gaussian Hilbert L2; their pointwise proper-cut action is bounded by their positive mass <=1.

Use a positive finite dyadic rule with scalar-Hermite multiplier error delta0 at each clock and mass <=2. Telescoping actual clocks, retaining one HS error and pointwise cuts on the other factors, gives

  ||κn,Q(H_tau)-κn(H_tau)||_(L2(Y);HS)
   <= D_n delta0 A^n tau^(-(n-2)) sqrt(D),
  D_n=(2n-1) 2^(2n-1) S_n B_n.                          (T1)

It is an integrated standard-Y statement, not a uniform arbitrary-caller quadrature estimate. For a requested relative tensor allowance epsilon, choosing delta0=epsilon tau^(n-2)/D_n is finite and explicit. Each clock then has O_n(log²(1/delta0)) nodes under the admitted rule. No inverse-tau Monte Carlo replication is implied.

The Appell result in the companion note gives the exact coefficient heat mismatch

  ||κn(H_tau)-κn(H)||_(L2(Y);HS) <= n! A^n tau sqrt(D).   (T2)

With alpha=q A/sqrt(u), q<=1, and a positive reference with buffer comparable to u, Hermite coupling turns (T1)-(T2) into normalized packet allowances C_n alpha^n(epsilon+tau)sqrt(D). Physical LAW multiplies by sqrt(u), and a terminal A-Lipschitz force multiplies again by A. Thus every tau/u/A loss is explicit. This does NOT yet bound native graph feedback or first paths.

One can select, for example, tau=alpha^gamma_n with any specified gamma_n>0 and epsilon=alpha^beta_n. Then the two analytical packet terms are C_n sqrt(u)[alpha^(n+beta_n)+alpha^(n+gamma_n)]sqrt(D). The choice is computable but must be reconciled with actual native guards, every offspring and the required target grade; it is not free accuracy.

## 3. Why the rank-four amplitude prescription does not repeat at rank five

The blind extension puts every nonroot amplitude at alpha and moves the entire normalization to a marked leaf:

  rho_root=d_h omega alpha/(R_h product_v sigma_v^(d_v-1)). (T3)

For a five-star, the product is sigma_center³. On a center endpoint panel with mass Delta and sigma_center² comparable to Delta+tau², ordinary positive weights give the certified scales

  max root radius <= Lambda alpha/tau,
  root-first sum <= Lambda alpha/tau,
  center caller-path sum <= Lambda alpha²/tau².           (T4)

The last expression comes from alpha² sum omega/sigma_center^4. In contrast, the quartic star had sigma_center² and center caller alpha²/tau. Setting tau=alpha, as in the sealed quartic packet, makes all three displayed budgets order one, not O(alpha). Endpoint panels with Delta comparable to tau² attain these powers up to their actual positive readout/node factors. An upper endpoint-mass estimate by itself is not a lower bound for an arbitrary quadrature node; the point is failure of the claimed uniform power certificate.

The guard problem is not caused by a rough non-gradient source. In the scalar stationary ideal graph with g(x)=Ax, the C3 center is identically zero while a normalized C0 leaf is identically one. The center output is then a fresh standard Z_center, and the root output is

  Y_root=rho_root Z_center+sqrt(1-rho_root²) Z_root.

Relative to its prescribed source-zero carrier Z_root, the residual derivative in Z_center is exactly rho_root. A fixed positive root readout multiplies it by that readout. This exposes the unattenuated root first even though the intended fifth cumulant is zero. It is a countertest of (T3)'s stationary/native graph extrapolation, not a claim that no alternative finite compiler can choose a different legal realization.

Moving the weight to a nonroot center, or redistributing amplitudes, can avoid this particular root problem. For example keeping rho_root=alpha and putting d_h omega alpha/(R_h sigma_center³) at the center moves its inverse shield behind a strict ancestor. But this changes the actual native radii, pair-order bills and all side-spine feedback. It is not covered by the quartic theorem. Choosing tau much larger than alpha improves these guards but worsens the substantive smoothing error (T2).

## 4. Exact star offspring, including the lower-rank current

Fix all original coarse roots and use the ideal multilinear target ONLY for finite coefficient identification. For an all-marked n-star, its distinguished root spine has three force occurrences. Let lambda be that spine's amplitude product including known readout factors and clock/shield normalization. The own center public is P; the n-3 side-output noises are independent S_j. Completing each leaf reverse shift and physically tilting the readout gives the formal symbol

  V(theta)=lambda theta² product_(j=1)^(n-2)(G_j+a_j theta). (T5)

Here G_1=P with a_1 the center's physical readout, which has zero force degree; for j>=2, a_j contains one side-leaf force amplitude. All G_j are independent standard Gaussians at this conditional stage. In multiple dimensions, (T5) is a multilinear tensor contraction, with its full coarse-root coefficient retained. The scalar formula is sufficient to refute a missing-current census.

The coefficient is nontrivial in the authorized original-gradient class. In scalar dimension use g(x)=A x/2+(A/(4k)) sin(kx), k>0. Then g(0)=0, A/4<=g′<=3A/4 globally, and its privately smoothed C3 target is (sigma³/A) E g⁽⁴⁾(q+sigma Z)=(sigma k)³ exp(-(sigma k)²/2) sin(kq)/4, which is nonzero at ordinary captured centers. The C0 leaves remain strictly positive. Taking k=1/sigma makes this a numerical-size normalized coefficient while preserving the Hessian sandwich. Thus the offspring census is not based on an inadmissible arbitrary tensor.

No exponential-integrability assertion is made for the unbounded polynomial target. For n>=5 its ordinary real MGF can fail to exist. Every coefficient below is the finite formal moment-cumulant jet, and bounded finite native fields have their own separately calibrated finite moments.

The intended first coefficient is

  E V = lambda product_j a_j theta^n.

The COMPLETE second coefficient is

  κ2(V)/2 = lambda² theta^4/2
      [product_j(1+a_j²theta²)-product_j a_j²theta^(2n-4)]. (T6)

For n=5, with (a_1,a_2,a_3)=(b,u,v), this is exactly

  lambda²/2 [theta^4 +(b²+u²+v²)theta^6
       +(b²u²+b²v²+u²v²)theta^8].                       (T7)

In particular a SIX-FORCE FOURTH current lambda²theta^4/2 remains with u=v=0. The six-force sixth current lambda²b²theta^6/2 remains too. Turning off side crosses does not turn off their Gaussian probes. Reversing the root sign does not reverse either even current. Conditional centering of a correction alone does not remove its nonzero law feedback.

For the first rank-four offspring, take two copies of the three-force root spine. Each retains its root and terminal physical mark. The own-center-public contraction and the two side-noise contractions join their centers by THREE parallel edges. Together with the four within-spine edges this graph has V=6, E=7, beta1=2, and M=4 physical marks. At each center the physical mark has been consumed. Its primitive derivative order remains four. The check

  sum_v(j_v-1)=6=M-2+2 beta1

is exact. Thus even this first extension needs a two-cycle current type, not only an all-marked tree or the rank-six tree already studied elsewhere.

Each intact three-force spine still retains two physical slots, so all fixed conditional Wick coefficients generated by (T5) have a restricted one-Hilbert graph certificate: group the intact spines; each Gaussian color occurs at most once in each initial group, so Wick pairings introduce no initial self-loop. Contract all shared edges between groups at each merge. This does not extend to arbitrary self-traces or to whole-old-bank covariance test norms.

## 5. A finite generator for ALL conditional star descendants

For k>=0 define the polynomial

  mu_k(theta)=product_j sum_(l=0)^floor(k/2)
      binom(k,2l)(2l-1)!! (a_j theta)^(k-2l),
  c_1=mu_1,
  c_k=mu_k-sum_(i=1)^(k-1) binom(k-1,i-1)c_i mu_(k-i).

Then the full formal conditional return through any fixed spine-amplitude cutoff K is

  sum_(k=1)^K lambda^k theta^(2k)c_k(theta)/k!.            (T8)

This lists every rank, multiplicity, sign and side-force occurrence. The checker produces the five-star output through k=6, verifies (T7), and verifies the arbitrary-n k=2 formula through n=9. Forces, clock multiplicities, old-root labels and physical ranks must accompany each monomial. Old-bank and auxiliary-root averaging are further connected cumulant operations on the SAME complete bank, not independent coefficient resampling.

With (T3), lambda scales as Lambda alpha³ omega/sigma_center³. The squared positive-clock endpoint budget is

  sum lambda² <= Lambda alpha^6/tau².                    (T9)

Consequently tau=alpha yields only an alpha^4 conservative feedback allowance. Even tau=sqrt(alpha) leaves an alpha^5 allowance, at the intended fifth-force grade. This is a certified allowance, not a lower bound on W2 for every original source. A slower heat exponent and amplitude redistribution may admit a bounded improvement; no such choice by itself proves all-order closure.

For example balancing ONLY the rank-five heat alpha^5 tau against the allowance alpha^6/tau² gives tau=alpha^(1/3) and exponent 16/3. This calculation is diagnostic, not a sealed rank-five producer: other histories, source first, priors, retained callers and graph returns still require checks.

## 6. Minimal next interface and the complexity answer

The next interface is now specific. It must return a programmable positive original-VALUE realization of the six-force/four-mark two-cycle coefficient in (T7), with its SAME old coefficient roots, C3 occurrence shields and squared original clock weights, and must enumerate all conditional, auxiliary and whole-old-bank descendants. Cutting two parallel edges and sharing two auxiliary Gaussian callers is a candidate opening into a six-vertex tree; both auxiliaries have zero physical readout. This is not automatically admitted: first/radius/clock and every auxiliary/old-root derivative-hit history need their actual certificates.

The imported eight-force one-cycle construction proves that auxiliary edge closure can work for a particular graph. It explicitly leaves public-tilt and whole-old-bank interfaces open. It cannot be cited as an arbitrary graph compiler. Likewise the refuted shield-uniform six-tree test cannot be used to remove the inverse-shield factors in (T9).

A full fixed-target generator would have to carry, for EVERY emitted history:

- original force occurrence count, physical rank, cycle count and distinguished-spine history;
- actual shield/variance/amplitude powers and numerical constant;
- original caller, shared-root and observer read sets;
- finite C_k source, selected-pair order, positive variance and covariance gaps;
- direct actual residual/private/caller first, source-zero carrier and energy certificate;
- a terminal genuine-gradient baseline/curl proof after the COMPLETE join;
- conditional/auxiliary/whole-bank feedback, moving covariance, regression and symmetrization;
- the exact complete query/root/replay recurrence and propagated absolute floors.

No law estimate may be differentiated to fill the actual-first/curl fields. No exact ideal mean-zero identity may erase a finite native mean floor. A rank-only queue is disproved by (T7). A force-only queue is not certified after the width costs in (T9). A computable well-founded effective-grade/history order, or an explicit finite acyclic expanded graph, is still required.

Therefore the present result DOES provide a finite computable all-rank analytical target/current generator, coefficient-only C2 smoothing with explicit rank-dependent choices, all-rank one-mark and proper-cut bounds, and a decisive rank-five admission/offspring test. It does NOT yet prove a finite arbitrary-order original-gradient positive algorithm, a computable all-order A threshold, or c(P)=o(P). Eventual sublinearity remains feasible in the logical sense that these countertests do not rule it out; it is not established, and no numerical constants/thresholds for that claim can honestly be supplied from the sealed interfaces.
