# An all-rank active VALUE source and a smoothing-grade admission ledger

2026-10-05. This note supplies a bounded, dimension-sharp, curl-free active source at every primitive rank. It does **not** claim that the resulting arbitrary decorated graph queue is a sealed native positive sampler. The precise remaining consumer contract is identified below.

## 1. Exact active source, including its origin

Let g=grad U on R^D, g(0)=0, and 0<=Dg<=A I, with A>0. Only g VALUE is executed. For k>=0, set a=(k+1)^(-1/2). Let x be one private vector, P_1,...,P_k active physical probes, z a captured caller, and t>0 a private coefficient width. Define

 f_k(x;P,z,t) = (A t)^(-1) 2^(-k) sum_eps (product_j eps_j)
   [g(z+a t(x+sum_j eps_j P_j))-g(z+a t sum_j eps_j P_j)].       (1)

Every sum is on its **same** common bank. There are exactly 2^(k+1) original VALUE calls before any exact duplicate caching. The two sites of a summand have the same captured center. No tensor, derivative, conditional mean or cumulant is an input.

The following hold pointwise:

- f_k(0;P,z,t)=0, exactly.
- |f_k(x)|<=a|x|.
- ||D_x f_k||op<=a, and D_x f_k is symmetric.
- ||D_(x,P_1,...,P_k) f_k||op<=1.
- ||D_z f_k||op<=1/t.

The improvement from 2/t to 1/t uses the Hessian sandwich: for symmetric 0<=H,H0<=A I, ||H-H0||op<=A. For one signed summand, the private derivative is aH/A and the j-th probe derivative is a eps_j(H-H0)/A. The block operator norm is at most a sqrt(k+1)=1. The triangle inequality over the normalized signs preserves this bound. With a caller z=z(y,Q), insert its actual derivative rows; no law estimate is differentiated. A joint caller/private bound is sqrt(1+||D z||op^2/t^2) when the caller coordinates and active coordinates are disjoint.

A scalar potential, used only to prove curl and not queried, is

 Phi_k(x) = [A a t^2]^(-1) 2^(-k) sum_eps product eps_j
   [U(q_eps+a t x)-U(q_eps)-a t x.dot g(q_eps)],
 q_eps=z+a t sum_j eps_j P_j.

Then grad_x Phi_k=f_k. Signs need not preserve convexity, but do preserve the genuine-gradient property. For X standard Gaussian independent of the probes, E|f_k(X)|^2<=D/(k+1), uniformly in the actual caller and probes.

For k>=1 the field is odd in **each** P_j separately. Thus averaging any one independent symmetric probe gives exactly zero. This is not a conditional-zero assertion after that probe has been exposed or tilted. For k=0, the private mean need not vanish. Exact origin zero is distinct from exact Gaussian mean zero. Finite arithmetic/calibration floors remain real floors.

## 2. Exact Gaussian response and all proper cuts

Let X,P_1,...,P_k be independent standard Gaussian D-vectors. Keep z fixed during this integration. Define the response tensor only analytically, in the displayed slot order:

 T_k(z,t) = E[f_k(X;P,z,t) tensor X tensor P_1 tensor ... tensor P_k].

Gaussian integration by parts, first in X and then in the P_j, gives

 T_k(z,t) = a^(k+1) t^k/A * E[D^(k+1) g(z+t Z)].             (2)

Here Z is standard D-Gaussian. Every active row has coefficient a t, so their combined variance is (k+1)a^2 t^2=t^2. The same-center anchor contributes zero to the X response. The sign from each probe differentiation cancels its Rademacher sign. Formula (2) is valid in weak Gaussian derivatives for the original C2 class. The tensor is fully symmetric because g is a gradient.

This is a response of an executed bounded-first vector field after native Gaussian integration. The algorithm does **not** execute f tensor X tensor P_1 ... P_k as a coefficient estimator. That distinction is exactly why the finite-VALUE first-slot-rank obstruction does not apply.

The products X_i P_(1,j1)...P_(k,jk) form an orthonormal family. Bessel's inequality therefore gives the one-energy estimate

 ||T_k||HS <= a sqrt(D).                                      (3)

For every p|q proper flattening, p+q=k+2 and p,q>=1,

 ||T_k||_(p|q,op)
  <= a^(k+1) 2^(k/2) sqrt((p-1)! (q-1)!).                    (4)

To prove (4), split tZ=(t/sqrt(2))(Z1+Z2). Leave one index on each side on the matrix Dg. Integrate the remaining p-1 derivatives against Z1 and q-1 against Z2. Contract a unit HS test tensor on each side to vector-valued Hermite scores. Their squared L2 norms are at most (p-1)! and (q-1)!, respectively. They are functions of independent banks. Pointwise ||Dg||op<=A, followed by Cauchy-Schwarz separately on the banks, gives (4). In particular the all-cut constant

 c_k = a^(k+1) 2^(k/2) sqrt(k!)

is explicit, dimension-free, and computable. This proof does not estimate a self-trace by a proper-cut bound.

For k>=1 and linear g(x)=H x, the entire executed source f_k is identically zero. The cancellation occurs before any tensor formation or norm. It is stronger than a zero target with a noisy tensor estimator.

## 3. Globally bounded preparation preserving exact private curl

A naive scalar cutoff multiplying f need not preserve curl. Use a known radial pullback instead. Fix R>0. Write r=|x| and psi_R(x)=h_R(r)x/r, with psi_R(0)=0. Define

 h_R(r)=r,                                      0<=r<=R,
 h_R(r)=R+R[v-v^3+v^4/2], v=(r-R)/R,            R<r<2R,
 h_R(r)=3R/2,                                  r>=2R.

Its derivative on the middle interval is 1-3v^2+2v^3, in [0,1]. Thus ||D psi_R||op<=1, |psi_R|<=3R/2, and a conservative direct radial calculation gives ||D^2 psi_R||op<=8/R. The map is C2 across the seams.

Execute

 F_(k,R)(x;P,z,t) = D psi_R(x)^T f_k(psi_R(x);P,z,t).          (5)

Only the same 2^(k+1) original VALUE calls and known radial arithmetic are needed. Execute the radial Jacobian as c I+d uu^T acting on a vector, in O(D) work/storage; no dense D-by-D matrix is needed. It is exactly grad_x[Phi_k(psi_R(x))], and its exact origin remains zero. Its uniform radius and complete first bounds are

 |F_(k,R)| <= 3aR/2,
 ||D_x F_(k,R)||op <= 13a,
 ||D_(x,P) F_(k,R)||op <=13,
 ||D_z F_(k,R)||op <=1/t.                                   (6)

Divide by 13 if a unit complete private/probe first is required, and multiply the known selected-response normalization by 13. A radius of natural vector size sqrt(D) is explicit; it is not hidden as a dimension-free uniform norm.

The source equals f for |x|<=R. For R=sqrt(D)+s, s>=0, Gaussian concentration and Cauchy-Schwarz give

 ||F_(k,R)-f_k||L2
 <= (5a/2) [D(D+2)]^(1/4) exp(-s^2/8) =: e_(k,D,R).         (7)

Indeed outside the ball the difference is at most (5a/2)|x|. The bounded response T_(k,R) obeys

 ||T_(k,R)-T_k||HS <= e_(k,D,R),
 every proper cut <= corresponding bound (4)+e_(k,D,R).     (8)

Bessel's inequality gives the first assertion; a flattening operator norm is at most HS. Taking

 s^2 >= 8 log(max(1, (5a/2)[D(D+2)]^(1/4)/epsilon))

makes each absolute response floor <=epsilon. This logarithmic dimension/tolerance choice is a real preparation bill. The exact oddness in every active probe survives (5), but approximate finite values can leave a mean floor.

## 4. VALUE precision, caller precision, retention and replay

If each executed original g VALUE has Euclidean absolute error <=delta_g, the source error before the optional /13 normalization is <=2 delta_g/(A t). The normalized sign average has l1 weight one, so there is no hidden 2^k precision amplification. The pullback has operator norm at most one. Choose delta_g<=epsilon_value A t/2. For A=0 the source is identically zero and no division is executed.

Perturbing captured z by delta_z costs at most |delta_z|/t; perturbing the concatenated private/probe inputs costs at most 13 times their Euclidean error. For the width, a direct differentiation with all roots frozen gives |partial_t F| <= a[2|psi_R(x)|+(sum_j |P_j|^2)^(1/2)]/t (use x in place of psi_R(x) for the unbounded field). This is an actual scalar geometry bound; finite interval propagation can use its supremum over the width interval. Radial and remaining scalar arithmetic require their own certified error propagation, including the actual query magnitudes. Width is a frozen program parameter in (6). Differentiating a Gaussian encoding of sampled clocks is a different first ledger and is not free.

A retained caller is preserved by evaluating the source at that same caller. Its Gaussian response averages X and the P_j. That identity is not a joint-law replacement if one of those variables is exposed as a retained observer. Storing a tape internally for exact replay does not expose it in the law. The all-rank consumer must declare which probes remain physical readouts, which are private, and which observer currents it retains.

For one smoothed tree, an explicit response-floor telescope is also available. Write b_k=a^(k+1), the non-width response factor in (2). Preparing each local bounded response to HS floor epsilon_v<=delta b_(k_v) makes its unnormalized derivative-tensor error <=A delta t_v^(-k_v). All other local cuts are at most A t_v^(-k_v)[2^(k_v/2)sqrt(k_v!)+delta]. For 0<delta<=1, telescoping the tree and its positive history sum gives the conservative total source-preparation error <= n S_n 2^(2n)(n-2)! delta A^n tau^(-(n-2)); there is no extra dimension factor in this absolute floor. Thus choose delta<=min(1,epsilon_prep tau^(n-2)/[n S_n 2^(2n)(n-2)!]) for normalized allowance epsilon_prep A^n. Split each epsilon_v across the actual preparation errors: for example radial response floor <=epsilon_v/3 using (7), original VALUE-response floor <=epsilon_v/3 by delta_g<=epsilon_v A t_v/6, and scalar-arithmetic response floor <=epsilon_v/3. The sum, not each separate floor, is what enters the telescope. Optional /13 source normalization is restored in the selected coefficient and does not change this telescope. All corresponding old calls must be rerun at these actual floors. This controls source-response preparation, separately from finite pair/filter errors. The exact first/curl claims concern the literal smooth VALUE field (1)/(5). A mere absolute VALUE rounding bound does not certify a Sobolev-first or curl error for an arbitrary discontinuous numerical oracle. An actual finite native implementation must retain its separately certified source/filter Sobolev and arithmetic-first floors; none is inferred from the response bias above.

One literal source call costs 2^(k+1) original VALUE calls, k+1 D-dimensional private/probe rows if not supplied by the caller, and O(k 2^k D+D) scalar/vector arithmetic with O((k+1)D) live storage, excluding the original oracle itself. The radial Jacobian action is linear in D. If a native pair/filter wrapper repeats it M_(k,b,epsilon) times, its bill is 2^(k+1)M_(k,b,epsilon), plus that wrapper's exact captured-call/readout work. We do not replace this by the cost of one ideal tensor.

A changed center, probe, private root, width, source order, requested floor, or original-g version requires replay of every affected original call and anchor. Caching has zero incremental cost only for the exact same complete version. Original HVPs may be used in a separately requested derivative sweep at recorded VALUE sites; none is used to execute (1) or (5).

## 5. Full structural endpoints are finite under genuine coefficient smoothing

This source gives a constructive alternative to moving high derivatives onto the deepest structural clock. Smooth every original force at width tau>0, keeping the literal ancestry. Split its total primitive heat so that every executed source width obeys

 t_i^2=(1-r_i^2+tau^2)/2 >=tau^2/2.

The additional coarse heat is included explicitly in the same bank. Every structural clock retains its original positive density r^(j-1)dr on the full interval (0,1). No factor (1-rho^2)^(-(n-2)/2) is introduced. Thus the sibling Hermite endpoint problem is absent from this representation, at the declared price of coefficient heat mismatch.

For n>=2 and the exact generated rank-n cumulant trees, let S_n be their summed integer weights:

 S_1=1;  S_n=sum_(i=1)^(n-1) binom(n,i)i(n-i)S_i S_(n-i).

The imported marked-tree heat bound is B_n=2^(n-2)(n-2)!. Each tree has n forces, 2n-1 positive clocks and total inverse private-width degree n-2. Put C_n=S_n B_n. Then

 all cuts <= C_n A^n tau^(-(n-2)),
 HS <= C_n A^n tau^(-(n-2)) sqrt(D).                         (9)

A union of structural endpoint sectors r>1-delta has absolute coefficient mass at most (2n-1)delta times the right side of (9). This is an optional truncation certificate, not a claim that the tail is zero. The full endpoint already has finite mass. An identical uniform envelope applies to primitive endpoints because of tau.

For the original conditional cumulant, the proved Appell estimate is

 ||kappa_n(H_tau)-kappa_n(H)||L2(HS)
 <= n! A^n tau sqrt(D).                                    (10)

A same-caller dyadic positive resolvent rule with scalar multiplier tolerance delta_0 has the imported complete quadrature allowance in L2 over the original standard Gaussian retained caller Y; this is not a uniform arbitrary-caller quadrature assertion:

 D_n delta_0 A^n tau^(-(n-2)) sqrt(D),
 D_n=(2n-1)2^(2n-1)C_n.                                    (11)

Thus this section settles full clock integrability of the **smoothed main coefficients**. It does not erase the substantive heat error (10), return feedback, or observer costs.

For n=8: there are 794,880 unmerged histories, S_8=32,049,561,600, B_8=46,080, and fifteen clocks per history. The source at a degree-d force vertex is k=d-1; the rank-eight star requires k=6 and 128 VALUE calls per raw source occurrence. Constants are large rank constants, not inverse-alpha exponents.

## 6. A computable root-heavy grade ledger, with its missing hypothesis explicit

Here is a useful sufficient allocation test. It is **conditional** on an unsealed finite marked-tree native comparison port described below. It must not be reported as a proved native rank-eight return.

Normalize widths against a fixed positive physical buffer and write alpha for the small force parameter. Choose

 tau=alpha^(1/n), beta=1/[2(n-1)],
 rho_nonroot=alpha^beta,
 rho_root=K_h alpha^(n-1/2).

K_h is the complete known scalar history/clock/readout/source-normalization coefficient, including product t_i^(-(d_i-1)), actual signs and positive clock weights. Its worst width exponent is n-2. Consequently the root's effective alpha grade is

 r_n=n-1/2-(n-2)/n=n-3/2+2/n.                              (12)

The exact product of amplitudes remains K_h alpha^n. The root's own caller may pay one more inverse tau, giving grade r_n-1/n. All nonroot caller paths have the root and their own positive beta exponent before their inverse tau. Complete path constants and node counts are included, rather than inferred from an ideal main coefficient.

The main coefficient first-hit envelope is C'_n alpha^n tau^(-(n-1)). An old endpoint with source-zero bank row zero and complete first O(alpha) therefore has the prospective one-hit allowance

 alpha^(n+1)tau^(-(n-1)) = alpha^(n+1/n).                  (13)

The heat mismatch (10) has exactly the same grade alpha^(n+1/n). If a conditional native return is quadratic in the actual root module **and has a width-uniform prefactor after all source/caller/frame costs have been charged**, it would begin at alpha^(2r_n); for every n>=4,

 2r_n > n+1/n.

At rank eight beta=1/14, tau=alpha^(1/8), r_8=27/4; the potential own-return grade is 27/2, root-caller grade 53/8, and prospective joined/heat grade 65/8. These are exact rational arithmetic facts, not an assertion that every descendant is root-quadratic or width-uniform. If the comparison port contributes an additional tau^(-ell), its actual own-return grade is 2r_n-ell/n. At rank eight that is 27/2-ell/8; using two full root-caller firsts already corresponds to ell=2 and grade 53/4. Admission at grade n+1/n requires ell<=n^2-3n+3 (ell<=43 at rank eight), with every further inverse-width factor included. The missing port must output that complete ell and its constant, not merely the word quadratic.

The missing hypothesis is a full arbitrary C_k marked-tree finite comparison port: it must show that after the selected main response all conditional/public/auxiliary descendants are either the exposed intended current or are quadratic in that root module, with dimension-sharp extended derivative frames, an explicit complete inverse-width exponent ell satisfying the stated admission test, and no uncharged linear old-bank/observer term. It must also provide the same-endpoint current extraction, exact covariance gaps, all finite means, and terminal signed-twin/curl return of the **completed native sampler**. The private source curl theorem (6) is not that terminal theorem. Fixed-order local C_k existence alone is insufficient.

Once this port and its inverse-width grade test are supplied, (9)-(13), a fixed untouched Gaussian keep, the literal first-hit rule and full replay give a bounded rank-n admission with error Lambda_n(alpha,D) alpha^(n+1/n)sqrt(D), plus separate numerical floors. Here Lambda_n must retain the actually computed rank/readout/source constants and public-log factors from the alpha-dependent finite clock list; no uniform constant is asserted before that computation. The current paper supplies the active-source part and its explicit geometry, but does not assume the missing port true.

## 7. What an induction-ready cost record must contain

The source definition and all-cut constants are an all-rank rule. The full native return is still a typed obligation. For every emitted module record:

- original force count, primitive k at every vertex, force degree, physical rank, surviving marked outputs, and cycle/ancestry labels;
- the exact common-bank read set, active probes, retained observers, source-zero carrier rows, positive reserve and actual covariance gaps;
- all scalar clock weights, t powers, amplitude exponents, numerical constants and first/caller path sums;
- complete selected-pair/filter order and every precision floor, including finite means and root/clock arithmetic;
- all first conditional, public-tilt, auxiliary and old-bank descendants at the same endpoint;
- source-level and completed-terminal gradient/curl certificates separately;
- complete VALUE queries, roots, arithmetic and nonincremental replays.

If M_n(delta_0) is the maximum nodes per clock, the unmerged main-source census obeys the explicit bound

 Q_main(n) <= T_n M_n(delta_0)^(2n-1)
   max_h sum_(v in h) 2^(d_v) M_native(d_v-1,b_v,epsilon_v)
   + Q_captures+Q_scalar_roots+Q_readouts+Q_replay,

where T_1=1 and T_n=sum_i i(n-i)T_i T_(n-i). The degree identity sum_v(d_v-1)=n-2 gives the sharper raw VALUE factor sum_v 2^(d_v)<=2^(n-1)+2(n-1) per history, attained by a star. At rank eight this is 142 calls for one complete history before wrapper replication, captures and replay. No M_native value is invented here. The original all-rank finite filter/pair implementation must supply it at the actually needed b and epsilon.

The exact summed raw source census over unmerged histories has the independent recurrence V_1=1 and

 V_n=sum_(i=1)^(n-1) [(n-i)(i+1)V_i T_(n-i)
                         +i(n-i+1)T_i V_(n-i)].

The convention V_1=1 is a combinatorial starting leaf of degree zero, not a claimed C_(-1) primitive. When a derivative picks a vertex, 2^(d_v) doubles; summing those increments over every left/right hit gives the recurrence. At n=8, V_8=27,020,800 raw original-g VALUE calls over all 794,880 histories at one clock tuple per history, before any pair/filter replication. This large finite rank constant is distinguished from the alpha-dependent clock grid and wrapper order.

For a genuine return queue, define Q(v)=Q_main(v)+Q_old_complete(required version)+sum_children Q(child)+Q_actual_replay. This becomes an algorithm only after an effective-grade invariant makes the child list finite and all widths/constants explicit. Exponential or factorial rank constants do not by themselves affect the inverse-alpha exponent. Conversely an alpha-dependent b, inverse smoothing width, minimum native radius or Monte Carlo replication can affect that exponent and must be retained.

Nothing above proves c(P)=o(P), arbitrary-order cancellation, or closure of all active-probe descendants. It isolates a constructive, all-rank, dimension-sharp source layer and an exact rational admission test that a stronger return theorem can consume without repeating the old raw-tensor or sibling-history obstruction.

## 8. Compositional all-cut theorem for open marked trees

This strengthening was checked jointly with the shared-variable induction work. It is a tensor-network theorem, independent of the missing native comparison port.

Let a finite connected tensor tree have one tensor A_v at each vertex. Every graph-theoretic leaf vertex must retain at least one external slot. Suppose every proper flattening of A_v has operator norm at most c_v. Then every proper flattening of the contracted tree tensor has norm at most product_v c_v. If one selected A_w has HS norm h_w, the full tree has HS norm at most h_w product_(v!=w)c_v.

For the all-cut claim, fix a nontrivial input/output partition of the external slots. Give every input terminal positive source mass and every output terminal positive sink mass, balanced in total, chosen so no tree edge has zero net flow. Such a choice is explicit: label terminals j=0,...,m-1, put w_j=2^(2^j), and let S_in and S_out be their respective sums. Input terminal j has mass w_j S_out; output terminal j has mass -w_j S_in. Each proper component cut off by an edge contains at least one terminal, and its complement does too. The net flow polynomial is a signed sum of distinct powers 2^(2^i+2^j) over input/output pairs; its largest nonzero term cannot cancel. Thus no edge flow vanishes.

Orient each tree edge with its unique conserved flow. Include terminal legs, directed inward for inputs and outward for outputs. Conservation ensures that every force tensor has at least one input and one output, so every vertex is used through a proper cut. The directed underlying tree is acyclic. Its contraction can therefore be evaluated as composition and tensor products of local bounded operators, with identity wires for carried slots. This proves the product bound. It does not require positive entries or independent tensor coefficients.

For the HS statement, root the tree at w and attach vertices away from it. Every newly attached vertex has one edge to the built group and either a future tree edge or an external slot, so its contraction uses a proper cut. Repeated HS-times-operator multiplication preserves one HS factor. The result holds pointwise in any common old coefficient bank; subsequent integration must still obey its actual ownership/observer rules.

There is a useful extension of the **HS statement only**. A loopless connected multigraph with a spanning tree whose leaves are all marked admits the same bound: root at a marked tree leaf and attach along that spanning tree, contracting all edges to the already built group at once. The future tree child or marked output keeps the new local cut proper. This argument does not prove all proper cuts of the full cyclic graph. Nor does it permit a local self-trace: an opened edge whose two endpoints land in one initial selected-spine chunk needs a new certificate.

Consequently (2)-(4) supply a dimension-sharp all-cut input for arbitrary open marked-tree chunks, including bank-hit trees whose surviving leaves are C1 rather than C0. Their positive native return, derivative-frame propagation through all nonlinear finite programs, and cyclic re-closure are still separate obligations.
