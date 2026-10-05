# Optional finite positive implementation of anchored Gaussian smoothing

2026-10-05. An optional finite original-VALUE implementation and transfer theorem. This note does not modify the sealed first or second bridge packets and does not posit an exact smoothed-source oracle.

The strongest result in this packet does NOT execute smoothing. The finite-residual stability and analysis-only mollification argument in BASELINE-SMOOTHING-TRANSFER.md transfers the ideal smooth comparison directly back to the original stage-two graph. The one-node observation in Section 3 explains why no extra smoothing VALUE multiplier is then necessary. The detailed grids below remain useful as an explicit approximation option and an honest backup bill; they are not required for the packet's bounded-block, or fixed-D, original-C2 A^(8/3) theorem.

## 1. Assumptions and what is implemented

Let g=grad U:R^D->R^D, U in C2, g(0)=0, and 0<=Dg(x)<=A I for all x, with 0<A<=1/2. Thus g is C1 and A-Lipschitz. The ideal source, used only in analysis, is

    g_e(x)=E[g(x+e Z)-g(e Z)],  Z~N(0,I_D), e>0.

The executable substitute is a deterministic finite positive sum

    h(x)=sum_(z in Q) p_z [g(x+e z)-g(e z)],
    p_z>=0, sum p_z=1.

All displayed source evaluations are original g VALUES. No expectation is executed. The node/weight version and e are fixed numerical inputs for a requested sweep; if a caller actually varies them or the source/anchor, their true dependencies stay live.

Two distinct certificates are used:

1. The ideal g_e supplies the smoothness required in the stage-two target-bias analysis.
2. The actual h supplies its own exact first/curl/energy/caller ports for guarded native completion. The finite graph's mean is transferred to the ideal graph by a uniform VALUE error bound.

The finite h is generally only C1. It is NOT permissible to apply the smooth stage-two Taylor proof directly to h or to differentiate an HVP.

## 2. Ideal smoothness and structural preservation

Gaussian convolution of the linearly growing g is smooth. Its anchored Jacobian is

    Dg_e(x)=E Dg(x+e Z),  0<=Dg_e<=A I,  g_e(0)=0.

It remains a gradient, of E[U(x+e Z)-U(e Z)-g(e Z).x]. For unit directions u,v,w, centering the bounded Jacobian by (A/2)I in the differentiated Gaussian kernel gives

    ||D²g_e(x)[u,v]||
      <= (A/(2e)) E|Z.u| = A/(sqrt(2*pi)e),
    ||D³g_e(x)[u,v,w]||
      <= (A/(2e²)) E|(Z.u)(Z.v)-u.v|
      <= A/(sqrt(2)e²).

For the second inequality E[((Z.u)(Z.v)-u.v)²]=1+(u.v)²<=2. These are induced Euclidean multilinear bounds, independent of ambient dimension. Consequently valid stage-two declarations are

    B_e=1/(sqrt(2*pi)e),   C_e=1/(sqrt(2)e²).

The source assumptions alone imply

    sup_x ||g_e(x)-g(x)|| <= 2 A e sqrt(D).

If g acts blockwise in any fixed orthogonal decomposition E_1+...+E_m with d_j<=b, both g_e and EVERY fixed finite h act in the same blocks. This is algebraic: the j-th output depends only on P_j x, even if the smoothing nodes are not rotationally invariant. For g_e the projected Gaussian on E_j is standard, so the same B_e,C_e apply on each block. A general source is covered analytically by a single block b=D, retaining the stage-two dimension dependence in K_(D,B_e,C_e,A) and C_D.

## 3. A global, deterministic clipped-grid coupling

The following construction has stronger VALUE control than a bounded-query-region approximation. It requires no Gaussian query-tail truncation.

For a requested scalar RMS quantization tolerance tau>0, put

    tau0=min(tau,1),
    R=2 sqrt(log(2/tau0)),
    n=ceil(R/tau0), a=R/n,
    q(z)=nearest point to clip(z,[-R,R]) in {-na,...,na}.

Use the symmetric tie convention, which affects a Gaussian null set only. The weights of k a are exact standard-Gaussian cell masses, with the endpoint cells extending to infinity. For -n<k<n they are

    w_k=Phi((k+1/2)a)-Phi((k-1/2)a),

and w_n=1-Phi((n-1/2)a), w_-n=w_n. The D-dimensional list is their product, with

    K_D=(2n+1)^D
       <= [3+4 tau0^(-1) sqrt(log(2/tau0))]^D.

The node weights are positive and have mass one. Let Z have iid standard coordinates and QZ=(q(Z_1),...,q(Z_D)). The scalar clipping error satisfies

    E[(|Z_1|-R)_+²]
       =2[(1+R²)barPhi(R)-R phi(R)] <= exp(-R²/2).

For R>0, the first inequality follows from R barPhi(R)<=phi(R) and 2 barPhi(R)<=exp(-R²/2); R=0 follows directly. Minkowski with the rounding error gives

    sigma=(E|Z_1-q(Z_1)|²)^(1/2)
          <= a/2+exp(-R²/4)<=tau0<=tau.

The independent coordinate errors are centered and identically distributed. Therefore

    E[(Z-QZ)(Z-QZ)^T]=sigma² I_D.

For EVERY fixed rank-d orthogonal projection P, including an unknown block projection,

    W2(Law(PZ),Law(PQZ)) <= sigma sqrt(d) <= tau sqrt(d).

The same statement holds after identifying P's range isometrically with R^d. This uses isotropic ERROR COVARIANCE, not rotational invariance of the finite law. A nontrivial finite support in D>=2 cannot be invariant under every orthogonal rotation. The grid itself has signed-coordinate/permutation symmetry and scalar covariance; that is enough here.

For the actual h formed from this list,

    sup_x ||h(x)-g_e(x)|| <= 2 A e tau sqrt(D).             (3.1)

If the true source has blocks E_j, the stronger simultaneous bounds are

    sup_x ||P_j[h(x)-g_e(x)]|| <= 2 A e tau sqrt(d_j).     (3.2)

Proof: under the displayed coupling, subtract the two anchored integrands. Its j-th difference is at most 2Ae||P_j(Z-QZ)||; integrate, apply Cauchy--Schwarz, and then sum squared block bounds. Equation (3.2) avoids incorrectly paying an ambient sqrt(D) error separately in every block. Neither (3.1) nor (3.2) depends on the query x or its radius.

The construction is deliberately sufficient, not cardinality-optimal. In particular the one-node rule QZ=0 has centered error covariance I and tau=1, and gives h=g exactly. The companion BASELINE-SMOOTHING-TRANSFER.md proves the required source-level weak transfer, admitting analysis-only smoothing with no extra original-VALUE multiplier. One must not claim that the exponential-D grid is necessary for that theorem. The tensor-grid bill describes this explicit increasingly accurate smoothing implementation, not a lower bound on all implementations.

For a desired absolute vector error eta>0, take tau=eta/(2Ae sqrt(D)). For a desired error eta_norm sqrt(D), take tau=eta_norm/(2Ae). The explicit inverse-tolerance exponent in generic source VALUE count is D, up to the displayed logarithm. When constant tau suffices under the companion weak-transfer theorem, this explicit grid multiplier is constant in A, though exponential in D; the one-node option removes even that multiplier. This note imposes no inverse-A tolerance choice without specifying the transfer theorem. For an arbitrarily requested vanishing uniform VALUE tolerance it is not a polylogarithmic-cost smoothing oracle.

## 4. Finite-precision positive weights and nodes

The exact Phi expressions above specify a mathematical rule; they need not be hidden exact-real oracle leaves. Here is a finite rational implementation with the same type of certificate.

Build the preceding grid at tolerance tau/2. One may choose a rational R at or above the displayed threshold, and rational a=R/n with n large enough for the stated rounding bound. This preserves all tail estimates. Approximate its symmetric scalar probabilities by nonnegative rational weights w'_k with exact mass one and

    TV(w,w')<=rho,   rho<=tau²/(16R²).

For example obtain certified rational lower bounds l_k<=w_k for k=1,...,n, with total positive-side deficit at most rho/2. Set w'_k=w'_-k=l_k and w'_0=1-2 sum_(k=1)^n l_k. This has exact mass one, is symmetric, and has TV at most rho. Known one-dimensional Gaussian integrals can be enclosed to any prescribed absolute precision by finite quadrature with explicit Gaussian tails; these setup arithmetic and precision costs are charged. Zero weights may be dropped.

A symmetric maximal coupling between w and w' moves a scalar node by at most 2R with probability rho. Composing it with the original Gaussian/grid coupling, independently in each coordinate, gives

    (E|Z_1-Q'Z_1|²)^(1/2)
       <=tau/2+2R sqrt(rho)<=tau.

The composite coordinate coupling can be symmetrized under (Z,Q') -> (-Z,-Q'); hence its mean error is zero. Its independent product again has error covariance sigma_actual² I. Thus (3.1)--(3.2) hold for the ACTUAL rational positive product weights. Exact rational products retain total mass one; subsequent arithmetic roundoff is a separately restored absolute floor. Do not truncate products and then silently use the exact-mass ports.

This method does not need a source-dependent modulus of Dg. The constants depend only on A,e,D,tau and known scalar Gaussian tails.

## 5. Known and unknown bounded blocks have different execution bills

### 5.1 Unknown block basis

Even if the analysis declares block size at most b, without an available block basis/oracle the safe implementation in Section 3 uses K_D=(2n+1)^D original full-source shift nodes. Its blockwise error still scales as sqrt(d_j), by (3.2). Its COST exponent remains D, not b. Blockwise smoothness in the proof does not magically expose the blocks to the executor.

### 5.2 Known block basis

If an orthonormal basis within each block is supplied as an actual numerical input, one common b-dimensional product rule suffices, where b=max_j d_j. For each b-tuple z, insert its first d_j coordinates into block j and concatenate these block vectors into the full-D shift v(z). Use the b-dimensional product weight of z. Every block marginal is exactly its d_j-dimensional product rule, since the unused coordinates sum out. As the source acts blockwise, dependence between block shifts does not affect any component expectation. Therefore the per-macro-VALUE node count is

    K_known=(2n+1)^b,

using the original full-D source oracle at x+e v(z) and e v(z). The simultaneous block errors (3.2) still imply (3.1); the ideal Gaussian and error couplings may be considered separately on each block. This construction does not need, and does not have, independent random shifts across distinct blocks. Its global error covariance need not be isotropic; its correct block marginals are sufficient.

For different block-specific tolerances/rules, separate block evaluation instead costs sum_j K_j, or a common-uniform cumulative-mass coupling produces at most 1+sum_j(K_j-1) concatenated nodes. These are fallback variants, not an unavoidable factor m for the common-grid rule.

Dense rotations, projection/storage arithmetic and basis precision must be charged. This is not an unpriced discovery of a hidden basis. The known-block option changes the execution interface relative to the sealed basis-free bridge graph.

### 5.3 General source

Set b=D in the ideal stage-two analysis and use K_D in execution. All dimension-dependent ideal bias and pair-quadrature constants remain explicit. No dimension-independent source qualification or arithmetic complexity is claimed.

## 6. Exact actual ports; why derivative closeness is unnecessary

For fixed nodes, weights, and e,

    h(0)=0,
    Dh(x)=sum_z p_z Dg(x+e z),   0<=Dh<=A I.

The potential sum is convex. When the original g is identically zero, every original VALUE is zero. At x=0, h(0)=0 is an algebraic cancellation, but the original leaves g(e z) are generally NONZERO. Reuse the same original VALUE records or cancel identical complete argument keys to preserve literal macro zeros numerically; do not erase the anchor leaves or their live caller dependencies. The nonzero original anchor g(e z) is not discarded.

The stage-two first/curl/energy proof used only anchoredness, the PSD Jacobian interval, positive clock weights and the actual conditional Gaussian rows. It therefore applies directly to h without its possessing D²h. In particular the actual residual has private dimension 5D, half-variance normalized first <=9A, curl <25A², declarations ell_E=9A, a_seed=3A, padding mu=A, and the actual captured caller origin and caller bound from the sealed ports derivation. A<=1/36 only discharges the displayed first-radius condition; every other native compiler guard still must hold at the actual parameters and dimensions.

Do not use an HVP of g_e as the derivative of h. A finite sum of shifted C1 sources generally remains C1: a continuous non-Lipschitz Jacobian cusp survives at one or more shifted node locations. More strongly, no source-uniform closeness of the finite sum Dh to Dg_e follows from a uniform VALUE approximation under these assumptions. A bounded continuous one-dimensional Jacobian can have narrow bumps at every finitely queried location, making the finite derivative average large while its Gaussian average is arbitrarily small. Integrating such a Jacobian yields an admissible C2 potential. Thus ideal derivative bounds and finite VALUE error cannot justify transplanting a smoothness certificate to the actual h.

Instead complete the actual h graph against its own mean using its exact ports, then compare only means by (3.1). This avoids any stability claim for entire compiled sample paths or compiler HVPs.

## 7. Transfer through the literal stage-two graph

Write F2[f] for the finite stage-two raw graph with every original source occurrence replaced by f, keeping identical finite rules and Gaussian rows. The sealed absolute VALUE propagation estimate applies to any two A-Lipschitz sources f,h with global uniform difference eta:

    sup_roots ||F2[h]-F2[f]|| <= L_F(A) eta,
    L_F(A)=7/2+14A+(11/2)A².                              (7.1)

This includes shifted/nested descendants; it is not a comparison only at unchanged arguments. Positive outer averaging has mass one and does not enlarge it. Hence the own-mean difference is at most L_F(A)eta in every caller value and in integrated L2. For residual E=F2-B, the corresponding bound is (L_F+1)eta; separately executed caller-origin subtraction doubles that worst-case residual bound. The target-mean comparison after correctly restoring baseline and origins is still the direct (7.1) bound.

Take f=g_e, eta=2Ae tau sqrt(D). The ideal stage-two bias theorem gives its stated K_(b,B_e,C_e,A), delta1,deltaL,deltaQ,delta_out terms. Complete the ACTUAL h graph with the native programs at its actual ports. Its source energy is still O(A²)(|z|+sqrt(D)), so its completion allowance remains Lambda_comp A^4 sqrt(D), subject to all native guards. Combining own-mean completion, (7.1), and the ideal bias yields an entirely finite positive law whose error to N(m3[g_e](Z),I) is at most

    {K_(b,B_e,C_e,A) A^4
      +delta1 A²+sqrt(11/8)deltaL A³+(deltaQ/2)A²
      +delta_out A C_F2+Lambda_comp A^4
      +2 L_F(A) A e tau} sqrt(D)
      +restored absolute floors.                         (7.2)

For a source with unknown blocks, b is its declared block bound and Section 3 supplies simultaneous projected error; the generic implementation still pays K_D. With no block qualification b=D. The pair rule normalization C_b and its b-dependent node bound remain those of the sealed theorem. Positive smoothing weights preserve source structure; positivity of the final law additionally requires the actual native positive variance shares, fills and buffers. It does not follow from smoothing weights alone.

### 7.1 Naive whole-source comparison to the original m3

For two globally A-Lipschitz sources with uniform difference eta, couple their genuine-history constructions. Their first-history means differ by at most eta, their second-history terms by at most (1+A)eta, and their terminal forces by at most

    L_m(A)eta,  L_m(A)=1+A+A².

The positive outer R1 operator is a contraction, so this is also an m3 bound. Consequently a valid additional term for comparing (7.2) to m3[g] is

    2 L_m(A) A e sqrt(D).                                (7.3)

This particular whole-source transfer contains A e and A^4/e². For fixed b its best generic power balance is e proportional to A, resulting in an O(A²) certificate. It does not prove an improved uniform C2 power or order four. This observation is only about this naive whole-source comparison. The sharper baseline-preserving and finite-residual weak transfer in BASELINE-SMOOTHING-TRANSFER.md gains an extra A and supplies the analysis-only A^(8/3) theorem. Its conclusion is not a consequence of the crude bound (7.3), and the crude bound is not an obstruction to it.

## 8. Complete original-VALUE, HVP and replay bills

Let

    M_E=5+3J1+5JL+6JQ,
    L_macro=Q_captured,macro+N_out N_B+M_E N_out N_E
                +Q_additional_macro_replay.

These are the fully expanded native occurrences at the NEW actual normalized ports. They include all caller captures, filter/clock/pair/mark/mode banks and replays; known scalar arithmetic is accounted separately. Replacing each macro h VALUE by a literal 2K-original-VALUE anchored sum gives the always-safe bill

    Q_original_VALUE
       <= Q_rule/grid/certificate_setup
           +2K L_macro+Q_known/numerical,

where K=K_D for basis-free execution, or the stated known-block count if that interface is authorized. Setup here denotes all charged setup work; the source-independent grid probability construction itself needs no original g VALUE. Any source/basis fitting adds its actual original-VALUE cost.

If g,e and the complete node/weight version are fixed, one may compute

    c_Q=sum_z p_z g(e z)

once as a genuine reusable capture, then h(x)=sum_z p_z g(x+e z)-c_Q. The improved original-VALUE count is

    K N_anchor_versions+K L_macro+Q_other_original_VALUE,

with N_anchor_versions=1 only when every expanded occurrence really shares that complete key. Anchor capture storage and replay are charged. Changing a source, e, node, weight, exterior source anchor or coefficient version invalidates the relevant cache; a discarded primal pays full replay. A globally fixed anchor does not excuse the separate stage-two caller-origin capture, whose descendants include nonzero x-dependent arguments.

For one direction v, the exact x-first is

    Dh(x)v=sum_z p_z [Dg(x+e z)v].

It costs K original potential HVPs at the K recorded shifted VALUE sites; symmetry supplies the same adjoint cost. The fixed c_Q anchor has zero x-tangent. A safe fully live sweep through the unreduced anchored sum uses at most 2K original HVP sites per macro VALUE whenever both shifted and anchor arguments depend on the caller. Any original source-parameter derivative service remains separately charged. No HVP is differentiated. No unqueried Hessian action or smoothed expectation is treated as free.

For multiple requested directions, count every HVP action unless the original oracle actually supports and charges a batch. Changed inputs or versions replay the appropriate full primal graph before a derivative sweep. Native compiler repetitions multiply these original leaf costs exactly as they multiply macro occurrences.

Per-original-leaf vector arithmetic is O(D). Creating/storing the generic shift list takes O(DK) scalar data in the explicit-list implementation; streaming can reduce storage but not the number of original leaves. Product-weight arithmetic/bit lengths, known Gaussian mass enclosures, original VALUE/HVP costs, block rotations, buffers, filters and numerical certification all remain visible costs. There is no dimension-independent arithmetic assertion.

## 9. Query radii, original oracle precision, and tails

The smoothing VALUE error (3.1) is GLOBAL. It does not require bounding the Gaussian inputs to a compiled source, nor adding a missing bad-event term. This is why the coupling construction is preferable here to a net restricted to ||x||<=X.

Nevertheless execution can query the original source farther out. At a macro input ||x||<=X, every original smoothing VALUE/HVP lies within

    ||x+e z||<=X+e R sqrt(D),
    ||e z||<=e R sqrt(D).

For block-isolated implementation the block radius is X_j+e R sqrt(d_j); the concatenated known-block list again has full norm at most eR sqrt(D). Any radius-dependent source oracle cost or precision certificate must be evaluated at these enlarged radii.

For one raw stage-two occurrence, include Z,G,N,M,L1,L2 in the 6D-dimensional root vector H. Unit Gaussian rows and A<=1/2 imply that all macro source inputs in its stencil have norm at most

    (1+5A+A²)||H|| <= (15/4)||H|| <4||H||.

Indeed ||S||<=(1+A+A²)||H||, each d has norm <=2A||H||, and the quadratic terminal S+/-d1+/-d2 is the largest bound. Thus on ||H||<=sqrt(6D)+t all original leaves have radius at most

    4(sqrt(6D)+t)+eR sqrt(D).

The standard Gaussian radial tail is <=exp(-t²/2). This is only a raw-occurrence radius bound. It is not silently reused for the fully expanded native program, which has additional roots, banks and transformations. If original oracle accuracy is certified only on bounded sets, the actual complete compiler's root dimension and source-argument growth bound must be used, and its out-of-region growth/error must be bounded before any tail term can be asserted. Merely assigning probability to the good region does not control an arbitrary inaccurate oracle outside it.

For original VALUE errors bounded by nu at all actually requested shifted and anchor sites, positive mass-one averaging gives an h VALUE floor at most 2nu, plus weighted arithmetic/node errors. This is not Knu unless a particular arithmetic implementation introduces that accumulation. Feed the resulting macro floor through (7.1), through the residual/caller capture floor when appropriate, and through every native numerical allowance. Analogous weighted HVP floors use the actual original HVP tolerance. Never divide a floor by smoothing error, A, measured residual energy, or a quantity that may vanish.

## 10. Exact scope and references

Established under the explicit global C2-potential/PSD-Jacobian assumptions:

- A finite positive anchored original-VALUE approximation to Gaussian smoothing with explicit global and unknown-block-projected errors.
- An honest dimension/tolerance/radius and original-VALUE/HVP/replay bill, including a rational positive-weight implementation.
- Exact C1-level actual source ports and a value-only transfer to the analytically smoothed stage-two target theorem.
- A finite guarded own-mean completion conditional on the imported native services and ALL their actual guards, not a new audit of those services.

Not established:

- A D²-Lipschitz certificate for the finite approximation.
- A free exact smoothing, original Hessian-modulus or hidden-block oracle.
- An arbitrary-accuracy dimension-free quadrature bill for an unknown block basis.
- Uniform order-four accuracy for general C2 sources, or arbitrary-order closure.
- This optional note does not reprove the sharpened baseline/residual weak transfer established in BASELINE-SMOOTHING-TRANSFER.md. That companion compares coherent source substitutions with an extra A factor and permits the one-node analysis-only route. The crude absolute floor (7.1) is not a lower bound and does not require tau to vanish in the companion proof.

Primary internal inputs read for this derivation:

- ../nonlinear-bridge-stage-two/STAGE-TWO-POSITIVE-VALUE-THEOREM.md, especially Sections 5--8.
- ../nonlinear-bridge-stage-two/SHIFTED-BIAS-AND-PORTS-DERIVATION.md, especially Sections 3--9.
- ../nonlinear-bridge-repair/NONLINEAR-ANCESTRY-BRIDGE-REPAIR.md, especially its complete bill and qualitative C2 boundary.
- The prior full resummed-mean bill and original-source/HVP conventions in /workspace/shared/log-concave-research-github-20261005-v12/research/cost/endpoint-stein-quadrature-20261004/positive-law-repair/higher-cumulant-gate/same-carrier-feedback/gaussian-rms-backbone/GAUSSIAN-RMS-RESIDUAL-RESUMMED-M3.md (the exact snapshot read for this note).
