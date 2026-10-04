# Affine native VALUE reverse: reconstructed source/audit

Publication copy: nonmathematical context and/or local paths were sanitized. Original and public SHA-256 values are recorded in `INVENTORY.json`; source/audit pins below identify their historical versions, not these edited bytes.

NEW RECONSTRUCTION. Last verified P audit d3cb3223957228945d7358dbc0dac580340a81f876a3c009586171190dd026bc bound R’s source f1f524e9ea878d551f9ac02846dd221f34e612bd53a50389722b363a7af25b93. Physical wrapper a445101f82e4539ecbd81be41e901038e473b7344043296603315cfb973f11df. These historical pins do not certify byte identity of this reconstruction. General DAG checker previously passed5,400 cases; hostile finite-twin audit e181c2640f20eca1b50f3f88c2c55c9c97a9b09650770793ba4260c02cd69e62 passed84 exact cases.

## Structural input

The COMPLETE fixed finite source graph consists of known affine maps and original gradients g_i(x)=H_i x+b_i, H_i symmetric. Its quadratic status is a supplied premise, not inferred from finite samples. Unroll each actual finite Picard iteration with its saved labels and shared ancestors. No implicit-limit source is substituted. Linear graph coefficients and Hessians are fixed in θ; offsets may move with θ.

For E:R^n→R^d,

    E(W)=J W+e0(θ),  L=P*J,  PP*=I,
    ||J||op≤A, ||L−L*||op≤κ, e=||J||HS.

P is the recorded coisometry, not the proxy B. All actual auxiliary columns remain.

## Exact VALUE reverse

At any fixed original anchor x0,

    H_i v=g_i(x0+v)−g_i(x0).

The same two VALUES give the primitive tangent and adjoint. Cache an identical zero/anchor once. A forward tangent pass uses this identity at every node; a reverse pass accumulates all shared output covectors, applies it once per visited node, then uses the TRANSPOSES of the known edges and input rows. This computes Jv and J*v exactly on the ACTUAL finite DAG. No HVP is evaluated in the sampler and no HVP output is differentiated.

Moving offsets cancel from the centered actions. Under the fixed-coefficient premise J is independent of θ, so these actions have zero caller derivative. Varying Hessians/edges require their separate derivatives.

## Five-pass physical-public orientation

At the same d-dimensional public G used by a covariance core, set x=P*G and execute

    u=Lx, v=L*x,
    T(G)=P[L v−(L u+L* v)/2].

This uses three forward and two reverse actions. It equals O_EG, where

    O_E=JJ*−Sym[(JP*)²]
       =P[LL*−Sym(L²)]P*.

With K=L−L*, LL*−Sym(L²)=−Sym(LK), so

    ||T||L2≤κe,  ||T||Lp≤C_pκe,  Lip_G T≤Aκ.

There are no coefficient-private roots, no mean bias and no caller path in the stated affine class. An n-dimensional Gaussian complement is unnecessary; all reverse n-vectors are actual work vectors whose storage is charged.

Adding−αT to the same-public forward core yields the exact extra covariance words

    −qO_E,
    α²(MO_E+O_E M*),
    αc(M²O_E+O_E(M*)²),
    α²O_E²,

withα=q/(2b), c=α²/(2b). When||M||≤A² their one-energy prices are A²κe, A4κe, Aκ²e. The complete original E mark is used, not termwise energies of F/H.

## Closed polynomial covariance family

Apply K_phys=JJ* by J* then J, without assembling a matrix. For fixed v,q and qA²/v≤1/2, define

    p_m(x)=Σ_(j=0)^m binom(1/2,j)(−x)^j,
    Y_m=√v p_m((q/v)K_phys)G.

Every nonconstant Taylor coefficient is negative, so p_m(x)≥√(1−x)≥1/√2 on[0,1/2]. The program is exactly Gaussian with a numerical covariance gap. A same-G coupling gives

    W2(Y_m,N(0,vI−qJJ*))≤C_(m,q,v)A^(2m+1)e.

Its actual carrier-subtracted energy is C_p A e, its first is C A² and its caller is zero. Every higher cumulant is exactly zero. Increasing m repairs the next known matrix word. Work is O(m Q_E) original VALUES on the full actual graph, plus known arithmetic. Five orientation actions cost at most six Q_E including one primitive-anchor cache pass; discarded caches are replayed.

For coefficient rounding require Σ_j |δc_j|A^(2j)≤ε/(√v√d), which bounds the added linear-Gaussian L2 error byε. Original source evaluation errors pass through the actual finite path/readout majorant. This is not a matrix-learning or inverse-A replica requirement.

The hostile actual twin has J=c e3*, P=e2*, L²=0, true covariance c², and T_phys=c²G. Its complete finite K matters: the audited fixture gives E=0 at K=2 and cZ at K=3. The reverse must preserve that version. Its conservative exact counts were E=14 VALUES, reverse28, cached origin14, T112 and degree-m polynomial14+42m.

## Nonlinear boundary

For nonlinear g_i, g_i(q_i+v)−g_i(q_i) is not H_i(q_i)v. Reversing that VALUE increment changes both coefficient and queries. Neither derivative convergence nor a source-aware joint calibration is supplied here. The closed affine family is genuine but restricted.
