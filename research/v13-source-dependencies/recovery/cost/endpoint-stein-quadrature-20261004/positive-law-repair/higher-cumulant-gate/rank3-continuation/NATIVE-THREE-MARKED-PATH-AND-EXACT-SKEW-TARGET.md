# Native three-marked path and the exact skew target

2026-10-04. Constructive source and normalization specification. The original-VALUE source admission and leading tensor identity below are explicit. The complete one-energy feedback estimate in Section 6 remains a port to audit before calling this a finished positive cubic-law packet. This note does not treat an analytical coefficient as a producer input.

## 1. Exact compressed target for the continuous OU displacement

Let g=grad U, and let H=integral_0^1 g(X_r)dr for the stationary OU path conditioned on X_1=Z. Write v=E[H|Z], C=Cov(H|Z), and kappa_3 for its conditional third cumulant. Let L be the standard OU generator in Z, and B=Dv. The discounted-additive-functional log generating equation gives

    (1-L)v=g,
    (2-L)C=2 B B^T,
    (3-L)kappa_3=6 Sym_3(B D C).                       (1)

Here (B D C)_ijk=sum_a B_ia partial_a C_jk. For g gradient, B is symmetric. Expanding the log generating function as sum_k lambda^(tensor k):kappa_k/k! verifies the factors 2 and 6. These identities use the continuous target directly; they do not cube a finite path quadrature.

The resolvent rule (k-L)^(-1)h=integral_0^1 r^(k-1)P_r h dr therefore gives

    C=2 integral_0^1 s P_s(B^2) ds,
    kappa_3=24 integral_0^1 q^2 dq integral_0^1 s^2 ds
          P_q Sym_3[ B . P_s((DB).B) ].                (2)

In indices the tensor inside the last expression is

    sum_(a,b) B_ia(Y) E[(partial_a B_jb)(X) B_kb(X)|Y],
    Y=qZ+sqrt(1-q^2)G0,
    X=sY+sqrt(1-s^2)G1.

The symmetrization is the standard average over permutations. The factor 24 comes from 6 in (1), 2 in C, and the two product-rule terms in D(B^2), which coincide after physical symmetrization.

Moreover

    B=integral_0^1 r P_r Dg dr,
    DB=integral_0^1 r^2 P_r D^2g dr.                   (3)

The second line is a heat derivative at positive clock, not a producer query or a C3 assumption. Expanding the three slots in (2) using (3) produces five positive scalar clocks and a leaf-center-leaf source-index tree with three original g occurrences: Dg, D^2g, Dg. Each force has one physical mark. The D^2g slot is privately smoothed before interpretation.

Scalar normalization check: for the polynomial fixture g(x)=x^2, v=(x^2+2)/3, C=2(x^2+1)/9, and (2) gives

    kappa_3(x)=16x^2/45+32/135.

It solves (3-L)kappa_3=6v' C' exactly.

Equation (2) is an economical exact source for a native rank-three tree. Differentiating the executed square-action C_sq is unnecessary and would create a forbidden producer HVP path. The original gradient VALUES at the three tree vertices are sufficient source primitives, if the complete native-tree law/feedback port below is verified.

## 2. A directly derived double-Riesz version for an admitted square gradient

For clarity, let G:R^n->R^n be an actual square-gradient VALUE source, radius ell, and centered energy e. For a scalar test direction theta write h=theta dot(G-EG). Gaussian Riesz is R h=integral_0^1 P_s Dh ds. Then

    kappa_3(h)=2 E[Dh . R((Rh).Dh)].                   (4)

Expand the derivative inside the outer R only after replacing its unexpanded Hilbert operator by the admitted positive finite Riesz rule. There are exactly two histories:

    2 integral dr ds E[
       DG(V) -- D^2G(W_r) -- DG(U_(r,s))
       +s DG(V) -- DG(W_r) -- D^2G(U_(r,s))],          (5)

with physical symmetrization and

    W_r=rV+sqrt(1-r^2)Z1,
    U_(r,s)=s W_r+sqrt(1-s^2)Z2.

Dashes mean the two source-index contractions of the three-marked path, not matrix multiplication with an unqueried tensor. Gradient symmetry permits the physical mark to be placed at each appropriate original vertex.

Approximating either unexpanded Riesz operation at Hilbert-operator tolerance delta perturbs kappa_3 by at most C delta ell^2 e. This follows before differentiating products: ||R G||2<=e, multiplication by DG costs ell, and the second R and DG cost another ell. Choosing delta=ell yields the required order-four coefficient floor ell^3 e. There is no need to approximate the original Hessian pointwise.

## 3. Positive independent shields at every original query

For each interior (r,s), the three Gaussian query covariance is the scalar block matrix

    K(r,s)=[[1,r,rs],[r,1,s],[rs,s,1]] tensor I_n.

Let

    sigma^2=min(1-r^2,1-s^2)/12.

Then K-sigma^2 I is positive definite. For example the Markov innovation factorization bounds ||K^(-1)||op by 12/min(1-r^2,1-s^2). Thus draw correlated coarse query roots Q from the KNOWN covariance K-sigma^2 I, and add one independent sigma-Gaussian shield per original force occurrence. This preserves the exact joint Gaussian queries in (5). The derivative tensors in the conditional native target are now genuine positive-heat coefficients.

All coarse roots stay in one declared original-query record and are averaged after the whole native packet. Replacing them by independent roots would change the target.

For a node of finite quadrature weight w, the positive dyadic rules give, up to their declared logarithmic multiplicities,

    sum w/sigma <=Lambda,
    sum w/sigma^2 <=Lambda.                           (6)

The second bound is logarithmic in the finite Riesz heat cutoff. The individual sigma may be a small fixed power of ell; it affects precision, not an inverse-heat empirical bank count. The same inequalities hold for the history weighted by s. A more general five-clock instance from (2) must carry the corresponding exact Gaussian query covariance and its own shield decomposition, rather than assume that it has precisely K(r,s).

## 4. Literal VALUE vertex and amplitude normalization

At a shielded original query Q, use the anchored normalized square-gradient source

    F_Q(u)=[G(Q+sigma u)-G(Q)]/(ell sigma).            (7)

It has fresh radius at most one, literal zero at u=0, and coarse-Q first at most 2/sigma. G(Q) is a captured source call at the identical complete key. Every changed Q rebuilds that complete source. If G itself is a finite original-gradient graph, all original ancestors are charged.

For a leaf, the native C0 source is F_Q((x+z)/sqrt(2)); multiply it by sqrt(2)/s0. Its selected-pair target is exactly

    H_leaf=(E DG at that query)/ell.

For the center, use t=1/sqrt(3) and the literal native C1 source

    C1(x,z;p)=[F_Q(t x+t z+t p)-F_Q(t x+t z-t p)]/(2t).

Apply the admitted finite first-chaos filter in p, with its fixed external OU/Sobolev preparation, and multiply by 1/(s0 t). Its selected-pair target is

    H_center(p)=(sigma/ell) E D^2G[p].                 (8)

Only the actual finite bounded filtered field is fed to a pair kernel. The unbounded exactly linear target in (8) is used for coefficient calibration, never inside a sampled covariance square root.

Take leaf and center pair-source amplitudes

    rho_leaf=rho_center=sqrt(ell).

Allocate a positive packet variance nu. Choose readout coefficients c=a=b=sqrt(nu)/2, and an independent final Gaussian buffer of variance nu/4. The actual packet reads

    c Y_root+a P_leaf+b P_center+sqrt(nu)/2 Zbuffer.

At the zero source its covariance is nu I. The root-source amplitude for a history whose kappa_3 coefficient is 2w is

    rho_root= w ell^2/(3 c a b sigma).                 (9)

For the second history replace w by sw. Known native normalization constants are exactly those already included in (7)-(8).

Indeed the product amplitude is rho_root rho_center rho_leaf=w ell^3/(3cab sigma). The native reverse symbol's one-cluster term is cab times this product times the normalized tree. Since the normalized tree equals (sigma/ell^3) times the original tree, its logarithmic third coefficient is w/3 times that tree. Equations (4)-(5) require exactly kappa_3/6=(2/6)sum w tree. Thus (9) has the correct sign and exact factor.

This is an explicit source scaling, not a request to evaluate kappa_3. A negative desired third current is produced by reversing the known root source sign, never a probability weight.

The admissibility check is on the ACTUAL source radii

    Lambda sqrt(ell),
    Lambda |w| ell^2/(nu^(3/2) sigma),                 (10)

and on the total variance allocation. Dyadic weights become smaller at small sigma, so small clocks do not force the second quantity large. Equal positive variance shares across the finitely many clock/history packets produce only public-log inverse-share losses. Leave a fixed independent numerical buffer after all packets are added.

Each node's original-query caller first is its source amplitude times Lambda/sigma. Strict-ancestor amplitudes are then multiplied along its actual path to the root. For the root this gives w ell^2/sigma^2; weighted summation is logarithmic by (6). Original external caller rows receive their actual additional G-caller factors. No derivative of a completed pair law is used to infer these actual first bounds.

The chronological pair replacement error for order b is bounded by

    Lambda sqrt(n) [rho_root^b
                    +rho_root rho_center^b
                    +rho_root rho_center rho_leaf^b]. (11)

For b>=4, the last two terms sum to order ell^4 sqrt(n) or smaller using (6). A larger fixed b can make these separate absolute floors negligible. They must not be divided by a realized e; the final statement must retain the absolute prior/numerical floors when e is small.

## 5. Literal original-g lifts when the source is a path force

Do not assert that a continuum path integral is itself an admitted gradient source. A useful alternative is to integrate mixed third cumulants at three literal times. Conditional on Z, encode the time triple by

    X_i=r_i Z+s_i R_i W, s_i=sqrt(1-r_i^2), R_i R_i^T=I_D,

where W has the fixed Gaussian dimension needed for the three correlated path sites. For each original occurrence define

    G_i(W)=R_i^T[g(r_i Z+s_i R_i W)-g(r_i Z)]/s_i,
    P_i=s_i R_i.                                      (12)

For interior times s_i>0, this is a genuine square-gradient VALUE source. Its derivative is

    DG_i=R_i^T Dg(r_i Z+s_i R_i W)R_i,

so radius<=A, centered energy<=C A sqrt(D), zero at W=0, and captured-Z first<=C A r_i/s_i. Its physical readout obeys

    P_i G_i=g(X_i)-g(r_i Z), P_iP_i^T=s_i^2 I_D.

The explicit constants disappear from mixed cumulants. Thus there is no 1/s_i physical readout, and no claim that a completed law-only mean is a native gradient.

The marked-row extension of the native reverse identity is literal: read out P_i times the corresponding ambient Gaussian root or public. In the reverse calculation its public shift is b P_i^T theta. The leading physical tensor has P_1,P_2,P_3 on its marks. Its known Gaussian variance is scalar; fill it to nu because each P_iP_i^T<=I. Original source blocks and actual native queries remain full square gradients throughout.

For the continuous H target, integrate over ordered triples rather than cube a finite H_Q. With r3=a,r2=ab,r1=abc the ordered time measure is 6a^2 b da db dc. The mixed-cumulant integrand is the true conditional Gaussian Markov triple. Real-node OU interpolation on dyadic panels gives a positive polylogarithmic clock rule for each nested moment/product. A conservative L6-product estimate has size A^3 D^(3/2); to obtain delta A^3 sqrt(D) accuracy, assign internal operator accuracy delta/D. This dimension loss is paid by log(D/delta) nodes, not silently called a one-energy estimate. Anchors in (12) are for the actual source only; use the OU moment identities when proving clock approximation, since the anchor terms cancel exactly from cumulants.

This direct mixed-time route and the five-clock resolvent route in (2) are alternative target representations. They should not be combined with a finite-path cube that restores a diagonal defect.

## 6. First port still requiring a complete proof

The positive finite native program, literal source guards, exact leading third coefficient, retained physical caller and prior replacement estimate have been specified above. A complete law theorem additionally requires the native packet's own feedback and owned-coarse-root averaging estimate in these actual nonuniform clock/amplitude scales.

The expected one-energy bound has the form

    Lambda sum_(histories,nodes) w^2 ell^5 e/sigma^4
       + absolute calibration/prior/numerical floors,

which is logarithmic at the dyadic clock boundary and is much smaller than ell^3 e. This is a target bound, not a proved inequality in this note. It must be established from the exact native reverse symbol and same-endpoint Hilbert-Riesz proof, including the finite first-chaos filter errors and all retained public/owned-coarse-root rows. Counting only the intended cubic coefficient does not establish it.

The old-root derivative, same-source one-energy mark and positive variance allocation must all survive that argument. In particular a globally Lipschitz generic f does not automatically admit (7): the native mean/pair source class requires a square gradient or the stated symmetric partial-gradient blocks. The explicit lifts (12), or an already audited full-gradient lift, are essential source qualifications.

After this port is closed, the independent-buffer positive quadratic-reference interpolation can compare the packet to the original buffered skew law through order four. That comparison needs the mean/covariance services at the same absolute grade. It still would not prove repeatability: the completed cubic packet is an actual VALUE program with its own residual first/curl, not a new original gradient. Re-entry must use its literal source qualification or a new valid lift.

## Source pin

The finite native primitives used here are LOW30, sections `Finite value filters and unexpanded Riesz clocks`, `A literal native vertex and complete tree program`, `Exact reverse-forest coefficients and finite feedback`, and `The marked whole-law remainder at the actual endpoint`, in research-source/High Acc Ideas/ai-bucket/30_low_acc.tex. The displayed primitives have no producer HVPs. HVPs occur only in requested first/adjoint sweeps at their recorded original-gradient VALUE sites.

## 7. Alternative all-ell vertex radii for a conservative re-entry guard

The square-root-ell allocation in Section 4 is optional. Taking

    rho_leaf=rho_center=ell,
    rho_root=w ell/(3cab sigma)

leaves the exact product and leading third coefficient unchanged. Dyadic positive weights satisfy sup w/sigma<=Lambda, so all original native source radii are Lambda ell under the same positive-share guard. The root coarse-caller sum is Lambda ell sum w/sigma^2<=Lambda ell. A fixed pair order b>=4 places the root prior floor at order ell^4 sqrt(n), and all descendant priors are attenuated further.

Only the root pair output is read in the final packet. The leaf and center marked publics are raw independent Gaussians and belong to the known source-zero carrier. Under the square-root-ell allocation, their pair residuals still reach the output through all strict ancestors, so a bare square-root-ell residual first cannot be inferred merely from a nonroot source radius. Under either allocation the complete radius must be proved by the actual path sums. The all-ell choice gives the simpler conservative O(Lambda ell) residual/caller return and avoids relying on a sharper attenuation claim before audit.

No small relative curl at that new residual's own first scale is asserted. Re-entry at an unchanged parent radius and re-entry as a genuine original gradient are different claims.
