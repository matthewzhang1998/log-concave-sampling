# The midpoint opposite source has no universal fixed-known-readout full-gradient lift

Publication copy: nonmathematical provenance wording and/or local paths were sanitized. Original/public SHA-256 values are recorded in `INVENTORY.json`; historical source/audit pins refer to the original versions.

New scoped theorem, 2026-10-04. This strengthens the actual-V-frame result in OPPOSITE-BASELINE-FULL-TAPE-AND-LINEAR-RETURN-GATES.md. The proof below includes its same-potential local-jet and padding gates.

## Exact scope and conclusion

Fix d>=2, M=I, r,epsilon>0, 0<a<=1/16, 0<h<=1/4, c^2+s^2=1 with s>0, and the complete Euclidean tape (W,V,P), where W=(S,U,Z), V is the d-dimensional auxiliary Gaussian block, and P is any finite unread padding. The known numerical geometry may determine a fixed matrix B, but B does not depend on the primitive, the root, or primitive oracle responses. Require BB*=I.

For each smooth anchored SAME primitive g=grad V0 with .4I<=Dg<=.6I, define

    x0=cS+sU,
    Delta=g(S)-g(S+epsilon Z),
    x1=x0+a Delta,
    tbar=S+(a/2)[g(x0)+g(x1)],
    N_h(W,V)=(kappa0/(2h))[g(tbar+hV)-g(tbar-hV)],
    kappa0=ra.                                                   (1)

N_h ignores P. There is no such universally fixed B for which every (1) has an exact representation

    N_h(W,V)=B grad Psi_g(W,V,P)

with a full first bound depending only on the displayed numerical geometry, dimension/padding, and original Hessian bounds. This remains false if each Psi_g is allowed arbitrary additional root/padding-dependent terms. In fact the proof forces B=[0, beta I, Q], after which beta!=0 gives an arbitrarily large full first over the same-potential class, and beta=0 permits no globally bounded full first for one explicit member of the class.

This is **not** a theorem against adaptive or source-specific frames, a law-only replacement, variable frame fields, a non-Gaussian graph selector, or a new weak joint constrained-law generator. An affine ambient lift evaluated only on its nonlinear feature graph also lies outside the claimed exact standard-tape representation. The theorem addresses the proposed universal fixed-known-coisometry entry to the existing small-full-gradient service.

It also does not automatically exclude a finite sum of different projected-gradient offspring, none individually equal to N_h, or a decomposition only of its required covariance current. The single directional-symmetry condition cannot be applied to such a sum without first proving a valid single common-gradient/readout representation for that sum.

## 1. Necessary directional symmetry

If N_h=B grad Psi_g for fixed B, then

    (D N_h) B* = B (D^2 Psi_g) B*

is symmetric. The same necessary condition holds when grad Psi_g is merely Lipschitz: it follows by distributional equality of the directional mixed derivatives, and the left side is smooth.

At one actual record put

    Tplus=Dg(tbar+hV), Tminus=Dg(tbar-hV),
    A0=Dg(x0), A1=Dg(x1), F=Dg(S), Fplus=Dg(S+epsilon Z).

Write B=[B_S,B_U,B_Z,B_V,Q] and define d-by-d matrices

    X=B_S*,  Y=cB_S*+sB_U*,  Z0=B_Z*,
    K=(D_W tbar) B_W*
      =X+(a/2)(A0+A1)Y
         +(a^2/2)A1[FX-Fplus(X+epsilon Z0)],
    L=h B_V*.                                                   (2)

The symmetry condition is exactly

    Tplus(K+L)-Tminus(K-L) is symmetric.                         (3)

There is no derivative of a Hessian in (2)-(3). These are original first/adjoint analytical identities, not executed HVP-valued programs.

## 2. Independent Hessian variations really are available for one potential

Start with g_base(y)=y/2. At a generic fixed W,V, the six points

    S, S+epsilon Z, x0, x1, tbar+hV, tbar-hV

are nonzero and pairwise distinct. For this linear base source, they are known linear maps of W,V. No two of these maps agree identically when a,epsilon,h,s>0: the terminal maps have their distinct V coefficients, the ancestor maps have an sU coefficient, and the two feedback/ancestor twins differ by nonzero Z coefficients. A finite union of proper linear subspaces cannot exhaust the root space. Thus a generic record exists for every nondegenerate numerical geometry.

Around each point q_i choose a disjoint small ball avoiding the origin. Add to the **same scalar potential** an independently chosen compact perturbation

    phi_i(y)=(1/2)(y-q_i)* E_i (y-q_i) rho((y-q_i)/delta_i),     (4)

where E_i is symmetric and rho is a fixed smooth compact cutoff equal to one near zero. Then

    grad phi_i(q_i)=0,  D^2 phi_i(q_i)=E_i.

Every original gradient VALUE at every one of the six points remains its base value. Consequently Delta, x1, tbar and both terminal queries remain **exactly fixed**, not merely fixed to first order. All six Hessians change independently by E_i.

There is a finite constant C_rho,d, independent of the ball radii, with

    sup_y ||D^2 phi_i(y)|| <= C_rho,d ||E_i||.

Choose the E_i in a sufficiently small open symmetric-matrix neighborhood of zero. Disjoint supports and the base Hessian I/2 ensure .4I<=Dg<=.6I globally. The supports avoid zero, so anchoring is preserved. Therefore (3) holds with all six Hessians varying independently in open neighborhoods while its literal same-potential record and B remain fixed. Open neighborhoods suffice for every linear-algebra step below.

## 3. Symmetry forces all original-root blocks of B to vanish

Use two elementary facts for d>=2:

- If H R is symmetric for every symmetric H, then R is a scalar multiple of I. Indeed H=I makes R symmetric and the remaining identities say R commutes with every symmetric matrix.
- If H R is a scalar multiple of I for every symmetric H, then R=0. Set H=I and then choose a nonscalar symmetric H.

First independently vary Tplus and Tminus in (3). Differences show that H(K+L) and H(K-L) are symmetric for all symmetric H. Hence both K+L and K-L, and therefore K and L, are scalar multiples of I. In particular

    B_V=beta I                                                   (5)

for a fixed real beta.

Now fix the other Hessians and vary A0. Formula (2) shows that H Y is scalar I for every symmetric H, so Y=0.

With Y=0, vary A1. The second fact gives

    FX-Fplus(X+epsilon Z0)=0                                    (6)

for every independently available F,Fplus. Vary F in (6); it gives H X=0 for all symmetric H, hence X=0. Equation (6) then gives epsilon Fplus Z0=0. Strong positivity and epsilon>0 imply Z0=0. Since Y=0 and s>0,

    B_S=B_Z=B_U=0.                                              (7)

Thus the only universally possible fixed-readout form is

    B=[0, beta I, Q],  QQ*=(1-beta^2)I,  |beta|<=1.              (8)

The classification uses d>=2. It does not silently apply the scalar commutant argument in d=1.

## 4. A midpoint SAME-potential fixture rules out the remaining frames

Use the scalar smooth family from the companion note, with 0<a<=1/16 and 0<h<=1/4:

    chi(u)=exp(1-1/(1-u^2)) on |u|<1, zero outside,
    eta=1/80, zeta=h/480,
    g_k,1(y)=y/2+(eta/k)(1-cos(ky))chi(8y)
                       +zeta chi(4(y-1-h)/h),  k>=48.

Embed it into dimension d by

    g_k(y)=(g_k,1(y1),y2/2,...,yd/2).

It is the gradient of one anchored smooth potential, with .45I<=Dg_k<=.55I globally. Fix original roots

    S=e1, U=-(c/s)e1, Z=0.

Here x0=x1=0 and, crucially for the midpoint source, tbar=e1. Along varying U1 while S,Z remain fixed, x0=x1 still coincide and

    tbar_1=S1+a g_k,1(cS1+sU1).

Thus this is exactly the earlier one-ancestor potential-difference test on the literal midpoint source, not a different query genealogy. At the displayed root,

    tbar_U1=a s/2,  tbar_U1U1=a s^2 eta k,
    g_k,1(1+h)+g_k,1(1-h)-2g_k,1(1)=h/480,
    g_k,1'(1+h)+g_k,1'(1-h)-2g_k,1'(1)=0.                      (9)

### beta != 0: padding cannot hide the forced potential difference

Let q=B(V,P) be the known coisometric row coordinate and complete it orthogonally by coordinates u. In the original variables, moving q by t e1 at fixed u moves V by beta t e1 and P by t Q*e1. The equation B grad Psi=N_h becomes grad_q Psi=N_h.

Choose fixed u,q0 giving V1=0. Move q from q0 to q0+e1/beta, so V1 moves from zero to one while the original roots W stay fixed. Since N_h does not read padding,

    Psi(W,q0+e1/beta,u)-Psi(W,q0,u)
       =(1/beta) integral_0^1 N_h,1(W,v e1) dv
       =(kappa0/(2 beta h^2))
          [P_k(tbar_1+h)+P_k(tbar_1-h)-2P_k(tbar_1)],          (10)

where P_k'=g_k,1 analytically. No potential VALUE is queried.

If grad Psi has full first at most L, the root gradient of the difference on the left of (10) is 2L-Lipschitz. Its exact U1/U1 derivative at the fixture is, by (9),

    kappa0 a s^2 eta k/(960 beta h).

Consequently

    L >= kappa0 a s^2 eta k/(1920 |beta| h).                  (11)

This diverges with k for every fixed nonzero beta. Since |beta|<=1, padding offers no small multiplicative rescue. Arbitrary potential terms constant along the q flow cancel in (10), even if they depend on all original roots and complementary coordinates.

### beta = 0: the unread-padding port is impossible for a nonaffine source

Now q=QP consists entirely of unread padding. Since N_h is independent of q, integration forces

    Psi(W,V,q+t e1,u)-Psi(W,V,q,u)=t N_h,1(W,V)               (12)

for all real t. If grad Psi has a globally bounded first L, the V1 gradient of the left side is 2L-Lipschitz for every t. Hence (12) would require D_V1V1 N_h,1=0 wherever that derivative exists.

But at the displayed root and V=e1, the terminal bump center has

    g_k,1''(1+h)=zeta(4/h)^2 chi''(0)=-1/(15h),
    g_k,1''(1-h)=0,
    D_V1V1 N_h,1=(kappa0 h/2)[g_k,1''(1+h)-g_k,1''(1-h)]
                      =-kappa0/30 !=0.                       (13)

Letting |t| grow in (12) contradicts every finite L. Thus this frame fails even for one fixed smooth member of the family. This argument concerns globally bounded full first, as required by the source interface; it is not a claim about a separately weighted local or Lp-only first contract.

Equations (8), (11) and (13) prove the stated fixed-known-readout theorem.

## 5. Consequences, costs and escape clauses

The midpoint change can improve other joint-current bookkeeping, but it does not supply a universal fixed-readout small-full-gradient realization of its exact opposite source on the original padded Gaussian tape. An arbitrary finite number of unread Gaussian blocks has already been included. Finite known orthogonal changes are included by the coisometry calculation. Heat/amplitude rescaling cannot remove the unbounded k in (11) when the source contract depends only on the available primitive first bounds.

The proof makes no demand to execute a Hessian derivative. The local Hessians and third jets only establish that a proposed exact VALUE-source type cannot have the claimed uniform bound. First/adjoint HVP access does not change that contradiction. Each numerical approximation or weaker law substitute would require its own absolute errors, true caller/readout, heat and source count, and exact retained covariance-field contract; none is obtained by discarding the equality premise above.

The theorem does not show that the entire K3 endpoint lacks a fixed-rank, zero-extra-copy repair. It isolates one natural purported entry into the existing full-gradient mixed kernel and rules it out uniformly. In particular these routes remain logically open:

- a source-specific/adaptive frame with its actual construction and derivative costs;
- a finite decomposition into genuinely admitted offspring with different fixed readouts, or a decomposition of the current rather than this exact source;
- an admitted source whose physical readout need not equal N_h pointwise on the full padded tape, coupled by a proved root-retained comparison;
- a new joint constrained-law generator or a complete endpoint covariance-field identity using other offspring;
- absorption of a retained baseline body through a marked parent after genuinely satisfying that parent's input-source conditions.

The mere fact that the affine ambient lift is small does not establish its nonlinear graph restriction, and this no-go does not prohibit using that lift in a new valid constrained-law theorem.

## 6. Diagnostics

check_midpoint_fixed_frame_obstruction.py builds the exact linear symmetry conditions for 600 independently varied Hessian tuples in each of d=2,3,4. In each case the kernel among the original W/V blocks is numerically one-dimensional, exactly the scalar V block; the remaining padding columns are unread and unconstrained. It also records a generic six-distinct-query same-potential base record and checks the beta-dependent first lower bound and the exact nonzero coefficient in (13). Results are in midpoint_fixed_frame_obstruction_checks.json.

The realizability of the local Hessian variations is proved in section 2, rather than inferred from these numerical matrix tests. The symmetry classification and no-go are analytic.
