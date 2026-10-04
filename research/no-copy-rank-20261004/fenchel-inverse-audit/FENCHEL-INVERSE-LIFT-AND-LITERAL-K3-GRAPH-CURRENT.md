# Fenchel inverse-gradient lift: a valid new block and its remaining K3 graph current

2026-10-04. Restricted source class: anchored original gradients g0=grad V0 with .4 I<=Dg0<=.6 I. The original producer oracle remains gradient VALUES; HVPs occur only in declared first/adjoint sweeps. This note constructs an exact gradient REFERENCE, gives its finite original-VALUE implementation, and tests its actual literal K3 graph. It does not prove a new complete covariance/mean producer or an all-rank recurrence.

## 1. Existing pins and what changes

The prior ambient lift is `p-native-twin/JOINT-GRADIENT-LIFTS-AND-CONSTRAINT-DESCENDANTS.md`, SHA-256 `64ed43694ae1497d50c3db93a9732beb00a7a3836bbf652fde6b5976981ecd56`. It has a six-VALUE genuine ambient gradient on eight independent d-blocks, but its K3 readout lives on a nonlinear graph. Its graph p-companion sum is a M*E+r M*d, and a natural full pullback has HVP-valued companions and unbounded original third-jet firsts.

The literal K3 source is `p-native-twin/THIRD-ITERATE-TWO-NODE-WEAK-TWIN-PORT.md`, SHA-256 `c8b20bf4eff6986f5e1f2d749e740beb935491a1043153ec62b39ef05bf4bca3`. We retain its same records, original-query aliases, and finite K=3 version:

    x0=cS+sU, c^2+s^2=1,
    Delta=gt(S)-gt(S+epsilon Z), d=a Delta,
    x1=x0+M*d, p_i=g0(x_i), t_i=S+aM p_i,
    E=r[gt(t0)-gt(t1)].                                 (1)

The useful change is exact: adding the signed dual potentials V0*(p0)-V0*(p1) removes the unmarked r M*d term from the COMMON p-companion sum. The opposite companion and nonlinear graph law remain.

## 2. Inverse gradient reference and finite original-VALUE implementation

Write m=.4, L=.6, q=(L-m)/(L+m)=1/5. Strong monotonicity gives a globally defined C1 inverse

    h(p)=g0^(-1)(p)=grad V0*(p),
    Dh(p)=[Dg0(h(p))]^(-1),  (5/3)I<=Dh<= (5/2)I.       (2)

These derivative matrices are analytical identities, not executed inverse-Hessian queries. Define a finite VALUE circuit, with a declared fixed integer N,

    u_0=0,
    u_(j+1)=u_j+2[p-g0(u_j)],  0<=j<N,
    h_N(p)=u_N.                                         (3)

Each step has derivative I-2Dg0 with norm <=q. Therefore

    |h_N(p)-h(p)| <= (5/2) q^N |p|,
    ||D_p h_N|| <= 2 sum_(j=0)^(N-1)q^j <=5/2.           (4)

For a varying captured caller theta, with direct primitive caller bound L_g and a caller-dependent p of bound L_p,

    ||D_theta h_N|| <= (5/2)(L_p+L_g),                   (5)

apart from any separately declared initializer/anchor path. For nonanchored primitives, replace |p| in (4) by |p-g0(0)| and include the actual stored g0(0) computation and caller derivative.

The actual first recurrence is

    J_(j+1)=(I-2H_j)J_j+2I, H_j=Dg0(u_j).

It uses one original HVP per stored iteration for a requested direction; reverse/adjoint replay uses the corresponding transposed recurrence. No derivative of H_j is needed. A discarded primal state costs its paid VALUE replay. First claims concern the finite real-arithmetic circuit, or coherent numerical source versions with their own first certificate; a bare rounding-error bound does not imply differentiability of a rounded implementation.

Crucially, h_N is generally NOT itself a gradient in dimension >=2. Starting at zero and using g0(0)=0,

    J_3-J_3* = 8[H_2,H_1].                              (6)

The independent audit gives a smooth strongly monotone two-dimensional example where this is nonzero. The admitted genuine-gradient object is h, with coherent finite VALUE restoration to h_N. We infer no derivative convergence from (4), and no uniform derivative convergence rate follows from only the Hessian sandwich.

If each g0 call has Euclidean absolute VALUE error at most delta_g, the inverse error due to these calls is at most (5/2)delta_g. Known finite arithmetic errors enter the same contraction recurrence. For an input profile ||p||_Lp<=P_p and a requested absolute inverse-output tolerance delta_h>0, it is enough to choose

    N >= max(0,ceil(log(5P_p/delta_h)/log 5)),
    delta_g <= delta_h/5,

so the truncation and provider errors each cost at most delta_h/2, with the actual later readout multiplier inserted BEFORE choosing delta_h. Known-arithmetic errors get a separate allocated part of this budget (or the two displayed allowances are reduced accordingly). Different source occurrences use their actual declared profiles/widths. This is logarithmic fixed-conditioning VALUE work, not an inverse-heat replication bank. An unknown or zero E mark does not select an absolute precision floor.

## 3. The exact two-block Fenchel reference

Set

    Q(x,p)=V0(x)+V0*(p)-x dot p.

Its gradient and Hessian are

    grad Q=(g0(x)-p, h(p)-x),
    D^2 Q=[[Dg0(x),-I],[-I,Dh(p)]].                      (7)

This is a globally bounded-first gradient: the block row bound gives ||D^2 Q||<=7/2. It is nonnegative and vanishes exactly on p=g0(x). On that graph its gradient is zero and its Hessian is positive semidefinite, with tangent kernel (v,Dg0(x)v).

It need not be jointly convex OFF the graph. In one dimension, if g0'(x)=.4 and g0'(h(p))=.6, the Hessian determinant is .4/.6-1=-1/3. Thus a penalty based on Q is not automatically a uniformly log-concave joint posterior in (x,p).

The finite executed pair

    (g0(x)-p, h_N(p)-x)

uses N+1 original gradient VALUES before exact anchors/caches. Its actual first is bounded by the same numerical block estimate, but it is generally only an approximate gradient reference. Its VALUE discrepancy from (7) is precisely (0,h_N(p)-h(p)), with (4)'s absolute floor.

## 4. Fenchel-corrected ambient K3 reference and finite program

Use independent coordinates Xi=(S,Z,x0,x1,p0,p1,d,lambda). Let 0<omega<=1 be a fixed numerical Fenchel weight; omega=1 recovers the direct modification of the previous ambient lift. Define the analytical potential

    L_F(Xi)/r = Vt(S+aM p0)-Vt(S+aM p1)
              + omega[Q(x0,p0)-Q(x1,p1)]
              + a[Vt(S)-Vt(S+epsilon Z)]-S dot d
              + lambda dot (x1-x0-M*d).                (8)

This is a SIGNED SADDLE reference, not a jointly convex sampling potential: it contains Q0-Q1 and the unconfined bilinear constraint multipliers. In particular the S,d Hessian block has a nonzero off-diagonal and zero d,d block. We do not sample exp(-L_F), assert its normalizability, or treat this reference as a positive Gaussian packet.

Its genuine full gradient reference has blocks

    G_S/r = gt(t0)-gt(t1)+a Delta-d,
    G_Z/r = -a epsilon gt(S+epsilon Z),
    G_x0/r = omega[g0(x0)-p0]-lambda,
    G_x1/r = -omega[g0(x1)-p1]+lambda,
    G_p0/r = a M*gt(t0)+omega[h(p0)-x0],
    G_p1/r = -a M*gt(t1)-omega[h(p1)-x1],
    G_d/r = -S-M lambda,
    G_lambda/r = x1-x0-M*d.                            (9)

No potential VALUE is queried. Execute G_N by replacing the two h values in (9) with the actual finite circuits h_N. Before identical-query caching this costs six original VALUES plus 2N inverse-iteration VALUES. The h_N trajectories are complete new query chains. They are not cached across independent random packets unless exact captured-caller dependency permits it.

The previous six-VALUE ambient source plus the two inverse blocks gives the conservative common bound

    Lip G <=20r,  Lip G_N<=20r.                         (10)

All actual callers are the literal original-gradient caller paths, including each inverse path (5), not derivatives of (4). At a standard independent 8d Gaussian input, the anchored ordinary profile is O_p(r sqrt(d)). The complete source VALUE restoration is

    ||G_N-G||_Lp <= (5/2)r omega q^N
                  ||(|p0|^2+|p1|^2)^(1/2)||_Lp.        (11)

Each first/adjoint sweep uses only original HVPs at the 6+2N recorded query sites. Finite approximation need not be curl-small because (11) is a VALUE restoration, not a derivative comparison.

One can omit the final feedback-constraint/Delta terms and use only the five blocks (S,x0,x1,p0,p1) for a smaller terminal/ancestor component. It then has four direct VALUES plus 2N inverse VALUES. That smaller component does not enforce x1=x0+aM*Delta; its precise remaining graph is still (1).

## 5. Exact restriction to the literal native graph

At the actual graph Gamma specified in (1), with lambda=0, the exact reference satisfies

    G_S(Gamma)=E,
    G_x0(Gamma)=G_x1(Gamma)=G_lambda(Gamma)=0,
    G_Z(Gamma)=-ra epsilon gt(S+epsilon Z),
    G_p0(Gamma)=ra M*gt(t0),
    G_p1(Gamma)=-ra M*gt(t1),
    G_d(Gamma)=-rS.                                    (12)

Thus

    G_p0+G_p1 = a M*E,                                 (13)
    G_p0-G_p1 = ra M*[gt(t0)+gt(t1)].                   (14)

Equation (13) is a genuine improvement over the old aM*E+rM*d sum. It holds for every strongly monotone g0 under the stated contract, not merely a quadratic fixture.

For the finite inverse circuit, let e_N(p)=h_N(p)-h(p). The common-mode identity becomes EXACTLY

    (G_N)_p0+(G_N)_p1 = aM*E
                         +r omega[e_N(p0)-e_N(p1)].    (15)

That displayed residual is an executable finite-program discrepancy. Its Lp bound is at most sqrt(2) times the right side of (11), due to the unnormalized common-sum readout; equivalently use (5/2)r omega q^N || |p0|+|p1| ||_Lp. It vanishes literally when p0=p1 if both inverse paths use the same finite version and identical-query arithmetic. Under merely bounded Hessians, no q^N-small derivative or extra relative E-grade is inferred from its small VALUE norm. Its absolute budget is selected before execution.

### One-energy test

Use the same scalar quadratic g0=gt=m id, M=I, with .4<=m<=.6, c,s>0. Then

    E=b Z, b=ra^2 m^3 epsilon,
    (G_p0-G_p1)/sqrt(2)
      =sqrt(2)ra m[S+a m(cS+sU)-(a^2 m^2 epsilon/2)Z].  (16)

Consequently its L2 norm is at least sqrt(2)ra m sqrt(d), whereas e_E=ra^2 m^3 epsilon sqrt(d). The ratio is at least sqrt(2)/(a m^2 epsilon). At epsilon=0 the source E vanishes exactly while the opposite companion remains nonzero. The full eight-block source also contains -rS, so its ordinary energy is larger still.

The new full ambient lift therefore does NOT have one-E energy. The common companion does; the opposite baseline and other channels must be retained as actual unmarked fields or treated by a new complete positive packet. Averaging their mean to zero does not delete their field energy or their correlations with the physical branch. A finite inverse does not change this conclusion even if its absolute error is arbitrarily small.

## 6. The exact tangent/current target remains

For a graph pair (x,p=g0(x)), the tangent derivatives of the two Fenchel constraints are

    D[g0(x)-p]=Dg0(x) Dx-Dp=0,
    D[h(p)-x]=Dh(p) Dp-Dx=0.                           (17)

These are true identities on the SAME query graph. They do not make p an independent standard Gaussian. The pushforward of a Gaussian x by g0 has density

    density_p(p)=phi(h(p)) det Dh(p),                   (18)

and jointly x=h(p) holds on a singular nonlinear graph. An inverse reference transports the input law; it does not restore the original Gaussian law without the Jacobian and all joint labels. The needed density derivatives can involve derivatives of Dg0 and are not admitted original C2 VALUE operations.

At a literal whole-input Gaussian clock the actual physical ancestor word is still

    D_U E = ra s[K0-K1],
    K_i=Ht(t_i) M H0(x_i),                             (19)

with all factors evaluated at their original same-record queries. The leading skew target is therefore unchanged:

    C_lead=P_S* B C_anc-C_anc* B*P_S,
    B=ra E_fine[K0-K1], C_anc=[cI,sI,0].                (20)

The inverse identity H0(x)^(-1)H0(x)=I in (17) is an exact tangent identity; using it to replace (19) by independently averaged inverse/terminal factors changes the joint query law.

Taking the full gradient pullback of (8) does not bypass this. The Fenchel gaps vanish identically along the graph, so their entire pullback contributes zero. The remaining terminal potential pulls back to r[Vt(t0)-Vt(t1)] plus the previously explicit feedback terms. In particular its ancestor companion contains

    ra[H0(x0)M*gt(t0)-H0(x1)M*gt(t1)].                 (21)

Executing (21) as a VALUE would require HVP outputs, and differentiating it would expose uncontrolled third jets. This is exactly the old natural-pullback problem, despite the bounded-first ambient reference (9). The correct finite alternative is to keep (9) ambient and solve its constrained-input current; that current has not been supplied here.

## 7. Positive penalties and approximate graphs: exact remaining terms

For a fixed x the Fenchel gap is strongly convex in p, so the conditional probability

    K_tau(dp|x)=Z_tau(x)^(-1) exp[-Q(x,p)/tau] dp        (22)

exists for every tau>0. But globally normalizing phi(x) exp[-Q/tau] does not preserve x: its marginal is proportional to phi(x)Z_tau(x). As tau tends to zero,

    Z_tau(x)/(2pi tau)^(d/2) -> sqrt(det Dg0(x)).        (23)

The independent audit gives an elementary affine example showing that an independent ambient Gaussian p prior changes even the limiting x covariance. Strong convexity of a fixed-p or fixed-x fiber is not joint convexity of Q.

If one instead normalizes (22) at EACH x, the exact joint x-score is

    -x + [p-E_(K_tau(.|x))p]/tau,                     (24)
    grad_x log Z_tau(x)=[E_(K_tau(.|x))p-g0(x)]/tau.

For the literal graph (1), preserving the old Gaussian W=(S,U,Z) and conditionally drawing the two fibers gives the exact density

    phi(W) product_(i=0,1) K_tau(dp_i|x_i(W)).

Its old-root score is

    -W + sum_(i=0,1) (D_W x_i)*[p_i-m_tau(x_i)]/tau,
    m_tau(x)=E_(K_tau(.|x))p,

and its p_i score is -(h(p_i)-x_i)/tau. These are the actual extra constraint-current fields, with the SAME nonlinear feedback derivative in D_W x1. Without individual conditional normalization, the old-root score instead contains p_i-g0(x_i), and the W marginal is reweighted by Z_tau(x0(W))Z_tau(x1(W)).

The conditional means and normalizers in these formulas are new analytical targets, not available oracles. Independently sampling the two p fibers destroys the exact paired mark: for quadratic g=m id, epsilon=0, native E=0, but p0,p1 sampled independently have variance m tau in each fiber and

    Var(rm a(p0-p1))=2r^2 a^2 m^3 tau I !=0.           (25)

For m=1/2 this is r^2 a^2 tau I/4. Thus independent positive penalty banks do not preserve the original one-energy graph.

A simpler EXECUTABLE approximate graph uses one SHARED Gaussian V:

    p_i^sigma=g0(x_i)+sigma V, i=0,1.                 (26)

It keeps all old aliases and adds two genuine defects

    g0(x_i)-p_i^sigma=-sigma V,
    h(p_i^sigma)-x_i,

the latter bounded by (5/2)sigma|V|. It induces exactly the already considered common terminal displacement a sigma M V:

    E_sigma=r[gt(t0+a sigma MV)-gt(t1+a sigma MV)].      (27)

It vanishes when the paired terminal queries coincide, and obeys

    |E_sigma-E|<=2ra sigma ||M|| |V|.                  (28)

The common p-companion now also contains the actual inverse-constraint difference

    r omega{h(g0(x0)+sigma V)-x0
                    -h(g0(x1)+sigma V)+x1}.           (29)

It is bounded by r omega min{5 sigma|V|, .5|x0-x1|}. Indeed the derivative of x -> h(g0(x)+delta)-x is B^(-1)(A-B), where A,B both lie between .4 I and .6 I, so its norm is at most .5. No product of the two small factors follows: the independent audit gives a same-potential scalar high-frequency fixture with this difference of order 1/k while |delta||x0-x1| is of order 1/k^2. Finite inverse errors add the actual common-sum version of (11) at the changed inputs.

On the nondegenerate M=I subclass with BOTH g0 and gt having Hessians in [.4 I,.6 I], the uncentered norm obeys ||E||_2 >= .4^3 ra^2 epsilon sqrt(d), by three applications of strong monotonicity. A centered bound also holds here: D_Z E=ra^2 epsilon Ht(t1)H0(x1)Ht(S+epsilon Z). Writing each Hessian as .5I+K with ||K||<=.1 gives this product as .125I+R with ||R||<=.6^3-.5^3=.091. Its symmetric part is therefore at least .034I. Conditional Gaussian first-chaos Bessel in Z, followed by averaging S,U, gives e_E=||E-E E||_2 >= .034 ra^2 epsilon sqrt(d). This prices (28) by C sigma/(a epsilon) e_E. It also prices the common-mode inverse floor by C omega q^N/(a^2 epsilon) e_E when the actual p input profile is O(sqrt(d)): choosing q^N below the desired relative grade times a^2 epsilon/omega costs only a logarithmic N. Neither estimate bounds the derivative of the small residual. Thus polynomially narrow sigma gives a legitimate relative restoration at fixed target rank and does not itself require inverse-sigma replicas. This is a possible approximation component, not closure: (29), the opposite baseline, and the old-root joint/current law remain actual fields. On degenerate/general M, a nominal sqrt(d) bound cannot be divided by an unknown actual e_E without a separately proved lower bound.

## 8. Cost, retained return, and precise outcome

- The raw K3 graph still uses six original VALUES. Its inverse-augmented ambient component uses at most 6+2N VALUES per complete occurrence. On the constrained graph original sites can be shared within that occurrence, but each independent new bank rebuilds its full private ancestors and inverse chains.
- Exact inverse reference firsts are bounded by 5/2; actual finite inverse firsts and callers have (4)–(5). Finite implementation uses N declared first/adjoint actions per inverse call when requested. No HVP is emitted as a new VALUE source, and none is differentiated.
- If a declared output/floor target is A^b times a known profile at fixed b, N=O_b(log(1/A)+log(profile/tolerance)). Therefore the inverse primitive introduces b_extra=0 in the inverse-heat copy exponent under fixed conditioning. All finite query tolerances, Gaussian dimensions, known matrix work, storage, finite clocks, shares and later readout norms remain priced.
- The reference (9) is genuine gradient; the executed G_N is a coherent finite VALUE restoration, with absolute error (11). A theorem requiring an ACTUAL gradient cannot silently receive G_N without that restoration step.
- On standard ambient inputs the source has ordinary O(r sqrt(d)) energy and first O(r); on the actual graph the physical block is E, the common p-block is marked, and the opposite/other blocks remain unmarked. Retaining all of them is allowed only with their actual profiles and full same-root first/caller paths.
- A global Fenchel penalty is positive as a weight but is not generally jointly log-concave and changes the root law. A conditionally normalized penalty has a new normalization/current requirement; independently sampled fibers destroy the paired E mark.
- The actual constrained K3 target is still (19)–(20). No stationary joint producer for it, no full one-E return, and no general c_P improvement follows from the inverse-gradient lift alone.

The finite positive result is the dual-completed ambient reference and its exact marked common-mode improvement (13), with a logarithmic-cost original-VALUE inverse implementation. The precise remaining requirements are the unmarked opposite mode and the nonlinear constrained-input/current law, not an unbounded first of the ambient inverse reference.

## 9. Verification

`check_fenchel_k3_lift.py` has 1,880 passing diagnostic checks covering finite inverse VALUE/first bounds, the explicit nonzero finite-inverse skew, ambient gradient blocks, the exact native graph and common-mode identity, shared terminal heat, opposite baseline energy, and independent-fiber false energy. Its reference inverse is used only inside the checker to inspect mathematical identities. The separate independent audit and its 1,849 checks cover the normalizer/current and off-graph convexity gates. Neither numerical check claims a general constrained-law producer.
