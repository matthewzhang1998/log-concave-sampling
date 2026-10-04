# Exact third-iterate two-node weak-twin port

New source-qualified boundary note, 2026-10-04. It identifies the smallest
literal signed-twin callback not covered by the proved scalar-shear repair.
It is an exact finite VALUE construction, not a counterexample to a possible
future repair. No new joint-heat theorem is asserted.

Source: surviving LOW30 SHA-256
`7767ac389d6b677eb9517b2585f7a24cb9f85ede47d099a653c635b7be750fb8`,
`t30:eq:twin`, `t30:eq:proxy-readout`, `t30:lem:finite-columns`.
The completed weak-family remainder is F_actual-H_fin on the same record;
this note takes the retained input F itself to be the explicit two-node graph,
so there is no separate old actual-to-retained tail in this example.

## 1. Dimensions, original sources, and protected Gaussian rows

Let X,Y,Z be independent standard Gaussians in R^d and W=(X,Y). Choose fixed
c0,s0>0 with c0^2+s0^2=1; for example c0=s0=1/sqrt(2). Let

    C=[I_d,0],  P=[c0 I_d,s0 I_d],
    C C*=P P*=I_d,  P C*=c0 I_d.

Let Pi select the Y coordinates. Then C Pi=0 and P Pi P*=s0^2 I_d, an exact
protected terminal row. All these maps are fixed under Gaussian/caller
 differentiation. Let M be a fixed known d-by-d contraction. Let 0<a<=1/16,
0<epsilon<=1, r>0, and set A_t=a M.

The original maps g0,g_t:R^d->R^d are gradients with Jacobian operator norm
at most one. They may be two occurrences of the SAME original potential at
different fixed/moving anchors. Their displayed anchored versions satisfy
g_i(0)=0. Every anchor and its actual caller path is kept. The only primitives
are these original gradient VALUES and known affine maps. FIRST/adjoint
sweeps use their actual original HVPs at the recorded queries.

The two-node directed force is

    p0(W)=g0(CW),
    p_t(W)=g_t(PW+a M p0(W)),
    F(W)=r p_t(W).                                      (1)

Its ancestor does not read the protected Y increment. Its original terminal
predecessor has physical coefficient a M. No rank-one or orthogonality
condition on M is imposed.

## 2. The exact finite simultaneous twin

Use the LOW30 signed twin with auxiliary perturbation epsilon Z only in the
MINUS terminal. Starting with all four p^[0]=0, the literal iteration is

    p0,+^[k+1] = g0(CW+a M* [p_t,+^[k]+p_t,-^[k]]),
    p_t,+^[k+1] = g_t(PW+a M p0,+^[k]),
    p0,-^[k+1] = -g0(CW),
    p_t,-^[k+1] = -g_t(PW+a M p0,+^[k]+epsilon Z).         (2)

Every occurrence uses the SAME W,Z and same finite versions. There is no
refresh between the plus/minus terminal difference or its later ancestor.
Define H_[K]=r p_t,+^[K] and E_[K]=F-H_[K].

At K=1, H_[1]=r g_t(PW). At K=2, p0,+^[1]=g0(CW), hence H_[2]=F exactly.
At K=3 the first nontrivial feedback enters:

    Delta(S,Z)=g_t(S)-g_t(S+epsilon Z),  S=PW,
    q0,fb=CW+a M* Delta(S,Z),
    H_[3]=r g_t(S+a M g0(q0,fb)),
    E_[3]=r{g_t(S+a M g0(CW))
                 -g_t(S+a M g0(CW+a M* Delta(S,Z)))}.    (3)

This identity fixes both signs and all aliases. It agrees with the affine
regression: if g0(x)=H0 x and g_t(x)=Ht x with symmetric fixed matrices, then

    E_[3]=r epsilon Ht (aM) H0 (aM)* Ht Z.

In the scalar equal-H case this is r a^2 epsilon h^3 Z. The physical forward
square can therefore miss a real auxiliary covariance even in this native
finite twin; the complete source record includes Z.

## 3. Genuine-gradient reference and the distinct coisometries

Stack the two directed rows as C_dir=[C;P], and write

    A_dir=[0,0;aM,0],
    C_tw=[C_dir,0; C_dir,epsilon e_t],
    S_tw=[A_dir+A_dir*,A_dir*; A_dir,0].

The fixed-point twin is the symmetric implicit system with original blocks
(g0,g_t,-g0,-g_t). Its exact reference
G_*=(r/beta) C_tw* p_* is a genuine full gradient. A finite
G_K=(r/beta) C_tw* p^[K] is an actual VALUE implementation and need not itself
be a gradient. Here

    beta=(s0^(-2)+epsilon^(-2))^(-1/2),
    B=beta[0,s0^(-1)I_d,-epsilon^(-1)I_d]

acts on the COMPLETE master record (X,Y,Z). Direct multiplication gives

    BB*=I_d,  B C_tw*=beta e_(t,+)*,
    B G_K=H_[K].

The recorded physical row is P_rec=[c0 I_d,s0 I_d,0], so

    P_rec P_rec*=I_d,   P_rec B*=beta I_d.

They are not the same row. All curl claims below use P_rec* E_[K], not B*E_[K].
The finite-column proof in LOW30 applies to every K, independently of
convergence of derivatives. Exact-gradient law comparisons use G_* and
VALUE-only restoration back to the SAME chosen G_K.

The norm of S_tw is bounded by3a<3/8, so the finite iterations and all their
source/zero versions can be chosen from a geometric VALUE-error budget.
For a requested gradient-mean numerical tolerance, K is the assigned finite
integer for the WHOLE batch. The named K=3 callback isolates the first new
shape; a final high-accuracy implementation usually has larger K.

## 4. Actual energy and complete first

Pointwise, |Delta(S,Z)|<=epsilon|Z|, so

    |E_[3](W,Z)|<=r a^2 epsilon |Z|.                    (4)

Thus ||E_[3]||_Lp<=C_p sqrt(d) r a^2 epsilon. It is an uncentered bound at the
actual full record. In fact for every K>=3, (2) gives
|p0,+^[K-1]-g0(CW)|<=a epsilon|Z| and the SAME bound(4) for E_[K]. This does
not require independent terminal copies or a differentiated small error.

For epsilon<=1, the feedback query has complete first

    ||D_(W,Z) q0,fb|| <=1+a(1+sqrt(1+epsilon^2))<=1+3a.

Therefore

    Lip F<=r(1+a),
    Lip H_[3]<=r(1+a+3a^2),
    Lip E_[3]<=r(2+2a+3a^2)<=3r.

The special auxiliary derivative is genuinely small:

    ||D_Z E_[3]|| <= r a^2 epsilon.

Its exact product, including sign, is

    D_Z E_[3]
      =r a^2 epsilon Ht(q_t,fb) M H0(q0,fb) M* Ht(S+epsilon Z),
    q_t,fb=S+a M g0(q0,fb).

This ordered product is not generally symmetric. It is not replaced by a
Hessian oracle or by independent copies of either terminal query.
For larger K the literal contraction recurrence gives the corresponding
uniform complete first and an auxiliary first bounded by a fixed multiple of
r a^2 epsilon. The actual full W first need not inherit the epsilon factor.

The direct P_rec term in each physical terminal Jacobian is symmetric after
P_rec* is applied. Bounding only the remaining ancestor terms gives

    ||Curl(P_rec* E_[3])||op
       <= 2ra+2ra(1+3a) <=5ra.                         (5)

The same finite terminal/resolvent proof supplies O(ra) at the later chosen K.
No VALUE convergence is used to infer (5).

At the weak canonical scales r=A, a=A, epsilon=A^.9, beta~A^.9:

    e_E,p<=C_p sqrt(d) A^3.9,
    full first E_K=O(A), full recorded curl=O(A^2),
    full G_K radius=O(A^.1), both physical DG_K sides=O(A).

Thus this is a literal weak-mark2.9 proxy-tail family, with actual finite K
and its separate numerical restoration budget. It is not merely an abstract
small-energy field fitting those inequalities.

## 5. Actual caller paths

Keep r,a,epsilon,M,P,C fixed under theta differentiation. Suppose the shifted
original-gradient occurrences have uniform direct caller bounds L0,Lt at
fixed displayed arguments. Then at K=3,

    ||D_theta F||<=r(Lt+a L0),
    ||D_theta H_[3]||<=r[Lt+a L0+2a^2 Lt],
    ||D_theta E_[3]||<=r[2Lt+2a L0+2a^2 Lt].             (6)

These bounds include BOTH terminal occurrences inside Delta and the shifted
ancestor occurrence. Any additional actual affine query/caller maps have their
known paths composed into(6). At canonical physical caller L_i=O(A^-1/2),
r=A, the return is O(sqrt(A)). Later iterates use the same bounded finite
resolvent. The finite proxy's full caller retains its beta division, whereas
its physical readout has the physical bound. Stored zero/anchor values remain
actual caller computations.

## 6. Smallest missing source interface

The proved scalar-shear repair uses one ancestor g0 at a linear Gaussian
coordinate and a terminal correction with a known rank-one orthogonal edge.
In(3), the ancestor query contains the nonlinear feedback

    a M*{g_t(S)-g_t(S+epsilon Z)}

on the SAME terminal/twin pair used by the remaining graph. Even one ancestor
and one terminal now produce mixed nonlinear ancestry. It is not legal to
freeze this query and replace its gradient by an independently heated
primitive while preserving the original coefficient. The issue persists for
d=1 through the shared terminal/twin records, and for d>1 also involves full
matrix products rather than one scalar derivative.

The original P,C rows already have a genuine protected gap. The finite graph
and its physical/column/caller bounds are explicitly admitted by LOW30. Thus
the next missing port is joint common-input heat/coefficient extraction for
the copied nonlinear feedback, with its intact E_K mark. It is not merely
absence of a protected terminal column or lack of a generic small-first
retained graph.

## 7. Exact conditional coordinates for the next joint-heat test

Define S=c0 X+s0 Y and U0=s0 X-c0 Y. Then S,U0,Z are independent standard and
CW=c0 S+s0 U0. Formula(3) becomes

 E_[3]=r{g_t(S+aM g0(c0S+s0U0))
      -g_t(S+aM g0(c0S+s0U0+aM*Delta(S,Z)))}.             (7)

Conditioning on S leaves a genuine independent Gaussian ancestor base U0 and
twin innovation Z, while the pair(S,S+epsilon Z) must remain coupled.
The feedback has a legal conditional gradient VALUE source

    k_S(z)=a[g_t(S)-g_t(S+epsilon z)],
    D_z k_S=-a epsilon Ht(S+epsilon z),
    Lip_z k_S<=a epsilon,  Lip_S k_S<=2a,
    a M*Delta(S,Z)=M* k_S(Z).

This is a gradient in z, with the exact anchored VALUES shown, and has
pointwise energy |k_S(Z)|<=a epsilon|Z|. It does NOT license treating
c0S+s0U0+M*k_S(Z) as Gaussian or replacing k_S by its conditional mean inside
g0. Its full root S remains in both terminal/twin calls and the outer terminal.

At a whole-input clock, apply its SAME affine Mehler substitution to the
COMPLETE (X,Y,Z) before using(7). Conditional coordinates acquire their actual
Gaussian means/widths; preserve those means and all zero-anchor offsets.
A future joint-heat construction must retain both the terminal/twin pair and
the ancestor base, or price a proved joint replacement. Separate marginal
stationarity and a small VALUE mark do not supply that replacement.

All primitive calls, repeated finite iterations, source-zero arrays and
FIRST/adjoint replays remain in the query count. A reference expectation,
derivative matrix or fixed-point limit is not executed.
