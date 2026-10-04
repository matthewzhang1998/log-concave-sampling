# One-pair known-matrix commutator adapter

Publication copy: nonmathematical context and/or local paths were sanitized. Original and public SHA-256 values are recorded in `INVENTORY.json`; source/audit pins below identify their historical versions, not these edited bytes.

New bounded source construction, 2026-10-04. This is a refinement of the
historical stationary skew adapter at a special coefficient type: one unknown
symmetric heat matrix and one known symmetric matrix. It is not a source
extraction theorem for a general nonlinear ancestor.

## 1. A single block pair and two known contractions

Let H=H* be the selected mean Jacobian of an admitted genuine-gradient source
h, normalized so that the complete source radius is at most rho. Assume the
selected matrix of its pair is s0 rho H. Let D=D* be known, ||D||<=1.
The genuine block-gradient source (h(x1),h(x2)) has selected matrix

    A=s0 rho diag(H,H).

Both blocks use the same captured coefficient roots and finite source version;
their input coordinates and all pair innovations have the ordinary declared
joint law. Let

    V=[D;I]/sqrt(2),   U=[I,-D]/sqrt(2).

These are contractions. Given physical standard passive p, execute the full
block pair at Vp+J Z0, J J*=I-VV*, and return U times that output plus F Z1,
F F*=I-UU*. These roots depend only on known geometry. The actual source-zero
carrier and every fill are retained. The pair law, at fixed captured roots,
has conditional reference

    Q(p) ~ N(Kp, I-KK*),
    K=U A V=(s0 rho/2)(H D-D H).

Indeed its conditional covariance is
U[I-A V V* A*]U*+I-UU*=I-KK*. Thus the whole output has standard marginal
after standard p is integrated. K is skew. Only ONE complete block pair is
executed; the formula does not compose two selected-matrix pair outputs.

## 2. Known affine state pullback

For a known symmetric contraction D and B=I+aD, 0<a<1/2, consider
f(x)=r[g(Bx)-g(x)], where g is an original gradient with first at most one.
At original common-input heat cX+vZ put Hbar=E Dg(B(cX+vZ)). The genuine
gradient pullback h_B(x)=B g(Bx) supplies Hprime=B Hbar B at this SAME heat.
Use two copies of the anchored source

    (rho/v) B[g(B(cX+v x))-g(cB X)].

Its actual declared radius is L0 rho for a fixed L0>=||B||^2. Apply the
compiler with that declared radius and its corresponding body guard. Its
selected matrix is exactly s0 rho Hprime; L0 is a numerical norm constant,
not an omitted amplitude division.

For clarity write the selected block as s0 rho diag(Hprime,Hprime). Set

    U=B^(-1)[I,-D]/C,  V=[D;I]B^(-1)/C,
    C>=sqrt(2)||B^(-1)||.

Then U,V are known contractions and, because B commutes with D,

    U diag(Hprime,Hprime) V
      =B^(-1)(Hprime D-D Hprime)B^(-1)/C^2
      =(Hbar D-D Hbar)/C^2.

All inverse norms are numerical. The source query remains B(cX+vZ), so its
Gaussian covariance has not been changed to that of BY+vZ. The complete X
first of each anchored source is O(rho/v). The adapter preserves this row
up to its known numerical factors; it does not query a derivative of Hbar.

## 3. Positive cross-fork normalization and paths

Let a separate marked pair, sharing the SAME captured roots and passive p,
have selected matrix M=s0 J. Suppose the true orientation field is

    O=-tau Sym(J W),   W=Hbar D-D Hbar,

where tau is a known scalar (for the affine f above, tau=ra). The one-pair
channel has exactly K=(s0 rho/C^2)W with the preceding source normalization.
Since K*=-K,

    M K*+K M* = 2 s0^2 rho O/(tau C^2).

With independent complete private banks conditional p and roots, emit
b Q_E+c Q+sqrt(v0-b^2-c^2)Z, where

    bc=-q tau C^2/(2 s0^2 rho),
    b^2+c^2<=v0-v_keep,   v_keep>0.

Its conditional-root marginal covariance is exactly v0 I-q O. Its finite
comparison pays |b| delta_E+|c| delta_Q. As in the historical positive fork,
use a fresh passive independent of the old consumer conditional roots, or
retain the full conditional-passive kernel instead of invoking its marginal.
The covariance field is cancelled at the same roots before they are integrated.

Keep the unmarked stationary Q body and its known carriers. Delete only the
marked b(Q_E-Z_E), with actual energy O(|b|e). The kept unmarked residual first
is O(|c|rho)=O(|tau|/|b|), and its old-root/caller row is
O(|tau|/(|b|v)) times the actual physical center injection. Its ordinary energy
can be sqrt(d)|tau|/|b| and is not assigned the marked e profile.

Compared with a two-pair product channel, no inverse-rho factor remains in
this retained residual/caller row. The numerical pair floor still pays
|c| Lambda_B sqrt(n)rho^B, and positivity imposes a real lower bound on rho
relative to tau/b. With tau=A^(1+g), b=A^zeta, rho=A^gamma, fixed geometry
admits zeta<g and gamma<1+g-zeta, subject to the actual body guard.

For the positive whole-source orientation clock with weights w_j and widths
v_j, multiply the COMPLETE per-clock fork by sqrt(w_j) and scale its covariance
share and target together. The established sums

    sum w_j/v_j^2=O(log),  sum sqrt(w_j)/v_j=O(polylog)

price independent root firsts and shared callers. Thus the weighted row guard
uses sqrt(w_j)|tau|/(|b|v_j), not an unweighted v_min penalty. The full known
variance fill and the same-root comparison are still necessary.

## 4. Marked difference variant and limits

Using the same complete block-pair records and fills at a marked Y and its
known source-zero carrier Z_E gives the literal difference Q(Y)-Q(Z_E).
Its energy is O(rho ||Y-Z_E||), because the pair carrier is passive-independent.
Readout tau/rho therefore gives actual energy O(tau e), private first O(tau)
and public O(tau A); root first is O(tau/v). The finite conditional-mean proof
uses separate standard/near-standard marginal bounds and the pair's mean
defect d, whose Lipschitz bound is O(rho). Its restored error is
O(tau[rho^(-1) delta_Q+delta_E+eta_E]), including the separately supplied
marked conditional-mean error eta_E. No joint Gaussian coupling of (Y,Z_E)
is invoked. This remains a forward
marked coefficient. The positive cross fork in Section3 handles the required
adjoint ordering without dropping a two-curl term.

Every VALUE call of the one full block pair, both duplicated gradient blocks,
all input/output fills, their precision, and discarded-tape replays must be
counted. Dense known maps are not free; scalar-block/native geometry can give
a cheap explicit implementation. The one-pair reduction supplies neither a
new general source-radius theorem nor an all-order complexity conclusion.
