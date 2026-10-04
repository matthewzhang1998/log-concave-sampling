# Hidden gradient ancestor: variance-preserving weak heat comparison

New candidate source extension, 2026-10-04. This replaces the scalar additive
heat/BV step by a dimension-free Gaussian OU comparison when the ancestor
Gaussian is independent of a fully heated terminal Gaussian. The coefficient
comparison is proved below. Positive-fork, old-reference, complete cost and
retention interfaces must be bound to this exact source before a completed
module is asserted. General correlated or nonlinear ancestor graphs remain
outside this note.

## 1. Literal source and exact conditional coefficient

Let Y in R^n and Z in R^m be independent standard Gaussians. Let g=grad V on
R^n and g0=grad V0 on R^m be admitted original-gradient VALUE maps, each with
Jacobian operator norm at most one. A common original potential restriction
is permitted when it has these literal VALUE and first contracts. Let U be a
known contraction from R^m to R^n, let 0<a<1/2, and put

    q(Y,Z)=Y+a U g0(Z),
    E(Y,Z)=r[g(q(Y,Z))-g(Y)],
    f(Y,Z)=(E(Y,Z),0) in R^(n+m), P=[I_n,0].

All original gradient calls in one f occurrence share that same(Y,Z). The
complete first is at most C r; set A>=C r. Denote kappa0=ra and the actual
centered Gaussian energy of f by e. No separate terminal energy replaces e.

With H=Dg and H0=Dg0,

    Df = [ r(H(q)-H(Y)),  kappa0 H(q) U H0(Z); 0,0 ].

The first block is symmetric. Thus the only full square curl is

    Curl f = kappa0 [ 0, H(q) U H0; -H0 U*H(q),0 ].     (1)

At one whole-input OU clock c^2+s^2=1, capture X=(X_Y,X_Z). The marked
callback is the COMPLETE source [f(cX+s w)-f(cX)]/s, with its own entire w
record. Its selected coefficient is J_f(X)=E_w Df(cX+s w). It does not read
any subsequently exposed curl-side innovation.

In the independent curl factor expose z=cX_Z+s Z_0, then average only the
terminal Y innovation. Put

    Hbar_X(z)=E_T H(cX_Y+s T+a U g0(z)).

The conditional curl C_f(X)=E Curl f(cX+s w) is kappa0 times the skew block
built from E_Z0[Hbar_X(z) U H0(z)]. The product with J_f(X) is a product of
two conditional means at the SAME X. The marked w and curl-side Z_0 are
independent by this identity and by the callback read sets; the ancestor
inside the marked w is never refreshed separately.

## 2. A full Hilbert derivative cut with the ancestor amplitude

For fixed physical unit x define the m-vector

    F_x(z)=U* Hbar_X(z) x.

Then ||F_x||<=1 and, UNIFORMLY in z,X,

    ||D_z F_x||HS <= C a/s.                             (2)

To prove it, first let H_s(q)=E_T H(q+sT). By gradient symmetry and Gaussian
first-chaos Bessel, the derivative of q -> H_s(q)x has matrix-HS norm at most
C/s. Indeed it is s^(-1) E[(H(q+sT)x) tensor T], after exchanging the input
indices using symmetry of the third derivative in the distributional sense.
The vector H x has norm at most one, so Bessel gives C/s with no sqrt(n).
This identity is justified for C2 V by Gaussian convolution and a smooth
approximation, or directly by the weak derivative/IBP formula.

Now Hbar_X(z)=H_s(cX_Y+a U g0(z)). Chain rule composes that Hilbert matrix
with a U H0(z), whose operator norm is at most a, and with U* at the output.
The ideal property of the HS norm proves (2). Only Dg0 is used. There is no
g0'' or uncontrolled higher derivative of the terminal potential.

For the standardized ancestor xi=(z-mu)/s, mu=cX_Z, equation(2) reads

    ||D_xi F_x(mu+s xi)||HS <= C a.                    (3)

This is stronger than an operator-valued Lipschitz bound; its full Hilbert
norm is what removes any sqrt(m) cost.

## 3. Fixed relative OU smoothing of the ancestor Hessian

Fix a numerical d0 in(0,1), for example d0=1/2, and set

    eta=s sqrt(1-d0^2),  z_c=mu+d0(z-mu),
    H0,OU(z)=E_T0 H0(z_c+eta T0).                      (4)

This is the Gaussian OU operator in the conditional law N(mu,s^2 I), so that
z and z_c+eta T0 have exactly the SAME marginal. It is not additive variance
inflation at a frozen ancestor.

For physical unit x and ancestor unit y, OU self-adjointness gives exactly

 E_z x*Hbar_X(z)U[H0,OU(z)-H0(z)]y
   = E_z [(P_OU-I)F_x(z)] * H0(z)y.

The Gaussian Dirichlet estimate and (3) imply

    ||(P_OU-I)F_x||_2 <= C a,
    || E_z Hbar_X(z)U[H0,OU(z)-H0(z)] ||op <= C a.      (5)

The second factor H0(z)y has norm at most one. No dimension factor enters.
All expectations here are at the actual conditional ancestor law, retaining
its correlation with the terminal center a U g0(z).

Let C_f,OU(X) be the skew block from (1) with H0,OU in place of H0 but the
terminal center still using the ACTUAL g0(z). Then

    ||C_f,OU(X)-C_f(X)||op <= C kappa0 a.

The intact whole-source marked cut ||J_f||_(L2_X;HS)<=e/s gives

    ||Sym[J_f(C_f,OU-C_f)]||_(L2_X;HS)
       <= C kappa0 a e/s.                              (6)

Summing a positive orientation clock uses sum w_j/s_j=O(1), so the complete
constant-orientation calibration error is C kappa0 a e. In particular one
gets an extra ancestor amplitude while eta remains a fixed fraction of s.

## 4. Actual full-gradient sources and the selected skew channel

At captured(X,z), use the full (n+m)-square block source

 h(y,w)=rho( [g(cX_Y+aU g0(z)+s y)-g(cX_Y+aU g0(z))]/s,
             [g0(z_c+eta w)-g0(z_c)]/eta ).             (7)

Both blocks are gradients on their own active coordinates. The source's
complete active radius is at most rho and selected matrix is
D=diag(Hbar_X(z),H0,OU(z)). Its actual coefficient-root first is O(rho/s),
including the g0(z) center and the separate g0(z_c) anchor. Fixed d0 ensures
eta is comparable to s, so no second shrinking source width is introduced.

Set J=[0,U;-U*,0], a known skew contraction. Two complete independent genuine-
gradient pair banks on(7), with the SAME captured X,z, composed through J and
a known sqrt(I-JJ*) fill, give selected matrix

    K=c_sel^2 rho^2 D J D.

This is exactly the skew block of C_f,OU/kappa0. The marked pair uses the
complete callback in Section1, independent private banks conditional on the
same passive and X,z. For M=c_sel J_f and O_R=-Sym(J_f C_f,OU), the positive
fork normalizer is

    b c_side=-q kappa0/(2 c_sel^3 rho^2).

Project the whole fork by the fixed physical coisometry P after the full
square construction when a physical n-output kernel is required. Unknown
heat matrices are never evaluated. The target physical orientation before
approximation is positive: if J_f=[A_X,B_X;0,0], A_X=A_X*, then
P[-Sym(J_f Curl J_f)]P*=B_X B_X*. The fork proof itself only uses the full
skew identity, so it does not execute an unknown positive matrix root.

## 5. Candidate quantitative join and explicit scope

With kappa=A^(1+g), a comparable to A^g and b=A^zeta, choose 0<zeta<g and
rho=A^gamma with zeta+gamma<g. The side's weighted actual root/caller rows are
Lambda kappa/(b rho), because eta is comparable to the original clock width
and sum sqrt(w)/s is logarithmic. Its retained output mass is therefore at
most A after a strict/logarithmic guard. The finite side prior, its true
source dimension, inverse s floors and known fill arithmetic remain charged.
The whole source makes no inverse-width sample copies.

The weak heat bias is kappa a e, which is smaller than b kappa e for zeta<g.
Root coefficient frames should be the same A kappa/s and energy kappa e/s as
in the affine fork: each differentiated heated matrix has its input-to-HS
bound, while all undifferentiated matrices have bounded operator norm.
This claim is to be checked on(7)'s complete X,z record before its owned-root
Price theorem is imported.

An independent orientation-only clock can use delta comparable to b. The old
LOW30 forward-target mean at mu=A^(4/3) has allowance A^(7/3)e. To make that
smaller than b kappa e require 1+g+zeta<7/3, in addition to the path guards.
At g=1, any zeta<1/3 leaves room for gamma>0. The known covariance/mean
reference and old keep must be instantiated exactly as in the independent
negative-fork consumer, not supplied by a law that retained no such observer.

This source class has an independent full terminal Gaussian and a hidden
gradient ancestor. The known-correlated case in which the terminal affine
base itself depends strongly on the ancestor is not covered: then (2) need
not retain the extra a. A general nonlinear ancestor map is likewise not
licensed. The proposed gain is local and source-relative; it does not prove
closure under arbitrary subsequent RAW composition or an all-order cost rate.
