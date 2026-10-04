# Exact leading curl of the third finite weak twin

New analytical reduction, 2026-10-04. It binds P's literal VALUE source `THIRD-ITERATE-TWO-NODE-WEAK-TWIN-PORT.md`, SHA c8b20bf4eff6986f5e1f2d749e740beb935491a1043153ec62b39ef05bf4bca3, and keeps its complete finite record. It does not yet execute the remaining leading curl coefficient.

## 1. The same-record coordinates and exact derivatives

Use the source's independent standard coordinates (S,U,Z), and write c=c0, s=s0, c^2+s^2=1. These are the FIXED original graph coefficients, not the later OU clock. Put

    x=cS+sU,
    Delta=gt(S)-gt(S+epsilon Z),
    d=(a/s) M* Delta,
    x1=x+aM*Delta=cS+s(U+d),
    t0=S+aM g0(x),       t1=S+aM g0(x1),
    K0=Ht(t0) M H0(x),   K1=Ht(t1) M H0(x1),
    D=Ht(S)-Ht(S+epsilon Z),   Hplus=Ht(S+epsilon Z).

Here Hi are the ANALYTICAL original Jacobians at their literal recorded queries. The executor continues to use original gradient VALUES. All blocks have their actual moving anchors. Define E=r[gt(t0)-gt(t1)], and the full square lift f=(E,0,0). Direct differentiation of this finite expression gives

    J_S E = r[Ht(t0)-Ht(t1)] +ra c(K0-K1)-ra^2 K1 M* D,
    J_U E = ra s(K0-K1),
    J_Z E = ra^2 epsilon K1 M* Hplus.                    (1)

In particular the same terminal/twin pair occurs in D and Hplus. There is no independent replacement of either occurrence.

The first term in J_S is symmetric, hence drops out of Curl f. Set DeltaK=K0-K1. The exact skew matrix has the decomposition C=C_lead+C_fb, where

    C_lead=ra [ c(DeltaK-DeltaK*), s DeltaK, 0;
                         -s DeltaK*,        0, 0;
                                  0,        0, 0 ],     (2)

and

    C_fb=ra^2 [ -K1 M*D +D M K1*, 0, epsilon K1 M*Hplus;
                              0, 0,                    0;
                    -epsilon Hplus M K1*, 0,            0 ]. (3)

Since ||Ki||<=1, ||D||<=2 and epsilon<=1, ||C_fb||op<=5ra^2. The bound is uniform in every actual record and caller, with no third derivative. Formula(3) keeps its complete aliases; the norm bound alone is enough only for the constant-target estimate below.

## 2. Exact shifted-Jacobian form of the leading term

For fixed S define

    Phi_S(u)=gt(S+aM g0(cS+su)).

Then Phi_S has first in u at most as and

    D_u Phi_S(U)=as K0,
    D_u Phi_S(U+d)=as K1,
    E=r[Phi_S(U)-Phi_S(U+d)].                           (4)

Conditional on (S,Z), U is genuinely independent Gaussian and d is fixed. Phi_S is generally NOT a gradient in u. Its S first is O(1), so freezing S removes the full independent terminal Gaussian used by the already proved hidden-gradient-ancestor port. The remaining leading coefficient is therefore the full shifted Jacobian difference in (4), not an arbitrary independent terminal heat or a product of averaged primitive Hessians.

A direct Gaussian density shift estimate can cost |d| divided by the innovation width. As |d| has an ambient-dimensional Gaussian profile, that estimate by itself is not a dimension-free one-energy calibration. No extra epsilon or a gain for (2) is inferred from the small VALUE difference.

## 3. What can already be omitted in the constant orientation target

Let e=||f-Ef||_2 and A bound its complete first. At each positive whole-input clock, denote its innovation width by v_j, its common coarse root by X, and let

    J_j(X)=E_V Df(c_j X+v_j V),
    C_fb,j(X)=E_V C_fb(c_j X+v_j V).

The marked and curl-side copies of V may be independent CONDITIONAL ON X because this is a product of their exact conditional expectations. Every query within either copy keeps the source's full internal aliases. Gaussian Bessel gives ||J_j||_(L2(X);HS)<=e/v_j. Thus, for the held positive clock,

    O_fb=-Sym sum_j w_j E_X[J_j(X) C_fb,j(X)]

satisfies

    ||O_fb||HS <=5ra^2 e sum_j w_j/v_j <=Lambda ra^2 e,
    ||O_fb||op <=5A ra^2 sum_j w_j <=Lambda A ra^2.      (5)

The same argument applies to the exact continuous covariance clock. It uses the actual whole E mark, not the larger energy of Delta or an individual terminal. No differentiated calibration error enters.

Consequently, an independent constant-reference orientation repair at target b(ra)e can discard this constant subleading contribution whenever a<=C b. At a=A^g and b=A^zeta with 0<zeta<g, that inequality has strict power slack. A fixed-gap Gaussian covariance comparison prices the constant omission by (5).

This does NOT grant a conditional-root law for C_fb, nor permit deleting it from an already observed random endpoint without its current. It is the constant target reduction used before a NEW independent fork is compared to an old constant forward reference.

## 4. Precise remaining constructive interface

The requested next producer may concentrate on (2): an actual whole-input heat channel for the SAME-root shifted derivative D Phi_S(U)-D Phi_S(U+d), with d=(a/s)M*[gt(S)-gt(S+epsilon Z)], contracted on its left with the intact J_j of the COMPLETE E source. It must preserve the joint terminal/twin record, source-zero anchors, complete private/caller paths and full tape cost. A bounded stationary law for one feedback marginal does not furnish this channel.

The exact joint-gradient feedback lift in P's `JOINT-TERMINAL-TWIN-GRADIENT-FEEDBACK-LIFT.md`, SHA 4a4d7a9ebdad83d1552a8278b0ede491e4be1e90161769d17a04808e9113b201, is a valid source component: its two blocks are a(Delta,-epsilon gt(S+epsilon Z)), the gradient of a[Vt(S)-Vt(S+epsilon Z)]. Its full first is at most3a and its anchored Gaussian energy is O(a epsilon sqrt(d)). Under whole-input heat the two terminal queries have Gram v^2[[I,I],[I,(1+epsilon^2)I]]. The lower companion and its actual captured-S profile must remain. This lift still needs a joint observer theorem before insertion inside the nonlinear g0 of (4).

No claim is made that later K iterates have exactly this two-variable gradient lift. Their finite unrolled feedback genealogy and numerical version must be re-expressed explicitly.
