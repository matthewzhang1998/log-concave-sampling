# Relative energy under a finer auxiliary width: a nondegenerate native subclass

New bounded source lemma, 2026-10-04. This does not close the general matrix orientation producer. It supplies a precise relative-energy comparison for two versions of the ACTUAL finite two-node twin, under extra nondegeneracy assumptions. No gradient convergence or high-order derivative assumption is used.

## Assumptions and literal versions

Use the complete same-record simultaneous twin c8b20bf4, with known ||M||<=1. In addition assume

    <g0(u)-g0(v),u-v> >=m0 |u-v|^2,
    <gt(u)-gt(v),u-v> >=mt |u-v|^2,
    sigma_min(M)>=sigma>0,

where m0,mt,sigma are supplied dimension-uniform numerical positive constants and both original first radii are at most one. Thus both original gradients are globally inverse-Lipschitz with constants m0,mt. The statements allow the SAME noncommuting original potential at both nodes. A singular M or an arbitrary nonconvex original source is outside this lemma.

For any width epsilon>0, let u_k=p0,+^[k], with the literal zero initialization and all Gaussian/caller records unchanged. Put E_K=F-r pt,+^[K]. K is an actual finite integer, not an exact-reference limit. Two widths may use different K>=3, but each keeps its own complete finite version and zero record.

## Complete auxiliary first, uniform in the finite iteration

Write q_k(Z)=S+aM u_(k-1)(Z) for k>=1 and

    Delta_k(Z)=gt(q_k(Z))-gt(q_k(Z)+epsilon Z).

Then

    u_(k+1)(Z)=g0(x+aM*Delta_k(Z)).

The actual auxiliary Lipschitz constants L_k=Lip_Z u_k satisfy

    L_(k+1)<=a epsilon+2a^2 L_(k-1),
    L_0=L_1=0.

Hence, for 2a^2<1,

    L_k<=a epsilon/(1-2a^2),
    Lip_Z q_k<=a^2 epsilon/(1-2a^2).                  (1)

These are direct finite recurrence bounds, independent of convergence of Jacobians.

## The contrast cannot collapse the auxiliary Gaussian

For Z,Z', insert the same q_k(Z) in the second contrast. The inverse-Lipschitz bound for gt and the two changes of q_k give

    |Delta_k(Z)-Delta_k(Z')|
      >=[mt-2a^2/(1-2a^2)] epsilon |Z-Z'|.            (2)

Assume eta0:=mt-2a^2/(1-2a^2)>0; for fixed mt this is a finite small-a guard. Compose (2) successively through M*, g0, M and the outer gt. For every K>=3,

    |E_K(W,Z)-E_K(W,Z')|
      >=c_* r a^2 epsilon |Z-Z'|,
    c_*=sigma^2 m0 mt eta0>0.                        (3)

The original directed F is independent of Z and cancels in this difference. All other original records W and exterior callers are frozen. Neither commutation of Hessians nor pointwise symmetry of a partial derivative is assumed.

The original pointwise upper bound remains

    |E_K(W,Z)|<=r a^2 epsilon |Z|.                    (4)

Thus this subclass has both bounds uniformly over every actual finite K>=3.

## Gaussian energy and width comparison

Use two independent Z,Z' at the SAME W. The conditional variance identity and (3) imply

    Var_Z(E_K|W)=.5 E_(Z,Z')|E_K(W,Z)-E_K(W,Z')|^2
                   >=c_*^2 r^2 a^4 epsilon^2 d.

Consequently, for the actual centered full-source energy,

    e_K=||E_K-E E_K||2 >=c_* r a^2 epsilon sqrt(d).    (5)

This lower bound also holds when the Gaussian W law carries the original correlated C/P geometry, because (3) is uniform in W. A stored mean or origin does not remove it.

Let b,s denote two auxiliary widths epsilon_b,epsilon_s and any actual K_b,K_s>=3. From (4),(5),

    ||E_(K_s,epsilon_s)||2
       <=[epsilon_s/(c_* epsilon_b)] e_(K_b,epsilon_b). (6)

The left side is the UNcentered actual finer remainder. Higher fixed moments satisfy the same ratio times the explicit Gaussian norm factor ||Z||p/sqrt(d). There is no inverse tiny force coefficient and no unproved heat-grade interpolation.

At r=a=A, epsilon_b=A^.9 and epsilon_s=A^.95, equation(6) supplies a genuine C A^.05 relative mark in this subclass. This is stronger than a comparison of nominal sqrt(d) envelopes, but does not make the full common-gradient lift's ordinary body small.

## Scope for a possible bounded mean bypass

If the separately admitted common coisometry realizes D=H_s-H_b as the physical projection of the difference of two genuine-gradient reference lifts, then E_b=D+E_s exactly on the complete common tape. A high-order projected gradient mean for D and an independent near-gradient mean for E_s can exploit (6). The finite-versus-exact-gradient VALUE restoration, actual D-body retention, carrier/zero paths and all work must still be supplied by that constructor. A marginal gradient mean is not by itself a standalone negative Cov(D) reserve.

For the displayed exponents the common full-gradient radius is A^.05, so this can at most support a bounded below-one auxiliary exponent refinement under the existing small-radius guard. It is not an arbitrary-order closure mechanism. The general matrix K3 same-query producer, singular/ill-conditioned M, and original sources without the stated inverse-Lipschitz constants remain outside this lemma.

The adjacent checker passes9,720 finite-source tests, including K=2 as the
excluded exact-zero case and K=3 through15 for noncommuting original
Hessians and nonsymmetric well-conditioned M. These diagnostics check the
literal recurrence; the independent-copy variance identity supplies(5).
