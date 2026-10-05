# Independent audit: augmented-score regrouping closes the analytical mixed K gate

2026-10-05 UTC. **Analytical PASS.** No finite producer or full endpoint join is claimed.

Assume g=grad U, g(0)=0 and 0<=Dg<=A I under U in C2. Let the conditional OU history have X_t=e^(-t)x+L_t V, where L_t is its isonormal Brownian row, and let

    I=integral_0^infinity e^(-t)g(X_t)dt,
    C_H(x)=Cov(I|x), K(x)=Cov(I,g(x-I)|x),
    R2 f=integral_0^1 r P_r f dr.

Then the following concrete replacement is valid:

    ||R2[K+C_H Dg]||_(L2(gamma);HS)
       <= [3/8+pi(2-sqrt(2))/16] A^4 sqrt(D)
       <0.491 A^4 sqrt(D).                           (1)

Consequently 2 Sym R2 K can be replaced by -2 Sym R2(C_H Dg), with at most twice this error. This does not assert a pointwise K factorization. The estimate uses an augmented Gaussian score which annihilates both Hessian vertices; it never differentiates Dg, a saved Hessian, or the Stein kernel.

## 1. Exact private Stein representation and orientation

Let V' be an independent copy of V, rho in [0,1], and V_rho=rho V+sqrt(1-rho^2)V'. Write J_t^rho=Dg(e^(-t)x+L_t V_rho), J_u=Dg(e^(-u)x+L_u V), and

    c(t,u)=L_t L_u*=e^(-|t-u|)-e^(-t-u).

The unsymmetrized private Stein kernel has the representation

    tau(V)=integral_0^1 d rho integral_0^infinity dt du
       e^(-t-u)c(t,u) E_(V')[J_t^rho J_u].            (2)

Thus E_V tau=C_H and K=-E_V[tau Dg(x-I)]. The positive scalar weight w=e^(-t-u)c(t,u) has total integral 1/4 (including rho). In particular ||tau||op<=A^2/4. Its realized product need not be symmetric.

Set E=g(x-I)-g(x). Since D_x I is symmetric and bounded by A/2,

    D_x E^T=(I-D_x I)Dg(x-I)-Dg(x),
    K+C_H Dg=-E[tau D_x E^T]
                   -E[tau(D_x I)Dg(x-I)].          (3)

The order of the last product in (3) is important. The second term has HS norm at most A^4 sqrt(D)/8 before applying R2, and therefore at most A^4 sqrt(D)/16 afterward.

## 2. An invariant augmented Gaussian direction

At a fixed outer clock r in (0,1), write x=r z+q G, q=sqrt(1-r^2), with fresh standard G. For a pair t,u>0 put m=min(t,u) and take the scalar Cameron-Martin function

    h_m(a)=sqrt(2)e^a 1_[0,m](a)/(e^(2m)-1).

It satisfies ||h_m||^2=(e^(2m)-1)^(-1) and, for every v>=m, L_v h_m=e^(-v). For v<m the latter response is multiplied by (e^(2v)-1)/(e^(2m)-1), a number between zero and one.

Let alpha=sqrt((1-rho)/(1+rho)), and define one augmented directional derivative for each physical coordinate k:

    A_k=q^(-1)partial_(G_k)-D_(h_m e_k,V)
                              -alpha D_(h_m e_k,V').

Both Hessian vertex vectors Y=e^(-t)x+L_t V_rho and Z'=e^(-u)x+L_u V are invariant under every A_k, because rho+alpha sqrt(1-rho^2)=1. Consequently W=J_t^rho J_u is invariant too. This is an exact invariance of its arguments, and requires no derivative of Dg.

The Gaussian score vector of these directions is

    S=G/q-V(h_m)-alpha V'(h_m),
    Cov(S)=sigma^2 I,
    sigma^2=q^(-2)+kappa^2,
    kappa^2=2/[(1+rho)(e^(2m)-1)].                   (4)

Conditional on z, S has zero covariance with Y and Z'. They are jointly Gaussian, so S is independent of both vertex vectors jointly. Its covariance is unchanged by this conditioning. S is also independent of the exposed z.

## 3. Product-current regrouping and one-energy bound

Gaussian integration by parts in the A_k directions gives exactly

    E[W D_x E^T | z]
       =E[W S E^T | z]+E[W D_h E^T | z],           (5)

where D_h is the simultaneous physical-direction derivative in the unrotated V bank only. The V' derivative of E is zero. Formula (5) can be obtained directly along the invariant Gaussian fibers. There is no product-rule derivative of W at any point.

Condition on the pair of vertex vectors and z. W is then fixed, ||W||op<=A^2, and the coordinates S_k/sigma are an orthonormal family in the conditional scalar L2 space. Conditional Bessel inequality, summed over the output coordinate of E, gives

    ||E[S E^T | Y,Z',z]||HS
       <=sigma (E[|E|^2 | Y,Z',z])^(1/2).          (6)

Multiplication by W, conditional Jensen, and then integration in z imply

    ||E[W S E^T|z]||_(L2(z);HS)
       <= A^2 sigma ||E||2
       <= A^4 sigma sqrt(D).                       (7)

The last step uses |E|<=A|I| and ||I||2<=A sqrt(D), since every unconditional X_v is standard Gaussian. No product of two dimension-sized energies occurs.

For the second term in (5), the explicit response bound L_v h_m<=e^(-v) gives

    ||D_h I||op<=A integral_0^infinity e^(-2v)dv=A/2,
    D_h E=-Dg(x-I)D_h I,
    ||D_h E||op<=A^2/2.

Therefore its HS norm is at most A^4 sqrt(D)/2 for each W, before scalar weight integration.

## 4. Integrability and explicit constant

The scalar integrals are

    integral d rho dt du w=1/4,
    integral d rho dt du w kappa=pi(2-sqrt(2))/8,
    integral_0^1 r/q dr=1,
    integral_0^1 r dr=1/2.

For the second identity, split t<=u and u<=t. On t<=u, w=e^(-2u)(1-e^(-2t)); the double integral with (e^(2m)-1)^(-1/2) is pi/16. Multiply by integral_0^1 sqrt(2/(1+rho))d rho=4-2sqrt(2).

Use sigma<=q^(-1)+kappa. The score term contributes at most

    [1/4+pi(2-sqrt(2))/16] A^4 sqrt(D).

The D_h E term contributes A^4 sqrt(D)/16. The final direct product from (3) contributes another A^4 sqrt(D)/16. Summing proves (1).

Endpoints r=1, min(t,u)=0 and rho=1 can be handled by limits after the displayed integrable majorants. Brownian Malliavin/Stein identities hold by Gaussian Sobolev closure for the stated C1 Lipschitz g. Alternatively smooth anchored convex-gradient approximations retain the bounds and pass to the same limits. No executed path grid is introduced.

## 5. Audit scope and remaining constructive work

The earlier pointwise Stein factorization and a separate bound of its endpoint-score and private-score pieces would be invalid; either can lose a second sqrt(D). Their invariant-score regrouping is the decisive step. The outcome specifically admits -C_H Dg after the outer R2 at order four.

This proof supplies an analytical coefficient target. It does not supply original-gradient VALUE production of C_H Dg, its transpose orientation, a signed mixed action inside a positive law, or the necessary source/caller/variance/floor/count contracts. Those finite gates must still be proved. It also does not replace Dg by Dj without another estimate.

Sources reviewed: CONDITIONAL-INNOVATION-LOCALIZES-THE-COVARIANCE-GATE.md; RESOLVENT-COVARIANCE-STABILITY-AND-FINITE-J-KERNEL.md; SAME-CARRIER-P3-MEAN-AND-FULL-COVARIANCE-GATE.md; THIRD-ORDER-NESTED-FORCE-MEAN-WITH-C2-DECOUPLING.md. The bridge-score mechanism was derived independently and then reconciled with the parent worker's matching derivation.

## 6. Exact author-text review and reproducible checks

The complete author note `AUGMENTED-SCORE-REMOVES-THE-NONLINEAR-K-HISTORY.md`, SHA256 `66aa022a4eac38187d1f277b3eae7a899502fc397f769a6d07e57aaba82e5833`, was read independently after writing Sections 1–5. Its matrix ordering, Stein coupling, invariant direction, conditional Bessel conditioning, regularity closure, constants, and declared finite-target scope all PASS. In particular its all-affine remaining target -4 Sym R2[(R2[B^2])Dg] follows from C_H=2 R2[B^2] because g is a gradient and B=D R1 g is symmetric.

`check_augmented_score_independent.py` runs 771 passing assertions, including 108 exact-rational conditional Gaussian covariance configurations, the transpose-order identity plus a noncommuting wrong-order negative check, high-precision positive-weight/shield integral checks, and the explicit linear-history bridge response. The result is saved in `augmented_score_independent_checks.json` with source pins. The integral values are 1/4 and 0.2300377961276525287..., giving the displayed error constant 0.4900188980638262644.... These diagnostics supplement, and do not replace, the preceding proof. They do not execute any history grid or claim finite VALUE production.
