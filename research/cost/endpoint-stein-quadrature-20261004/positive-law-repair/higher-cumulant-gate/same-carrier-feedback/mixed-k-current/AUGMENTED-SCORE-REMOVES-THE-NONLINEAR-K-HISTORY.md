# An augmented Gaussian score removes the nonlinear K history

2026-10-05. Analytical reduction. Independently derived and checked by the mixed-current analyst; full written audit is separate. This is not a finite VALUE producer.

## Result

Let X_t^x=e^{-t}x+Y_t be the genuine conditional OU history, I_x=∫_0^∞e^{-t}g(X_t^x)dt, C_H(x)=Cov(I_x), and K(x)=Cov(I_x,g(x-I_x)). Assume only g=∇U, g(0)=0 and 0≤Dg≤AI. With R2 f=∫_0^1 rP_r f dr,

    ||R2[K+C_H Dg]||_{L2(γ);HS} ≤ C A^4 sqrt(D).       (1)

Consequently the audited covariance localization becomes

    Cov(F2|Z)=Cov(H_j|Z)−2 Sym R2[C_H Dg](Z)
                                      +O_{L2HS}(A^4 sqrt(D)).    (2)

The remaining finite mixed target is an all-affine cubic matrix word:

    −4 Sym R2[(R2[B²])Dg],  B=D R1 g.                (3)

This theorem does not make C_H, B or Dg executable coefficients. It does remove the true nonlinear conditional history from the target, while pricing rather than dropping its effect. In particular it does not replace the old history with a cheap common-G path.

## 1. A derivative-free representation of the conditional Stein matrix

Write Y_t=L_t V for the standard isonormal Brownian input V. Then

    L_t h = sqrt(2)∫_0^t e^{-(t-v)}h(v)dv,
    L_t L_u*=c(t,u)I,
    c(t,u)=e^{-|t-u|}−e^{-(t+u)}.

Let V' be independent and V_ρ=ρV+sqrt(1−ρ²)V'. Define J_t^ρ=Dg(e^{-t}x+L_tV_ρ), J_u=Dg(e^{-u}x+L_uV). The Gaussian Riesz formula gives the analytical Stein matrix

    τ=∫_0^1 dρ ∫_0^∞dt du e^{-(t+u)}c(t,u) J_t^ρ J_u,             (4)

where its V' dependence is averaged whenever used as a conditional Stein kernel. Thus

    K=−E[τ Dg(x−I)],  Eτ=C_H.                       (5)

All orderings in (4)-(5) are intentional. Neither each integrand nor τ is asserted symmetric. Formula (4) uses only the bounded first derivative of g in its analytic certificate. Its total positive weight is 1/4.

Set E(x,V)=g(x−I_x(V))−g(x), H*=Dg(x−I), Jx=D_x I. Since Jx is symmetric and ||Jx||op≤A/2,

    D_xE=(H*−Dg(x))−H*Jx,
    H*−Dg(x)=D_xE*+Jx H*.                          (6)

Therefore

    K+C_HDg=−E[τ D_xE*]−E[τ Jx H*].                (7)

The second term is immediately O(A^4 sqrt(D)) because ||τ||op≤A²/4. We must handle the first without differentiating τ.

## 2. The early-history Cameron shift

For m=min(t,u)>0 put

    h_m(v)=sqrt(2)e^v 1_{[0,m]}(v)/(e^{2m}−1),
    ||h_m||²=1/(e^{2m}−1),
    α_ρ=sqrt((1−ρ)/(1+ρ)).

For v≥m, L_v h_m=e^{-v}; for v<m,

    L_v h_m=e^{-v}(e^{2v}−1)/(e^{2m}−1),

which lies between zero and e^{-v}. Shift V in direction h_m a and V' in direction α_ρ h_m a. This shifts both original Hessian vertices in (4) by exactly e^{-t}a and e^{-u}a, respectively. The joint shift therefore reproduces a derivative in the endpoint x at both vertices, without ever evaluating that derivative.

Its action on the full unrotated history is bounded using only Dg:

    ||D_h I||op≤A∫_0^∞ e^{-2v}dv=A/2,
    D_h E=−H*D_h I,  ||D_h E||op≤A²/2.             (8)

The second-bank direction acts trivially on E. The crucial difference from endpoint differentiation is the order-A² bound in (8).

## 3. Include the outer Gaussian before integration by parts

Fix an outer clock r<1, exposed Z=z, and put q=sqrt(1−r²), x=rz+qG, where G is standard independent of V,V'. Let W=J_t^ρJ_u. Consider the Gaussian directional operator

    A_a=q^{-1}∂_{G,a}−D_{V,h_m a}−α_ρD_{V',h_m a}.

It annihilates BOTH Hessian-vertex arguments. Hence it annihilates W distributionally, with no derivative of Dg required. Its Gaussian score is the D-vector

    S=G/q−V(h_m)−α_ρV'(h_m),
    Cov(S)=σ²I,
    σ²=q^{-2}+2/[(1+ρ)(e^{2m}−1)].                 (9)

Conditioned on z, S has zero covariance with each of the two Hessian vertex vectors, exactly because their directional changes cancel. Joint Gaussianity therefore makes S independent of the pair of vertices and of W. This is independence of the augmented score only, not independence of the original history.

Gaussian integration by parts, with W constant in this direction, yields the exact matrix identity

    E[W D_xE*|z]=E[W S E*|z]+E[W(D_hE)*|z].         (10)

This can be proved for smooth g first. Alternatively, condition on the two vertex vectors; the directional Gaussian integration by parts then never differentiates W at all. That proves (10) directly for the admitted bounded continuous Hessian, by Gaussian Sobolev closure.

## 4. Conditional Bessel pays one physical energy

Condition on z and on the two Hessian vertices. The coordinates S_i/σ remain an orthonormal Gaussian family. Applying Gaussian Bessel separately to each physical component of E gives

    ||E[S E* |z,vertices]||HS
          ≤σ (E[|E|²|z,vertices])^{1/2}.             (11)

As W is now fixed and ||W||op≤A², Jensen gives

    ||E[W S E*|z]||HS≤A²σ ||E||_{L2|z}.             (12)

No norm of the D-vector S is taken; that would introduce a spurious sqrt(D). The conditional chord energy satisfies

    ||E||_{L2|z}≤A||I||_{L2|z}
                     ≤A²(|z|+sqrt(D)),              (13)

and after integrating z~γ the sharper unconditional bound ||E||2≤A²sqrt(D) follows from stationary Gaussian marginals and Minkowski. The last term in (10) is at most A^4sqrt(D)/2 in HS by (8).

It remains only to integrate the positive weight in (4) and the r dr outer weight. The elementary estimates

    ∫_0^1 r/sqrt(1−r²)dr=1,
    ∫_0^1dρ∫_0^∞dt du e^{-(t+u)}c(t,u)=1/4,
    ∫_0^1dρ∫_0^∞dt du e^{-(t+u)}c(t,u)
       sqrt(2/[(1+ρ)(e^{2min(t,u)}−1)])
                     =π(2−sqrt(2))/8               (14)

show that σ≤q^{-1}+sqrt(2/[(1+ρ)(e^{2m}−1)]) is integrable. Equations (7)-(14) prove (1), with an absolute constant. A conservative explicit value is 3/8+π(2−sqrt(2))/16<0.491. Constants may be loosened without changing the claim.

## 5. Admitted regularity, genealogy and exact remaining gate

All continuous histories in this proof are analytical comparison objects. Finite smooth positive path averages can establish Gaussian Sobolev identities and pass to the history limit with the displayed uniform first/energy majorants; they are not an executed path-grid algorithm. The bridge shift preserves the exact conditional OU covariance c(t,u). The ρ rotation is the genuine Gaussian Riesz coupling of the complete private path.

No Hessian is an executed producer, no saved HVP is differentiated, no small conditional-mean field is differentiated as if its derivative were small, and no strong conditional expectation is supplied. The coherent high-dimensional shift inside E is retained until its one-energy chord bound is paired with the augmented score.

The source-qualified finite task is now (3), to error A^4sqrt(D), inside a fixed positive variance reserve. Positive clocks for its nested resolvents and finite original-g VALUE actions must still include the unsmoothed original vertex Dg(x), the correct common centers, noncommutative matrix orientation, complete replay roots, actual caller/first/origin ports, absolute precision floors and mode restoration. This analytical reduction alone is not a full same-carrier order-four PASS.
