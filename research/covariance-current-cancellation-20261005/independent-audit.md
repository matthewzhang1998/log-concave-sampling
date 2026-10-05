# Independent audit: exact covariance/time-current cancellation

2026-10-05. Analytical statement for the actual positive path. The proposed identity is correct. It removes the claimed covariance/time contribution of order L⁴e at unit buffer. It does not remove the ordinary quadratic-in-cumulants feedback, whose leading grade is L⁵e.

## 1. Conventions and exact identity

Let C=C_t=(v/2)I+(1−t²)Σ, S=C^(1/2), x=SZ1, and let

- q_r,t(x)=K_r:H_(r−1)^C(x), preserving K_r's first physical output;
- p_t=Σ_(r=3)^m a_r(t)q_r,t, with finite m;
- J_ia=∂_a p_i and H=I+J;
- D_t=Σ_r a_r(t) ∂_t[q_r,t(S_t Z1)], with Z1 fixed;
- W=tB+x+p_t(x)+ηZ2, with B,Z1,Z2 independent.

The derivative in D_t does not differentiate a_r. Coefficient derivatives remain in the other path currents. No symmetry of J is assumed. The Hermite convention is the score convention E[H_n^C(X)F(X)]=E[D^nF(X)].

Write ρ_C for the centered Gaussian density. Since H_n^Cρ_C=(−1)^nD^nρ_C and C′=−2tΣ,

    ∂_t(q_r,t ρ_C)=−t Σ:D²(q_r,t ρ_C),
    ∂_tρ_C=−t Σ:D²ρ_C.

Expanding this weighted heat equation gives the fixed-x derivative

    (∂_t q_r,t)(x)=−t Σ:D²q_r,t(x)
                        +2t Dq_r,t(x) Σ C^−1x.

Here Σ and C commute, so S′S^−1=−tΣC^−1. Adding the moving-root derivative therefore gives

    ∂_t[q_r,t(S_tZ1)]
        =−t Σ:D²q_r,t(x)+t Dq_r,t(x)ΣC^−1x.

Consequently D_i=−tΣ_ab∂_abp_i+tJ_iaΣ_ab(C^−1x)_b. Condition on B,Z2 and use Gaussian integration by parts in x:

    E[(C^−1x)_b J_ia ∂_iφ(W)]
      =E[(∂_bJ_ia)∂_iφ(W)
                 +J_ia H_jb ∂_ijφ(W)].

The derivative-of-J term cancels the Hessian-of-p term exactly. Thus

    E[D_t·∇φ(W)]=t E[(JΣHᵀ):D²φ(W)].                  (1)

Because D²φ is symmetric, (JΣ):D²φ=(ΣJᵀ):D²φ, even when J is not symmetric. This proves

    E[D_t·∇φ(W)]−t E[(ΣJᵀ):D²φ(W)]
      =t E[(JΣJᵀ):D²φ(W)].                            (2)

This is equality of currents under expectation, not a pointwise identity for D_t. It holds for every finite collection of the indicated Hermite maps, with arbitrary differentiable scalar coefficients a_r. In the consumer, a_r=1−t^r.

The argument is valid for smooth compactly supported tests by polynomial Gaussian integrability; truncation/Sobolev approximation extends the same finite identities under the hypotheses in the source notes. No derivative of B beyond those required by the already imported Stein recurrence is introduced.

Matrix scope: the stated isotropic-buffer path works for every symmetric positive semidefinite Σ. Do not apply the intermediate principal-square-root formula S′S^−1=−tΣC^−1 to a noncommuting covariance path. Such a path requires its actual root-transport matrix and corresponding covariance-root current. No commutation between J and Σ or C is used here.

## 2. Explicit one-energy bound for the new rank-two current

Let h_r=||K_r||HS and let k_r bound every nonempty proper physical flattening of K_r. Put λ=λ_min(C)>0. Define J_r=Dq_r,t(x), so J=Σ_r a_rJ_r. With n=r−2,

    (J_r)_ia=(r−1) Σ_(j,u1,...,un)
        K_(r;i,j,u1,...,un)(C^−1)_ja H_n^C(x)_(u1,...,un).

For r,s≥3, set n=r−2, d=s−2, and

    β_rs=(r−1)(s−1)
       [Σ_(p=0)^min(n,d) {p! binom(n,p)binom(d,p)}²
                                     (n+d−2p)!]^(1/2).

Then the dimension-free estimate is

    ||J_r Σ J_sᵀ||_(L²;HS)
      ≤ β_rs ||Σ||op λ^(-(r+s)/2)
                           min(h_r k_s,h_s k_r).       (3)

Proof. The product formula for the two Gaussian Hermite factors contracts p pairs of Gaussian slots, leaving chaos degree n+d−2p. The deterministic graph has two K vertices, one derivative-slot edge C^−1ΣC^−1, p Gaussian edges C^−1, and the surviving Gaussian-slot transforms C^−1/2. Both vertices keep their physical output slot. Start either vertex in HS norm and add the other through all its existing edges. That attachment uses 1+p≥1 slots and leaves at least the physical output, hence is a proper cut. HS-times-operator multiplication pays one h and one k, with edge-transform bound

    ||Σ||op λ^−2 λ^−p λ^(-(n+d−2p)/2)
      =||Σ||op λ^(-(r+s)/2).

There is no self-trace or second Hilbert energy. Symmetrization of surviving Gaussian slots decreases HS norm. Hermite isometry supplies √((n+d−2p)!), and different p have orthogonal chaos degrees. Their squared sum gives β_rs. Choosing the better anchor proves the minimum in (3).

In particular β_33=4√3, β_34=6√10, and β_44=18√15. Summing (3) yields, at v=1 with bounded coefficients |a_r|≤1,

    ||JΣJᵀ||_(L²;HS)≤C_m ||Σ||op h s≤C_m L² h s,
    h=Σ_(r=3)^m h_r,  s=Σ_(r=3)^m k_r.               (4)

The reserve transfers this rank-two current to a velocity with factor η^−1. Since λ≥v/2 and η=√(v/2), an (r,s) covariance pair costs at most

    C_m L² min(h_rk_s,h_sk_r) / v^((r+s+1)/2).         (5)

The old linear bill L²h is replaced by the quadratic bill L²hs. Smallness is unnecessary for identities (1)–(3); it is used only to simplify the final finite-rank bound.

## 3. Corrected current generator and scaling

In the all-rank generator, replace A1=D_t and the rank-two summand −tΣJᵀ together by

    A1=0,
    A2 receives +tJΣJᵀ.

Every partition-generated term, every coefficient-derivative term, and the terminal Stein current remain unchanged. All are still evaluated at the same W_t, and JΣJᵀ is independent of the reserve Z2.

Thus the corrected version of the unit-buffer bound (A8) is

    W2(Law(W1),Law(W0))
      ≤C_m [L^m e + L² h s
           +Σ_(r=3)^m h_r Σ_(j=1)^(r−1) s^j].        (6)

For general v, set L′=L/√v, e′=e/√v, h′_r=h_r/v^(r/2), k′_r=k_r/v^(r/2), apply (6) to the primed quantities, and multiply by √v. This supplies every width loss, including (5).

The existing actual-cumulant estimates h_r≤L^(r−1)e/r and k_r≤c_rL^r/r! now give, under a fixed-rank computable small-radius guard,

    W2≤C′_m e [(L/√v)^m+(L/√v)^5+(L/√v)^7]
       ≤C″_m e (L/√v)^min(m,5).                       (7)

The leading ordinary nonlinear pair has bill h3k3/v^(5/2)=O(L⁵e/v^(5/2)). The leading covariance pair now has bill L²h3k3/v^(7/2)=O(L⁷e/v^(7/2)). More generally ordinary (r,s) pairs have grade L^(r+s−1)e, while covariance (r,s) pairs have grade L^(r+s+1)e, in unit-buffer normalization.

Consequences:

- m=3: leading bound L³e/v^(3/2), from the Stein remainder.
- m=4: leading bound L⁴e/v², still from the Stein remainder. The existing fourth-rank theorem remains valid, but its attribution of a leading L²h3 covariance/time bill is unnecessarily pessimistic.
- m=5: leading bound L⁵e/v^(5/2).
- m≥6: the same raw-map construction is bounded at grade L⁵e/v^(5/2), with no surviving L⁴e covariance/time ceiling.

For the fourth-rank consumer specifically, the new covariance bills are

    L²h3k3/v^(7/2),
    L²(h3k4+h4k3)/v⁴,
    L²h4k4/v^(9/2),

in addition to the unchanged ordinary nonlinear bills and Stein remainder.

## 4. A real raw-map obstruction remains at grade five

The improved grade-five ceiling is not solely an artifact of these norm estimates. At t=0, Gaussian-chaos orthogonality gives

    Cov(W0)=vI+Σ+E[p0(x0)p0(x0)ᵀ],

because every q_r has degree r−1≥2, is orthogonal to x0, and distinct r have orthogonal chaos degrees. For a scalar source,

    Var(W0)−Var(W1)
       =Σ_(r=3)^m (r−1)! K_r²/C0^(r−1),
    C0=v/2+Var(B).

Take B=εb(G) for a fixed centered smooth Lipschitz scalar b with κ3(b)≠0. Then the leading variance surplus is

    2K3²/C0²=κ3(B)²/(18C0²)=Θ(ε⁶/v²)

as ε/√v→0. The elementary second-moment lower bound gives

    W2≥|√Var(W0)−√Var(W1)|=Θ(ε⁶/v^(5/2)).

Since L and e are both proportional to ε, this is precisely the grade L⁵e/v^(5/2). For example b(G)=cos(G)−exp(−1/2) is a smooth Lipschitz centered choice with nonzero third cumulant. Higher raw Hermite ranks add nonnegative variance and cannot cancel this obstruction. Higher-order targets therefore require correcting the polynomial reference (or an equivalent constructive cancellation), rather than retaining the old linear covariance/time bill.

## Audited sources

- /workspace/shared/fourth-cumulant-return-20261005/comparison/POSITIVE-FOURTH-CUMULANT-BUFFERED-CONSUMER.md
- /workspace/shared/rank-indexed-positive-returns-20261005/ALL-RANK-APPELL-AND-POSITIVE-CONSUMER.md

No source file was modified. These conclusions are analytical and do not supply executable cumulant tensors, original-VALUE producers, or caller-path regularity.

## 5. Additional audit: leading skew-square correction target

This separate weak-order computation is also correct. Fix C>0, Q=C^−1, x∼N(0,C), an independent Gaussian reserve ηZ2, and U=x+ηZ2. For a symmetric rank-three K, put

    p_i=K_iab H2_ab^C(x),
    N_ij=K_iab Q_ac Q_bd K_jcd,
    M_ijbd=Σ_(a,c) K_iab Q_ac K_jcd.

Repeated a,b,c,d are summed in N; only a,c are summed in M. The exact Hermite product is

    H2_ab H2_cd=H4_abcd
       +Q_ac H2_bd+Q_ad H2_bc+Q_bc H2_ad+Q_bd H2_ac
       +Q_ac Q_bd+Q_ad Q_bc.

After contracting K⊗K and integrating by parts at U, its four single contractions agree, as do its two double contractions. Therefore

    (1/2)E[p_i p_j ∂_ijφ(U)]
      =(1/2)(K⊗K):E D⁶φ(U)
         +2M:E D⁴φ(U)+N:E D²φ(U).                    (8)

Set

    s_i(x)=−N_ij H1_j^C(x)−2M_ijbd H3_jbd^C(x).       (9)

Then E[s·∇φ(U)]=−N:ED²φ(U)−2M:ED⁴φ(U). More precisely, with a formal scalar amplitude δ and suitable smooth tests,

    Eφ(U+δp+δ²s)
      =Eφ(U)+δ K:E D³φ(U)
         +(δ²/2)(K⊗K):E D⁶φ(U)+O(δ³).                (10)

Thus (9) removes the unwanted connected rank-two/rank-four terms at second weak order and preserves the intended disconnected rank-six term. Only the relevant symmetrization of M contributes to (9) or (8); M need not be fully symmetric as written.

For λ=λ_min(C), the same marked two-vertex argument proves

    ||N||HS≤h k λ^−2,     ||M||HS≤h k λ^−1.

The linear and cubic Gaussian chaoses in (9) are orthogonal. Hermite isometry therefore gives the slightly sharper one-energy estimate

    ||s(x)||₂≤5 h k λ^(-5/2).                         (11)

In scalar C=1, p=aH2 indeed gives s=−a²H1−2a²H3, with ||s||₂=5a². The map x↦x+p+s followed by an independent Gaussian reserve defines a genuine positive law without an invertibility assumption.

Qualification: (8)–(10) use derivatives at the Gaussian base U, not at the nonlinear positive endpoint. They certify the leading analytical correction target; they do not establish an exact cancellation of all nonlinear path currents, a higher-order Wasserstein theorem, an all-order iteration, or a native original-VALUE port. Cross terms with other cumulant maps and higher Taylor terms still require their own complete consumer analysis.

## 6. General covariance/root gauge

The exact current reduction also has a gauge-independent formulation. Let C(t)=S(t)S(t)ᵀ>0 with differentiable invertible S, Ω=S′S^−1, and A=ΩC. Then C′=A+Aᵀ. With x=S Z and the same convention excluding derivatives of the scalar a_r,

    D=(1/2)C′:D²p +Dp(Ω−C′C^−1)x.

Since ΩC−C′=−Aᵀ, Gaussian integration by parts makes the coefficient of D²p equal to (A−Aᵀ)/2, whose contraction vanishes. Hence

    E[D·∇φ(W)]=−E[(J Aᵀ(I+J)ᵀ):D²φ(W)].

The actual moving-root current is A(I+J)ᵀ:D²φ. Its baseline is (C′/2):D²φ, and its feedback is AJᵀ:D²φ. Adding that feedback to the D-current leaves

    −(1/2)E[(J C′ Jᵀ):D²φ(W)].

Indeed AJᵀ−JAᵀ is antisymmetric, and only Sym(A)=C′/2 contributes in the remaining quadratic term. Substituting C′=−2tΣ recovers +tJΣJᵀ. This extension does not permit replacing the actual root-transport feedback by a commuting-root formula on a noncommuting covariance path.

## 7. Final source review and content pin

Reviewed source:

    /workspace/shared/covariance-current-cancellation-20261005/EXACT-COVARIANCE-CURRENT-REDUCTION.md

SHA-256:

    c86de158d9b7a1203c64d0563517a42803d7db4aeef979ac30a9b99adf45f32b

The pinned source passes the analytical review above, including its arbitrary-root version (3a), explicit beta_rs and reserve-width coefficient, all-rank small-radius rate, bounded-Hessian scalar gradient counterexample, and qualified weak-order skew-square correction. Its Gaussian-base wording has been corrected so that (10) is explicitly at U=x+reserve rather than the nonlinear W_t.

The deliberately loose numerical rate constant is sufficient. Indeed for L≤1/2, h≤(2/3)L²e; s≤M_m L³≤1/2; the nonlinear sum is at most 2hs; and L²hs≤(1/4)hs. Thus the source's 8 C_m(1+M_m)² dominates the final rate. The new covariance-pair constant is also smaller than the imported C_m: the crude bound beta_rs≤m³2^m m! sqrt((2m)!) together with its width factor is already far below 2^(10m²)((2m²)!)⁴ for m≥3. Doubling the imported constant safely covers the unchanged consumer and the new covariance term.

The six-force correction graph's stated mark/cycle counts are consistent: joining two trees with one edge gives a tree with four remaining legs; joining them with two edges gives one cycle and two remaining legs. No new executed program is introduced by the current regrouping itself. This review does not independently re-admit the imported native producers, root/caller interfaces, or pre-existing cycle gates; the source explicitly leaves those interfaces open and does not claim otherwise.

The symbolic checker's results are separate evidence; this source pin records the analytical review, not an unobserved checker outcome.
