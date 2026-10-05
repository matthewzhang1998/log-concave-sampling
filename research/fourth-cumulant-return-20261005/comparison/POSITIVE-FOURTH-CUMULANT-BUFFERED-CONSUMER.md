# Positive fourth-cumulant buffered consumer

2026-10-05. Analytical positive-law comparison under C1/Gaussian assumptions. Every residual current below is evaluated at the same actual positive endpoint. All tensors and polynomial maps are analytical references, not executable tensor-source ports.

## 1. Result and normalization

Let B=f(G)-Ef(G) be a centered R^D-valued map of a standard Gaussian bank, with f C1, complete Lipschitz constant L, and e=||B||₂. Write

    Σ=Cov(B),    K3=κ3(B)/6,    K4=κ4(B)/24.

Fix v>0 and split its independent Gaussian buffer equally. Put

    C_t=(v/2)I+(1-t²)Σ,    S_t=C_t^(1/2),    η=√(v/2).

For x in R^D let H_r^C(x) be the symmetric covariance-C Gaussian Hermite score, characterized by

    E[H_r^C(X) F(X)]=E[D_x^r F(X)],   X~N(0,C).

In particular H2^C(x)=(C^(-1)x)^⊗2-C^(-1), and H3^C is the cubic product minus its three covariance-linear contractions. Define

    q3,t(x)=K3:H2^(C_t)(x),
    q4,t(x)=K4:H3^(C_t)(x).

With independent Gaussian roots Z1,Z2, independent also of G, the positive reference is

    R_v=S_0 Z1+q3,0(S_0 Z1)+q4,0(S_0 Z1)+ηZ2.         (1)

There are numerical c,C>0 such that, whenever L/√v<=c,

    W2(Law(B+√v Z),Law(R_v)) <= C L⁴e/v².              (2)

A more informative bill, before absorbing the smaller terms, is

    C[L⁴e/v² + L⁵e/v^(5/2)
                 + L⁶e/v³ + L⁷e/v^(7/2)],             (3)

plus higher products controlled by these terms when L/√v<=c. The leading term includes the covariance/time-path correction L²h3; the second term includes L²h4 and the leading b3² feedback. The mixed b3b4 feedback has scale L⁶e/v³, and b4² has scale L⁷e/v^(7/2).

No logarithmic dimension loss is necessary for this fixed-degree analytical consumer. Section 5 proves its finite Wick-contraction bounds directly. No claim about native original-VALUE production, clipping, root ownership, or caller-path regularity follows from (2).

## 2. Tensor bounds available from C1/Gaussian assumptions

For r=3,4, let h_r=||K_r||HS and let k_r be the maximum operator norm over every nonempty proper bipartition of its physical slots. Full symmetry means that for K3 there is only the 1|2 cut, up to transposes/permutations; for K4 one needs the 1|3 and 2|2 cuts.

The centered Stein recurrence gives

    h3<=L²e/3,    h4<=L³e/4.                             (4)

The corresponding proper-cut bounds also follow from the same Gaussian hypotheses:

    k3<=L³/3,    k4<=L⁴/4.                              (5)

Here is a direct proof of (5). Gaussian Poincaré gives Σ<=L²I. The centered matrix Q=B⊗B-E[B⊗B] has covariance operator at most 4L⁴ on HS matrix space, as proved by differentiating B^T M B once. Consequently the 1|2 flattening of κ3=E[B⊗Q] has operator norm at most 2L³.

For κ4, define W=B^⊗3 minus its three covariance-linear contractions. The proof in FOURTH-CONDITIONAL-CUMULANT-COMPARISON.md gives Cov(W)<=36L⁶I. Thus the 1|3 flattening κ4=E[B⊗W] has operator norm at most 6L⁴. In the 2|2 flattening,

    κ4=Cov(Q)-Σ⊗Σ-(Σ⊗Σ) followed by input-slot exchange.

The three summands have operator norms at most 4L⁴, L⁴ and L⁴. Therefore this cut also has norm at most 6L⁴. Division by 6 or 24 proves (5). These estimates are analytical consequences, not extra tensor-source assumptions.

We prove (2) at v=1 first. Then η=1/√2 and C_t>=I/2. All coordinate transforms below have uniformly bounded operator norm. Scaling is restored in Section 7.

## 3. Exact positive path and the covariance feedback

Let

    a_t=1-t³,    d_t=1-t⁴,
    p_t(x)=a_t q3,t(x)+d_t q4,t(x),
    x_t=S_t Z1,
    W_t=tB+x_t+p_t(x_t)+ηZ2.                            (6)

The endpoints are W_0=R_1 and W_1=B+Z. Every W_t is a genuine probability law, and ηZ2 is untouched and independent.

At x=x_t define

    J_t=D_x p_t(x),    H_t=I+J_t,
    E_t=D_x²p_t(x),    F_t=D_x³p_t(x).

For time derivatives write

    dot b3,t(z)=∂_t[q3,t(S_t z)] with z held fixed,
    dot b4,t(z)=∂_t[q4,t(S_t z)] with z held fixed,
    D_t=a_t dot b3,t(Z1)+d_t dot b4,t(Z1).              (7)

The phrase “z held fixed” matters: these derivatives include the moving covariance and the moving visible Gaussian root. They are not derivatives with x held fixed.

Because C_t is an affine function of Σ,

    dot S_t=-t ΣS_t^(-1).

Gaussian integration by parts in Z1 gives exactly

    E[(dot S_t Z1)·∇φ(W_t)]
       =-t E[Σ:D²φ(W_t)]
        -t E[(ΣJ_t^T):D²φ(W_t)].                       (8)

The second term is the covariance-path sampler feedback; it cannot be discarded.

Let T_j be the all-rank centered Stein tensors for B. The imported recurrence gives

    E T2=κ3/2=3K3,    E T3=κ4/6=4K4,
    ||T4||_(L²;HS)<=L⁴e.

The Gaussian-bank Riesz identities, applied only to G while Z1,Z2 are held independent, then yield the exact finite identity

    d/dt Eφ(W_t)
      =3t² K3:E D³φ(W_t)+4t³ K4:E D⁴φ(W_t)
       +t⁴ E[T4:D⁵φ(W_t)]
       -3t² E[q3,t(x_t)·∇φ(W_t)]
       -4t³ E[q4,t(x_t)·∇φ(W_t)]
       +E[D_t·∇φ(W_t)]
       -t E[(ΣJ_t^T):D²φ(W_t)].                        (9)

In particular neither the test argument nor the source-bank law is Gaussianized along this calculation.

## 4. Hermite identities and all five current ranks

For clarity define the following contractions, preserving the first output slot of K:

    (K3[H,H])_ijk =sum_ab K3_iab H_ja H_kb,
    (K3[E])_ij =sum_ab K3_iab E_jab,

    (K4[H,H,H])_ijkl =sum_abc K4_iabc H_ja H_kb H_lc,
    (K4[E,H])_ijk =sum_abc K4_iabc E_jab H_kc,
    (K4[F])_ij =sum_abc K4_iabc F_jabc.

Symmetrizing their final physical output slots is optional, because every contraction below is against a symmetric derivative of φ.

Two Gaussian integrations in x_t give

    E[q3,t(x_t)·∇φ(W_t)]
       =E[K3[H_t,H_t]:D³φ(W_t)]
        +E[K3[E_t]:D²φ(W_t)].                          (10)

Three integrations give

    E[q4,t(x_t)·∇φ(W_t)]
       =E[K4[H_t,H_t,H_t]:D⁴φ(W_t)]
        +3E[K4[E_t,H_t]:D³φ(W_t)]
        +E[K4[F_t]:D²φ(W_t)].                          (11)

The factor three in (11) is the three partitions of three differentiations into a pair and a singleton. The third-derivative term is essential: q4 is cubic, so F_t is nonzero.

The leading constants in (10)-(11) cancel the first two terms of (9). Thus

    d/dt Eφ(W_t)=sum_(r=1)^5 E[A_r,t:D^rφ(W_t)],        (12)

with the following complete residual coefficients:

    A1,t = D_t,

    A2,t = -t ΣJ_t^T
           -3t² K3[E_t]
           -4t³ K4[F_t],

    A3,t = -3t²(K3[H_t,H_t]-K3)
           -12t³ K4[E_t,H_t],

    A4,t = -4t³(K4[H_t,H_t,H_t]-K4),

    A5,t = t⁴ T4.                                      (13)

Every coefficient is independent of Z2. Every test is evaluated at the same actual W_t. Equations (12)-(13) retain the moving-root feedback, the quadratic-map Hessian feedback, the cubic-map second-derivative feedback, and the cubic-map third-derivative feedback.

## 5. Dimension-free one-energy bounds for the feedback

Set s=k3+k4. Uniformly in t, the following estimates hold with numerical constants:

    ||D_t||₂+||ΣJ_t^T||_(L²;HS)
         <=C L²(h3+h4),

    ||K3[E_t]||_(L²;HS) <=C h3 s,
    ||K4[F_t]||HS <=C h4 k4,

    ||K3[H_t,H_t]-K3||_(L²;HS)
         <=C h3(s+s²),

    ||K4[E_t,H_t]||_(L²;HS)
         <=C h4(s+s²),

    ||K4[H_t,H_t,H_t]-K4||_(L²;HS)
         <=C h4(s+s²+s³).                               (14)

The first line follows from Hermite isometry. Indeed

    H_r^(C_t)(S_tz)=(S_t^(-1))^⊗r H_r^I(z),

and ||∂_t S_t^(-1)||op<=C L². Differentiating these explicit coefficient transforms bounds dot b3 and dot b4 by C L²h3 and C L²h4. First spatial derivatives of the Hermite polynomials similarly give ||J_t||_(L²;HS)<=C(h3+h4), and ||Σ||op<=L².

Here is a complete tensor-contraction justification of the remaining lines, avoiding dimension-dependent matrix-operator moments.

Each first, second or third spatial derivative of q3 or q4 is one coefficient tensor K3 or K4, transformed on its differentiated and Gaussian slots by bounded matrices, contracted against a centered Gaussian Hermite polynomial of degree at most two. Derivative multiplicities are bounded numerical constants. In particular:

- J_t has one free physical output and one chain-rule input, and Gaussian degree one or two.
- E_t has one free physical output and two chain-rule inputs, and Gaussian degree zero or one.
- F_t has one free physical output and three chain-rule inputs, and Gaussian degree zero; only K4 contributes.

Expand H=I+J in each displayed feedback. Each monomial has one anchor tensor K3 or K4 and at most three additional coefficient tensors. The deterministic chain-rule contractions connect every additional tensor directly to the anchor. Each tensor factor retains at least one final physical output slot. None of those output slots is a Gaussian slot or is summed away.

Apply the finite Hermite product formula. It expands each product into a finite sum over cross-factor Gaussian contractions and surviving Hermite tensors. There are no contractions within a single factor, because each individual factor is already Wick ordered. Every resulting coefficient network remains connected through its chain-rule edges. Every vertex still retains its physical output slot. The number of diagrams, their multiplicities and the remaining Hermite degrees are bounded numerical constants; the largest original Gaussian degree here is six.

For each resulting deterministic coefficient network, start with the anchor tensor in HS norm. Add the other vertices in any order connected to the vertices already added. Contract all its edges to existing vertices at that step. This uses a nonempty subset of its slots. At least its physical output slot remains uncontracted, so the selected flattening is a proper cut. The matrix inequality

    ||XY||HS <= ||X||HS ||Y||op

therefore bounds this addition by the appropriate k3 or k4. Future-edge slots and all surviving Gaussian/physical slots are retained in the current HS tensor. Bounded covariance transforms multiply the estimate by numerical constants. This argument also covers extra cross-factor Wick edges: when the latter endpoint is added, those edges are simply among the slots contracted at that step. It introduces no self-trace and no second Hilbert energy.

The final surviving Hermite contraction has L² norm at most the square root of its degree factorial times this kernel's HS norm. Summing the finite diagrams proves a bound of

    C × (anchor HS norm) × product(other proper-cut norms).

Applying this rule gives exactly the last five lines of (14). In particular it gives C h3k3 for b3², C(h3k4+h4k3) for b3b4, and C h4k4 for b4². Terms containing two or three Jacobian factors contribute the explicitly displayed higher powers of s.

The “every vertex retains a physical output” condition is important: this proof is for the displayed same-endpoint feedback graphs, not an assertion about arbitrary tensor networks with fully closed vertices or self-traces.

## 6. From exact currents to Wasserstein length

A rank-r current in (12) can transfer r-1 derivatives to the untouched ηZ2 buffer. If A_r,t has not been fully symmetrized, symmetrization can only reduce its HS norm. Hermite isometry and conditional Jensen yield an admissible continuity-equation velocity of L² norm at most

    √((r-1)!) η^(-(r-1)) ||A_r,t||_(L²;HS).

All coefficients are independent of Z2, as required. Summing the ranks and integrating t from zero to one gives

    W2(Law(W1),Law(W0))
      <=C[L⁴e+L²(h3+h4)
              +h3(s+s²)+h4(s+s²+s³)].                 (15)

This is a positive-law dynamic bound; it does not infer a sampler from a truncated cumulant series or differentiate a Wasserstein estimate.

Using (4)-(5), L<=c implies s<=c' and gives

    W2(Law(W1),Law(W0))
       <=C[L⁴e+L⁵e+L⁶e+L⁷e]
       <=C L⁴e.                                       (16)

The leading nonlinear-feedback contributions before absorption are

    h3k3=O(L⁵e),
    h3k4+h4k3=O(L⁶e),
    h4k4=O(L⁷e).

The L²h3 covariance/time-path contribution is O(L⁴e), so it remains part of the dominant bill. L²h4 is O(L⁵e).

The reserve Gaussian gives smooth positive endpoint densities. Polynomial Gaussian moments and the finite-bank C1/Sobolev Stein identities justify the differentiations for smooth tests; standard truncation and Sobolev approximation extend the result to the stated C1/Lipschitz setting. No derivative of f beyond Df is used. Uniformly Lipschitz strong finite-bank limits pass by the same moment/conditional-approximation argument as in the fourth-cumulant comparison note.

## 7. Exact buffer scaling

Apply (15) to B/√v. The normalized parameters are

    L'=L/√v,    e'=e/√v,
    h3'=h3/v^(3/2),    k3'=k3/v^(3/2),
    h4'=h4/v²,        k4'=k4/v².

Multiplying the output back by √v yields (2). The leading individual bills are

    L⁴e/v²,
    L²h3/v²,            L²h4/v^(5/2),
    h3k3/v^(5/2),
    (h3k4+h4k3)/v³,
    h4k4/v^(7/2).                                      (17)

The higher products in (15) are paid at their actual scaled powers, and are absorbed only under the explicit normalized small-radius guard L/√v<=c. With (4)-(5), (17) reduces to (3).

In particular this is a buffer-explicit bound, not a fixed-buffer constant silently used as v tends to zero.

## 8. What this result does and does not establish

This proves the analytical positive fourth-cumulant consumer for the actual centered source and its actual third and fourth connected cumulants. It pays every same-endpoint feedback of the quadratic and cubic polynomial reference maps with one physical Hilbert energy. The reference uses only a visible Gaussian carrier and an independent positive reserve; it requires no nonpositive signed law.

A native finite original-VALUE construction must still realize the proper cubic reference coefficient and its complete lower-order feedback, own and consume its Gaussian roots correctly, reconcile independent packet references, control first/caller paths, and pay all numerical and finite-clock errors. Formula (1) is not permission to execute K3 or K4 as arbitrary tensor-source inputs.

## Source context

- Existing rank-three positive consumer: /workspace/shared/recovery-20261004/cost/endpoint-stein-quadrature-20261004/positive-law-repair/higher-cumulant-gate/rank3-continuation/POSITIVE-QUADRATIC-REFERENCE-FOR-SKEW-GAUSSIANIZATION.md.
- Centered Stein recurrence: /workspace/shared/recovery-20261004/cost/endpoint-stein-quadrature-20261004/positive-law-repair/higher-cumulant-gate/order-reentry/ALL-RANK-CENTERED-STEIN-CURRENT-RECURRENCE.md.
- Fourth-cumulant comparison and mixed Wick-cubic covariance proof: /workspace/shared/fourth-cumulant-return-20261005/comparison/FOURTH-CONDITIONAL-CUMULANT-COMPARISON.md.

No imported source file is modified.
