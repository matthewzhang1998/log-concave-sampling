# Paired Gaussian rotation: all-cut covariance and cubic clock sum

2026-10-04. Exact coefficient theorem for the finite native node. This proves the delta-squared coarse-covariance estimate needed by the intermediate-scale B mark. It does not realize a negative tensor current as a new positive packet and does not supply an all-current recurrence.

## 1. An exact Price identity on the entire shared bank

Let F:R^m->(R^D)^(tensor 3) be continuously differentiable, with the bounds below, and fix an arbitrary retained affine center x. Let Q,P be independent standard m-dimensional Gaussians. For 0<=delta<=1/2 and s>=0 put

    q_+ = x+s[sqrt(1-delta^2)Q+delta P],
    q_- = x+s[sqrt(1-delta^2)Q-delta P],
    E_delta=F(q_+)-F(q_-),
    tau=1-2delta^2.

Both marginals are N(x,s^2 I), and their cross covariance is s^2 tau I. Thus E E_delta=0. Gaussian Price differentiation gives exactly

    E[E_delta tensor E_delta]
      =2s^2 integral_tau^1
          E[sum_(a=1)^m partial_a F(x+sG)
                    tensor partial_a F(x+s[tG+sqrt(1-t^2)H])] dt,       (1)

where G,H are independent standard Gaussians.

Indeed the left side is twice the difference between the diagonal tensor moment and the tensor moment at correlation tau; their common mean terms cancel. Differentiating the latter correlated Gaussian expectation in t differentiates each factor once and contracts their common input index, with coefficient s^2. Integrating gives (1). One may first smooth and cut off; the stated bounded derivative cuts and finite Gaussian moments justify passage to the limit. There is no second derivative of F in the right side.

Equation (1) retains the FULL m-dimensional bank and its existing query ancestry. It neither splits the bridge heat among original vertices nor treats q_+ and q_- as independent.

## 2. Tensor contraction lemma: every proper cut is controlled

Suppose every proper matrix flattening of the rank-four tensor D F(q), including the root-versus-all-three-physical cut, has operator norm <=K1 uniformly in q.

If A and B are two such derivative tensors, contract their root indices:

    T_(ijk,lmn)=sum_a A_(ijk,a) B_(lmn,a).

Then every proper cut of the resulting six-physical-slot tensor T is <=K1^2.

To prove this for an arbitrary external bipartition, view A and B as matrix families A_a and B_a, with rows the external slots on the left and columns those on the right. The global cut is sum_a A_a tensor B_a. Operator Cauchy-Schwarz bounds it by

    ||sum_a A_a A_a*||op^(1/2)
       ||sum_a B_a* B_a||op^(1/2),

or the reversed orientation. These two factors are squares of local flattening norms, with root index a assigned to the opposite side. Whenever the global cut is nontrivial, one orientation avoids an empty local cut in both tensors; the two factors are therefore bounded by K1. This includes cuts which put an entire old subtree on one side.

In the three-versus-three flattening, ordinary matrix multiplication also gives

    ||T||HS <= min(K1 ||A||HS, K1 ||B||HS).

The root-versus-three cut alone implies ||D F(q)||HS<=sqrt(m)K1. Consequently ||T||HS<=sqrt(m)K1^2. More refined source-qualified Hilbert energy can be substituted if separately available; it is not inferred merely from proper cuts.

## 3. Exact delta-squared estimates

The integration interval in (1) has length 2delta^2. Thus its covariance tensor V_delta obeys

    every proper cut(V_delta) <= 4s^2 delta^2 K1^2,
    ||V_delta||HS <= 4s^2 delta^2 sqrt(m) K1^2.                      (2)

These estimates hold at each retained x, not just after an additional average in x. They require precisely the all-proper-cut derivative input above. Generic proper cuts of F alone do not imply that input.

The same proof permits pointwise varying Hilbert energies: one marked factor in (1) carries its actual HS norm, while the other is controlled by K1. Gaussian marginality of each endpoint allows a single L2 Hilbert energy to be integrated. For the application below the explicit sqrt(5D) bound is enough to retain one dimension-sized Hilbert factor.

## 4. Substitute the audited native cubic-node cuts

For one exact five-clock cubic node, use the normalized analytic tensor F=A_Q from the native adapter. Its entire original coarse bank has m=5D, and the independently audited source-qualified derivative bounds are

    proper cuts(D_Q A_Q) <= C(sigma_1^(-1)+sigma_2^(-1)+sigma_3^(-1)).

The same sigma_i are the original independent shields. No sigma_i is enlarged. The old cubic amplitude is

    epsilon_node = 4 w A^3/(c a b sigma_2),

up to the already declared source/readout conventions. After restoring all known readout/variance-share factors, (2) bounds the node's new pure rank-six coarse-mixture coefficient by

    Lambda A^6 delta^2 sqrt(D)
          w^2 sigma_2^(-2)
             (sigma_1^(-1)+sigma_2^(-1)+sigma_3^(-1))^2.             (3)

This is a finite-coefficient statement about the actual declared node's analytic current. Finite bounded-field, pair/filter, clipping and numerical errors retain separate absolute floors; the analytic tensor is never sampled as a covariance matrix.

## 5. Positive dyadic clock sum

For the original-force clock r_i, sigma_i^2=(1-r_i^2)/2. On an endpoint dyadic cell, its positive quadrature mass Delta_i is at most a constant times sigma_i^2; the final cutoff cell satisfies the same bound. The bounded polynomial r_i factors do not worsen it. The total positive masses away from the endpoint are bounded.

Use (sum_i sigma_i^(-1))^2<=3 sum_i sigma_i^(-2). For i=2 the worst clock factor is

    Delta_2^2/sigma_2^4 <= C,

so its sum costs only the recorded logarithmic number of endpoint panels. For i=1 or 3, the corresponding factors are

    Delta_i^2/sigma_i^2 <= C Delta_i,
    Delta_2^2/sigma_2^2 <= C Delta_2,

and are summable. The other three clock factors have bounded sums of squared positive masses. Therefore

    sum_nodes w_node^2 sigma_2^(-2)
         (sum_i sigma_i^(-1))^2 <= Lambda.                         (4)

All inverse positive-share factors and the finite clock count are already included in the public logarithm Lambda.

If completed nodes use independent entire coarse banks conditional on the retained physical caller, their signed differences remain independent under the whole-bank rotation: the standard marginal P has independent blocks even when it is correlated with a physical readout R. Cross-node covariances therefore vanish, and (4) applies to the completed sum. This requires the literal independent-node read set. A shared outer random caller or any other shared coefficient root must either remain retained or be included as a distinct shared block with its cross-node currents; it cannot be silently declared independent.

## 6. Fractional B mark consequence

In the exact B-pair reference of the intermediate-scale corollary, P is still standard Gaussian marginally. Its cross matrix with retained physical input R does not change (1)-(4). Set delta=A^p, 0<p<1. Then the new pure rank-six covariance coefficient is bounded by

    Lambda A^(6+2p) sqrt(D).

For p=1/2 this is Lambda A^7 sqrt(D), with every proper cut bounded by Lambda A^7. The conditional quartic six-force coefficient remains at grade six and requires its own pointwise-coarse correction or full-current handling. Public tilting of both terms generates additional ranks, which must stay in the same endpoint current record.

The proof is a grouped coefficient estimate and respects the original Gaussian Gram. It does not claim an independent shield sqrt(sigma_i^2+u), does not differentiate a completed output, and does not replace an expectation by a producer. The executed producer remains the two whole-kernel copies plus the correctly oriented B pair. Its actual first and full VALUE count are those of the companion construction.

## 7. History-specific strict spine for the cubic node

Expand D_Q A only after keeping the entire Price bridge intact. Retain the actual coarse-root injection and the derivative-hit vertex v on the first three-vertex path, and w on the second. Contraction of the root index joins those vertices by one bridge edge. A center has eccentricity one and a leaf eccentricity two, so the joined diameter is

    ecc(v)+ecc(w)+1 in {3,4,5}.

The corresponding one-B-leaf attachment at v or w has diameter max(2,ecc(v)+1), in {2,3}. The bridge diameter is strictly larger than both in every cubic hit pair. No hit is pooled away, and no independent old shield is enlarged. The complete original query copies and all coarse injection weights remain in the record.

This verifies strict spine for these two cubic paths only. It is not a claim that arbitrary different older trees can be pooled and ordered by their largest diameter.

## Scope of the gain

The previously open delta-squared all-cut estimate is proved for the literal native cubic-node derivative contract and the stated independent-node completion. The companion B-mark corollary explains uniform sqrt(D) recalibration of the imported finite provider floors under changed coarse labels. Unrelated clock-target errors under new callers, full public-tilt/current export, and the sign-programmable conditional quartic correction remain separate obligations. The result does not assert a general all-order closed grammar or eventual-sublinear query exponent by itself.
