# Separate the outer cutoff from the tail delay

2026-10-05. Generalization of the reusable ledger, prompted by the independent short-prefix construction. This uses only the already pinned variable-buffer services; it claims no improved prefix or skew producer.

The main note's 10/3 prospective prefix-only ceiling refers ONLY to its old tied cutoff eta=w. An outer endpoint cutoff is a separate public parameter and need not equal the force-tail delay. The unmodified 16/5 optimum remains unchanged even after making that distinction.

## 1. Independent public scales

Keep q=1-w and sigma^2=1-q^2 in the delayed tail and its Gaussian conditional disintegration, but run whole covariance/mean services for

    0<t<=1-eta,

and the fresh original F_Q RAW branch for 1-eta<t<1. Assume

    A<=eta<=w<=1/2.

At bulk nodes the actual conditional variance is

    v(t)=(1-t^2)(1-q^2)/(1-q^2 t^2),
    1/v(t)=q^2/(1-q^2)+1/(1-t^2).                       (1)

Thus c eta<=v(t)<=2w. Every allocated service variance is at least c eta and hence at least c A. The exact variance allocations, independence/readset constraints, full private-bank integrations and common-carrier rotations are unchanged. The actual variable-buffer residual first

    Lambda[A+A^(3/2)/sqrt(u)]

therefore stays O(Lambda A). The mean and self-reserve normalized radii are at worst Lambda A^(1/2); the mixed-K radius is at worst Lambda A^(2/3). Check the literal constants in the imported guards as usual. All positive gaps remain actual fractions of u.

## 2. Integrated debts

For 0<eta<=w<=1/2, elementary powers of (1) give

    int_bulk v^-1 dt <= C[w^-1+log(1/eta)],
    int_bulk v^(-3/2) dt <= C[w^(-3/2)+eta^(-1/2)],
    int_bulk v^(-7/6) dt <= C[w^(-7/6)+eta^(-1/6)],
    int_bulk v^(-1/2) dt <= C w^(-1/2),
    int_near v^(-1/2) dt <= C[eta/sqrt(w)+sqrt(eta)].     (2)

For example (a+b)^p<=2^(p-1)(a^p+b^p) when p>=1. The only endpoint integrals needed are integral_eta^1 x^-1 dx, integral_eta^1 x^(-3/2) dx, and integral_eta^1 x^(-7/6) dx. On the near branch use sqrt(a+b)<=sqrt(a)+sqrt(b). Since eta<=w, eta/sqrt(w)<=sqrt(eta), but both terms are retained in the exact ledger.

The finite positive dyadic rule obeys the same inequalities up to a numerical constant: insert t=1-eta as a panel boundary; on a panel with endpoint gap h, total weight is O(h) and all strictly interior nodes have gaps comparable to h. Its contributions to these sums are O(1), O(h^-1/2), and O(h^-1/6), respectively. Summing panels down to eta gives the first three bulk bounds. The near square-root sums run toward the early positive midpoint and are bounded by C sqrt(eta). No service is executed at t=1 or below the minimum bulk variance. The already prescribed finite outer tolerance and early midpoint can be chosen much smaller than eta with public-log node cost.

## 3. Generalized finite-depth recurrence

The exact same finite-depth input/output source type and covariance restoration proofs from the main note now give

    e_(k+1) <= A e_k
      +C sqrt(D)[A^2 w^(3/2)+A^3 w
          +A^3(eta/sqrt(w)+sqrt(eta))
          +A^4(w^-1+log(1/eta))+A^4 eta+A^4]
      +Lambda_(k+1) sqrt(D)[
          A^5(w^(-3/2)+eta^(-1/2))
          +A^6(w^(-7/6)+eta^(-1/6))]
      +absolute floors.                                  (3)

Whole covariance target restoration contributes A^5 times the fourth integral in (2), dominated by the displayed mean/Gram term. The higher covariance-mixture term A^7 times the v^-3/2 integral is also dominated. The completed mean input is still used only in the bulk; the old F_Q supplies the fresh near branch with its numerical O(A) first.

Choose eta=A. The near branch is grade7/2, while the original prefix A^2 w^(3/2) and Gaussianization A^4/w still balance at w=A^(4/5), grade16/5. Thus the main unchanged-mechanism conclusion is robust.

If a new prefix producer improves the first term, the earlier 10/3 cap from eta=w no longer applies. With the residual A^3w still present, deleting ONLY A^2w^(3/2) from (3) would permit grade7/2 at w=A^(1/2), where A^3w and A^4/w balance; the near branch also has that grade. This is an upper envelope of a prospective ledger, not a proved prefix theorem. The short-prefix worker's actual new errors must be inserted rather than deleted.

Likewise a new third-cumulant service must supply its own variable-buffer law error and fourth-order weak remainder. Equation (3) does not certify those producers. It merely separates two parameters that were unnecessarily tied in the earlier scalar bookkeeping.

## 4. Costs

Replace the old bulk and near sets in the main note's literal count and root sums by the new sets. Reentries, origins and all original VALUES remain fully charged. There are still public-log many positive outer nodes; the smallest executed u is comparable to eta rather than w. Precision allocations use these actual smaller gaps and native normalizations. No inverse-eta replication count or hidden carrier subtraction is introduced. At any fixed depth, counts and actual dimensions/D remain public-log polynomials under the imported native guards.

## 5. Weighted source closure below eta=A

The simplifying eta>=A restriction in Sections 1-4 is sufficient, not necessary. For this graph one can instead take

    A<=w<=1/2, eta=A^beta<=w, 1<beta<3/2,

with beta fixed, and impose the literal native guards at every actual u>=c eta. The singular-sum estimates (2), recurrence (3), and full costs are unchanged.

The key distinction is between individual service residual first and the first of the final positive source. At one bulk node the mean-service contribution has actual carrier-subtracted private/caller first bounded by

    R_j <= Lambda[A+A^(3/2)/sqrt(v_j)].

The Gram and mixed-K contributions satisfy the same safe majorant. No downstream native consumer requires this INDIVIDUAL physical residual first to be O(A); the complete return is fed directly into the original F, and only the final positive sum is passed to the next own-mean compiler. With the literal common-carrier rotations, all Gaussian row operators have norm one. Hence weighted triangle and Minkowski give

    sum_bulk omega_j R_j
      <= Lambda[A+A^(3/2) sum_bulk omega_j v_j^(-1/2)]
      <= Lambda[A+A^(3/2)/sqrt(w)] <= C Lambda A.        (4)

The near branch still has numerical O(A) first. At each terminal node the potentially O(A) Hessian difference acts along a scalar multiple of the common carrier and is symmetric. Every nonsymmetric/off-carrier terminal path has an additional original A first times R_j (or the original wA^2 term). Summing those bounds with (4) proves aggregate square-lift curl O(Lambda A^2). Source-zero vanishing on the identical anchor record and the complete root dimensions give aggregate remainder energy O(Lambda A^2)(|z|+sqrt(D)). The total source first and retained-caller first also remain O(Lambda A). Thus the SAME actual reusable input/output class closes, without imposing O(A) on each unweighted boundary node.

Native admission must still be checked pointwise. At the smallest u comparable to eta=A^beta, the substantive powers are

    mean radius A/sqrt(u): A^(1-beta/2),
    raw self-reserve first A^(3/2)/u: A^(3/2-beta),
    mixed-K radius A/u^(1/3): A^(1-beta/3).

All are small for fixed beta<3/2, after substituting every actual public-log factor and the imported numerical thresholds. Also u>=A^3 for sufficiently small A. The concrete beta=7/5 gives powers 3/10, 1/10, and 8/15, respectively, and a near-RAW exponent 3+beta/2=37/10.

Do not include beta=3/2 in this strictly-power theorem: the self-reserve power would be zero. A boundary theorem with an explicit logarithmic multiplier in eta would need a separately checked quantitative native threshold. This note deliberately keeps strict power slack.

This extension alone still does not improve the original prefix/Gaussianization balance beyond16/5. Its purpose is to supply a valid aggregate source port and error ledger when independently proved prefix/cumulant improvements remove those older forcing terms.
