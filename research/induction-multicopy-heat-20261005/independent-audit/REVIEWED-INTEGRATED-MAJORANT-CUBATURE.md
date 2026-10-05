# Sharper positive cubature from the integrated n-cubic majorant

2026-10-05. A finite-cost addendum. The four main proof files remain unchanged. This improves the conservative original-width global-Lipschitz bill for n old cubic copies; it still has polynomial inverse-tolerance/cutoff cost and does not establish public-log cubature.

## 1. Pointwise positive derivative majorant

Fix n>=2 and all input finite cubic sums. Factor out the common alpha^(3n): F and epsilon below are normalized coefficient quantities, so the actual coefficient error has the additional alpha^(3n). Consider the full tensor-valued BKAR integrand F(t), summing old clocks and force-hit histories with their exact weights but not the new replica-tree parameters. The following can be applied one replica tree at a time or with a triangle sum over all trees.

Gaussian Price differentiation adds an edge between two DISTINCT replicas and one force hit at each end. It preserves the existing marked spanning tree, so the all-cut/one-Hilbert bound applies directly to the augmented graph. Do not first contract the two new open bank indices as a self-trace.

The local clock proof with two extra hits supplies a known nonnegative majorant

    B(t)=sum_lambda C_lambda eta^(-kappa_lambda)
                     product_i (delta_i(t)+eta)^(-a_i,lambda),
    sum_i a_i,lambda=n/2,                               (1.1)

for the sum of ALL possible Price pair/force-hit derivative terms, after the old-clock sums. Every C_lambda is a computable scalar from the finite rank census, actual input weights, row bounds and C_k constants. Neither F nor a tensor entry is evaluated to construct B. Physical/native readout constants, if included in F, remain included in C_lambda.

The all-n q=2 theorem implies

    integral_[0,1]^(n-1) B(t) dt
       <= C(n,2) L^(6n-1) eta^(-(n-2)).                 (1.2)

Enlarge the known C(n,2) by the exact finite Price pair multiplicity if its chosen q-hit convention did not already count it. This is a computable rank factor, not a width or tolerance power.

## 2. Comparability inside dyadic cells

Partition every edge gap s=1-t into dyadic panels down to eta and one final panel [0,eta]. In any resulting product cell, every s_e+eta varies by at most a factor two. Since

    delta_i+eta=min_(e incident i)(s_e+eta),

the same is true of each delta_i+eta. Each term in (1.1) therefore varies by at most 2^(n/2), and so does the whole positive sum:

    sup_cell B <=2^(n/2) inf_cell B.                     (2.1)

This includes the final [0,eta] panels. It does not require the raw gaps themselves to be bounded below by eta. Any subdivision of a dyadic product cell retains the same comparison.

## 3. Positive midpoint error and node count

Subdivide the panels into intervals of length at most h and tensor their positive midpoint rules. Every node is interior; node weights obey u<=2s, all weights are positive, and panel masses remain exact.

For a point t in a small product cell and its midpoint t0, follow the line segment from t0 to t. Every min-path covariance entry changes with speed at most ||t-t0||_infinity<=h/2, at each differentiability point of that path. Price differentiation and (1.1) then give

    ||F(t)-F(t0)||HS <=(h/2) sup_cell B sqrt(D).

A generic perturbation handles a segment in a min-path ordering wall, using continuity. Integrate over the cell, use (2.1), and sum the cells. The numerical coefficient error is at most

    (h/2) 2^(n/2) C(n,2) L^(6n-1) eta^(-(n-2)) sqrt(D).  (3.1)

Thus an explicit sufficient choice for error epsilon sqrt(D) is

    h=min(1, 2 epsilon /
             max(1,2^(n/2) C(n,2)L^(6n-1)eta^(-(n-2)))). (3.2)

One edge then has at most 2+ceil(1/h)+ceil(log2(1/eta)) nodes; the tensor rule has that count to the power n-1. In big-O notation, with all rank constants charged,

    Q_T <= C_n [1+epsilon^(-1)L^(6n-1)eta^(-(n-2))
                   +log(1/eta)]^(n-1).                 (3.3)

This is considerably smaller than using the maximum original-width Price derivative everywhere. It is still polynomial in inverse tolerance and original cutoff. All generated source instances, actual readout shares, native/filter floors and replays must be billed at this ACTUAL Q_T. No claim of source admission at every intended grade follows merely from finiteness.

There is no tensor oracle, Monte Carlo coefficient estimator, target reheating, or hidden alteration of repeated old-clock labels in this rule. It is a positive deterministic scalar quadrature chosen from a proved envelope.
