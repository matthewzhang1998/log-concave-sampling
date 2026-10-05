# Rank-explicit exact n-cubic heat and clock bounds

2026-10-05. This strengthens the general computable rule with a closed-form safe cutoff exponent for every finite n. It uses the SAME exact m_star allocation for all derivatives, preserves the original target, and states rather than hides the higher-rank loss.

## 1. Statement

Take n>=2 full finite old cubic coefficient sums sharing their complete original Gaussian bank. Every original primitive and relevant structural endpoint variance is >=eta>0, eta<=1/2. Let L>=1 dominate one plus all original dyadic panel counts and log(1/eta). The old positive mass assumptions and zero-readout/retained-variable contract are those of THREE-CUBIC-EXACT-HEAT.md.

At each replica use the exact n-copy Gaussian geometry and execute m_star=max_j m_j in its extra independent private heat. Let q>=0 be an additional analytical coefficient-bank/caller derivative order. For q=1 this is the complete captured-caller coefficient bound needed for the actual first/path ledger; q>1 is only an analytical heated-tensor statement, not an executed higher derivative of g.

Define

    e(n,q)=max(n+q-4,0).

For the full finite tree/hit/history census,

    all proper cuts <= C(n,q) L^(6n-1) alpha^(3n) eta^(-e(n,q)),
    HS <= the same bound times sqrt(D).                  (1.1)

C(n,q) includes the explicit tree/force-hit count, actual input masses and rows, local finite C_k constants, physical permutations and declared normalization. No cutoff power is inside C(n,q). Before physical Wick contractions, the main has

    a=3n, N=3n, K=3n-2, R=3n, M=3n-1.

For fixed original input histories the tree/bridge-hit count is

    n^(n-2) 3^(2n-2),

and q ordered extra hits multiply this by at most (3n)^q. The largest main adapter order is n, reached when a replica of degree n-1 receives all hits at its original C1 center. The actual raw source count is the sum of 2^(k_v+1), not one tensor query per vertex.

For node scalar factors b_h including all positive clock/bridge masses and inverse widths, and the q-hit scalar envelope b_h^(q), the same proof gives

    sup_h b_h^(q) <= C(n,q) eta^(-e(n,q)),
    sum_h b_h^(q) <= C(n,q)L^(6n-1)eta^(-e(n,q)),
    sum_h [b_h^(q)]^2 <= C(n,q)L^(6n-1)eta^(-2e(n,q)).    (1.2)

Here b_h^(1) includes the sum of all actual row-weighted 1/t caller factors. Fixed inverse readout products and native source/frame constants must be included explicitly before using this in a law comparison. A chosen quadrature must preserve the stated positive panel envelope and pay its own error/cost.

## 2. One cubic with arbitrarily many hits

Let d be its total bank/caller hit count and let p=(p_1,p_c,p_3) be its actual inverse-width powers, so p_c starts at one and sum p_v=1+d. Define

    a_d=max((d-1)/2,0),  kappa_d=max((d-3)/2,0).

The universal executed heat dominates every m_j, so the analytical choice of a largest-power force does not change the source. Its positive old-clock sum and its individual node contribution obey

    J_d(delta) <= C_d L^5 eta^(-kappa_d)(delta+eta)^(-a_d),       (2.1)

with L^5 omitted for the node bound.

Proof. If all p_v<=2, original private widths and endpoint masses give only logarithms. Otherwise choose a force j with p_j>2 and put a=p_j/2-1. Sum its clock first, gaining

    C m_j^(-a)(delta+eta)^(-a).

Expand m_j^(-a) by the sum of its two other primitive inverse powers and its relevant structural inverse power. For any one resulting term, the TOTAL exponent on the remaining endpoint variables is

    sum_(v!=j) p_v/2+a=(d-1)/2.

Each endpoint mass pays exponent one. The total unpaid exponent is at most max((d-1)/2-1,0)=kappa_d: if any exponent exceeds one, subtract at least one paid mass; if none does, there is no unpaid exponent. Bound unpaid powers by the actual eta floor. Also a<=a_d. Replacing its delta power by a_d changes only a fixed C_d because delta+eta<=3/2. At most five original clock sums contribute their logarithms. The pointwise proof uses the same endpoint powers and node mass bounds. This proves (2.1).

## 3. Ordered edge-gap sectors pay the remaining singularity

Fix one replica tree T. Its n-1 edge gaps are s_e=1-t_e. Replica i has d_i=degree_T(i)+q_i hits, where sum q_i=q. Since every tree degree is at least one,

    sum_i a_(d_i)=(n+q-2)/2.                            (3.1)

Also

    sum_i kappa_(d_i) <= max(n+q-4,0)/2.                 (3.2)

Indeed, if S={i:d_i>3} is nonempty, twice the left side is sum_(i in S)(d_i-1)-2|S| <= sum_i(d_i-1)-2=n+q-4. If S is empty it is zero.

Order the edge gaps s_pi(1)<=...<=s_pi(m), m=n-1. On this sector the independent variance delta_i is the smallest incident edge gap. Assign a_i to that edge; let b_e be the sum of assignments. Then the bridge integrand after original-clock summation is bounded by

    eta^(-sum kappa_i) product_(j=1)^m (s_pi(j)+eta)^(-b_pi(j)),
    sum_e b_e=(n+q-2)/2.                                (3.3)

For any nonnegative edge exponents b_1,...,b_m, the ordered dyadic integral has cutoff power

    kappa_bridge=max(0, max_(1<=k<=m)[sum_(j<=k)b_j-k])    (3.4)

and at most L^m logarithmic factors. A direct proof uses dyadic levels ell_1>=...>=ell_m>=0, capped at J comparable to log2(1/eta). A cell's mass times its integrand is proportional to 2^(sum_j (b_j-1)ell_j). Write r_k=ell_k-ell_(k+1)>=0. Its exponent is sum_k r_k B_k, where B_k=sum_(j<=k)(b_j-1) and sum r_k=ell_1<=J. This is at most J max(0,B_1,...,B_m); the number of integer level vectors is at most (J+1)^m. The final [0,eta] cell has the same bound using its mass eta and regularized integrand. This proves (3.4).

In our case every nonempty prefix has k>=1, so

    kappa_bridge <=max(sum_e b_e-1,0)
                  =max(n+q-4,0)/2.                      (3.5)

Adding (3.2) and (3.5) yields exactly e(n,q). The n old clock blocks contribute L^(5n), and the bridge cells contribute L^(n-1). Summing the finite order sectors and tree/hit census gives (1.1).

For nodewise bounds include each actual bridge weight u_e<=C s_e. The same dyadic exponent calculation with clipped gaps max(s_e,eta) proves the same eta power without counting cells. Max times sum gives the squared literal-weight bound in (1.2). It is not a variance assertion for independently replaced old clocks.

## 4. What the formula does and does not close

The main/first cutoff exponents are:

    n=2: 0,0;
    n=3: 0,0;
    n=4: 0,1;
    n=5: 1,2;
    general n: max(n-4,0), max(n-3,0).

All omitted borderline factors are the declared L powers. Thus three-copy main and first have no algebraic inverse-cutoff penalty; the fourth-copy main still has only logarithms. Higher main/caller terms have an explicit safe power which the actual amplitude/root and original truncation ledger must pay.

These powers are upper bounds. Sharper covariance allocation, cancellation or positive-return organization may improve them. No lower bound or impossibility claim follows. No original heat is replaced, and no new heat-mismatch error is introduced. Original finite-clock approximation debt and conservative new cubature cost remain real separate obligations.

At fixed n a native source with all the actual finite constants can be built from original VALUE calls by the marked-spanning-tree compiler. The width-zero conditional port and exact Riesz endpoint join can consume (1.1)-(1.2) only under their actual ownership, old-first, positive-reserve, readout and numerical prior guards. In particular a marginal pair LAW certificate does not authorize exposing its internal carrier as a retained variable. This theorem does not assume any new common-carrier amalgamation.
