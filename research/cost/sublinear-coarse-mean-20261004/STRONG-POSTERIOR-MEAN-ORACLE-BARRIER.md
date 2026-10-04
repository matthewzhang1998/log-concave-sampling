# A strong posterior-mean barrier for original-gradient control variates

2026-10-04. A uniform-class, target-specific information bound. This is not a lower bound for posterior sampling, Gaussian-buffer law comparison, or arbitrary outer algorithms.

## Result first

Let the target be the scalar canonical posterior force

    mu_A(V) = sqrt(A) E[ V'(X) ],
    density(X) proportional to exp[-x^2/(2A)-V(x)].

There is a universal c>0 such that every adaptive randomized program using at most N original gradient VALUE / scalar first-action query sites has, on some C2 potential with

    (1/4) <= V'' <= (3/4),   V'(0)=0,

strong mean-square error at least

    sup_V E|M-mu_A(V)|^2 >= c A^2/(N+1)^3,  0<A<=1.

Exact oracle responses, arbitrary other computation, private Gaussian generation, arbitrary adaptivity, and knowledge of the exact mode 0 are all granted. Thus the result also applies to a finite numerical implementation. It includes an arbitrary control variate, smoothing operation, or cubature whose complete original queries are included in N. It does not cover an additional free expectation, potential-value, or posterior-sample oracle.

Consequently a uniform strong target O(A^j), including a finite-source target whose bias from mu_A is O(A^(j+1)), requires

    N(A) >= c_j A^[-2(j-1)/3]

up to constants. For an admitted fixed-order cost bound N(A)<=Lambda_j(log(1/A)) A^(-c_j), necessarily

    c_j >= 2(j-1)/3 = (P-1)/3,   P=2j-1.

This rules out eventual c(P)/P -> 0 while retaining the existing uniform strong force-mean contract. It still leaves a possible smaller fixed slope; it does not assert that the displayed slope is attainable. In particular it cannot rule out the law-only Gaussian-buffer repair being studied separately.

## 1. Why the quadratic test is insufficient

For V(x)=lambda x^2/2 at caller 0, the exact posterior is N(0,A/(1+lambda A)) and mu_A=0. A centered affine or antithetic Gaussian control variate removes the entire mean fluctuation at constant original-query cost. At a general known caller y its exact canonical mean is sqrt(A) lambda y/(1+lambda A), also analytically tractable once lambda is identified. Hence this fixture cannot establish a barrier for an improved control variate.

The following genuine C2 family keeps the mode, gradient and Hessian at the original caller fixed, but hides independent nonlinear information away from it. The normalization and posterior dependence on V are included exactly.

## 2. Explicit C2 family

Define

    phi(t) = t^2(1-t)^2 on [0,1], and 0 otherwise,
    Phi(t) = integral from -infinity to t of phi(s) ds,
    B = integral phi = 1/30.

phi is nonnegative and C1 on the real line, Phi is C2, and |phi'|<=1. For integer n>=2 put w=1/n, a_i=1+(i-1)/n, and signs theta_i in {-1,1}. Fix lambda=1/2 and eta=1/4. Set

    V_theta(x) = lambda x^2/2
       + eta A w^2 sum_(i=1)^n theta_i Phi((x/sqrt(A)-a_i)/w).

Its original gradient and Hessian are

    V_theta'(x) = lambda x
       + eta sqrt(A) w sum_i theta_i phi((x/sqrt(A)-a_i)/w),
    V_theta''(x) = lambda
       + eta sum_i theta_i phi'((x/sqrt(A)-a_i)/w).

The open supports are disjoint. Thus 1/4<=V_theta''<=3/4 everywhere, including the cell boundaries by continuity. The potential is genuinely C2. All perturbations vanish with all displayed derivatives at x=0; its exact mode for caller 0 is 0. The same construction may use a nonnegative C-infinity bump instead, with eta rescaled by its finite first-derivative bound.

A query of V' or V'' at any x reveals at most the one sign whose open support contains x. Outside the supports it reveals none. Granting that this sign is revealed perfectly, even at a point where the numerical coefficient vanishes, only strengthens the algorithm. Several VALUE/first actions at a common site still reveal at most one sign; counting them as separate queries only weakens the lower bound. The premise is the stipulated original gradient/first-action oracle, not an oracle for V itself.

## 3. Exact posterior sensitivity, including normalization

Write U=X/sqrt(A). Its density is

    p_theta(u) = Z_theta^(-1) exp[-(1+lambda A)u^2/2-R_theta(u)],
    R_theta(u) = eta A w^2 sum_i theta_i Phi((u-a_i)/w).

For signs or for their continuous interpolation in [-1,1],

    |R_theta(u)| <= eta A B/n <= eta B/2.

Consequently p_theta/q_A lies between exp(-2 eta B/n) and exp(2 eta B/n), where q_A=N(0,(1+lambda A)^(-1)); using A<=1 only enlarges those bounds. This accounts for the actual normalizing constant.

Gaussian tails justify integration by parts and give the exact identity

    mu_A(theta) = -E_theta U.

Differentiating the normalized finite-dimensional exponential family gives

    partial_(theta_i) mu_A(theta)
      = eta A w^2 Cov_theta(U, Phi_i(U)),
    Phi_i(u)=Phi((u-a_i)/w).

This is an analysis of the target, not an executable derivative oracle. The covariance includes both the changing numerator and changing normalizer.

There is a universal c0>0 such that

    Cov_theta(U,Phi_i(U)) >= c0

for every A<=1, n>=2, i, and interpolated theta. Indeed, for independent U,U', monotonicity yields

    Cov(U,Phi_i(U))
      = (1/2) E[(U-U')(Phi_i(U)-Phi_i(U'))] >= 0.

On U in [2,5/2] and U' in [-1,0], Phi_i(U)=B, Phi_i(U')=0 and U-U'>=2. Include the reversed event to cancel the factor 1/2. The density-ratio bounds and q_A's variance in [2/3,1] give a uniform positive lower bound for the product of the two interval probabilities. Therefore

    partial_(theta_i) mu_A(theta) >= k A/n^2

with one universal k=eta c0>0. The sign of this derivative is the same for every other sign configuration. In particular, flipping one coordinate from -1 to +1 changes mu_A by at least 2k A/n^2. No approximate Gaussian-mean substitution has been made.

## 4. Adaptive randomized query lower bound

Put the independent uniform Rademacher prior on all n signs, and condition also on the algorithm's private randomness. Exposing each queried sign is a stronger oracle than the actual VALUE/first responses. After any adaptive transcript, every unexposed sign remains independent uniform Rademacher: query selection and stopping only use previously exposed signs and private randomness.

Condition on a transcript exposing m signs. Let I be the unexposed indices and f(theta_I)=mu_A(theta). For each i in I, its first Walsh coefficient is

    b_i = E[f(theta_I) theta_i | transcript]
        = (1/2) E_(theta_(I\{i})) [f(theta_i=+1)-f(theta_i=-1)]
        >= k A/n^2.

The coordinate Rademachers are orthonormal in conditional L2. Projection/Bessel therefore gives

    Var(mu_A | transcript) >= sum_(i in I) b_i^2
                                >= (n-m) k^2 A^2/n^4.

Any transcript-measurable output has conditional squared error at least this conditional variance, even when granted the optimal conditional mean. With a deterministic cap N, take n=2 max(1,N); then m<=N and

    E|M-mu_A|^2 >= k^2 A^2/(16 max(1,N)^3).

An average over the finite sign prior lower-bounds the worst-case risk, proving the claim. The same argument covers a program whose worst-case expected number of query sites is at most N: E m<=N suffices after taking expectations, with n=2 max(1,ceil(N)). Queries after all signs are exposed can only increase work.

The result allows the exact posterior mode as side information because that mode is identically zero throughout this family. The original caller, original gradient there and caller Hessian are also identical. A numerical caller profile cannot hide the risk: at y=0 the original force is zero, so no unbounded-caller issue arises.

## 5. Application to the complete direct-mean family

The admitted outer proof requires its actual executable force statistic M_j to satisfy

    ||M_j - sqrt(A) m_A(y)||_2 <= C_j sqrt(D) A^j

up to a separately selected numerical allowance. At D=1,y=0 the target here is exactly mu_A. If instead the statistic first targets the named finite source with strong error C_j A^j, the finite-source law/force bias C_j A^(j+1) only changes the constant by triangle inequality. Thus the lower bound applies to the complete mean program, not merely to its coarse bank. Its count N must include all source descendants, mode iterations, changed-query ancestors and every original first-action site.

This caller is also present in the baseline's actual selected outer chronology: at the Chebyshev-Lobatto endpoint c_i=0, the integral row A_i,l is zero and y_i^[d]=x for every depth and every captured refresh. Taking x=0 realizes the hard known-center target without declaring any other hidden Picard caller deterministic. The proof does not require conditioning on a probability-zero refreshed momentum.

Writing N<=C_j L^(d_j) A^(-c_j), the lower bound implies

    k C_j^(-3/2) L^(-3d_j/2) A^(1+3c_j/2)
        <= C'_j A^j.

For fixed j and A tending to zero this is impossible if 1+3c_j/2<j. Constants, logarithmic degrees, and the positive small-heat threshold may depend arbitrarily on j, but must be uniform over the chosen Hessian-sandwich class. That is the uniform class used by the baseline theorem. A theorem with uncontrolled potential-specific constants or thresholds, or a supplied quantitative Hessian modulus, is a different claim.

This invalidates no admitted upper bound. The existing exponent P-1 is consistent with the necessary exponent (P-1)/3. It says a successful eventual-sublinear route must weaken/replace this strong-mean consumer, rather than only improve its control variate, exact Gaussian integration, smoothing, sample allocation, or ancestor reuse.

## 6. Genuine noncommuting two-dimensional extension

Let

    T = [[1/2,1/16],[1/16,1/2]],   E11=diag(1,0),
    V_theta(x) = (1/2) x^T T x
       + eta A w^2 sum_i theta_i Phi((x1/sqrt(A)-a_i)/w),
    eta=1/8.

The Hessian is T+eta theta_i phi'(t) E11 at a point in cell i. Its spectrum stays in [5/16,11/16] by Weyl's elementary bound. At two points with different perturbation coefficients, their commutator is a nonzero scalar multiple of [T,E11], since T12=1/16. Thus this is a same-potential genuinely noncommuting fixture.

Every vector gradient query or first/adjoint action still reveals at most one sign, because only x1 selects a cell. Under the standardizing x=sqrt(A)u, integrating out u2 leaves a centered one-dimensional Gaussian in u1 with precision

    s_A = 1+A/2 - (A/16)^2/(1+A/2),

times exactly the same exp(-R_theta(u1)) factor. For A<=1 its variance lies in a fixed compact positive interval. The first canonical-force component is -E U1 and has the same covariance sensitivity and query lower bound. The second component and all cross-Hessian actions cannot reveal an additional hidden sign.

This extension is not needed logically, because a uniform multidimensional theorem must already handle D=1. It verifies that the obstruction is compatible with the noncommuting class, rather than being an artifact of commuting quadratics.

## 7. What a constructive replacement must change

A viable weaker contract can target the buffered law

    Law[host + B Z + physical_row * M]

directly without demanding an executable M be strongly close to the deterministic mean. The final independent Gaussian, nonlinear decoder and full retained caller must be part of that theorem. Estimating the mean after the fact from repeated outputs would incur new inverse-accuracy work, so there is no contradiction with the barrier above.

In particular the lower bound does not prove that weak buffered laws require a linear heat exponent; it does not price higher-order law counterpackets or their complete offspring. Those are separate constructive gates. It does show that a supposed no-copy mean repair returning the existing strong L2 mean port at sublinear heat cost cannot be closed under the current oracle and uniform C2 assumptions.

## 8. Proof and diagnostic boundary

The proof is Sections 2-5, using exact posterior identities, a finite independent-sign prior, and elementary orthogonal projection. Numerical checks only test formulas, normalization-sensitive finite differences, the fixed-Hessian quadratic case, the noncommuting extension, and exponent bookkeeping. They are not evidence of a general oracle lower bound by themselves.

Baseline sources: NEXT-COMPLETE-MEAN-COST-GATE.md; ANY-ORDER-DIRECT-MEAN-OUTER.md; ELEMENTARY-SEED-AND-COMPLETE-SOURCE.md; LINEAR-MILESTONE-AND-SUBLINEAR-COST-MAP.md. No old Fenchel/opposite-current module is used or re-tested.
