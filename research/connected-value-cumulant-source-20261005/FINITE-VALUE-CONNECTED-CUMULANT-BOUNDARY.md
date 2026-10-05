# Connected VALUE cumulants: exact replica construction and a dimension-uniform source obstruction

2026-10-05. Independent construction/test. The result is negative for a precisely stated coefficient-estimator class, not for all native VALUE sources.

## Result first

1. There is a completely explicit, derivative-free, conditional-mean-oracle-free replica kernel for every conditional cumulant. All partition subtractions can be performed on one common finite Gaussian record before any further integration.
2. That construction fails the requested one-Hilbert-energy source contract already at rank four. An exact all-distinct-coordinate calculation gives energy squared of order D^n, including at the critical rank eight. This is true even for the original admissible linear gradient and its literal additive conditional OU history.
3. More generally, any executed random coefficient tensor whose first physical mode lies in the span of M evaluated VALUE vectors needs M at least proportional to D to have one-HS energy while approximating the desired rank-four or rank-eight cumulant. The proof is a rank/nuclear-norm argument. It survives arbitrary coupling, antithetic or Gaussian rotation groupings, data-dependent scalar coefficients, clipping, and approximation error smaller than a fixed fraction of the target. No cutoff or frequency assumption is involved.
4. The remaining possible route is genuinely different: an active bounded vector/gradient source whose cumulant tensor appears only after native integration. Such a source is not an executed random tensor estimator and need not obey the pointwise rank bound. The present note neither constructs nor excludes it. The all-rank analytical Appell/Riesz bounds do not themselves provide it.

## 1. Literal conditional OU target

Fix the retained caller y and write the conditional OU process as

    X_t = exp(-t)y + sqrt(2) integral_0^t exp(-(t-s)) dW_s,
    H_y = integral_0^infinity exp(-t) g(X_t) dt.

This is the same additive history as the correlation-clock integral after r=exp(-t). Each physical coordinate has an independent Gaussian future when g is coordinate-separable. A finite positive clock rule gives

    H_Q(y,Z) = sum_j w_j g(d_j y + L_j Z),

on its exact joint Gaussian history record, with the same future shared among its time nodes. The following replica formula is exact for H_Q. It is also exact as an analytical formula for H_y. Replacing H_y by H_Q has its own quadrature allowance; the formula does not erase it.

## 2. Exact unknown-mean-free connected kernel

Let B_1,...,B_n be independent CONDITIONAL copies of H_Q at the same y. All replicas use the same frozen time rule and their own complete future record. For a partition pi of [n] with k blocks, let

    mu(pi) = (-1)^(k-1) (k-1)!,
    (n)_k = n!/(n-k)!.

Assign the blocks injectively to n replica labels and define

    K_n = sum_pi [mu(pi)/(n)_|pi|]
                 sum_(injective a:pi->[n]) tensor_(slot i) B_(a(block_pi(i))).

The slot order is the physical tensor order. The expression is symmetric after the whole partition sum. Every occurrence with the same replica label is the SAME B vector, not a separately resampled occurrence. Taking the conditional expectation gives exactly the moment-cumulant partition formula:

    E[K_n | y] = kappa_n(H_Q | y).

No E[H_Q|y], covariance, or cumulant is queried. The kernel uses exactly n complete H_Q copies, hence n times its original-g query count, plus finite rank-dependent tensor arithmetic. The unexpanded product count can be very large in n but has no hidden derivative oracle. This gives the exact algebraic construction requested by a replica/Mobius proposal.

### 2.1 Exact off-diagonal energy after the COMPLETE partition grouping

Suppose the coordinates of B are independent and identically distributed, centered, with variance sigma^2. This includes y=0 and g(x)=a x in the original class.

For n DISTINCT physical coordinates i_1,...,i_n, each map a:[n]->[n] selects one distinct monomial. Its coefficient in the grouped kernel is

    (-1)^(k-1)(k-1)!/(n)_k,  k=|image(a)|.

Different maps are orthogonal in L2: some physical coordinate has two different replica labels, whose centered values are independent. There are S(n,k)(n)_k maps with k labels. Consequently

    E |(K_n)_(i_1,...,i_n)|^2 = c_n sigma^(2n),
    c_n = sum_(k=1)^n S(n,k) ((k-1)!)^2/(n)_k.

The exact constants are

    c_4 = 10/3,
    c_8 = 46212/35.

There are (D)_n distinct ordered coordinate tuples, so

    E ||K_n||_HS^2 >= (D)_n c_n sigma^(2n).            (1)

This is a bound on the ALREADY GROUPED random source. It is not a sum of absolute values of partition terms, nor an independence assumption between partitions. On these off-diagonal components the target cumulant is zero, so (1) is also a source-variance lower bound.

For the literal linear conditional OU example, y=0 and g(x)=a x,

    Y=integral_0^infinity exp(-t)X_t dt,
    Cov(Y,X_t)=t exp(-t),  Var(Y)=1/4,
    H=aY, sigma^2=a^2/4.

All kappa_n with n>=3 vanish, yet (1) remains positive with the displayed dimension power. A finite positive history rule has a different fixed positive sigma^2 and the same conclusion.

Even using a known zero mean does not fix rank four. With independent centered B,C,

    K_4^centered = (B^tensor4 + C^tensor4)/2
                         -3 Sym(B^tensor2 tensor C^tensor2)

is unbiased. Its all-distinct component contains eight orthogonal monomials of coefficient magnitude 1/2, so its variance is exactly 2 sigma^8 and its full HS energy squared is at least 2(D)_4 sigma^8.

## 3. A stronger rank obstruction, including arbitrary coupling and clipping

Here is a useful broader source test that does not depend on Gaussian independence of the queried values, on the time grid, or on a particular Mobius representation.

Let T be an executed random rank-n tensor. Assume that at each sample its 1|(n-1) flattening has rank at most M. A sufficient condition is that its first physical mode lies in the span of at most M available VALUE vectors:

    T = sum_(a=1)^M V_a tensor W_a.

The tensors W_a may depend on ALL the same roots, values, scalar inner products, finite antithetic rotations, and clipping decisions. They need not be independent of V_a. Known transformed query vectors and any additional Gaussian vector occupying the first physical slot must simply be counted in M. Thus an arbitrary finite homogeneous outer-product kernel built from M queried vectors is included, as are its scalar-coefficient regroupings.

Suppose the target is

    K = c sum_(i=1)^D e_i^tensor n,   c != 0.

Every proper flattening of K has D singular values equal to |c|. In particular, at 1|(n-1),

    ||K||_* = D|c|,  ||K||_HS=|c|sqrt(D),  ||K||op=|c|.

If T is unbiased, convexity of nuclear norm and the rank bound imply

    D|c| = ||E T||_* <= E||T||_*
                         <= sqrt(M) E||T||_HS
                         <= sqrt(M E||T||_HS^2).

Therefore

    E||T||_HS^2 >= D^2 c^2/M.                         (2)

This conclusion is uniform in every Gaussian query row and every scalar coefficient. Coalescing sites, antithetic sites, arbitrary original widths and taking very large coefficients cannot change it.

### 3.1 Approximate unbiasedness and bounded preparation

Suppose instead the full tensor bias satisfies

    ||E T-K||_HS <= eta sqrt(D).

A D-by-D^(n-1) matrix has rank at most D, so the nuclear norm of this bias is at most eta D. Consequently

    E||T||_HS^2 >= D^2 (|c|-eta)_+^2/M.               (3)

For eta<=|c|/2, an O_n(A^n sqrt(D)) source-energy certificate with c=c_0 A^n, c_0!=0 fixed, requires M>=c'_n D. Rank-only, fixed-clock, and public-log M do not suffice.

Clipping or regularizing scalar coefficients preserves the span assumption. It therefore cannot improve the energy to the desired scale while retaining the stated bias allowance. The bound does not require a globally polynomial source or unbounded field.

If the same random tensor is required to have an L2 proper-cut bound, rank(T)<=M also gives

    E||T||op^2 >= E||T||_HS^2/M
                    >=D^2 (|c|-eta)_+^2/M^2          (4)

for that 1|(n-1) cut. Thus the executed coefficients fail both the one-energy and bounded-cut ports when M is sublinear in D.

### 3.2 Integrating the retained caller does not remove the obstruction

The same conclusion holds if the contract only controls an L2 average over a standard Gaussian caller. For a coordinate-separable source the conditional tensor is

    K(y)=sum_i c(y_i)e_i^tensor n.

The scalar c(y) is continuous, and the witness in Section 4 has c(0)!=0, so m=E|c(Y)|>0 for standard scalar Gaussian Y. Conditional unbiasedness and the rank bound give

    E_(Y,Q)||T||_HS^2 >= E_Y (sum_i |c(Y_i)|)^2/M
                              >= D^2 m^2/M.

If the integrated HS bias is at most eta sqrt(D), apply the nuclear triangle inequality and Cauchy-Schwarz after averaging Y. The lower bound becomes D^2(m-eta)_+^2/M. Thus the obstruction is not confined to the zero-probability caller y=0.

### 3.3 Gaussian-probe vector readouts do not hide the energy

Let P_1,...,P_(n-1) be fresh independent standard physical Gaussian probes, independent of the coefficient roots, and put

    F = T[P_1,...,P_(n-1)].

Successive Gaussian isometries give exactly

    E_P |F|^2 = ||T||_HS^2.

Using one public with its Wick/Hermite polynomial gives the corresponding fixed factorial multiple for symmetric coefficient slots. Thus merely avoiding materialization of T, and computing its probe contraction instead, does not evade (2)-(3). The independence here matters: an active source whose queried g sites depend on these probes is a different class.

## 4. A bounded-Hessian original-gradient witness with nonzero target at ranks four and eight

It remains to verify c!=0 in the ORIGINAL source class, rather than postulating a diagonal abstract tensor.

In one dimension take a=1/2 and

    g_epsilon(x)=a x+epsilon sin(x)+epsilon^2(1-cos(x)).

It is the derivative of

    U_epsilon(x)=a x^2/2+epsilon(1-cos(x))
                                +epsilon^2(x-sin(x)).

The exact anchor is g_epsilon(0)=0 and

    g'_epsilon(x)=a+epsilon cos(x)+epsilon^2 sin(x).

For 0<epsilon<=1/8, this lies strictly inside [0,1]. Hence the Hessian sandwich is uniform. Scaling g_epsilon by A gives the original 0<=Dg<=A class and scales every rank-n cumulant by A^n.

At y=0 write the scalar additive history as

    H_epsilon=aY+epsilon S+epsilon^2 C,
    S=integral_0^infinity exp(-t)sin(X_t)dt.

Cumulants are finite polynomials in epsilon because H_epsilon has that exact polynomial dependence and every required moment is finite. At epsilon=0, H_0 is Gaussian. Multilinearity and Gaussian integration by parts yield, for even n>=4,

    d/d epsilon kappa_n(H_epsilon)|_(epsilon=0)
      = n a^(n-1) (-1)^(n/2-1) I_n,

    I_n = integral_0^infinity t^(n-1) exp(-nt)
                     exp(-(1-exp(-2t))/2) dt.

Indeed kappa(Y,...,Y,sin(X_t)) equals

    Cov(Y,X_t)^(n-1) E[(d/dx)^(n-1)sin(X_t)],

with n-1 copies of Y. This identity is obtained by differentiating the joint Gaussian exponential tilt, or by repeated one-dimensional Gaussian integration by parts. The remaining slots from the epsilon^2 term contribute no first derivative.

The integral has explicit positive bounds

    exp(-1/2)(n-1)!/n^n <= I_n <= (n-1)!/n^n.

Thus the derivative is nonzero at BOTH n=4 and n=8 (negative at both). There is a fixed sufficiently small positive epsilon for which the scalar cumulant c is nonzero. No claim that every epsilon<=1/8 has the same sign is needed.

Apply this same scalar gradient independently in all D coordinates. Their conditional future histories are independent, so their joint cumulant tensor is exactly

    kappa_n(H_epsilon|y=0)=c sum_i e_i^tensor n.

This is the target in Section 3. The obstruction therefore concerns an actual anchored smooth convex gradient, with a fixed Hessian bound, and the genuine conditional additive OU history.

For a nonempty finite positive clock rule, replace I_n by the finite positive sum

    sum_j w_j Cov(Y_Q,X_(t_j))^(n-1)
                    exp(-Var(X_(t_j))/2).

With positive nonzero times and weights its covariance factors are positive. The same witness proves nonzero kappa_n for the finite history target, separately from the continuous clock approximation.

## 5. Fixed-geometry polynomial groupings have the sharper D^n failure

Section 3 supplies the robust lower bound D^2/M. For a fixed scalar-row query dictionary the standard outer-product representation gives a sharper, exact identity.

Let q_a=L_a Z, a=1,...,M, and V_a=g(q_a) be coordinatewise evaluations. Fix y=0. Let v=(g(L_1 z),...,g(L_M z)) denote the vector of values for ONE physical coordinate, and let G=E[v v^T]. If

    T = sum_(a_1,...,a_n) C_(a_1,...,a_n)
                              V_(a_1) tensor ... tensor V_(a_n),

where C is a fixed fully grouped scalar tensor, then for distinct physical coordinates

    E|T_(i_1,...,i_n)|^2
       = <C,G^tensor n C>
       = ||(G^(1/2))^tensor n C||^2.                  (5)

This is a Gram identity, including ALL cross terms. It is not a triangle estimate. Consequently, if the right side is positive,

    E||T||_HS^2 >=(D)_n <C,G^tensor n C>.             (6)

For g(x)=a x+b sin(x)+c(1-cos(x)), with b,c nonzero, and distinct nonzero scalar Gaussian rows L_a, G is positive definite. To see this, a relation sum u_a g(L_a z)=0 consists of a linear polynomial, a constant, and a finite sum of exponentials exp(+/-i L_a z). Distinct frequencies are linearly independent. Within each antipodal pair, the nonzero sine and cosine coefficients force both u values to vanish. Thus all coefficients vanish. Exact duplicate and zero query rows should first be removed from the dictionary.

So deterministic antithetic/rotation groupings can make selected linear-source coefficients cancel, but a nonzero grouped C still has strictly positive (5) for this bounded-Hessian nonlinear source. Fixed geometry then gives the D^n energy failure. The coefficient in (6) can depend badly on the finite rows; unlike Section 3 this sharper assertion is not claimed uniformly as the dictionary changes.

## 6. Exact boundary of the result and next admissible construction question

These results DO NOT exclude:

- An active source f(P,Q) whose g query sites depend on the physical Gaussian publics P, where the desired tensor is obtained only after native integration/IBP. There may be no pointwise low-rank random tensor T(Q) to which Section 3 applies.
- A native gradient pair or filter implementing a high-rank matrix/tensor response after averaging, without strongly estimating that response.
- A source with executed high-rank deterministic identity/control tensors not contained in the finite-VALUE first-slot span. Such a proposal must supply its own source-qualified construction and one-energy/cut proof; it is not covered merely by calling it a Mobius kernel.
- A dimension-linear number of VALUE vectors, although that does not meet the desired rank-only/public-log count.
- A tensor estimator whose energy is allowed to grow as a power of D. It is exact but unsuitable for the requested source port.

There is also no all-source impossibility claim about cumulants. The imported analytical all-rank recurrence already proves that the TARGET has one-HS energy and all proper cuts. Equations (1)-(6) show why taking a finite unbiased connected tensor estimator as its VALUE producer does not inherit those target bounds.

A constructive successor therefore needs to state an active bounded vector/gradient source, list its exact original-g query rows and all caller/private firsts, and prove that its native integration produces the cumulant. Applying a Gaussian score to a raw tensor kernel or adding independent probes does not by itself do this. The low-first/caller and all-cuts contracts are not filled by an unbiasedness identity. No positive-consumer or all-order cost conclusion is licensed by the present replica construction.

## 7. Checks and source qualification

`check_connected_value.py` enumerates all set partitions through rank eight, independently checks the Stirling counts and exact variance factors, directly enumerates every replica assignment through rank five, checks the centered rank-four monomial sum, evaluates an exact antithetic Gaussian-value Gram matrix, checks the witness integrals against their proved bounds, and checks finite rank/nuclear inequalities on random small tensors. `checks.json` and `check-run.txt` record the run.

The analytic proofs above, rather than the numerical diagnostics, establish the result. No native pair compiler was executed, no source production or cutoff-uniform all-order theorem is claimed, and no external files were published.
