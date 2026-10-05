# Odd-cumulant cancellation: exact capability and target identity

2026-10-04. This answers whether coefficients and independent averaging make the new order gain repeatable. It separates cancellation useful for a Gaussian mean service from preserving the third current of a non-Gaussian target law.

## 1. What coefficient averaging does exactly

At an exposed retained caller u, let F_i be independent COMPLETE copies of the SAME conditional source law. For deterministic real a_i and every existing cumulant tensor,

    kappa_k(sum_i a_i F_i | u) = (sum_i a_i^k) kappa_k(F|u).

The source's full internal descendants must be copied; shared hidden ancestors between copies generally invalidate this formula. Coefficients/schedules must be frozen before differentiation. A negative coefficient is ordinary deterministic arithmetic in a positive random-variable construction, not a signed probability.

The existing nine-copy near-gradient mean uses eight coefficients 1/6 and one coefficient -1/3. Exactly,

    sum a_i=1, sum a_i^2=1/3, sum a_i^3=0,
    sum a_i^4=1/54, sum a_i^5=-1/324.

Thus it preserves the mean and kills the third cumulant, not all odd cumulants. A fixed finite list can cancel a prescribed finite list of odd ranks with computable order-dependent coefficients/copies. Its fourth-power sum is positive for any nonzero real list; additive independent Gaussian noise cannot remove that fourth cumulant.

No FINITE coefficient list preserves arbitrary means while cancelling every odd rank k>=3. To prove this, group its nonzero coefficients by distinct absolute magnitude b_j>0 and let d_j be positive-minus-negative multiplicities. All higher odd cancellations imply

    sum_j d_j b_j^3 (b_j^2)^m = 0,   m=0,1,... .

The first number-of-distinct-magnitudes equations form an invertible Vandermonde matrix, so d_j=0 for every j. Consequently sum_i a_i=0, contradicting mean preservation. This obstruction is only to that finite linear iid-copy construction.

## 2. Exact all-odd symmetrization after separating the mean

Let m=E(F|u), X=F-m, Sigma=Cov(F|u), and take complete iid F_1,F_2. The executable source

    D=(F_1-F_2)/sqrt(2)

is exactly symmetric, centered, and has covariance Sigma. Every existing odd cumulant vanishes. If F is L-Lipschitz in its COMPLETE standard Gaussian tape, then D is L-Lipschitz in the combined independent tapes: squared block first bounds sum to L^2. Its centered energy is exactly e=(E|X|^2)^(1/2), since Cov(D)=Sigma. No uncomputed m is used to execute D.

To preserve a nonzero target mean one still needs an executed mean service and must pay its actual covariance and errors. Simply adding m is not an oracle operation. A unit-buffer mean law supplies an additional independent Gaussian variance share, which must fit the bridge reserve. One cannot use the mean service one is trying to construct as an unproved premise of its own construction.

A separately derived iterated Stein identity yields a fourth-order matching-Gaussian comparison when kappa_3=0, under the same fixed independent buffer. This improves the Gaussianization of D. It is a theorem about the law of D, not a claim that D has the law of X.

## 3. Why this does not raise the raw posterior bridge by itself

The order-three endpoint currently compares a compiled law to

    h[Z-F_cont]+sqrt(1-h^2)N,

where F_cont is the SAME continuous second-substitution force on the Markov OU path. At a fixed retained Z, the private centered F_cont generally has a nonzero third cumulant of order A^3. That cumulant belongs to the analytical law reference; it is not only Monte Carlo noise introduced by a mean compiler.

Replacing centered F_cont by D removes this true current. The new D-plus-Gaussian law is more Gaussian, but its distance to the old target can still be order A^3. Thus it cannot establish an A^4 endpoint bound merely by applying the improved Gaussianization theorem. One must restore the correct third current by an actual positive source, or use a different exact reference and prove its full target identity.

By contrast, the desired target of a Gaussian mean service N(m,Sigma_fill) has zero third cumulant. There the nine-copy cancellation is useful and is ALREADY in the imported near-gradient mean constructor. The re-entry mean's remaining A^4 allowance includes its covariance/curl orientation term e_E*kappa_E, not an uncancelled third cumulant. Further odd cancellation alone does not remove that term.

## 4. Explicit C2 Markov-path separator

This section shows that the third current in Section 3 is genuinely present within the allowed convex C2 class, not merely possible for an unrelated source.

Work in D=1, condition the stationary Markov OU path on Z=X_1=0, and write

    Cov(X_r,X_s)=min(r,s)/max(r,s)-rs,    0<r,s<1.

Choose epsilon=1/20 and the smooth anchored gradient

    g_0(x)=[x+epsilon log cosh(x)]/(1+epsilon),
    g_A(x)=A g_0(x).

Its Hessian satisfies

    0 < (1-epsilon)/(1+epsilon) <= g_0'(x) <= 1,

so g_A obeys the required 0<=g_A'<=A and g_A(0)=0. Define the analytical leading force

    H_0=int_0^1 g_0(X_r)dr=(L+epsilon J)/(1+epsilon),
    L=int_0^1 X_r dr,   J=int_0^1 log cosh(X_r)dr.

The entire path is invariant under sign reversal. L is odd and J is even, so after centering J_c=J-EJ,

    kappa_3(H_0)=[3epsilon E(L^2 J_c)+epsilon^3 E(J_c^3)]
                    /(1+epsilon)^3.                    (4.1)

The exact covariance identity is

    Cov(L,X_r)=int_0^1 Cov(X_s,X_r)ds=-r log r.

Gaussian integration by parts twice therefore gives

    E(L^2 J_c)=int_0^1 r^2(log r)^2 E[sech^2(X_r)]dr.

Since Var(X_r)<=1, P(|X_r|<=1)>1/2 and sech^2(1)>1/4. Thus

    E(L^2 J_c) >= (1/8)int_0^1 r^2(log r)^2 dr =1/108.

Also 0<=log cosh x<=x^2/2, and Minkowski gives

    ||J||_3 <= (15)^(1/3)/3,
    |E J_c^3| <= ||J_c||_3^3 <= 40/9.

At epsilon=1/20, the numerator in (4.1) is at least

    3/(20*108)-(40/9)/20^3 = 1/1200.

Consequently

    kappa_3(H_0) >= 20/27783 > 0.                        (4.2)

The exact second-substitution force has

    F_cont=A int_0^1 g_0(X_r-A int_0^1 g_0(X_(r tau))d tau)dr
          =A H_0+O_Lp(A^2),

for every fixed p, by Lip(g_0)<=1 and Gaussian moment/Minkowski bounds. Therefore

    kappa_3(F_cont|Z=0)=A^3 kappa_3(H_0)+O(A^4).

In particular it is bounded below by a positive computable constant times A^3 for sufficiently small A. A deterministic finite path approximation is used only to justify the Gaussian calculations and passed to its Lp limit; no continuum is executed as a query oracle.

For a direct W2 separation, center the conditional force and add any fixed independent scalar Gaussian buffer sigma N. The symmetrized force D has a symmetric law, hence E sin(hD+sigma N)=0. On the original source,

    E sin(h[F_cont-EF_cont]+sigma N)
       =-exp(-sigma^2/2) h^3 kappa_3(F_cont)/6+O(A^5).

The linear term is zero; the fifth-moment remainder follows from the fixed-p bounds. Since sin is 1-Lipschitz, this is a lower bound on W1 and hence W2. For fixed h>0 and small A,

    W2(original buffered law, symmetrized buffered law) >= c_h A^3.

Adding the SAME deterministic mean does not change this lower bound. This separator is conditional at an allowed retained caller; it does not by itself assert an unconditional lower bound after mixing a different caller distribution. It is sufficient to refute a claimed uniform conditional A^4 replacement at this port.

## 5. Complete-cost implication and scope

At fixed target odd ranks, a finite coefficient pattern only multiplies the complete source cost by its literal number of independent copies. This can be an order-dependent constant and introduces no inverse-A exponent by itself. It still pays all original-gradient ancestors, anchors, conditional modes, private roots, covariance reserve actions, precision and changed-source replays. It must preserve the actual conditional target and variance budget.

Thus odd cancellation is real and already part of the mean toolbox. It is not by itself a repeatable all-order posterior-law constructor. The first remaining positive-law task is to reproduce the required third current (or bypass it with a proved alternate full-law reference), followed by all even-current and covariance/orientation repairs. This note proves no impossibility for other positive law constructions and no c(P)/P conclusion.
