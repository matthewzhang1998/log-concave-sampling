# Single-history mean: exact centered/resummed covariance current

2026-10-05. Exact analytical identities under C2. The proposed fourth-order reduction to an unshifted covariance current remains OPEN; this note neither proves nor refutes it. No finite path-grid producer is supplied.

## Result and scope

Let g=grad U, g(0)=0, 0<=Dg<=A I. Write

    I_H(x,V)=int_0^infinity e^(-t) g(e^(-t)x+L_t V) dt,
    I_Q(x,H)=sum_j v_j g(tau_j x+d_j H),
    d_j=sqrt(1-tau_j^2), sum_j v_j=1,
    mu_a(x)=E I_a(x), xi_a=I_a-mu_a, a in {H,Q}.

Here V is the genuine conditional OU Brownian input, L_t L_u*=c(t,u)I with c(t,u)=e^(-|t-u|)-e^(-t-u). H is one common independent standard root for the entire cheap packet. The positive clock rule has exact first moment 1/2. All histories in this note are analytical comparison objects.

There is an exact, coherent-mean-preserving identity comparing the TRUE continuous I_H directly to the FINITE cheap I_Q. It uses a convolution homotopy of the actual centered laws, not equality of their moments. Put, with the two private banks independent conditional on x,

    mu_theta=(1-theta)mu_Q+theta mu_H,
    F_theta=mu_theta+sqrt(theta)xi_H+sqrt(1-theta)xi_Q,
    T_theta=x-F_theta, Delta_mu=mu_H-mu_Q.

Then

    j_H-j_Q = -int_0^1 E[Dg(T_theta)] Delta_mu dtheta
              +(1/2) int_0^1 E[D2g(T_theta):(tau_H-tau_Q)] dtheta.   (1)

The Stein matrices tau_H and tau_Q are the complete Gaussian Riesz coefficients defined below, each with its own independent auxiliary rotation bank. Their expectations are exactly the corresponding conditional covariances. Formula (1) is a WHOLE weak current under C2. D2g is never a bounded assumption, an executable coefficient, or a separate pointwise factor in the nonsmooth limit.

If the mean quadrature guarantee is ||Delta_mu||_(L2 gamma)<=delta A sqrt(D), the first term after R1 costs at most delta A^2 sqrt(D). There is no additional strong path-quadrature assumption in (1).

The missing fourth-order theorem is stated precisely in Section 5. The exact identity does not close that theorem, and a generic Hilbert/Bessel estimate does not close its contracted vector trace.

## 1. Exact Riesz coefficients with their complete genealogy

For a finite-dimensional standard Gaussian source B and a square-integrable C1 vector function F(B), take an independent B' and B_rho=rho B+sqrt(1-rho^2)B'. Its analytical Stein coefficient is

    tau_F(B)=int_0^1 E_[B'][ DF(B_rho) DF(B)* ] drho.             (2)

For a smooth scalar phi,

    E[(F_i-EF_i) phi(F)] = sum_j E[(tau_F)_ij partial_j phi(F)],
    E tau_F=Cov(F).                                             (3)

This follows from the Gaussian covariance/Riesz formula. It is NOT a claim that F is Gaussian, that tau_F is deterministic, or that tau_F is pointwise symmetric. Infinite Gaussian histories follow by finite analytical approximations and Sobolev closure.

For I_H, let V_rho=rho V+sqrt(1-rho^2)V' and

    J_t^rho=Dg(e^(-t)x+L_t V_rho), J_u=Dg(e^(-u)x+L_u V).

Then

    tau_H=int_0^1 drho int_0^infinity dt du e^(-t-u)c(t,u)
                         E_[V'][J_t^rho J_u | V].               (4)

For I_Q, with H_rho=rho H+sqrt(1-rho^2)H' and J_j^rho=Dg(tau_j x+d_j H_rho),

    tau_Q=int_0^1 drho sum_(i,j) v_i v_j d_i d_j
                         E_[H'][J_i^rho J_j | H].               (5)

The order of each matrix product is intentional; the final D2g contraction only sees its symmetric part. All expectation identities keep the entire relevant private bank. In particular, tau_H is correlated with xi_H, and tau_Q with xi_Q, in (1). They are not replaced by their expectations.

The bounded first derivatives give

    ||tau_H||op <= A^2/4,
    ||tau_Q||op <= A^2 beta_Q^2, beta_Q=sum_j v_j d_j<=1.          (6)

The first constant uses int int e^(-t-u)c(t,u)dt du=1/4. These are analytic bounds, not producer declarations. The corresponding unconditional physical energies satisfy ||I_a||2<=A sqrt(D), ||xi_a||2<=A sqrt(D), with harmless improvements available for H.

## 2. Derivation and C2 definition of the exact current

Initially take a smooth regularization with the same Hessian bound and zero anchor. Set f_theta(x)=E g(T_theta). Direct differentiation on 0<theta<1 gives

    partial_theta f_theta
      =-E[Dg(T_theta)]Delta_mu
       -(1/(2sqrt(theta))) E[Dg(T_theta)xi_H]
       +(1/(2sqrt(1-theta))) E[Dg(T_theta)xi_Q].                  (7)

Apply (3) in the H-history bank while keeping the Q bank fixed. The derivative of T_theta with respect to xi_H is -sqrt(theta)I. Thus its term becomes +(1/2)E[D2g(T_theta):tau_H]. Applying (3) in the Q bank gives -(1/2)E[D2g(T_theta):tau_Q]. Integrating theta proves (1), with the sign TRUE covariance minus CHEAP covariance.

Equation (7), including all its correlations, is also the derivative-free definition of the combined current on the right of (1). Its endpoint singularities are integrable because Dg is bounded and xi_a has finite L2 energy. More explicitly, its non-mean portion is bounded in L2(x~gamma) by

    (A/2)[ ||xi_H||2/sqrt(theta) + ||xi_Q||2/sqrt(1-theta) ].

Mollify U, subtract the resulting constant gradient at zero, and let the smoothing width decrease. The regularized gradients and Hessians converge locally to g and Dg, with uniform Lipschitz/linear-growth bounds. All terms of (7), its theta integral, and the endpoint j_H-j_Q converge by dominated convergence. This defines the whole D2g current uniformly under the original C2 hypothesis. Neither isolated D2g factors nor their absolute expectations are asserted to converge.

The outer R1 is an L2(gamma) contraction. Therefore

    ||R1 int E[Dg(T_theta)]Delta_mu dtheta||2
                    <= A ||Delta_mu||2 <=delta A^2 sqrt(D).     (8)

This is the precise mean-quadrature floor. It is valid although the mean vector may have O(1) norm when D=A^(-2).

The existing independently admitted second-order decoupling theorem, applied separately to I_H and I_Q, also gives the already-known total bound

    ||R1(j_H-j_Q)||2 <=C A^3 sqrt(D)+delta A^2 sqrt(D).           (9)

Combining (8)-(9) bounds the WHOLE current in (1) at C A^3 sqrt(D)+2delta A^2 sqrt(D). This is a one-physical-energy bound on the combined target, not on every absolute integrand. It does not improve the old mean order.

## 3. A second exact identity when the clock measures agree

This formulation is useful for seeing precisely what the proposed unshifted current discards. First use any common finite positive clock rule for both a Markov packet I_M,Q and the common-root packet I_Q. In tau coordinates define

    k_M(tau,sigma)=min(tau,sigma)/max(tau,sigma)-tau sigma,
    k_C(tau,sigma)=sqrt(1-tau^2)sqrt(1-sigma^2),
    Delta_k=k_M-k_C.

Let the centered Gaussian clock process Y^lambda have covariance lambda k_M+(1-lambda)k_C. It exists as sqrt(lambda)Y^M+sqrt(1-lambda)Y^C with independent whole banks. Set

    I_lambda=sum_j v_j g(tau_j x+Y_j^lambda),
    J_j^lambda=Dg(tau_j x+Y_j^lambda).

The marginal at each individual clock is unchanged. Consequently E I_lambda=mu_Q for every lambda; this identity retains the common coherent conditional mean exactly.

Gaussian covariance interpolation gives

    j_M,Q-j_Q
      =(1/2)int_0^1 dlambda sum_(i,j) v_i v_j Delta_k(i,j)
             E[D2g(x-I_lambda):(J_i^lambda J_j^lambda)].         (10)

There are no inner D2g terms: all diagonal covariance increments Delta_k(i,i) are zero, and those are exactly the terms that could differentiate one inner VALUE twice. Mixed clock derivatives hit two different inner sites and give (10). The expression remains valid if clock values coincide, because their covariance increment is then zero.

The same interpolation applied to E[I_lambda I_lambda*] yields

    C_M,Q-C_Q=int_0^1 dlambda sum_(i,j) v_i v_j Delta_k(i,j)
                                    E[J_i^lambda J_j^lambda].  (11)

The full sum on the right is symmetric. Formula (11) uses constant means, not moment-to-law inference.

Equations (10)-(11) also hold with the common Lebesgue clock measure, yielding the true continuous history versus a continuous common-root comparison. This passage is a finite-dimensional analytical approximation, with strong L2 limits and uniform first/energy majorants. It is NOT a polylogarithmic executed Markov quadrature. In particular, (10) with finite Q does not by itself compare I_M,Q to I_H at order four. The direct finite-Q identity (1) avoids making that assertion.

For Lebesgue clocks, k_M<=k_C and

    int int k_M=1/4, int int k_C=pi^2/16,
    int int |Delta_k|=(pi^2-4)/16.                               (12)

This sign statement concerns scalar clock covariances; it does not assert a Loewner ordering between arbitrary vector-force covariance matrices.

## 4. Exact relation to an unshifted or a centered leading term

Write Delta_C=C_H-C_Q. From (1), add and subtract the formal unshifted whole current to obtain

    R1(j_H-j_Q) = (1/2)R1[Delta_C:D2g(x)] + E_mean + E_trace,

    E_mean=-R1 int_0^1 E[Dg(T_theta)]Delta_mu dtheta,
    E_trace=(1/2)R1 int_0^1 E[(tau_H-tau_Q):
                                    (D2g(T_theta)-D2g(x))]dtheta. (13)

Equations (13) are literal for smooth regularizations. Under C2 each proposed separated current must itself be admitted as a whole weak limit; equation (1), or the derivative-free form (7), is the definition that is already proved. One cannot assume separate convergence merely because their sum converges.

If mu_H=mu_Q=mu, a mean-retaining alternative is

    (1/2)R1[Delta_C:D2g(x-mu(x))],                              (14)

with exact residual obtained by subtracting D2g(x-mu) inside (13), while T_theta=x-mu-sqrt(theta)xi_H-sqrt(1-theta)xi_Q. With unequal means the analogous leading current is the theta integral at x-mu_theta. Neither choice erases the random Stein-coefficient correlation. Centering is an exact bookkeeping improvement, not a proof that its residual is order four.

The older coherent-shift covariance counterexample does not refute (13): its single-I rank-one contribution at D=A^(-2) is at the allowed A^4 sqrt(D)=A^3 scale. No contrary claim is made here.

## 5. Precise missing analytic and finite contracts

To admit the advertised unshifted candidate, one must prove the following for this source class, uniformly over smooth regularizations and over the physical dimension:

    ||E_trace||_(L2 gamma) <= Lambda A^4 sqrt(D),               (15)

together with existence and consistent convergence of both separated whole currents in (13). Lambda may contain the explicitly admitted public logarithms and priced analytic cuts, but no extra physical sqrt(D) or inverse-A power. A centered/resummed replacement can instead specify its actual leading target and prove the corresponding residual and convergence statements.

An equivalent same-clock version asks for (15) with the weighted integrand

    Delta_k(i,j) [D2g(x-I_lambda)-D2g(x)]:(J_i J_j).

The augmented outer-Gaussian/early-bridge score can preserve the affine Hessian vertices. One score provides a sound one-energy matrix Bessel estimate in the mixed-K theorem. That theorem does NOT automatically bound the present contracted vector trace. Two unstructured Gaussian scores produce the Hilbert norm of the coefficient, or the norm of a quadratic score, and can introduce an additional sqrt(D). Likewise an arbitrary small-energy, small-curl vector field is not automatically covered by the gradient-input Hermite contraction used for a one-Hilbert connected block.

A sufficient new lemma must use the actual gradient/chord structure of g(x-F)-g(x), every live ancestor, and an explicit trace/Hodge or equivalent cancellation. It must prove the contraction norm, not cite generic tensor Bessel as if an operator norm were a Hilbert-Schmidt norm. No such lemma is proved in this note.

Even if (15) passes, a separate finite native consumer is needed for the retained whole D2g response current. It must use only original g VALUES, preserve the correct covariance genealogy and captured Gaussian endpoint, replay all nonlinear ancestors at changed arguments, and price positive clock compression, endpoint cuts, absolute floors, actual first/caller/origin ports, and restored mode versions. Neither tau_H nor D2g nor a saved HVP is a producer. The admitted finite covariance service alone does not supply this mean-current consumer.

## Status

PROVED ANALYTICAL IDENTITIES: (1)-(8), (10)-(12), with the C2 whole-current interpretation. Equation (9) imports the audited old decoupling bound. The order-four estimate (15), the separate limiting currents in (13), and the finite VALUE current consumer remain OPEN. Independent review is required before any broader PASS.
