# Connected-cumulant resolvent recurrence for the exact discounted OU force

2026-10-04. Analytical target construction, not a finite original-VALUE compiler. This note isolates an arbitrary-order connected clock recurrence without assuming a finite path quadrature has the right law.

## 1. Exact target and generating equation

Let X_t^z be standard OU with generator L=Delta-x dot grad, started at X_0=z. Let g:R^D->R^D be globally A-Lipschitz, with g(0)=0, and define the analytical discounted additive functional

    H_z=int_0^infinity e^(-t) g(X_t^z)dt
       =int_0^1 g(X_r)dr,

where the second notation is the conditional Markov path ending at X_1=z used by the existing endpoint proof. This is not the cheap same-root packet. Let

    M(z,lambda)=E exp(lambda dot H_z),
    K(z,lambda)=log M(z,lambda).

The Gaussian-Lipschitz first bound and finite Gaussian energy of H imply all finite exponential moments of each scalar projection. The conditional OU Markov property gives, for every t>=0,

    M(z,lambda)=E_z[exp(lambda dot int_0^t e^(-s)g(X_s)ds)
                                  M(X_t,e^(-t)lambda)].

Differentiating this dynamic identity at zero gives

    (L_z-lambda dot partial_lambda)M+(lambda dot g)M=0,
    (L_z-lambda dot partial_lambda)K+|grad_z K|^2
                                              +lambda dot g=0.   (1)

For a smooth regularization these are classical local identities. Under the stated C1/Lipschitz source, finite-order differentiation at lambda=0 is justified by the finite Gaussian moments; the spatial identities can be understood weakly and passed from smoothing. No third derivative of the potential is required or executed.

## 2. Exact connected hierarchy

Write the tensor cumulants at the SAME fixed z as

    K(z,lambda)=sum_(k>=1) kappa_k(z):lambda^(tensor k)/k!.

Only a finite jet at lambda=0 is needed for any fixed k; no global convergence of the series is asserted. With standard averaging symmetrization over the k output slots,

    (1-L)kappa_1=g,
    (k-L)kappa_k
       =sum_(i=1)^(k-1) binom(k,i)
                     Sym(D kappa_i dot D kappa_(k-i)),  k>=2.      (2)

The dot contracts the spatial derivative/root index, not a physical tensor slot. All connected mean subtractions are already included because K is the log moment-generating function. Raw products of moments cannot replace these kappa_i.

The resolvent has the positive OU-clock formula

    R_k F=(k-L)^(-1)F=int_0^1 q^(k-1) P_q F dq.                    (3)

The correct solution is selected by the growth inherited from the original conditional moments, or directly by the Markov representation; arbitrary homogeneous solutions of the elliptic equation are not introduced.

Let v=kappa_1=R_1 g, B=Dv, and C=kappa_2. Since g is a gradient, B is symmetric. Equations (2)-(3) give

    C=2 R_2(B B*)=2 int_0^1 q P_q(B^2)dq,                          (4)
    kappa_3=6 R_3 Sym(B dot DC)
            =6 int_0^1 q^2 P_q[Sym(B dot DC)]dq.                  (5)

Equation (4) is exactly the independently audited positive square-clock covariance identity. Equation (5) is its connected third-rank successor. It has the correct original common coarse argument inside B and DC; replacing them by factors at unrelated centers changes the target.

For arbitrary non-gradient g use BB* in (4). No commutation of matrix-valued factors is implicit in (2), (4) or (5).

## 3. Independent martingale check of the third coefficient

The conditional martingale for H is

    M_t=int_0^t e^(-s)g(X_s)ds+e^(-t)v(X_t),
    dM_t=sqrt(2)e^(-t)B(X_t)dB_t.

Here M_0=v(z), and M_infinity=H in Lp. Its quadratic variation proves (4). Ito applied to the centered cubic tensor gives

    kappa_3=6 Sym int_0^infinity e^(-2t)
                  E_z[(M_t-v(z)) tensor B(X_t)B(X_t)*]dt.

Pair the martingale M_t-v(z) with the Doob martingale of the terminal matrix B(X_t)B(X_t)*. The covariance identity gives

    E_z[(M_t-v(z)) tensor A(X_t)]
       =2 int_0^t e^(-s)
          P_s[B dot D P_(exp(-(t-s))) A](z)ds,

where the displayed P_s on this line denotes time-s OU evolution, equivalently P_(exp(-s)) in correlation notation. Substituting A=BB*, then changing variables q=exp(-s), u=exp(-(t-s)), gives

    kappa_3=12 Sym int_0^1 q^2 dq int_0^1 u du
                              P_q[B dot D P_u(BB*)].

Since DC=2 int_0^1 u D P_u(BB*)du, this is exactly (5). Gaussian smoothing regularization justifies any rough matrix coefficient before passing to the conditional moments. This is an analytical check, not permission to query DB or Dg as producer values.

## 4. Useful dimension-safe bounds at the first successor

At a fixed z, the private Gaussian-Lipschitz first of H is at most A, so Cov(H|z)<=A^2 I. Its actual z derivative has operator norm at most A/2, because D_zX_t=e^(-t)I and int_0^infinity e^(-2t)dt=1/2.

For any unit z-direction u, differentiate the covariance using ONLY this first derivative:

    D_u C=Cov(D_u H,H)+Cov(H,D_u H).

The matrix covariance inequality and the bounded vector D_u H yield

    ||D_u C||HS <=2 sqrt(||Cov(H)||op) ||D_uH-E D_uH||2
                   <= A^2.                                      (6)

Thus DC has operator norm R^D->HS at most A^2, and full HS norm at most A^2 sqrt(D). Also ||B||op<=A/2. Therefore the source inside (5) has full HS norm at most (A^3/2)sqrt(D), giving the conservative bound

    ||kappa_3(z)||HS <= A^3 sqrt(D).                               (7)

This is a one-physical-energy estimate uniform in z. The general Gaussian-image Stein argument gives another independent bound of this type. None of these norm bounds exposes kappa_3 as an executable tensor oracle.

## 5. What positive clock quadrature alone supplies

For any Hilbert-valued F in Gaussian L2, R_k is a scalar Hermite multiplier 1/(k+n) on chaos n. Thus any positive finite rule with

    sup_(n>=0)|sum_j w_j q_j^(k-1+n)-1/(k+n)|<=delta

satisfies ||R_k F-Q_k F||_2<=delta||F||_2. The already admitted positive dyadic rule, with shifted chaos degree, supplies this at fixed k using O_k(log^2(1/delta)) nodes. This applies to the WHOLE connected coefficient field F at its actual Gaussian carrier law.

It does not justify multiplying separately approximated high-order tensor fields without their actual product/first bounds. It also does not provide DC or the products in (2) as VALUE sources. Those require finite native response/gradient-lift producers, with their actual common-root laws and complete banks. A quadrature proof for each scalar expectation is insufficient for the final full-law action.

## 6. Consequence for the finite construction problem

The hierarchy (2) is a repeatable analytical target with explicit coefficients and no direct finite-path diagonal. It suggests a native tree construction from the same initial genuine-gradient g and positive OU resolvents. A successful finite compiler must still supply, at every order:

- a tensor-free original-VALUE action for the coefficient at the SAME captured caller;
- the proper-cut/one-energy norms needed by products and response feedback;
- complete conditional mean, covariance, even/odd current and numerical descendants;
- positive variance reserve and an exact same-endpoint law comparison;
- old-root ownership, gradient reference/restoration, caller/first/zero ports and complete original-query costs.

The hierarchy only describes cumulants of H. The true posterior reference at order p additionally has nonlinear resolvent/Picard feedback. Those feedback currents cannot be ignored or replaced by the additive H cumulants at every rank. No finite all-order positive posterior program and no c(P)/P theorem is established here.
