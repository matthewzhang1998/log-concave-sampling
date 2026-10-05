# All-fixed-rank centered Stein current with one energy

2026-10-04. Analytical recurrence, pending independent audit. This separates the fully repeatable law-comparison step from the still-required positive finite coefficient producers.

## 1. Definition and one-energy norm

Let X=f(V)-Ef(V), where V is standard Gaussian in any finite dimension, f is C1 and L-Lipschitz, and e=||X||2. Let R be the Hilbert-valued Gaussian Riesz operator grad(-L_OU)^(-1). Set

    T_0=X,
    T_j=Sym_(j+1)[R(T_(j-1)-E T_(j-1)) dot Df], j>=1,
    c_j=E T_j.

For j=1, T_0 is already centered and c_1=Cov(X)=Sigma. The dot contracts the old Gaussian derivative/root index and leaves j+1 physical slots. Every symmetrization is the standard average. No derivative of T_j is taken. Hilbert Riesz contraction, variance contraction and ||Df||op<=L give inductively

    ||T_j||_(L2;HS)<=L^j e.                           (1)

All these are analytical fields. They are not original-gradient producer queries and are not asserted to satisfy a native gradient source contract.

## 2. The constants are the actual connected cumulants

For every fixed j>=1,

    c_j=kappa_(j+1)(X)/j!.                            (2)

One proof is scalar polarization. Contract every output slot against the same theta. The first Stein identity and repeated centered Riesz identities give the finite integration-by-parts expansion against a polynomial test phi:

    E[(theta dot X) phi(theta dot X)]
      =sum_(j>=1) c_j[theta^(tensor(j+1))]
                            E phi^(j)(theta dot X),

where only finitely many terms survive for any polynomial. Comparing the resulting moment recursion with the ordinary moment-cumulant recursion gives c_j=scalar kappa_(j+1)/j!. Polarization identifies the full symmetric tensors. This uses finite Gaussian moments only, not convergence of an infinite cumulant series.

In particular

    ||kappa_(j+1)||HS <=j! L^j e.                     (3)

The recurrence remains valid when earlier non-Gaussian cumulants are nonzero: each previous T is centered explicitly before the next Riesz operation. None of those constants is replaced by a moment tensor without its cumulant subtraction.

## 3. Exact finite expansion at the same buffered endpoint

For independent standard Gaussian Z and fixed sigma>0, define

    Y_t=tX+(sigma^2 I+(1-t^2)Sigma)^(1/2)Z.

The start is the matching Gaussian and the end is X+sigma Z. The first covariance-preserving derivative is

    d/dt E phi(Y_t)=t E[(T_1-Sigma):D^2 phi(Y_t)].

For every fixed integer m>=2, iterate the centered Riesz identity m-1 times to obtain the EXACT finite current expansion

    d/dt E phi(Y_t)
      =sum_(j=2)^(m-1) [t^j/j!] kappa_(j+1):E D^(j+1)phi(Y_t)
           +t^m E[T_m:D^(m+1)phi(Y_t)].              (4)

The sum is empty at m=2. To isolate the rank-(m+1) cumulant as well, write T_m=c_m+(T_m-c_m) or take one more step. Every test in (4) is evaluated at the SAME actual positive Y_t. There is no Gaussian replacement of the test argument and no discarded sampler feedback.

If kappa_3,...,kappa_m all vanish, then (4) has only its final term. Integrating m derivatives in the independent Gaussian buffer and using Hermite isometry gives a velocity of L2 norm at most

    sqrt(m!) t^m sigma^(-m) L^m e.

Hence

    W2(Law(X+sigma Z),N(0,sigma^2I+Sigma))
       <=sqrt(m!)/[(m+1)sigma^m] L^m e.              (5)

At m=2 this is the cubic theorem; at m=3 it is the fourth-order theorem with kappa_3=0. An anisotropic buffer Q>=qI has the same bound with sigma^m replaced by q^(m/2), proved by covariance differentiation and buffer Hermite integration without commuting matrices improperly.

The fixed buffer, finite moments and integrable displayed velocity justify the Wasserstein dynamic bound; smoothing/Sobolev approximation gives the C1/Lipschitz case. The input dimension does not occur in (1)-(5).

## 4. What this settles about repeatability

The analytical connected-current expansion is genuinely available at arbitrary FIXED rank with explicit computable factorial constants and one physical energy. Thus lack of higher derivatives of f does not block THIS comparison step. Odd cancellation can remove selected constants in (4), but nonzero even cumulants remain and must also be realized/canceled by a true positive source.

Equations (4)-(5) do NOT give an executed coefficient field, a positive flow for an isolated odd/even current, a tensor oracle, or native gradient admission for the completed source. A repeatable finite algorithm still needs source-qualified positive packets for every retained constant/current and for all their mean/covariance/feedback descendants. Proper-cut, caller, owned-root, numerical and variance-budget bounds of those packets are separate requirements.

This is therefore an all-rank ANALYTICAL recurrence and a precise law-consumer target. It proves no all-rank finite original-VALUE compiler and no c(P)/P claim by itself.
