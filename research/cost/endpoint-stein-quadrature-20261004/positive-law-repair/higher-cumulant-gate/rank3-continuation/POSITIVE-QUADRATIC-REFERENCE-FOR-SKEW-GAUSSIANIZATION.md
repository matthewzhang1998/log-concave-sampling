# Positive quadratic reference for skew Gaussianization

2026-10-04. Analytical comparison, pending independent audit. The tensor and covariance in this note are references only; the companion native path supplies their VALUE-based leading-current realization.

Let X=f(V)-Ef(V), Sigma=Cov(X), e=||X||2, and Lip(f)<=L<=1. Let kappa=E X^(tensor 3), K=kappa/6. The earlier centered-Riesz calculation gives

    ||K||HS<=L^2 e/3=:h,
    every proper cut of K <=L^3/3=:k.

Fix sigma>0, eta=sigma/sqrt(2). Put

    C_t=eta^2 I+(1-t^2)Sigma, S_t=C_t^(1/2),
    b_t(z)=K:[(S_t^(-1)z) tensor (S_t^(-1)z)-C_t^(-1)],
    a_t=1-t^3,
    W_t=tX+S_t Z1+eta Z2+a_t b_t(Z1),                  (1)

where V,Z1,Z2 are independent standard Gaussian blocks. Every W_t is a genuine positive probability law and eta Z2 remains an untouched independent buffer. Its endpoints are

    W_1=X+sigma Z,
    W_0=S_0 Z1+b_0(Z1)+eta Z2.                        (2)

Consequently a positive Gaussian plus a quadratic correction realizes the original skew source through order four:

    W2(Law(W_1),Law(W_0))
       <=C_sigma [L^3 e+L^2 h+Lambda k h+Lambda^2 k^2 h]
       <=C_sigma Lambda^2 L^3 e.                     (3)

Here Lambda is a fixed Gaussian matrix-series logarithm in the dimension; sharper finite Wick contractions can remove it from this particular analytical lemma. The bound is one energy. The standard public-log small-radius guard used by finite source compilers is more than sufficient. No tensor or covariance in (1) is a producer query.

## Exact same-endpoint proof, including feedback

Use B=tau-Sigma and the symmetrized second Stein tensor T from the companion fourth-Gaussianization note. Then E T=kappa/2=3K. Define Q from R(T-E T) multiplied by Df; ||Q||_(L2 HS)<=L^3e.

Differentiate E phi(W_t). Gaussian covariance differentiation against Z1, now including the actual b_t dependence, and then two centered V-Riesz steps give

    d/dt E phi(W_t)
      =3t^2 K:E D^3phi(W_t)+t^3 E[Q:D^4phi(W_t)]
       -3t^2 E[b_t dot grad phi(W_t)]
       +a_t E[dot b_t dot grad phi(W_t)]
       -t a_t E[F2_t:D^2phi(W_t)],                    (4)

where

    F2_t=Sigma S_t^(-1)(D_z b_t)^T.

The last term is the covariance-path sampler feedback. Omitting it would be an incorrect use of the unperturbed covariance-preserving path.

Two Gaussian integrations in Z1 give

    E[b_t dot grad phi(W_t)]
      =E[K(H_t,H_t):D^3phi(W_t)]
       +a_t E[C2_t:D^2phi(W_t)],                      (5)

with

    H_t=I+a_t(D_z b_t)S_t^(-1),
    (C2_t)_ip=sum_(a,b,j,k) K_iab
          (S_t^(-1))_aj (S_t^(-1))_bk
          partial_j partial_k(b_t)_p.

The leading K term of (5) cancels exactly the first term of (4). Thus the remaining exact currents have ranks one through four and coefficients

    rank 1: a_t dot b_t,
    rank 2: -t a_t F2_t-3t^2 a_t C2_t,
    rank 3: -3t^2 [K(H_t,H_t)-K],
    rank 4: t^3 Q.                                    (6)

They are all independent of Z2, and all tests are evaluated at the same W_t.

Since Sigma<=L^2 I and S_t>=eta I, Gaussian Hermite isometry gives

    ||dot b_t||2+||F2_t||_(L2 HS)<=C_sigma L^2 h,
    ||C2_t||HS<=C_sigma k h.

The second line is the Hilbert norm of a Gram-type contraction of two K tensors, bounded by one proper cut times one Hilbert norm. D_z b_t is a Gaussian matrix series with both variance operators bounded by C_sigma k^2. Gaussian matrix-series moments therefore give

    ||K(H_t,H_t)-K||_(L2 HS)
       <=C_sigma[Lambda k h+Lambda^2 k^2 h].          (7)

No derivative of f beyond its bounded first has entered this argument. The b_t derivatives concern a fixed explicit quadratic analytical field. Finally integrate rank-minus-one derivatives of (6) in the untouched eta Z2 buffer, using the corresponding Hermite isometries, and use the Wasserstein dynamic length bound. This proves (3).

## Joining the native packet

The native three-marked-path feedback lemma compares its finite positive output with a positive quadratic reference having its analytically identified third tensor. Independent packet references can be combined by conditional Gaussian regression onto their total visible Gaussian carrier, with the private-regression feedback paid quadratically in their coefficient scales.

Mean and covariance Gaussian services may supply the remaining mean and covariance, on independent complete banks and with a fixed untouched reserve. If their errors and the native coefficient calibration are at most the right side of (3), the combined positive output approximates the original buffered skew law at that order. The different Gaussian-coordinate choices for the quadratic reference are reconciled by the same private-Gaussian regression lemma; a third-cumulant match alone is not the justification.

This closes the analytical law consumer once the finite services and their precise source-qualified clock/first/caller/floor returns have been verified. It is not an assertion of arbitrary-order mean/covariance services or of native admission for the completed cubic output.
