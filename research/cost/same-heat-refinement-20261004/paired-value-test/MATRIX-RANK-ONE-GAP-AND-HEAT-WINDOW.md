# Matrix test of the scalar paired covariance rule

2026-10-04. This screens the new empirical rank-one covariance fill, not native Jacobian/Gram reserves or all posterior samplers.

## Verdict

The direct matrix analogue does **not** supply a dimension-safe no-copy source. It normalizes one sampled outer product, whose size is controlled by the whole signal energy, rather than by the covariance operator norm. In the relevant early-heat/high-dimensional regime, capping loses the covariance target or the needed replica count grows with dimension and heat.

The scalar positive covariance companion remains valid in its stated endpoint scope. Its success does not change the admitted general c_P family.

## 1. Exact dimension test

Conditional on the common record, take the admissible affine vector source

    H=ell G_d,   H'=ell G_d',
    Cov(H)=ell^2 I_d,

with independent complete standards. Its full source first is ell and its centered energy is ell sqrt(d). The proposed matrix correction is

    S=(H-H')(H-H')*/2 = ell^2 U U*,
    U=(G_d-G_d')/sqrt(2) standard.

Thus

    E S=ell^2 I_d,
    ||S||_op=ell^2 |U|^2.

The covariance operator is small independently of d, but the executed sample matrix typically has size d ell^2. The scalar keep split allocates variance v_keep=sigma^2/2 to the corrected branch and preserves another sigma^2/2 as an untouched keep. A strict positive square-root gap for this matrix branch needs ||S||<v_keep, up to its stated fixed margin.

For a Gaussian source no finite cap-free gap is uniform on all records. If d ell^2 is much smaller than v_keep, one may cap a rare tail with its actual VALUE/mean/first price. If d ell^2 is larger than v_keep, the cap acts on typical records and cannot be called a numerical tail.

## 2. Trace capacity: why a cap cannot preserve the target for free

Let T be any random positive semidefinite correction satisfying

    rank(T)<=m,    0<=T<=v_keep I.

Then pointwise tr(T)<=m v_keep. Hence an unbiased covariance target E T=Cov(H) requires

    m v_keep >= tr Cov(H)=d ell^2.                     (1)

For one sampled outer product m=1. If d ell^2>v_keep, no rank-one positive cap meeting the gap can retain E T=Cov(H), regardless of its thresholding rule. This is a trace-capacity identity, not only a concentration estimate.

A randomized number of such directions obeys the corresponding expected-work requirement E m>=d ell^2/v_keep. A finite fixed-rank correction depending only on the desired order cannot hide this dimension-dependent demand in an order constant.

Using m independent paired samples in an empirical covariance matrix can raise its rank, but needs at least the count (1) before concentration and a strict gap are considered. That necessary count is not by itself sufficient. Every complete H call and its actual first/tape must be charged.

For an isotropic radial cap, the missing trace appears as a positive residual covariance. More generally, it must remain as an explicit covariance/current error or be corrected by a different source theorem. It cannot be discarded after matching a small operator covariance in expectation.

## 3. R=4 scalar schedule versus the early heat window

The scalar two-layer construction uses

    sigma=a^(5/4),
    N_final=a^(-2/3),
    ell approximately a N_final^(-3/2)=a^2,
    v_keep approximately a^(5/2).

Here the ell formula is deliberately optimistic for the matrix screen: it grants a dimension-safe covariance envelope. The actual shared-clock source can have worse dimension coherence; that does not improve the following requirement.

Equation (1) gives

    m >= c d a^(3/2).                                (2)

Thus a fixed-rank/no-copy version is confined to d a^(3/2) bounded, and a strict-gap/tail argument needs additional slack. That condition is not available in the general early-heat window.

There are two visible ways to pay the deficit within this particular rule:

1. Keep N_final=a^(-2/3), but use at least order d a^(3/2) complete covariance samples before further tail/gap overhead.
2. Increase the quadrature count so that a fixed number of samples has enough trace capacity. With ell=a N_final^(-3/2), this requires

       N_final^3 >= c d a^(-1/2),
       N_final >= c d^(1/3) a^(-1/6).                 (3)

At the usual illustrative accuracy balance sqrt(d) a^R comparable to one, d is of order a^(-2R). At R=4, (2) is m at least order a^(-13/2), and (3) is N_final at least order a^(-17/6). These are far above the scalar counts. This substitution is an explicit accuracy-balance illustration, not an assertion that d and a must satisfy that relation in every call.

The same calculation for a target R with the scalar keep choice sigma^2=a^(R-3/2) gives

    m N_final^3 >= c d a^(7/2-R).                     (4)

At d of order a^(-2R), a fixed-rank rule requires N_final at least order a^[-(R-7/6)]. Its heat exponent is linear in R. This is a limitation of the proposed empirical rank-one normalization certificate, not a universal complexity lower bound.

## 4. What existing matrix source machinery does differently

The original rectangular covariance constructions use an actual bounded Jacobian/adjoint action and a stationary Gaussian completion with an operator gap. Their full source, root, orientation and first contracts are substantive. They are not equivalent to replacing the matrix covariance by one rank-one output outer product.

The scalar rule cannot inherit those dimension-safe interfaces merely because E S has the right covariance. To replace (1), one needs a matched full-rank/operator-gap source construction, with its actual original-query count and all retained/caller currents. Invoking that construction is precisely returning to the native covariance-closure gate rather than bypassing it.

## 5. A useful exact VALUE calibration, with its scope

There is a separate way to remove the simplest quadratic random-clock fixture without singular importance weights. For the first Picard force integral,

    I=integral_0^(pi/2) cos(s)b(cos(s)z0+kappa sin(s)Z) ds,

the exact identity is

    I=(pi/4)b(z0)+(1/2)b(kappa Z)
      +integral_0^(pi/2) cos(s) E_s ds,
    E_s=b(cos(s)z0+kappa sin(s)Z)
          -cos(s)b(z0)-sin(s)b(kappa Z).               (5)

All three VALUE sites use the same actual path records. For every linear b, E_s vanishes pointwise, so the integral is deterministic in its quadrature clocks. The two anchor calls may be cached. There is no inverse-sin clock factor, and its direct first/caller is bounded using only Db.

Under general C2, E_s still has order-a energy and first, with its complete shared-root and higher law currents. Identity (5) does not establish a stronger posterior grade, a general dimension-safe covariance fill, or a completed two-layer feedback repair. It is a concrete better-calibrated source on which the missing coupled theorem could be tested.
