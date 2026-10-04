# A polylogarithmic scalar law refinement at the fluctuation scale

2026-10-04. New bounded-source extension of SCALAR-POSITIVE-HERMITE-MEAN-LAW-FACTORY.md. Independent derivation; no external manuscript used. The same missing vector and smooth first/paired-source ports remain.

## Statement

Let W be a complete scalar source, with a KNOWN deterministic/exposed anchor b0 and bound |W-b0|<=S. Let sigma>0 be a known target Gaussian standard deviation. For 0<delta<e^(-2), there is an explicit finite positive program with

    W2(output, N(E W,sigma^2)) <= C sigma delta,
    number of complete W records <= C (1+(S/sigma)^2)
                                  [log((e+S/sigma)/delta)]^3.

All records counted here are independent complete replays at the exposed original caller. Known scalar arithmetic has count O((1+(S/sigma)^2)L^3+L^2), with L the logarithm displayed above; it is polynomial-logarithmic when S/sigma is bounded by public logarithms. The program needs no original potential VALUE, expectation, or exact integration oracle. The prototype is in the same Gaussian/exact known scalar arithmetic model as the preceding factory; numerical threshold errors still receive a separately selected budget.

Thus a scalar Gaussian buffer comparable to the source's bounded fluctuation scale suffices for arbitrary law accuracy with only polylogarithmically many complete records. A buffer much smaller than that scale pays the explicit squared ratio (S/sigma)^2. There is no claim that this ratio is always harmless.

Unlike the fixed-K statement, the rank K below grows logarithmically with accuracy. This does NOT generate an exponential replay tree: all its K source products are made from one additive list of K independent batch averages. The literal count is displayed below.

## 1. Executed pilot, fresh averages, and known clipping

Write rho=S/sigma and

    L=log((e+rho)/delta), R=sqrt(64 L),
    epsilon=1/(16R), a=sigma epsilon,
    K=ceil(8L), B=ceil(8L),
    M=ceil(C0 (1+rho^2) R^2 L),

where one universal sufficiently large C0 is fixed in advance; C0=2^20 is more than sufficient for the elementary tail bounds used below.

Draw M complete pilot records and let b be their average. This b is an actually executed, known random anchor. It is not E W. Retain the complete pilot if desired. For each i=1,...,K, draw a new independent batch of M complete W records, form its average A_i, and execute

    W_i^clip = min(b+a, max(b-a, A_i)),
    X_i=(W_i^clip-b)/sigma.

Conditional on the pilot, the X_i are independent and satisfy |X_i|<=epsilon. Feed them, this actual anchor b, and sigma into the explicit positive Hermite factory, using cutoff R, rank K, and rejection cap B. The guard epsilon R+epsilon^2/2<log(5/4) holds.

The total complete provider count is exactly (K+1)M, before any optional source-specific anchor producer. It is O((1+rho^2)L^3). There is no product of M through K nested levels, and no resampling of W during Gaussian rejection proposals.

## 2. Pilot and clipping error

Let mu=E W. The pilot and each fresh average lie in [b0-S,b0+S]. Hoeffding's elementary bounded-sum estimate gives

    Pr(|b-mu|>a/4) <=2 exp[-M a^2/(32S^2)] <=2 exp(-8L)

when S>0 and C0 is as selected. For S=0 the mean is already known and one outputs its Gaussian translate exactly. The event G={|b-mu|<=a/4} is the good pilot event.

Conditional on any good pilot, clipping can occur only when |A_i-mu|>3a/4. Since both b and A_i lie in the original bounded interval, the clipping displacement is at most 2S. Therefore, writing mu_clip(b)=E[W_i^clip | pilot],

    |mu_clip(b)-mu|
       <=4 S exp[-9M a^2/(32S^2)]
       <=4 S exp(-8L).

This is an actual bias bound at the random pilot anchor. No unknown centering occurs in the program.

On the bad pilot event, |b-mu|<=2S always. The factory's standardized density is between 1/2 and 3/2 times the standard Gaussian, with a Gaussian fallback of the same moment scale. Thus its squared distance from an independent N(mu,sigma^2) output is bounded by C(S^2+sigma^2) in conditional expectation. Multiplying by the bad-event probability costs at most C sigma^2(1+rho^2)exp(-8L), which is O(sigma^2 delta^2) by the definition of L.

The entire pilot may be retained in this comparison: fresh batch means are integrated conditional on that retained pilot, and the comparison Gaussian has mean mu independent of the pilot. This does not permit retaining the K fresh batches that the Hermite averaging integrates out.

## 3. Growing-rank Hermite and tail bounds with explicit constants

The previous factory's finite-K C_K notation is not sufficient when K grows. Instead use bounds uniform in K.

For every batch and every z,

    |sum_(k=1)^K P_k H_k(z)/k!|
       <=exp(epsilon |z|+epsilon^2/2)-1.

The same bound applies to the marginal polynomial with P_k replaced by m^k, where m=(mu_clip(b)-b)/sigma and |m|<=epsilon. Hence, uniformly in K,

    ||1_(|Z|>R) S_m(Z)||_2 <= C exp(-R^2/8).

For example square the positive exponential bound, complete the square against the Gaussian density, and use R>=1 and epsilon<=1/(16R). This pays the entire growing polynomial tail; no C_K is hidden.

Hermite orthogonality gives the remainder bound

    ||exp(mZ-m^2/2)-sum_(k=0)^K m^k H_k(Z)/k!||_2
       <=exp(epsilon^2/2) epsilon^(K+1)/sqrt((K+1)!)
       <=C 16^(-K).

The explicit cutoff constant is again bounded by the L2 tail. Thus the density-ratio error is at most

    C[16^(-K)+exp(-R^2/8)] <= C exp(-8L).

The elementary maximal-common-density W2 conversion in the previous note gives conditional W2 <=C sigma exp(-4L) from the unclipped target N(mu_clip(b),sigma^2). The finite rejection cap contributes at most sqrt(5)sigma 3^(-B/2), also O(sigma exp(-4L)). These hold for every realized pilot, not just good pilots.

On the good event add the clipping mean bias. On the bad event use the moment bound from Section 2. A same-pilot coupling and then integration prove the stated C sigma delta bound.

## 4. Query, arithmetic, and numerical budget

The provider count is (K+1)M=O((1+rho^2)L^3). Complete private Gaussian dimension is the provider dimension times this literal count, plus at most B+1 Gaussian proposals and their uniforms. Every old root, changed caller, mode iteration and original VALUE/first site in a W record is paid inside its complete replay. Pilot records are not reused as fresh factors.

Compute prefix products P_k and the K Hermite recurrence values per proposal. This is O(KB) known scalar operations, in addition to the O(KM) additions for averages and clipping. The known factorial/polynomial intermediate scales have logarithmic size O(K log K + K log R), polynomial in L. At a selected absolute accuracy, smaller products may be rounded away; arbitrary exact real source values do not literally have this finite bit length. Known Gaussian densities, finite Hermite integrals and acceptance thresholds can then be assigned O(poly(L)) bits for those scales. Provider/anchor precision, Gaussian generation, and propagated rounding require separate certified budgets. This is not a completed bit-complexity or numerical caller port, and any exact-proposal or numerical-oracle convention must remain explicit in a complete integration theorem.

For rho bounded by public logarithms and delta=A^q at fixed q, the new complete original-query heat exponent is zero. More generally its explicit extra exponent is that of rho^2. This refinement does not solve an expensive provider's posterior bias; it only transforms its mean into a buffered law.

The hard clipping and accept/reject steps still do not supply a smooth first-action map. A smooth quantile implementation is a separate prospective repair, not included in this theorem. No general-D positivity theorem or complete hidden-source induction is asserted.

## 5. Scalar compatibility with the new endpoint Stein packet

For the packet H_Q(Z,G) of ENDPOINT-STEIN-PACKET-AND-POLYLOG-QUADRATURE.md in D=1, retain Z and evaluate the actual known anchor b(Z)=-H_Q(Z,0), paying its full node list. Replace each fresh G by its executable clipping to [-L0,L0] and set W=-H_Q(Z,clip G). Its global G-Lipschitz bound A gives

    |W-b(Z)|<=A L0,
    |E_G W + v_Q(Z)| <=A E[|G|1_(|G|>L0)],

uniformly in the retained Z. The latter bound is O(A exp(-L0^2/2)). Taking L0 proportional to sqrt(log(1/A)) makes this bias arbitrarily high order at fixed target. All conditional record independence is literal.

The present refinement therefore approximates the buffered conditional law

    Z-v_Q(Z)+sigma N

with no new inverse-heat count even when sigma is comparable to A: now rho=A L0/sigma is only logarithmic, and the pilot plus growing-rank construction repairs the unbatched small-epsilon guard. This is a valid retained-Z scalar law-only compatibility result.

The independent buffer variance sigma^2 remains REAL. For a scalar quadratic g(z)=Bz, v_Q(Z)=BZ/2, so the target of this compatibility step has exact variance

    (1-B/2)^2+sigma^2,

whereas the posterior variance is (1+B)^(-1). These differ at order B^2 and/or sigma^2 unless an additional paid correction is made. The construction does not silently remove that covariance debt, identify the nonlinear Stein target with the posterior, or close higher law currents.
