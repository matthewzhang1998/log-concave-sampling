# Independent audit of the bounded scalar Gaussian mean law factory

Audit date: 2026-10-04 UTC.

## Conclusion

PASS for the isolated scalar marginal-law lemma in the stated ideal Gaussian/uniform sampling and real-arithmetic model. No fatal gap was found in normalization, positivity, Hermite averaging, error control, or the capped rejection law. The construction needs exactly K fresh complete provider executions, an executable known anchor, and a certified conditional bound. It does not produce a strong mean statistic, a smooth first/adjoint port, a paired-caller contract, or a complete source recurrence.

The final source audited is `../SCALAR-POSITIVE-HERMITE-MEAN-LAW-FACTORY.md`, SHA-256:

    728e26f5c7de342c6cd9b4d10e8dde6e6514ee33afea20759649d2d75a90c65f

The source was re-read after its author incorporated the audit's initial scaling and proposal-count clarifications. This review used that source and independent calculations only. No external peer manuscript was read or used, and no external communication occurred. The source itself was not modified by this audit.

## Exact conditional normalization and positivity

For every batch x, C_x is precisely the integral of the cutoff polynomial against phi. Therefore integral q_x = 1 + C_x - C_x = 1, with no division and no unknown normalizer. The displayed boundary formula for I_k follows from (phi H_(k-1))' = -phi H_k and has the correct sign. Odd k yield I_k = 0.

The absolute-coefficient majorant is valid. Expanding the probabilists' Hermite polynomial and dropping signs gives

    sum_(k>=0) epsilon^k |H_k(z)|/k!
        <= exp(epsilon |z|) exp(epsilon^2/2).

Every batch satisfies |P_k| <= epsilon^k, so the guard gives |S_x| <= 1/4 inside the cutoff and |C_x| <= 1/4. Hence globally 1/2 <= r_x <= 3/2. There is no bad-batch exception. In particular the cutoff and correction create a positive conditional density before the fresh bank is averaged out.

The hypothesis should be read as 0 <= epsilon, sigma > 0, and R >= 1. The guard itself implies epsilon < 1, so the epsilon <= 1 assumption used in the tail estimate is automatic.

## Hermite expectation and error

All fresh records must have the same conditional mean and be jointly independent conditional on the exposed caller. Repeated independent executions of the same complete provider meet this condition. It is not sufficient to specify only pairwise independence for a general K. With the intended conditional independence,

    E[P_k | caller] = mu^k,  mu = E[X | caller].

This yields the asserted exact marginal density. The only use of mu is analytical; no execution step evaluates it.

Hermite orthogonality gives the stated exact squared remainder. Its bound follows from

    (K+1+j)! >= (K+1)! j!,

and hence the remainder norm is at most exp(epsilon^2/2) epsilon^(K+1)/sqrt((K+1)!). For fixed K and epsilon <= 1, |S_mu(z)| <= C_K epsilon (1+|z|^K). Repeated integration by parts gives

    integral_(|z|>R) (1+|z|^(2K)) phi(z) dz
        <= C_K (1+R^(2K)) exp(-R^2/2),  R >= 1.

Taking square roots proves the stated polynomial-tail bound. Since the full polynomial has zero integral, C_mu equals minus its removed-tail integral and is bounded by that tail's L2 norm. Thus the exact triangle decomposition is

    L2 profile error <= remainder norm + tail norm + |C_mu|
                     <= remainder norm + 2 tail norm.

This establishes the advertised T_K with constants depending on K. Arbitrarily high fixed order is an asymptotic statement in a shrinking epsilon, with R and B chosen accordingly. Merely increasing K while keeping epsilon and R fixed does not remove the cutoff error.

## Elementary W2 conversion

Let p and q be the two compared densities and h = min(p,q). Match mass h identically. If the residual mass is m > 0, couple (p-h)/m and (q-h)/m independently. The residual cost obeys

    m E|U-V|^2 <= 2 integral z^2 [(p-h)+(q-h)] dz
                = 2 integral z^2 |p-q| dz.

If m = 0, the claim is immediate. Applying Cauchy-Schwarz against phi then gives

    W2(p,q)^2 <= 2 sqrt(E_phi Z^4) ||(p-q)/phi||_L2(phi)
               = 2 sqrt(3) ||(p-q)/phi||_L2(phi).

The factor, use of the unnormalized residuals, and square-root loss in the final W2 bound are all correct. The result does not require an unproved high-order transport theorem. Translation by b and scaling by sigma multiply W2 by sigma.

## Rejection cap and provider cost

For every batch, the acceptance probability is exactly (2/3) integral q_x = 2/3. Fresh Gaussian/uniform proposal pairs therefore give conditional failure probability delta = 3^(-B), independent of the batch. Conditional on success at any attempt, the output law is q_x. Consequently the capped law, after averaging the bank, is exactly

    (1-delta) qbar + delta phi.

The proposed coupling with qbar has squared cost at most 5 delta because qbar has second moment at most 3/2 and phi has second moment 1. In fact, the independent phi variable has zero mean, so this specific coupling gives the stronger bound (5/2) delta. The source's constant 5 is a valid harmless relaxation.

The uncapped expected number of attempts is 3/2. With stopping upon acceptance, the actual capped number is at most B and its expectation is (3/2)(1-3^(-B)); an all-failed run then needs one fresh fallback Gaussian. A literal fixed-length implementation may generate unused dummy proposals. The final source now says "at most B" and records the fixed-tape option.

The original-query bill is K complete provider calls plus any nonfree anchor producer. No rejection attempt invokes the provider again. The O(KB) arithmetic estimate is valid as an operation count; it is not a bit-complexity estimate and does not erase Gaussian generation, high-precision arithmetic, or the provider's own internal cost. The word "complete" must continue to include all provider ancestors and repeated calls.

## Explicit asymptotic choices

Write L = log(1/A), with A tending to zero. Suppose epsilon <= C A^a L^u and sigma <= C A^(-s) L^v, where a > 0 and s >= 0 are fixed. Set q' = q+s. It is enough to choose fixed K and positive constants c_R, c_B so that

    a(K+1) > 2q',
    a + c_R/4 > 2q',
    c_B log(3)/2 > q',

and use R = sqrt(c_R L) and B = ceil(c_B L). Eventually R >= 1 and epsilon R + epsilon^2/2 <= log(5/4); the strict inequalities absorb all fixed logarithmic powers, including those from R^K and sigma. For a sigma bounded only by fixed logarithms, s = 0. Constants may depend on the fixed exponents and target order. This verifies the final source's corrected high-order statement.

If sigma(theta) is caller-dependent, the bound is pointwise conditional. A uniform joint or host error bound still requires the corresponding uniform sigma bound, an integrable squared conditional-error bound, or an appropriate host estimate. No unbounded random caller scale is free.

## Known anchor and retained information

The algorithm uses b, sigma, the realized W_i, and explicitly computable constants only. The center of the density is the supplied b, not E[W]. A random anchor can be part of the exposed caller, provided the provider hypotheses hold conditional on it and its production cost is charged. Reusing source randomness in the anchor does not automatically preserve the required conditional source law or independence.

The refusal to retain the private bank is essential. For a concrete counterexample, let the X_i be independent symmetric signs times epsilon. Then mu = 0 and the unretained output is exactly phi for every K. If X_1 is retained but the other records are averaged out, then

    q_(given X_1)(z) = phi(z)[1 + X_1 z 1_(|z|<=R)].

Its conditional mean is

    X_1 m_2(R),
    m_2(R) = E[Z^2 1_(|Z|<=R)]
           = 1 - 2[R phi(R) + Phi(-R)] > 0.

The conditional W2 distance from the declared zero-mean target is at least epsilon m_2(R), regardless of K. At R=2 and epsilon=0.1 this lower bound is approximately 0.073853587. Thus a high-order unretained law theorem cannot be reused as the same high-order theorem conditional on the hidden bank. An exposed caller supplied before the bank is a different, permissible conditioning variable.

Similarly, sharing one symmetric sign across X_1 and X_2 gives E[X_1 X_2] = epsilon^2 while E[X_1]E[X_2] = 0. This explains the complete independent-execution requirement without relying on a source-specific example.

## Unclosed gates and numerical implementation

The construction truthfully supplies no smooth first/adjoint source representation or declared two-caller coupling. The accept/reject program has discontinuous decisions, and the cutoff adds further nonsmoothness. A different representation might have useful regularity, but none is proved here. A marginal W2 statement alone does not supply those ports. No dimension-free vector or shared-coordinate-bank theorem follows from the scalar calculation.

Likewise, a certified source-specific boundedness or clipping argument around the executable anchor remains a prerequisite. Clipping error in the provider mean would contribute directly to the Gaussian target error and must be charged separately. No posterior bias theorem or eventual-sublinear recurrence is established by this factory.

The source explicitly reserves finite-precision error and moment/profile costs. This is an honest scope restriction. Exact normalization belongs to the ideal real-arithmetic density; finite-precision acceptance probabilities do not preserve exact 2/3 success or the exact mixture identity automatically.

A small supplemental observation shows that moment control is available for a uniformly controlled threshold approximation, though it does not close all implementation details. Suppose exact Gaussian/uniform proposals are retained and the implemented acceptance threshold a_tilde_x(z), clipped to [0,1], differs uniformly from a_x(z) = (2/3)r_x(z) by at most eta <= 1/6. Its total success probability is at least 1/2, so every accepted-proposal density is bounded by 2 phi; mixtures with the Gaussian fallback have second moment at most 2. Couple the two B-capped programs using the same proposals and uniforms. At the first differing decision j, the mismatch probability conditional on its proposal is at most eta. One program outputs that proposal, while the other's remaining output has conditional second moment at most 2. Summing over j gives

    TV(exact capped, approximate capped) <= B eta,
    W2(exact capped, approximate capped)^2 <= 6 B eta.

This observation is conditional on the batch and survives averaging it; physical scaling adds sigma^2. It still requires a certified uniform threshold error, and it does not account for approximating Gaussian generation, source records, or the anchor. It grants no smooth first-action port.

## Independent diagnostic results

`check_scalar_factory.py` is a deterministic numerical sanity check. It enumerates all sign vertices of bounded batches for K in {1,2,4,6,8} and R in {1,2,4,6}, with epsilon at 97 percent of the guard boundary. For a fixed z the batch ratio is multi-affine in the X_i, so batch extrema are attained at sign vertices; for each vertex the script checks the cutoff endpoints and real derivative roots, plus the constant exterior ratio. Floating-point roots make this a diagnostic rather than a certified interval proof.

Across 20 parameter cases:

- Maximum numerical normalization discrepancy: 4.44e-16.
- Maximum discrepancy between weighted complete-bank averaging and the mu polynomial: 1.11e-15.
- Observed global ratio range: [0.75996965, 1.24003035].
- The profile triangle bound, weighted second-moment Cauchy-Schwarz bound, and second-moment envelope passed in every case.
- The retained-record and shared-ancestor counterexamples agree with their exact formulas.

The full results are in `diagnostics.json`, and `diagnostics.out` contains the concise summary. The analytic proof above establishes the general result; numerical checks do not replace it.
