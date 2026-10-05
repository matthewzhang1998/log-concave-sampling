# A simultaneous positive family join, with rank-eight grade 57/7

2026-10-05. This is a one-stage complete mixed-bank estimate, with all higher old/new mixing bounded at once. It is not a LAW-to-RAW re-entry theorem. No unproved multi-copy heat estimate is used in the rank-eight application.

## 1. Exact hypotheses and complete old-first certificate

Let B be the complete standard-Gaussian coefficient innovation bank, independent of the retained caller Y and all other primitive tapes. Original repeated bank rows remain repeated. Every root shared with a coefficient caller is included in B with its actual derivative paths. Captured query means may depend on fixed Y.

Let W(B,U,Y) be the old additive endpoint and all retained observer coordinates. Its source-zero B row is zero. The FULL stacked B-Jacobian has the certificate

    sup_(B,Y) ||D_B W||_(L4(U);op) <= L alpha^s,   s>0.      (1)

A uniform L2(U) bound suffices for the old-only term below when the new reference row is independent of U; L4 is retained to match the supplied ports. The certificate is a direct first of the literal old finite program. A law-error derivative is not a substitute. A merely joint-average L4 bound in B,Y is also not a substitute for (1) under the stated L2 coefficient hypothesis.

All new groups own their dynamic roots. They use positive covariance/readout shares fixed in advance and a single untouched physical keep `sqrt(kappa) G`, kappa>0. The keep is independent of every coefficient and readout field. Retained B-independent observers contribute zero to (1); an observer reading B must be included in the full stacked derivative. An unattenuated copy of B fails s>0.

After the frozen-coefficient native comparisons, let the sum of the new main reference fields be

    F0(B,H,Y)=sum_i F0_i(B,H,Y),
    barF(H,Y)=E_B F0(B,H,Y),
    F=F0-barF.

H is the complete new reference Gaussian row or collection of rows, with fixed B-independent geometry, independent of U and B. Each F0_i is the actual finite Hermite transport of its coefficient, not a tensor oracle executed by the algorithm. When an aggregate tensor heat estimate is used, all its terms must be recoupled to the SAME reference H with the same inverse-covariance Wick convention before that estimate is applied. Separate unrecoupled H_i references do not support cancellation in an aggregate tensor estimate. Its total physical polynomial degree is at most m. Assume

    ||F||_(L2(B,H,Y)) <= E sqrt(D),
    sup_(B,Y) ||D_B F||_(L4(H);op) <= J.                    (2)

The same bounds at fixed Y with a uniform coefficient energy are sufficient; an integrated-Y energy also suffices when the Jacobian constants are uniform in Y. All physical inverse-variance/Hermite factors belong in E and J. The original currents may have different physical ranks; use their complete summed field, or sum its actual one-Hilbert energy bounds before applying the lemma.

## 2. A single exact identity sees every mixed current

Consider the common positive path

    X_t = W(B,U,Y) + H_physical + barF(H,Y)
                       + t F(B,H,Y) + sqrt(kappa)G,         (3)

with the B-independent source-zero Gaussian rows included once. Embed F only in the physical-output coordinates if observers are stacked. Every X_t is a real positive pushforward.

Let `R_B=D_B N_B^(-1)` be the Gaussian Riesz operator on centered B-functions. Conditional on H,Y,U,

    E[F . grad phi(X_t)]
      = E[(R_B F)(D_B W+t D_B F)^T : Hess phi(X_t)].         (4)

This identity is exact, not a second-cumulant truncation. In particular all higher mixed cumulants between old and new packets, and between distinct new packets sharing B, are already included in its left and right sides. No equality conditional on a retained B is asserted.

Set `A_t=(R_B F)(D_B W+t D_B F)^T`. One index of A_t is a physical output direction. Transfer that physical test derivative through the independent keep. The resulting vector/current is `A_t^T G/sqrt(kappa)` in the full stacked endpoint coordinates, with the natural restriction of G to the physical index. Thus its L2 norm is bounded by

    kappa^(-1/2) ||A_t||_(L2;HS).                           (5)

The keep is independent of A_t, so (5) uses exactly one Hilbert factor.

For each fixed B,Y, the coefficient of R_B F remains a Gaussian polynomial in H of degree at most m. Hypercontractivity gives

    ||R_B F||_(L4(H);HS) <= 3^(m/2)
                             ||R_B F||_(L2(H);HS).

Riesz is an L2 contraction in B. Conditional Hölder in H,U and the uniform Jacobian bounds therefore imply

    ||A_t||_(L2;HS)
       <= 3^(m/2) E [L alpha^s+t J] sqrt(D).                (6)

For the old part one may remove the hypercontractive factor by independence of H and U, but the displayed common factor is safe. Integrate the L2 continuity velocity, or use the associated same-endpoint transport-current inequality, to obtain

    W2(Law(X_1),Law(X_0))
       <= C_(m,kappa) E [L alpha^s+J/2] sqrt(D),
    C_(m,kappa)=3^(m/2)/sqrt(kappa).                        (7)

The same proof applies to the appropriate physical conditional comparison with retained B-independent variables. Equation (7) is the full bank-mean join. It is stronger than multiplying a coefficient first-hit envelope by the old first: its old term requires coefficient **energy**, not one extra original heat derivative.

## 3. Aggregate heat, own return, and floors

Suppose a complete target family has nominal main grade a, total original inverse-width degree at most K, and one common `tau=alpha^gamma`. Assume its finite original-VALUE source generator and positive clock rule certify

    E <= C_E alpha^a tau^(-K),
    J <= C_J alpha^a tau^(-(K+1)),
    target-heat error <= C_H alpha^a tau sqrt(D).            (8)

The last inequality is required for the **complete** target family. An aggregate cumulant heat bound must not be allocated independently to selected histories without proof.

At each graph with N original forces, put all known coefficient factors `K_h`, including clocks, signs, original widths and inverse response/readout normalizations, at the root. Choose a fixed rho exponent beta at each nonroot, and

    rho_root,h = K_h alpha^(a-beta(N-1)),
    |K_h| <= Lambda_h tau^(-K).

The frozen-coefficient native port gives

    own <= C_own alpha^(2a-2beta(N-1)) tau^(-2K) sqrt(D),    (9)
    C_own = Gamma_mix [sum_h Gamma_h |tau^K K_h|]^2,
      <= Gamma_mix N_groups sum_h Gamma_h^2 |tau^K K_h|^2.

Here Gamma_h bounds the actual packet energy/first/main-regression constants, and Gamma_mix is the complete joint current/keep-transfer constant. This finite square explicitly includes cross-group branches h != g. It includes all actual node/readout/gap/moment constants and cannot be replaced by a sum of uncoupled individual-law errors or an independent variance saving. The supplied common-H native path first sums every main tensor at one positive endpoint; every unmatched intrinsic branch then has two possibly different root factors. This is the same joint extension now stated in the rank-eight input, Section 4.

Combining the positive conditional native comparison, (7), ordinary same-root coupling for target heat/quadrature replacement, and all propagated numerical floors gives

    error / sqrt(D) <=
       C_H alpha^(a+gamma)
      +C_(m,kappa) C_E L alpha^(a+s-gamma K)
      +(C_(m,kappa)/2) C_E C_J alpha^(2a-gamma(2K+1))
      +C_own alpha^(2a-2beta(N-1)-2gamma K)
      +complete propagated finite floors / sqrt(D).        (10)

This is a quantitative multitype one-stage recurrence. The four grades respectively belong to heat restoration, all old/new mixing, self/new-new bank mixing, and intrinsic native own return. All original first/curl/Sobolev obligations stay attached to the literal programs.

Balance the first two terms at

    gamma = s/(K+1),      P = a+s/(K+1).                    (11)

A sufficient admission test is

    2a-gamma(2K+1) >= P,
    2a-2beta(N-1)-2gamma K >= P,
    a-beta(N-1)-gamma(K+1) >= s,                            (12)

plus all real source/first/gap/readout/floor inequalities. The last line is the conservative root-caller first. All other new paths contain its strict root ancestor and their actual positive nonroot amplitudes, but every path and injection row is still summed. Strict exponent margins absorb fixed public logarithms when necessary.

For a finite collection with different (a_i,K_i), use each actual E_i,J_i,own_i and aggregate E=sum E_i, J=sum J_i in (7). Cross terms are included by the product. Using different floors for different source families additionally requires their exact common-bank construction and each family's complete-target heat estimate; it is not granted by (10) alone.

## 4. The complete rank-eight family now reaches 57/7

Use the full rank-eight cumulant generator and the ownership/old-first/retention setting of

`/workspace/shared/induction-native-rank-generator-20261005/RANK8-SMOOTHED-NATIVE-RETURN.md`.

No histories are removed: it has 794,880 generated histories, exact integer weight sum `S_8=32,049,561,600`, inverse-heat constant `B_8=46,080`, fifteen original clocks per history, and 27,020,800 raw original VALUE calls over all histories at one clock tuple per history. Original repeated clocks and bank roots remain common. The existing same-caller positive dyadic rule and all its actual nodes are retained at the changed tolerance.

Here a=N=8, K=6, m=7 and s=1. Choose

    tau=alpha^(1/7),
    rho_nonroot=alpha^(1/14),
    rho_root,h=K_h alpha^(15/2).                            (13)

The four error grades and root-caller grade are

    target heat       : 8+1/7                  = 57/7,
    whole old mixing  : 9-6/7                  = 57/7,
    new-bank self     : 16-13/7                = 99/7,
    conditional own   : 15-12/7                = 93/7,
    root caller first : 15/2-7/7               = 13/2.       (14)

The actual root envelope is alpha^(93/14), and its square gives 93/7. Thus the complete rank-eight positive original-VALUE return reaches grade **57/7**, improving the conservative 65/8 balance. Their difference is exactly 1/56. Nothing about the old full B-first improves: it is still the independently supplied O(alpha) certificate.

Explicit conservative coefficient constants are

    C_E = H_7 2^15 S_8 B_8 * actual current normalization,
    H_7 = sqrt(7!) v_0^(-7/2),
    C_H = 8! * actual Hermite/current transfer,
    D_8 = 15 * 2^15 S_8 B_8.

For C_J include the complete finite matrix-frame/readout constants, actual affine bank injections, and the safe one-hit scalar envelope

    2^15 * 2 * 8 * S_8 * 2^7 * 7!.

The actual physical reference variance lower bound v_0 and keep kappa are separately reserved. If a growing finite node census makes inverse readout shares larger, insert those computed factors in C_E,C_J,C_own and verify that they are only the declared public logarithms.

The positive clock allowance is `D_8 delta_0 alpha^8 tau^-6 sqrt(D)`. A sufficient changed clock multiplier tolerance is now `delta_0<=c alpha`; the old `alpha^(7/8)` tolerance does not automatically pass. The source-response telescope likewise needs local response scale `delta<=c alpha`, with the actual normalization constants divided out. At the minimum width, the original VALUE floor has scale `c A alpha^(8/7)`, or smaller if the complete native/Sobolev floors demand it. Every affected old capture, source origin, same-center anchor and ancestor is replayed at its required version.

A sufficient starting native order is `b>14(P+J_amp)`, followed by evaluating the complete prior and all actual floors. Here `J_amp` is the computed downstream inverse-alpha amplification; it is not the coefficient Jacobian constant C_J. Changed tolerances and widths require a freshly frozen program.

No numerical terminal-curl bound follows merely from response precision. The canonical literal-gradient terminal baseline, the new actual first/path sum, and all numerical Sobolev floors remain required. No individually consumed physical public may be retained externally.

## 5. All primitive ranks n>=3 have the same one-stage improvement

Where the supplied all-rank generator and fixed-graph native ports have their stated complete finite constants, set

    a=N=n, K=n-2, m=n-1, s=1,
    beta=1/[2(n-1)], gamma=1/(n-1).

Then the complete **single rank-n family** has the prospective admitted grade

    P_n=n+1/(n-1),
    own=2n-3+2/(n-1),
    self=2n-2+1/(n-1),
    root caller=n-3/2.                                     (15)

For n>=3 both own and self are strictly above P_n and root caller is strictly above one. This is a fixed-rank source-return theorem under the same full old-first and reserve/precision assumptions. It does not say that correcting rank n repairs every smaller old current or that the law residual can be re-entered.

The raw source count at one clock tuple is V_n from the supplied exact generator recurrence. Its unmerged total is

    sum_(history,node) sum_v 2^(k_v+1)
       M_native(k_v,b_v,epsilon_v)
       * Q_actual_original-caller-replay(v)
    + all captures, readouts and versioned replays.          (16)

The exact max-plus inverse-alpha cost rule is in `WEIGHTED-MIXED-QUEUE.md`, Section 9. Its audited LOW30 serial recurrence certifies local `q_native=0` for these eligible fixed-order sources with fixed-power/public-log floors. The original positive resolvent rule likewise has exponent zero. Hence the new fixed-rank stage adds no inverse-alpha query exponent beyond the actually requested old complete versions and caller/ancestor replays. If those also use the admitted fixed-order public-log compiler, the complete fixed-rank port has query exponent zero relative to unit-cost original VALUE calls. In particular (15) by itself establishes neither `c(P)=o(P)` nor a globally rank-n target sampler.

## 6. What the result says and what it does not

The rank-eight conclusion consumes a complete generated family, its exact smoothing mismatch, the actual full old-first defect, every new/old mixed effect through an exact nonperturbative identity, and the native root-square remainder at a common positive endpoint. This is more than a conditional feedback calculation.

Its endpoint is nevertheless a positive reference carrying a specified mean cumulant field, with the rest of the old ledger unchanged. To make it an arbitrary-order recursion one still needs source-qualified repair of the residuals (especially heat restoration), preservation of all terminal/source ports, and complete native/clock/replay query exponents. The weighted heated queue identifies a sufficient potential; it does not supply those missing ports by renaming an error.
