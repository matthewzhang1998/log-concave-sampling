# A positive shared-passive skew fork with a retained stationary baseline

## 0. Status and exact input

This is a conditional construction for an ADMITTED common-root heat word. Its cross covariance has the true-orientation ordering automatically. It uses whole stationary channels and keeps an unmarked nonlinear side in the retained graph. It is not an energy-kappa-e counterpacket, and it does not yet identify the heat word with the actual nonlinear native covariance defect.

The finite pair/chain source is `astra-raw-interface/true-orientation/STATIONARY-MARKED-PAIR-SKEW-CHANNEL.md`, SHA-256 `ace5e761c38177388cf809cf27ea4c30e8f6c943a1b92753843e4bf11f2ff7bf`. Its independent R audit is `INDEPENDENT-STATIONARY-MARKED-SKEW-AUDIT.md`, SHA-256 `40fe1d02fa503f4b68a3e48c3f71c61a5f1299cdc93f46db47cb25a3ea915e74`.

Fix theta and captured coefficient records R. The new passive p is standard conditional on R. The intended insertion uses a NEW p independent of the old consumer conditional on R. All source versions are fixed. Known scalar coefficients, widths, row maps and variance allocations may depend on theta and fixed clock indices, but not on p or the owned Gaussian R records. Any theta derivatives of these choices must be charged separately; the displayed caller rows hold with them fixed.

Use the complete marked pair Q_E(p), with conditional reference

    Q_E^*(p) = M_R p + sqrt(I-M_R M_R*) Z_E^*,
    M_R = c_sel J_R.

Its known source-zero carrier Z_E is standard and independent of (p,R). The ACTUAL residual Q_E-Z_E has Lp energy at most Lambda_p e_p and declared residual first/caller. Its conditional-p comparison error is delta_E(R). The noise in its analytical reference is not identified with the actual Z_E.

Independently of that pair, conditional on (p,R), run the WHOLE two-pair skew chain Q_C(p) from the source above. Put

    D_R = diag(H_t(R),H_0(R)),
    L = [P*,C*],      J_0 = [[0,I],[-I,0]],
    W_R = L D_R J_0 D_R L*,
    C_R = tau W_R,
    K_R = c_sel^2 rho^2 W_R / 2.

Here H_t and H_0 are the SAME specified primitive heat Hessians across both chain banks, tau is the actual known signed word coefficient, and |tau| <= kappa. Thus W_R and K_R are skew. No division by tau is executed, including when tau=0. The whole chain has joint conditional-p Gaussian reference with selected coefficient K_R, marginal standard at fixed R, and error

    delta_C(R) <= epsilon_outer(R) + Lambda rho epsilon_inner(R).

The two finite errors are the actual high-order gradient-pair priors/restoration errors. Known Gaussian completions and all primitive replays remain part of each channel.

## 1. Literal positive program and exact conditional kernel

Choose b>0 and set the KNOWN signed coefficient

    c = -q tau / (b c_sel^3 rho^2),

where q is the desired covariance coefficient. Require a fixed positive reserve

    v - b^2 - c^2 >= v_keep > 0.

Execute

    Y = b Q_E(p) + c Q_C(p) + sqrt(v-b^2-c^2) Z_keep.       (1)

Z_keep is fresh, independent and unobserved by both channel comparisons. All marked-pair private records are independent of all skew-chain private records conditional on (p,R). The two skew pair banks share the SAME R and source version, but have independent private tapes as specified in the source.

At the reference level, conditional on p,R,

    mean(Y^* | p,R) = (b M_R + c K_R) p,
    Cov(Y^* | p,R) = v I - b^2 M_R M_R* - c^2 K_R K_R*.    (2)

After integrating the new p at fixed R,

    Cov(Y^* | R) = v I + bc(M_R K_R* + K_R M_R*)
                = v I - q O_R,
    O_R = -Sym(J_R C_R).                                (3)

Indeed M K*+K M* = -c_sel^3 rho^2 Sym(J_R W_R), and the chosen bc gives (3). This is an exact noncommuting matrix identity. It does not replace J_R* by J_R. No two-curl term is discarded. The term J_R C_R has the orientation required in the supplied target O_R.

Positivity follows directly from (1) and its independent keep. In particular the conditional reference covariance is at least v_keep I. No unknown matrix square root is executed.

If p is observed outside this module, retain the FULL joint kernel (2). Formula (3) alone is then insufficient. Using a new private p avoids that observer problem without changing the common coefficient roots R.

## 2. Finite conditional comparison and endpoint chronology

At fixed R, couple the two actual channels to their conditional-p references with independent coupling innovations given p. The triangle inequality gives the safe joint-p error

    |b| delta_E(R) + |c| delta_C(R).                       (4)

There is no centering assumption and no unproved root-mean-square improvement. The skew chain's conditional-p comparison follows by first comparing its inner pair, transporting that error through the outer incoming first O(rho), and then comparing the outer pair at the resulting standard-marginal intermediate passive. The conditional pair kernel permits retaining p through this step; there is no extra observer of the replaced private bank.

All old consumer records may be captured during this comparison EXCEPT the independent keep used by any subsequent buffer argument. Independence conditional on R and the new p makes (4) valid with those old spectators. Only after both channels have their joint reference may the new p be integrated.

To cancel an old +q O_R covariance field, the old and new references must carry the IDENTICAL coefficient records R through the covariance comparison. A separately averaged O_R cannot be substituted. If the sum of the two conditional covariance fields cancels exactly at fixed R, no O_R root mixture remains from that canceled field. If R is later owned while an uncanceled coefficient remains, its same-endpoint root/Price currents must be supplied. This note does not infer such a current from a coefficient HS bound.

For the admitted near-gradient marker with first A and kappa >= A^2, one available bill is

    delta_E <= Lambda [kappa + A^(7/3)] e + epsilon_E.

Thus (4) costs at most

    Lambda |b| [kappa + A^(7/3)] e
      + |b| epsilon_E + |c| delta_C.                       (5)

The pair priors are numerical, not automatically relative to a vanishing actual e. Set their absolute precision and source-version restoration before execution, including the readout |c|.

## 3. The actual retained baseline and its energy

The useful retained decomposition of the ACTUAL program is

    Y = B_ret + b(Q_E-Z_E),
    B_ret = b Z_E + c Q_C(p) + sqrt(v-b^2-c^2) Z_keep.       (6)

Only the marked residual is deleted in (6). Its actual deletion energy is |b| e_p. Z_E is independent of the entire skew chain and of its p,R; it is the literal source-zero carrier, not an analytical coupled noise.

Replacing Q_C by its stationary reference shows that B_ret has standard marginal variance v at each fixed R, with error |c| delta_C(R). Its conditional-p mean is c K_R p. Consequently (6) is a retained STATIONARY-channel baseline, not a bare Gaussian row conditional on every public label.

Let Z_C be the skew chain's known source-zero carrier. Its source-dependent residual has ordinary Lp size at most Lambda_p rho sqrt(d), with the actual finite source-profile/numerical envelope included. The two gradient-source blocks have dimension 2d; the number of private coordinates used to evaluate them is not substituted for d in this ordinary source bound. The primitive block source is anchored and has first rho. At a changed inner passive, its propagated residual is controlled using the declared incoming first and complete ordinary profile, not by borrowing an Lp coupling from an L2 theorem.

Hence the unmarked body retained in (6) has ordinary energy

    ||c(Q_C-Z_C)||p <= Lambda_p sqrt(d) |c| rho             (7)

plus assigned numerical errors. It is generally larger than kappa e and need not vanish when the marked source energy vanishes. It is therefore not eligible for deletion at a kappa-e mark. The known total carrier b Z_E+c Z_C+sqrt(v-b^2-c^2)Z_keep has exact covariance v I.

The full execution still includes Q_E and its complete tape. Only after the literal deletion in (6), and only if no other retained field reads those private records, can their surviving known carrier be represented by its exact Gaussian coisometry. This does not compress the actual execution or its future numerical prior.

## 4. Complete first, actual ancestor substitution and retained graph premise

At primitive widths s_min, the safe skew-channel source-dependent full first is O(rho), while a saved q_i root/caller contributes O(rho/s_min) times its actual q_i derivative. The residual public p-first is O(rho^2). Thus the side in (6) returns

    new private first:  Lambda |c| rho,
    root first:         Lambda |c| rho/s_min times the saved-root row,
    p-first:            Lambda |c| rho^2,
    physical caller:   Lambda |c| rho/s_min times the actual q_i caller. (8)

Add the marked residual first/caller multiplied by b. Whole Gaussian carriers have their known O(1) first and are not included in a claim of small residual first.

For retention one needs more than (8). Require each actual primitive center q_i(R) to have an ELIGIBLE kept-center replacement qhat_i on a specified active Gaussian selector, with its actual state/profile error Delta_i. Replace q_i by qhat_i simultaneously in EVERY complete skew-chain source call, including the anchoring value g_i(q_i), with all moving-origin paths retained. The Lipschitz comparison gives

    || c Q_C(q)-c Q_C(qhat) ||p
      <= Lambda |c| rho/s_min sum_i ||Delta_i||p.           (9)

Under |c| rho/s_min <= A this is an A-weighted old-center deletion row. The retained implementation is the SAME gradient-pair chain at those kept centers, not a new covariance oracle.

A controlled-block retention certificate must additionally verify the actual source syntax: each weighted pair source, its known Gaussian rows, positive-power internal edges, final original-gradient terminals, anchors and protected independent keep. Its output original-gradient mass is bounded by Lambda |c| rho/s_min; its old ancestor graphs are kept through their actual selectors. The standard weighted whole-block proof can then be applied to this explicit graph. This is a premise to be checked for the named native q_i, not a consequence of an arbitrary saved query or of small first norms alone.

After substitutions the active read set includes all nonlinear records still used by c Q_C(qhat), plus the known carrier row from the deleted marked arm. Only truly unread records may be pruned or replaced by an exact known coisometry. The full actual tape, all dropped evaluations and all direct/nested source calls remain charged.

## 5. A nonempty fractional-gain regime

For fixed numerical constants, write

    kappa=A^(1+g),   b=A^z,   rho=A^gamma,   s_min=A^w.

A sufficient derivative/caller guard is

    z>0, gamma>0, w>=0,       z+gamma+w <= g.              (10)

For canonical saved centers with root first O(1) and physical caller O(A^-1/2), (10) makes the retained side's full root first at most O(A) and physical caller at most O(sqrt(A)). Other actual center rows must be inserted in (8), rather than reset to this canonical type.

The side coefficient has exponent 1+g-z-2gamma. Under (10), with gamma<1, that exponent is positive; b^2+c^2 therefore fits a numerical variance share for small enough A. Strict inequalities leave logarithmic and fixed-constant margins. The marked law term improves kappa e by A^z, while the deleted marked energy improves e by A^z. The retained unmarked body (7) has exponent 1+g-z-gamma and cannot be assigned either marked bound.

For example, g=1/2, z=1/8, gamma=1/8, w=1/8 leaves a strict derivative margin and gives a positive 1/8 gain on the marked orientation law. This only illustrates the supplied-word port. It does not assert that the actual native common-query word can be represented at those widths with an admissible matching error.

If delta_C has the supplied fixed-order floor Lambda sqrt(d) rho^B, its amplified contribution is Lambda sqrt(d) A^[1+g-z+gamma(B-2)]. For any fixed desired nominal accuracy this can be assigned a finite B subject to the ACTUAL B-dependent radius guard. Arbitrary-order work or a uniform B-independent complexity is not claimed.

## 6. Literal work and remaining construction gate

One fork executes one complete marked pair and TWO complete block-gradient pairs, plus the displayed known completions. Its literal bill is

    Q_marked_pair + 2 Q_block_pair + Q_fill,

including source-center construction/replays, actual original VALUES, all first/caller sweeps and storage. Common centers may be cached only on identical arrays. General dense known fill matrices carry their real arithmetic cost; the cheap native specialization requires actual scalar-block rows or another certified completion.

The exact supplied-word covariance identity, finite comparison and retained decomposition are established above. A native true-orientation repair still requires:

- the SAME-root/common-input coefficient identity connecting J_R and C_R to the actual covariance-minus-forward defect;
- admissible heat/finite-source calibration without independent primitive smoothing of a correlated nonlinear genealogy;
- the named saved-center retained-graph certificate and its actual profile/caller/deletion errors;
- a same-endpoint joining proof for any old observer of the coefficient roots, or exact conditional cancellation before they are averaged.

These are concrete input contracts. The positive fork removes the forward-adjoint ordering obstacle once those inputs are supplied; it does not remove their source or cost obligations.

## 7. Explicit controlled-block import and bounded checks

The retention premise in section 4 is checked by re-running the following literal source table, rather than by declaring a new normalizer output to be a primitive gradient. In LOW30 (`research-source/High Acc Ideas/ai-bucket/30_low_acc.tex`), the relevant statements are `b27:compiler:paths`, `b27:nc:raw-paths`, `b27:raw:graph-certificate`, `t30:eq:normalized-return` and `t30:eq:retained-graph`.

For each block h_i=rho[g_i(qhat_i+s_i x)-g_i(qhat_i)]/s_i, the original VALUE multiplier rho/s_i and the SAME incoming query width s_i give internal edge rho. The qhat_i affine/old-block row is never divided by s_i. The second stationary pair consumes the whole output of the first; its inherited source path gains the second rho factor. J_0 is only a known swap/sign, while B_in and B_out retain their actual block-row bounds. After the external multiplier c the total outgoing original-gradient mass is at most Lambda |c|rho/s_min. The old kept qhat_i graphs remain whole admitted blocks, and every original centering query is included. A later native original-force terminal and its protected carrier must still be attached in the host's declared physical units; a sum of completed pair outputs is not retyped as a single original gradient.

The correct LOW30 pin is `7767ac389d6b677eb9517b2585f7a24cb9f85ede47d099a653c635b7be750fb8`. At fixed compilation orders, the positive internal rho powers and the host's small output mass permit the usual absolute-block guard. The active selector, state error and protected-row checks for the actual qhat_i are still required.

The adjacent `check_positive_skew_fork.py` passed 2,728 exact/reference and rational-guard checks. They verify the conditional kernel, its positive keep, the noncommuting cross sign, the retained baseline, an affine case whose entire projected orientation is the two-curl contribution, and the nonempty width regime. Maximum matrix discrepancy was below 3.9e-16. These tests do not establish the remaining native heat identity or saved-center eligibility.
