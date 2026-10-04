# Stationary marked-pair skew channel: reconstructed mathematical source

Publication copy: nonmathematical context and/or local paths were sanitized. Original and public SHA-256 values are recorded in `INVENTORY.json`; source/audit pins below identify their historical versions, not these edited bytes.

RECONSTRUCTED NEW VERSION, 2026-10-04. This is a new mathematical reconstruction. It is not asserted byte-identical to the earlier source. Last verified original source pin: ace5e761c38177388cf809cf27ea4c30e8f6c943a1b92753843e4bf11f2ff7bf, formerly STATIONARY-MARKED-PAIR-SKEW-CHANNEL.md. Original scoped audits: C dc59da051022264bd597c326c42b78ad888c6328fa80c64e638b124a28633bbe and R 40fe1d02fa503f4b68a3e48c3f71c61a5f1299cdc93f46db47cb25a3ea915e74. Those audits bind the old bytes; this reconstructed version needs its own binding check.

## Scope

This is a finite VALUE response to a specified skew heat-matrix word times a FORWARD marked coefficient. It preserves one actual incoming energy mark and avoids shifting a standard passive by an arbitrary nonlinear mark. It does not identify separately heated primitive Hessians with the original correlated nonlinear-query word. The marked adjoint/two-curl term in true orientation is a separate obligation.

All sources and finite versions are fixed. The fixed numerical pair-selection constant is c_sel (called s0 in LOW30); primitive widths s_t,s_0 are different parameters. Choose an ACTUAL known signed word coefficient τ with |τ|≤κ. Program readouts use τ, while bounds use κ. A norm bound is not substituted for the actual coefficient.

## 1. Marked pair

Let f_R be a complete square marked source with first A, curl at most κ, conditional centered energy e_R, and its actual label/caller profiles. The public G is standard CONDITIONAL on R,θ. In the intended construction R is independent of G conditional on θ. This is separate from later independence from an old consumer.

Use the admitted near-gradient VALUE pair

    Y=Pv(f_R;G,Ω_E)=Z_E+R_E^val,
    M_E(R)=c_sel E Df_R.

Z_E is its KNOWN standard Gaussian source-zero carrier, independent of G,R. The actual carrier-subtracted residual has integrated Lp energy at most Λ_p e_R and residual first O(A) under the declared padding. The conditional-G comparison has error δ_E(R) to

    M_E(R)G+sqrt(I−M_E(R)M_E(R)*) Z.

This reference is marginally standard at fixed R after G is integrated. Its innovation is not identified with the actual Z_E. No joint law of (Y,Z_E) is assumed.

For κ≥A² and μ=A^(4/3), the existing pair supplies

    ||δ_E||L2(R)≤Λ[κ+A^(7/3)]e+ε_E,
    e²=E_R e_R².

Every original source width, finite floor and caller profile of this imported pair remains charged. Its conditional mean bias η_E(R,G)=E[Y−Z_E|R,G]−M_E(R)G is separately available and is bounded by δ_E; the numerical first-chaos filter can give a sharper bound.

## 2. Two gradient blocks and known Gaussian rows

At the SAME R, save actual finite labels q_t(R),q_0(R). Let g_t,g_0 be exact original gradient sources with first at most one. Define the square block-gradient source on R^(2d)

    h_R(x_t,x_0)=ρ( [g_t(q_t+s_t x_t)−g_t(q_t)]/s_t,
                    [g_0(q_0+s_0 x_0)−g_0(q_0)]/s_0 ).

It has full private radius at most ρ. Write

    D_R=diag(H_t(R),H_0(R)),
    H_i(R)=E Dg_i(q_i(R)+s_i Z_i).

A fixed-order genuine-gradient pair on h_R has selected coefficient c_sel ρ D_R, incoming first O(ρ), and residual/private first O(ρ). Its conditional-input error envelope has Gaussian L2 bound ε_h(R), with the actual √(2d)ρ^B floor, finite source restoration, and B-dependent radius guard.

The two pair banks use identical source versions, q_t,q_0,s_t,s_0, so their analytical D_R is the SAME matrix. Their complete private tapes are independent conditional R.

For fixed known contractions P,C:R^n→R^d, let

    L=[P*,C*]:R^(2d)→R^n,
    J=[[0,I],[-I,0]],
    B_in=L*/sqrt(2),      B_out=L/sqrt(2).

Both B rows are contractions; J is a known orthogonal Gaussian rotation. Known rows/scalars are fixed in the exterior caller and coefficient roots unless their derivatives are separately charged.

## 3. Complete stationary chain

Execute, at input x∈R^n,

    x0=B_in x+sqrt(I−B_in B_in*) Z_in,
    y1=PB(h_R;x0,Ω_1),
    y2=PB(h_R;J y1,Ω_2),
    K_R(x)=B_out y2+sqrt(I−B_out B_out*) Z_out.        (1)

The square roots are known constant covariance completions. No selected matrix is evaluated. Each pair bank and each fill is complete. For paired evaluation, the same pair tapes and BOTH same fill tapes are reused.

The exact Gaussian-reference chain is stationary with selected matrix

    N_R=(c_sel ρ)² L D_R J D_R L*/2.                 (2)

The matrix W_R=L D_R J D_R L* is skew and satisfies exactly

    W_R=P*H_t H_0 C−C*H_0 H_t P.                    (3)

This is the specified heat word. Symmetry of each H_i gives the reverse product; no completed VALUE output is declared a new gradient.

## 4. Actual marked difference

Use the SAME complete chain tapes at x=Y and x=Z_E and emit

    T(G)=2τ/(c_sel²ρ²)[K_R(Y)−K_R(Z_E)].            (4)

The outer pair carrier and output fill cancel pathwise. The incoming first of the whole chain is O(ρ²), hence

    |T|≤Λκ|Y−Z_E|,
    ||T||Lp≤Λ_p κ e_p.                              (5)

This is an actual executed one-energy bound. It preserves all original incoming and coefficient aliases.

The complete chain-private first is O(κ/ρ): differences of private derivatives use their ordinary O(ρ) bounds. Marked-pair common carrier paths cost O(κ), and its residual paths cost O(κ A). The public first is O(κ A), since Z_E has no G path. An old coefficient label in q_i contributes at most

    Λ κ/(ρ s_min) times its actual q_i first,

plus the κ multiple of the marked-pair caller. If a coefficient root is later owned, this row belongs to the complete private first; conditioning on it temporarily does not remove its cost.

## 5. Mean calibration at proper marginals

Let Ψ_R(x)=E_chain K_R(x), d_R(x)=Ψ_R(x)−N_Rx. At a standard input, the two conditional pair comparisons give

    ||d_R(G_std)||2≤δ_C(R),
    δ_C≤ε_h,outer+Λρ ε_h,inner≤Λε_h.

The actual incoming first and ||N_R|| imply Lip(d_R)=O(ρ²).

Y can be coupled to its stationary reference conditional on G,R with error δ_E. That reference is marginally standard conditional R after G is integrated. Therefore

    ||d_R(Y)||L2(G,Ω_E)≤δ_C+Λρ²δ_E.

Separately Z_E is exactly standard conditional R, so ||d_R(Z_E)||2≤δ_C. These are MARGINAL defect estimates. They do not replace the actual joint pair (Y,Z_E) by a fictitious joint Gaussian law. The exact linear difference N_R(Y−Z_E) is retained until its conditional expectation.

Consequently

    E[T|G,R]=τ W_R M_E(R)G+error,

with integrated error bounded by

    Λκ[ρ^−2 δ_C+δ_E+η_E].                           (6)

η_E≤δ_E is allowed, or use the actual smaller first-filter mean bias. All fixed c_sel factors are included in Λ. The gradient-prior error is genuinely amplified by κ/ρ²; it is not an unproved source-sensitive relative prior.

This fixes the earlier shifted-passive issue: an arbitrary p+v need not have a standard marginal, while both current inputs have the correct exact or nearby standard marginal. Actual same-tape coupling is used only for the energy/first proof.

## 6. Positive coefficient use and work

If every new R record is independent of the old endpoint conditional only on G,θ, put all R records inside the new private comparison and leave a fixed independent Gaussian keep untouched. At numerical primitive widths with complete owned-root first O(κ/ρ), the independent-buffer Riesz bill is

    Λκ²e/ρ.

The general root/caller bill includes the displayed 1/s_min and actual old q_i path. At ρ≥κ/A the numerical-width bill is at most Aκe. Add(6), including the marked pair's κδ_E row and all finite priors. Actual residual energy is κe, complete fresh first at most A, public first κA, and exterior caller is the stated κ/(ρs_min) row.

If R is retained as an external observer, the all-private comparison is unavailable. Its conditional coefficient field and later root/current obligations must be returned instead.

The literal bill is

    Q_marked-pair+4Q_block-pair+Q_fill,

including original-gradient VALUES, state/center construction, FIRST/adjoint sweeps, all cache replays and storage. There are two independent pair banks, each evaluated at both changed passives. Common centers are cached only on identical arrays. Each block source call expands into its two original-gradient differences.

Known fill matrices are not automatically cheap. General dense P,C require their real arithmetic/storage cost. The cheap native specialization needs the actual scalar-block Gaussian rows or another supplied bounded-cost completion; scalar-block row Gram and fills are finite known coefficient matrices tensored with I_d, including singular aliases.

With ρ=A^γ, e_bound≈sqrt(Dim)A^3.9 and actual source dimension n_h≤ΛDim A^−κ_n, the amplified side prior is

    Λκ sqrt(Dim) A^[γ(B−2)−κ_n/2].

A strict fixed B>2+(4.9+κ_n/2)/γ pays the nominal Aκe_bound budget, subject to its actual radius guard and source tolerance. This is an assigned numerical budget, not a relative promise when the actual e is smaller.

An optional odd-mark version uses Y^±=Pv(f_R;±G,Ω_E), same Ω_E, and readout τ/(c_sel²ρ²)[K_R(Y^+)−K_R(Y^−)]. It has the same mean target but executes both marked endpoints unless actual common expressions are cached. The primary implementation is(4).

## 7. Arithmetic and boundaries

Cancel identical known carriers algebraically where possible and allocate all remaining absolute evaluation/rounding error through 2|τ|/(c_sel²ρ²) and the expanded finite path norm. Fixed floating-point cancellation is not a proof of(5). Precision affects the declared bit/logarithmic bill, not a hidden replica factor.

The original adjacent checker passed4,486 known-row/reference/mark/grade cases: contraction pads, noncommuting block identity, covariance I−N_RN_R*, factor2, same-tape marked difference and prior inequalities. C independently checked350 further cases. These historical checks bind the original source; they must be rerun or reconstructed for this new version.

Two main native obligations remain:

1. Identify the H_i(R) heat matrices with the actual correlated primitive-query word of the SAME finite E. Independent primitive smoothing is not justified by marginal query legality or whole-source OU stability.
2. The coefficient in(6) is skew-word times FORWARD M_E. The true orientation −Sym E[R_E Curl_f] equals Sym E[Curl_f R_E*]. Replacing R_E* by R_E leaves a genuine marked two-curl term, indispensable in the hostile affine fixture. A separate positive-fork use can recover the correct cross ordering, but its unmarked side must then be retained and certified.

Thus this is an executable supplied-heat-word generator, not an all-source true-covariance theorem.
