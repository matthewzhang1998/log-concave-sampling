# Aggregate source ports with a smaller endpoint cutoff

2026-10-05. This is separate from the completed standalone grade-17/5 audit. It checks the proposed extension to cutoff η=K A^(3/2), where K is a sufficiently large declared public-log factor. It does not audit an independent skew-correction target comparison.

The full written `AGGREGATE-PORT-EXTENSION.md`, SHA256 `0be5daf191b9618d36350e02458168fd0b1ef14466abb4d10d294d5d3aab971b`, was subsequently read and agrees with the passed argument below, including its actual-guard requirement and its explicit limit on claiming a new skew service.

## Conclusion

The aggregate full-first, curl, and energy argument is valid. There is no hidden maximum-node amplification from sharing the aligned carrier, embedding node tapes in the complete tape, or replaying the entire source in a final compiler. The result uses the positive outer weights and the actual node path certificates; it cannot be obtained merely by differentiating a LAW error.

Let the actual preterminal residual at node t have uniform private-and-caller first

    L_t ≤ Λ[A+A^(3/2)/√v_t].

The same positive finite dyadic rule satisfies

    Σ_LAW ω_t v_t^(−1/2) ≤ C[w^(−1/2)+1],

uniformly in the cutoff η. This is the integrable endpoint case: a panel of gap d contributes O(√d). Consequently, when w≥A,

    L_bar := Σ_LAW ω_t L_t + Σ_RAW ω_t O(ΛA) ≤ CΛA.

Unlike the old pointwise argument, this remains true if some L_t are larger than A.

## Why the carrier rotation does not lose this gain

For each node rotate the entire known terminal carrier rX_t+hG_h into the same D-dimensional coordinate G. Keep that node's perpendicular coordinates. The map from the global tape to any one full node tape is a scalar-block isometric embedding followed by an orthogonal map, so its operator norm is one. Its coefficients are deterministic and frozen before differentiation.

Write the node terminal as g(S_t−Δ_t), with

    S_t=rtz+κ_tG,  κ_t=√(1−r²t²),
    Δ_t=r[d_t+P_Q(X_t−d_t,Y,N)].

The prefix contributes at most C Aw(1+L_t) to the full first of Δ_t. Therefore

    Σω_t Lip(Δ_t) ≤ C[L_bar+Aw(1+L_bar)] ≤ CΛA.

For E_t=g(S_t−Δ_t)−g(S_t), its common-G derivative contains

    κ_t[Dg(S_t−Δ_t)−Dg(S_t)],

which is symmetric. Every other common-G or off-G derivative contains Dg times a derivative of Δ_t. Thus the square-lift curl of Σω_tE_t is bounded by C A Σω_tLip(Δ_t)≤CΛA². Positivity and the triangle inequality give this directly; no maximum over t is introduced.

The full first of the terminal sum is bounded by

    A Σω_t[κ_t+Lip(Δ_t)] ≤ CΛA,

and the same bound applies to its caller first. The baseline Σω_tg(S_t) is a genuine gradient with first at most A.

For energy, use the literal zero-anchor property and the actual complete dimension d_t. After rotation,

    ||Δ_t||_p ≤ C_p Lip(Δ_t)(|z|+√d_t).

All √(d_t/D) are declared public-log factors. Weighted Minkowski therefore gives ||Σω_tE_t||_p≤Λ_p A²(|z|+√D). Independence across nodes is unnecessary for this upper bound. Shared G causes no additional factor.

Final native reentry substitutes a new caller and a fresh complete global tape into the same deterministic positive sum. All node certificates are uniform in the caller, and the scalar rotations remain frozen. Consequently the weighted bounds apply to every reentry; replay changes the count, not the norm argument. Numerical error allocation must still retain every node's actual individual path weight and buffer rather than replace them by an averaged execution tolerance.

## Native scales and limitations

On t≤1−η, v_t≥η/2. At η=K A^(3/2), the worst normalized mean/Gram radius is O(ΛK^(−1/2)A^(1/4)); the raw self-reserve radius is at most a declared public-log factor divided by K; and the mixed-K radius is O(K^(−1/3)A^(1/2)).

Choose K to satisfy the actual numerical self-reserve window, not merely to obtain a nonnegative A exponent. The native mean and pair windows must still be checked. Also impose η≤w. Calling K a public-log polynomial requires the imported guard constants to be covered by that declared polynomial; the power count alone is not a numerical admission certificate.

For proposed combined parameters w=A^(5/6), h=A^(7/4), η=K A^(3/2), the new prefix balance is A^(15/4), the near term is O(√K A^(15/4)+K A^(49/12)), and the intrinsic mean/Gram contribution includes **ΛA^(15/4)**. Accordingly, any combined grade-15/4 theorem must retain a leading public-log prefactor. That prefactor cannot be absorbed into a numerical constant using a small-A window, because these terms have exactly the leading exponent. A separate valid skew-correction comparison is required to replace the old bulk A⁴/w error, which would otherwise remain dominant.
