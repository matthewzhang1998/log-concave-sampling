# Independent audit: all-layer currents, retained path, and reserve allocation

Publication copy: nonmathematical provenance wording and/or local paths were sanitized. Original/public SHA-256 values are recorded in `INVENTORY.json`; historical source/audit pins refer to the original versions.

Date: 2026-10-04. Separate extension of the pinned scalar terminal audit.

Reviewed notes and exact SHA256 pins:

- `ALL-LAYER-MARTINGALE-CURRENTS-AND-PATH-GAP.md`: `dbf2b9b04c72cd83970273cb487b2ac800ff9358dec945e1ef489f054b829d6b`.
- `INTERMEDIATE-GAUSSIAN-RESERVE-ALLOCATION-TEST.md`: `02db073eb1dca83c37ac8a8924cc34798550bb02f572e427a2053d053b2037c1`.

Exact reviewed snapshots are retained beside this audit. The separate `ALL-LAYER-MANIFEST.json` pins those snapshots and the independent checks. It does not modify the main scalar audit's manifest or verdict.

## Verdict

**Scoped PASS for the paired-current identity, the retained three-time covariance obstruction, and the stated conditional-C2 reserve allocation test.** Two source qualifications were corrected during review before these versions were pinned: the impulse clock is now smoothly Gaussian-generated, with its unbounded amplitude-dependent first explicitly exposed; and the current identity now requires the future map to have no unaccounted direct access to the selected clock bank.

The reported gap is between an actual finite quadrature-source path and its own conditional clock mean. It is not an exact-flow or posterior-law lower bound. The reserve exponents describe the current worst-case conditional certificate, not a lower bound for all reserve-based or randomized algorithms. No new general `c_P` family is proved.

The primary checker passes on the reviewed files. A separate independent script passes **82 checks**, including a vector-valued paired-current fixture, smooth capped Gaussian-clock path firsts, a seven-count `N`-stratum path-gap test, exact reserve exponent arithmetic, and snapshot integrity.

## 1. The paired rank-one/rank-two current is exact

Condition on a record `R` independent of the selected fresh bank except for the entering records already fixed before that bank is sampled. The gathered input vector is `q+D`, with `q` fixed and `E[D|R]=0`. Every remaining effect of that bank must pass through this vector. Every other independent primitive can be frozen into `R`, including the final keep, or averaged separately.

For a fixed `C1` remaining map `F_R`, put `J_t=DF_R(q+tD)` and `Y_t=F_R(q+tD)`. One derivative of the composite test gives

    E[phi(Y_1)-phi(Y_0)]
      = integral_0^1 E[(J_t D)·grad phi(Y_t)] dt.

Split `J_tD=(J_t-J_0)D+J_0D`. Conditional centering annihilates `E[(J_0D)·grad phi(Y_0)]`. Integrating the difference of test gradients then gives

    integral_0^1 (1-s)
      E[(J_0D)·Hess phi(Y_s)(J_sD)] ds.

Thus the report's `V_t`, tensor orientation in `Q_s`, and weight `1-s` are correct. Only the symmetric part of that tensor matters because the test Hessian is symmetric. No derivative of `DF_R` is taken. The derivation therefore uses first actions of the finite graph, rather than an unavailable original third derivative or a derivative of an HVP child.

The independent vector fixture uses a map from three input coordinates to two output coordinates, mixed derivatives, and a nondiagonal test Hessian. Its residual is approximately `1.1e-17`; this checks the tensor orientation beyond a scalar-only example.

If an observer reads the selected raw clock bank outside this gathered vector, `F_R` is no longer such a fixed map. Conditioning on that observer may destroy `E[D|R]=0`. The repaired text correctly requires an additional compared coordinate/current instead. Applying this identity separately to successive replacement graphs can telescope all layers, but does not eliminate any of their regenerated currents.

## 2. Nonlinear feedback is a real conditional first-order cost

For `F(q)=B-w b(q)`, the rank-one field is exactly `-w[b'(q+tD)-b'(q)]D`. A bounded continuous `b'` has no uniform quantitative modulus under the stated `C2` class. The elementary bound is consequently `O(|w| a E|D|)`.

For the smooth near-kink force in the report and a uniform centered residual on `[-epsilon,epsilon]`, the exact mean is

    a c { 1/2 [sqrt(epsilon²+eta²)
                +(eta²/epsilon) asinh(epsilon/eta)] - eta }.

After dividing by `a epsilon`, this tends to `c/2`. Its derivative lies strictly between `a(1/2-c)` and `a(1/2+c)`. The uniform residual is smoothly realizable using `2Phi(G)-1`, so the example uses neither a discontinuous clock nor an inadmissible Hessian.

This example is conditional at a fixed caller. It does not establish the same lower bound after averaging an unobserved nondegenerate Gaussian caller with its actual genealogy. The notes maintain that distinction correctly.

When the future map is linear, the rank-one field vanishes. The rank-two field does not vanish; the earlier scalar note's common-clock random-scale fixture explains why conditional centering alone is not a dimension-safe weak theorem.

## 3. The transverse retained-path obstruction survives legal smooth clocks

The reviewed clock is `S=s_-+(s_+-s_-)Phi(G)` in an interval strictly between the first and second observation times. Hence every displayed kernel read stays away from its kink; a small fixed `C1` smoothing can leave every such read unchanged.

For each `S`, write the three-position output as `A(S)X+B(S)Z`. The conditional mean given `X,Z` is `bar A X+bar B Z`; its support lies in the two-dimensional coefficient plane. For a unit normal `w` to that plane,

    w·q_bar=0,
    E[(w·q_S)²]=E[(w·D)²].

The cross covariance between `q_bar` and `D` is zero because `E[D|X,Z]=0`. If the last displayed energy is positive, then the mean-path covariance has zero quadratic form in direction `w`, while the clock-error covariance has a positive one. Therefore subtracting that covariance from the baseline Gram is impossible. Projection onto `w` also gives the exact lower bound

    W2(Law(q_S),Law(q_bar)) >= sqrt(E[(w·D)²]).

An independent final Gaussian unused by these three path positions contributes no covariance to this projection.

The reviewed unguarded clock derivative grows like `|X|+|Z|`, as the note now states. Its conditional/moment first exists, but this is not a global uniform source-first admission. To verify that this is not merely a discontinuous or unbounded-first artifact, the independent script also replaces `X,Z` by explicit odd `C1` caps of independent Gaussian roots. Their positive variances replace the two Gaussian variances in the same Gram calculation. The conditional mean still lies in the same plane, the transverse variance remains positive, and the complete clock/source first is globally bounded. The finite directional first matches direct finite differences, including clock motion and cap derivatives.

That capped enhancement is itself just a source-geometry witness. It does not make its reference path an exact posterior flow or establish all native host contracts.

## 4. The `N`-stratum transverse scale is genuinely `a² N^-3`

The reserve note correctly distinguishes an upper residual-energy scale from a necessary reserve size. Here the latter has a nondegenerate fixture.

Replace the single impulse by `N` independent uniform times, one per stratum of a fixed interval `[s_-,s_+]`, with weights `h/N`. The conditional mean coefficient vectors are unchanged with `N`. Let

    f_X(s)=w·sin_+(times-s) cos(s),
    f_Z(s)=w·sin_+(times-s) sin(s).

Both are smooth throughout that fixed interval. At the specified three ordered times, their weighted derivative energy is nonzero. Independence and conditional centering give

    Var(w·D_N)
      =(a lambda h/N)² sum_j
        [Var(X) Var_j(f_X(S_j))+Var(Z) Var_j(f_Z(S_j))].

For a continuously differentiable nonconstant coefficient, the stratum-variance expansion gives

    sum_j Var_j(f(S_j))
      = (s_+-s_-)/(12N) integral_(s_-)^(s_+) |f'(s)|² ds
        + o(N^-1).

Thus the transverse variance is `Theta(a² N^-3)`. The normal plane varies continuously with small `a`; the nonzero derivative energy at the limiting harmonic plane keeps the constants uniform at sufficiently small heat.

For the independent capped fixture at `a=0.03`, the quantity `N³ Var(w·D_N)/a²` ranges from about `1.8644e-5` at `N=1` to `1.8744e-5` at `N=64`. This checks the advertised scale while preserving the same three retained observers.

This lower scale is attached to the specified path-source comparison. It does not turn the scalar quadrature upper bound into an information lower bound for exact nonlinear flow or arbitrary posterior sampling.

## 5. Intermediate reserve allocation and the actual mean price

A covariance subtraction in a zero-Gram direction requires a new positive reserve in that direction. For the nondegenerate fixture, strict covariance domination requires

    tau_j >= c a N_j^-3/2.

A fully positive matrix reserve may need additional directional/eigenvalue and dimension accounting. The scalar necessary direction alone suffices for the displayed allocation test.

Under the current conditional-C2 future-map contract, a perturbation of normalized size `tau_j` passing through `M-j` force levels costs at most

    C sqrt(a) a^(M-j) tau_j.

This first-power conditional profile cannot generally be replaced by a quadratic one: at `q=0`, the near-kink Gaussian fixture has mean

    E b(tau G)
      = a c [E sqrt(tau²G²+eta²)-eta]
      ~ a c tau sqrt(2/pi).

An independent centered final Gaussian cannot erase that mean. Further affine force levels propagate it by their actual factors of `a`.

Requiring this worst-case certificate to fit `a^R`, while preserving the displayed covariance gap, gives exactly the report's count exponent

    max(0, (2/3)[R-(M-j+3/2)]).

In particular the nearest earlier layer gives `max(0,(2R-5)/3)`. These are constraints of the stated safe allocation/certificate. The transverse quadratic fixture and the nonlinear mean fixture do not by themselves prove that one potential simultaneously saturates every stage of every possible algorithm. The note expressly excludes such a general lower-bound claim.

A complementary scalar buffer calculation gives the same arithmetic: if the residual has both energy and Lipschitz scale `e`, adding an independent width `tau` costs `e²/tau` for its conditional buffered comparison, while the current strong propagation prices the added reserve by `tau`. Their sum is minimized at `tau=e` and remains order `e`. This is another certificate calculation, not a claim that those upper bounds equal every method's actual error.

If a stronger observer-qualified theorem instead supplied a `tau^p` future price, the formal count exponent is

    max(0, (2/(3p))[R-(M-j+1/2)-p]).

The algebra in the report is correct. The last-only posterior keep theorem's `p=2` does not transfer to arbitrary intermediate conditional nodes.

## 6. Exact value cancellation also cancels the reserve's unread status

The value identity

    b(Q+tau G)+[b(Q)-b(Q+tau G)]=b(Q)

is legal with two distinct original VALUE queries and reuse of the shared first query. Its firsts require only the corresponding Hessian actions.

But the corrected observer sees the unbuffered `Q`, either directly or through `(Q+tau G)-tau G`. For a vector path, the joint reference `(Q+tau G,G)` still has a deterministic reconstructed zero-Gram direction when `Q` does. A linear observer recovers that direction. It is therefore invalid to count the same Gaussian as an unread reserve after making this exact correction.

Likewise, a frozen linear HVP correction has zero Gaussian mean and does not cancel the near-kink mean price. Differentiating that saved HVP would require an unavailable higher original derivative. A new finite VALUE or genuine marginal-law counterpacket remains possible, but would be additional research, not a consequence of the exact telescope.

## Reproduction and limits

Run `OPENBLAS_NUM_THREADS=1 python check_all_layer_independent.py`. The script verifies the two note snapshots before testing the currents, source geometry, stratified scale, and formal reserve arithmetic.

The pinned notes correctly identify a live reusable all-layer gate. Their conclusions are conditional/current/source constraints and an explicit finite retained-path witness. They do not establish a posterior-law impossibility result, a dimension-uniform high-order positive compiler, or a new family recurrence.
