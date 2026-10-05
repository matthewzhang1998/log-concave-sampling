# Independent audit of the retuned mean/Gram compiler

2026-10-05. This is a mathematical and finite-interface audit of the proposed scalar retuning, not a numerical execution of LOW30. It uses the literal imported raw-split, self-reserve, native path and serial-count contracts. No sealed source file was modified.

## Verdict

PASS for the following guarded and scoped claims:

1. For a physical positive variance share `u`, changing the normalized raw-split padding to the literal declared radius `mu=r`, and using gradient order seven, gives the proposed six physical error powers. For `0<A,u<=1` and `u>=c A^(3/2)`, they are bounded by `Lambda A^4/u` times the inherited one-energy/retained-label profile.
2. The actual carrier-subtracted private and retained-caller first is `Lambda[A+A^(3/2)u^(-1/4)]`. In particular it is `Lambda A` for `u>=c A^2`; the error window above is stronger for small A.
3. The ideal positive Gaussian covariance mixture with the documented Gaussian-to-HS Lipschitz coefficient `C A^2` has an all-real-parameter conditional subGaussian linear-functional bound. With a fresh independent keep Gaussian of variance comparable to `u`, its conditional KL divergence from the matched Gaussian is `C D A^8/u^4`, uniformly in the exposed label, without dimension-smallness.
4. The completed Gram W2 law still has its existing mixture term `C A^4 u^(-3/2) sqrt(D)`. Retuning the native piece does not remove or improve that term.

The conclusion is conditional on actual native admission, all complete banks and retained labels, and all absolute numerical floors. It does not establish an unbuffered RAW approximation, a zero-buffer endpoint theorem, a dimension-free numerical KL value, or uniformly numerical compiler constants independent of public logarithms.

## 1. Exact substitution, with logarithms made distinguishable

Freeze the public scalar buffer `u` and the pre-enumerated public logarithmic factors. Write

    r=rho=kappa A/sqrt(u), delta=alpha A, a=beta A,
    mu=r, b=7, z=M rho.

Here `kappa,alpha,beta,M` are separate known public-log factors. Their distinction prevents a changing use of `Lambda` from obscuring a literal padding choice. The final `Lambda` can absorb their fixed powers and the fixed-order-seven compiler constants. The actual padding is the full number `r`, including `kappa`; it is not the bare power `A/sqrt(u)`.

The LOW30 `b27:raw:split` error is

    Lambda sqrt(n) [rho z^(b-1)+r^3 delta+r^2 delta mu
                    +r^2 a delta+r^4 delta mu^(-1/2)+r^4 delta].

After the physical readout `sqrt(u)`, the six terms are exactly

    M^6 kappa^7 A^7/u^3,
    kappa^3 alpha A^4/u,
    kappa^3 alpha A^4/u,
    kappa^2 beta alpha A^4/sqrt(u),
    kappa^(7/2) alpha A^(9/2)/u^(5/4),
    kappa^4 alpha A^5/u^(3/2).

After suppressing only the displayed public factors, their ratios to `A^4/u` are

    A^3/u^2, 1, 1, sqrt(u), A^(1/2)/u^(1/4), A/sqrt(u).

For `u>=c A^(3/2)` and `A,u<=1`, the ratios are at most

    c^(-2), 1, 1, 1, c^(-1/4) A^(1/8), c^(-1/2) A^(1/4).

Thus the asserted domination is correct. The first term makes `u~A^(3/2)` the binding algebraic scale. For general gradient order b, its ratio at `u=A^(3/2)` is `A^((b-7)/4)`, so seven is the smallest integer order that attains this endpoint uniformly as A tends to zero. Orders lie on the native `2+(1/16) N` grid, so seven is a constructed fixed order.

There is also a useful normalized interpretation. The intrinsic endpoint is `r^3 delta`. The reserve terms remain at most this scale when

    r^2<=mu<=r, a<=r, r<=1.

Choosing `mu=r` is the largest allowed padding at that scale and minimizes `r^2/sqrt(mu)`, the self-reserve first bound. Public constants need not obey `a<=r` literally if their ratios are carried explicitly in `Lambda`; the power comparison above is the unconditional stated calculation.

## 2. Actual first ports and exact guards

The raw-split actual residual first, apart from its explicit carrier, is `Lambda(r+r^2 mu^(-1/2))`. With the same readout and padding this is

    Lambda [kappa A+kappa^(3/2) A^(3/2)u^(-1/4)].

For `u>=c A^2`, the second bare power is at most `c^(-1/4) A`. For `u>=c A^(3/2)` it is even at most `c^(-1/4) A^(9/8)`.

This is a direct graph/path certificate. The LOW30 raw-split proof explicitly extends its last-source physical-path argument to external-center derivatives and fixed zero-tape source constants. It is not a derivative of a W2 claim. The existing actual source must retain the original full-source label bounds; a private first bound alone does not create a retained-label bound.

The following guards must remain literal:

- `u>0` is a public scalar frozen before all first/adjoint sweeps. No derivative is silently taken through a data-dependent variance or padding selection.
- `0<mu=r<=1` is required by the self-reserve lemma.
- The normalized source and its genuine-gradient component have their actual recorded complete first/curl/energy and final-physical-path bounds.
- The order-seven native call in raw split is applied to `(z/rho)G=M G`, not to the bare unamplified source. Its actual source radius, and every internal vertex after finite variance/readout allocations, must satisfy the imported `r_*(7,1/16,Lambda_7)`-type guards. Merely checking the bare monomial `A/sqrt(u)` is insufficient.
- The self-reserve first is `Lambda r^(3/2)` before physical readout, and all its recorded numerical smallness, filter and covariance allocation conditions remain live.
- All log factors are fixed from a sufficient public budget containing the original scales, `|log A|`, `|log u|`, requested absolute precision, dimension and the order-seven finite graph. At `u~A^(3/2)`, the bare source and reserve first powers are `A^(1/4)` and `A^(3/8)`; these powers allow, but do not numerically prove, the required log-dependent smallness.

There is no universal numerical A cutoff certified here. If the caller uses `u~A^(4/5)` as in the existing delayed-tail join, the bare normalized source and reserve first powers become `A^(3/5)` and `A^(9/10)`, respectively. The physical reserve first is `A^(13/10)`.

## 3. Retained variables, original VALUES, full roots and cost

This change is a scalar and fixed-order substitution in the existing finite source interfaces. The native estimate targets the normalized source's own mean. The inherited `F_Q` force-mean discrepancy `C A^3 sqrt(D)` is unchanged: normalization and physical readout cancel. The Gram first-coefficient absolute calibration floor may still be set to any fixed grade by logarithmic filter degree, and its matching `sqrt(u)` and `1/sqrt(u)` factors cancel upon restoration. Covariance target discrepancies of size `C A^4 sqrt(D)` continue to pay their actual positive-gap factor `u^(-1/2)`. None of these target comparisons becomes a freely adjustable numerical error.

The change must preserve the following exact composition:

- Each rectangular coefficient mean call keeps its genuine retained `(X,p)` records until its own complete conditional law comparison.
- Every private-H-dependent anchor remains inside each actual source occurrence. Captured origins are restored on their exact retained caller records, using identical-anchor numerical reuse.
- Complete mean banks are fresh conditional on their retained labels. `p` is integrated only after the conditional coefficient-mean comparison, producing the exact `B_Q B_Q*` orientation.
- The positive coarse-clock Gram assembly uses the same aggregate `1/sqrt(v_mean)` normalization, not an inverse individual node weight.
- Exposed labels are not confused with owned coarse roots. The latter are integrated in the covariance-mixture law comparison. No reference Gaussian is appended as an observer of already integrated private roots.
- Known source-zero carriers retain all perpendicular Gaussian coordinates. The number of roots is not numerically unchanged by raising the native order: replace the native mean graph and its complete dimensions by the order-seven graph everywhere in the old bookkeeping formulas.

The LOW30 raw-path and serial propositions apply at every fixed constructed order. Their constants therefore remain public-log polynomials at order seven. Counts for each raw gradient/remainder occurrence must be re-expanded using that graph and the new `mu`; these include all original-g leaves, source/mean/pair/filter/response/origin copies, Gaussian tapes and replays. The old structural count formulas still apply with the updated counts. This audit does not replace those counts by unit-cost actions or assert a favorable polynomial exponent.

Changing `mu` changes finite filter degree and precision through `log(1/mu)`. It adds no inverse-mu replication count under the imported finite-filter and serial contracts. Absolute primitive errors are allocated only after the complete new downstream paths are enumerated, so that their weighted sum is the requested absolute floor. Original scalar rows, rotations, roots, modes, anchors and caller restorations stay in that ledger.

## 4. Uniform conditional concentration of the actual covariance target

Fix the exposed label `Y=y`. Let the owned complete coarse tape G be standard Gaussian and let

    C=C_mix(G;y), F=C-E_G C,
    Lip_(G -> HS)(C)<=L_C<=C_0 A^2,
    0<=C<=R I, R<=C_1 A^2.

These are the actual coefficient-target hypotheses established in the existing one-energy mixture note. The PSD/operator bound also gives

    e_C^2=E||F||_HS^2<=E||C||_HS^2<=D R^2<=C D A^4.

For every fixed symmetric T, the scalar function `Tr(T C(G;y))` is Gaussian-Lipschitz with constant `L_C ||T||_HS`. Gaussian concentration gives, for every real t,

    log E_G exp(t Tr(T F)) <= (t^2 L_C^2/2)||T||_HS^2.

All constants are uniform in y. The same statement holds with T measurable in the frozen exposed label. It does not hold on this argument for T chosen after observing G. The centered operator is the true conditional mean `E_G C`, not an empirical or finite-compiler approximation.

The dimension factor in e_C is legitimate and appears once. A coherent random scalar covariance times `I_D` has Gaussian-to-HS Lipschitz constant proportional to `sqrt(D)` and cannot be substituted while keeping the hypothesis `L_C=O(A^2)`.

## 5. Posterior entropy plus Stein: a dimension-safe proof

For `nu>0`, introduce the analytical conditional mixture

    V|C ~ N(0,S_C), S_C=nu I+C,
    Sigma=nu I+E C.

Let `pi_v` be the posterior law of C given V=v and `pi` its prior. The entropy variational inequality and the preceding all-t concentration yield, for every test matrix T,

    Tr[T E(F|V=v)] <= KL(pi_v||pi)+(L_C^2/2)||T||_HS^2.

Optimizing T in the finite Hilbert space of symmetric matrices gives

    ||E(F|V=v)||_HS^2 <= 2 L_C^2 KL(pi_v||pi).

After integration,

    E||E(F|V)||_HS^2 <= 2 L_C^2 I(C;V).

The mutual information is bounded without a small-dimensional expansion:

    I(C;V) <= E_C KL(N(0,S_C)||N(0,Sigma))
            <= e_C^2/(4 nu^2).

For the second inequality, write Gaussian KL as one-half the Bregman divergence of `-log det` at Sigma. Its Hessian on a symmetric increment E is `Tr(S^(-1) E S^(-1) E)<=nu^(-2)||E||_HS^2` along the entire segment, since both endpoints have floor nu. The Taylor integral contributes the factor one-half, and Gaussian KL contributes another one-half. No series requiring `D A^4/nu^2` to be small is used.

Conditional Gaussian integration by parts gives V the Stein matrix

    T_V(v)=E(S_C|V=v).

The existing Gaussian Stein-to-W2 bound then gives

    W2^2(Law(V),N(0,Sigma))
       <=nu^(-1) E||E(F|V)||_HS^2
       <=L_C^2 e_C^2/(2 nu^3).

This matches the already proved one-energy mixture theorem. In fact that theorem alone is sufficient for the next step; the entropy route is an independent derivation of the same intermediate estimate.

## 6. Exact keep-buffer KL conversion

Let `h>0` and add `sqrt(h) Z` with a genuinely independent standard Z. For any two laws P,Q of finite second moment,

    KL(P*N(0,hI) || Q*N(0,hI)) <= W2^2(P,Q)/(2h).

To prove it, choose any coupling `(V,V_ref)`, mix the two Gaussian kernels `N(V,hI)` and `N(V_ref,hI)` using that same coupling, apply joint convexity/data processing of relative entropy, and use their exact equal-covariance Gaussian KL `|V-V_ref|^2/(2h)`. Minimize the coupling cost.

Consequently the final ideal mixture of total base variance `u=nu+h` satisfies

    KL(Law(V+sqrt(h)Z) || N(0,uI+E C))
       <=L_C^2 e_C^2/(4h nu^3).

If both `nu>=c_1 u` and `h>=c_2 u` for fixed positive constants, this is

    C D A^8/u^4.

For the equal split `nu=h=u/2`, the exact displayed coefficient is `4 L_C^2 e_C^2/u^4`. This is the required orientation of KL: mixture relative to the matched Gaussian. The opposite orientation is not supplied. No assumption that the right side is less than one is needed.

An existing independent keep root can supply h. The same total u must be split; an additional h cannot be appended while pretending the final covariance remains `uI+EC`. For the ideal mixture, one can also analytically split its isotropic u floor. For an actual finite return, the final keep must be an untouched independent executed root in its literal tape if this is to be a producer-level KL certificate.

## 7. Important boundary: ideal mixture versus actual finite return

The uniform clean `C D A^8/u^4` bound above concerns the ideal covariance mixture. KL has no triangle inequality. Native W2 closeness to that mixture therefore does not, by itself, give a KL bound for the finite native return.

A valid finite-return extension is available when the independent keep is preserved. Before adding it, establish a conditional W2 bound to the matched Gaussian:

    e_pre(y) <= e_native(y)+e_mix(y)+e_restore(y).

Then add the same independent keep to both laws and obtain

    KL(actual final law|y || matched final Gaussian|y)
       <= e_pre(y)^2/(2h).

The native conditional retained-label profile and its public-log factors must be restored. An L2-in-y native W2 certificate supplies an integrated-in-y KL certificate by this formula; it does not become a uniform-in-y certificate. At `u<=1`, the retuned native power `A^4/u` is below the covariance-mixture W2 power `A^4/u^(3/2)`, but constants/logs and all covariance target-restoration floors still count.

## 8. Diagnostic scope

The accompanying independent Python check verifies the exact monomial substitution, dominance thresholds, first-port powers, order-seven minimality at the proposed endpoint, numerical Gaussian concentration samples, and a one-dimensional smooth PSD covariance-mixture entropy/posterior calculation. Its tensor-product example extends to arbitrary D by exact independence, including dimensions with a large nominal KL bound.

Numerical quadrature is a diagnostic only. It does not prove the all-t inequality, the Hilbert posterior argument, or native admission, and it does not execute the LOW30 graph. The arguments above prove those mathematical implications from their explicit hypotheses.

## Local source pins

- `/workspace/shared/law-only-variance-join-20261005/LAW-ONLY-SHRINKING-BUFFER-JOIN.md`, especially Sections 4, 6, 9–11.
- `/workspace/shared/v9-curation-work/frozen/prerequisites/research-source/High Acc Ideas/ai-bucket/30_low_acc.tex`, labels `b27:raw:self-reserve`, `b27:raw:split`, `b27:compiler:mean`, `b27:compiler:paths`, `b27:nc:raw-paths`, `b27:compiler:serial`.
- `/workspace/shared/recovery-20261004/cost/endpoint-stein-quadrature-20261004/positive-law-repair/higher-cumulant-gate/same-carrier-feedback/GAUSSIAN-COVARIANCE-MIXTURE-ONE-ENERGY-LAW.md`.

The diagnostic results record the exact SHA256 of each reviewed local source and of this audit.
