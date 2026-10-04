# Intrinsic physical dimension for the fixed marked cap

Date: 2026-10-04. This is a proof extension of the imported finite constructor contracts, not an application of the theorem-as-stated and not a new marginal-Wasserstein hypothesis.

## Result

The fixed marked mean/pair cap used by Exact32 admits a **fixed-output-subspace variant** with the physical prefactor `sqrt(d)` in every substantive current, fourth-channel, and genuine direct allowance. Here `d=Dim`, the complete Gaussian input dimension is `n`, `P:R^n -> R^d` is the recorded coisometry, `Q=P^T P`, and the original square VALUE source is `D=P^T F`, so `QD=D`.

Use exactly the native original-source, first/adjoint, independent-shield, proper-cut, retained-label, and one-energy-profile contracts. In particular an arbitrary low-rank map supplied only through its marginal law is not an admissible substitute. Assume the required fixed moments of centered original VALUE replays, including the prescribed same-tape differences, are `L_p sqrt(d) rho^(1+z)`. The original radius and curl bounds remain `L rho` and `L rho^(1+k)`.

For `c*=199/200`, starting `k>=8/9`, and the Exact32 working cap `B=52/5`, define

    F_c(z,k)=min{4+z+c min(k,2), B, B−2+2k},
    G_c(z,k)=F_c(z,k)+c min(k,1).

The variant targets exactly `N(E D,I_n)` and, with the incoming `x` retained,

    K_M(x)=Mx+(I_n−MM^T)^(1/2)Z,   M=s0 E JD.

Its substantive bounds are `L sqrt(d) rho^F` at currents of rank at least four and separately at fourth Gaussian channels, and `L sqrt(d) rho^G` at genuine direct type. Numerical direct and rank-two allowances are selected separately. The ordinary incoming slope stays `L rho`. There are the same 10849 levels, 2877 main replays, force-ten queue, and polynomial-logarithmic multiplier of the complete original-source query bill. No new inverse-heat query exponent occurs.

A convenient explicit implementation change is to replace the output of **each complete marked mean or pair call** by

    T_Q(Y,Z_perp)=QY+(I_n−Q)Z_perp,

with one fresh standard Gaussian independent of the whole call and its retained inputs. The incoming pair variable is untouched. Do not perform this projection separately on the self-square and curl-orientation covariance summands before combining them.

## 1. Why the complete-return projection is legal

Since `QD=D`, one has `QM=M`. The mean target has mean in `ran Q`; the pair target has mean in `ran Q` and covariance `I−MM^T`, which is block diagonal with an identity block on `ker Q`. Its perpendicular Gaussian is independent of the physical output and the retained incoming variable. Thus `T_Q` fixes both exact targets, including the conditional pair law at every incoming `x`.

This operation is a common Markov postprocessing. For currents it replaces `T_j` by `Q^(tensor j)T_j` on output indices only; derivatives never hit the incoming variable. For a fourth channel it replaces the physical channel map `R` by `RQ`, a contraction, before its ordinary Gaussian fill. For direct entries, couple the fresh perpendicular Gaussian identically and use `|Q(y−y*)|<=|y−y*|`. Consequently it preserves all comparison types and all retained-label restrictions. An enclosing ideal Gaussian chain still contributes its actual tensor, fourth, or first power. It adds no new endpoint observer.

The exact self helper may have target covariance `c^2 I−Sym(M^2)`, and its orientation helper may have cross-subspace blocks. Those intermediate targets are not individually fixed by `T_Q`. Their exact covariance identity is first retained unchanged; only the complete call, whose total target is `I−MM^T`, is projected. Old pair projections are harmless because each old pair's own conditional target is fixed exactly.

The projector is known from the actual recorded coisometry. It is not a projection onto a sampled Jacobian range or an unknown averaged range. It requires known linear arithmetic and Gaussian generation, no original source query or derivative. It does not differentiate a completed output.

## 2. The dimension bookkeeping lemma

For an original replay with output in `S=ran Q`,

    rank(JD)<=d,   ||JD||_HS<=sqrt(d)||JD||op,
    rank(JD−JD^T)<=2d,   ||JD−JD^T||_HS<=sqrt(2d)||JD−JD^T||op.

Scalar-affine/sign/OU VALUE filters preserve the fixed output space. Conditional centered Gaussian moments of an ordinary radius-`r` S-valued source satisfy `||D−E_private D||_p<=C_p sqrt(d) r`. This follows by vector Gaussian Poincare plus scalar Gaussian concentration of its norm (or Hilbert-valued Gaussian Sobolev estimates). It applies uniformly in frozen translates. Improved marked energies use their original integrated profile instead of this ordinary estimate.

More importantly, if a tensor has its force-output slot in S and its output-versus-all-other-indices operator cut is at most `a`, then

    ||T||_HS^2=trace(T T*)<=d a^2.

This estimate also holds conditionally on old roots. A permutation or transpose of tensor slots preserves its full Hilbert norm. Tensor contractions in a tree carry one full-Hilbert factor and proper operator cuts on all other factors. Scalar Gaussian Riesz operators, Hermite projections, positive averaging, and conditional Jensen are dimension-free Hilbert operators. They do not replace `sqrt(d)` by `sqrt(n)`.

This statement does not assert that products, transposes, or mixtures remain in one common low-rank range. Only a single Hilbert estimate is carried. In particular it is not legitimate to estimate each factor in Hilbert norm or to use an independent `sqrt(n)` estimate for an identity Gaussian baseline.

## 3. Exhaustive local ledger

LOW31, `mcd:thm:cap-ten`, lines 1901–1932, lists all local types. Apart from numerical errors, the only local classes without the improved VALUE mark are:

1. `rho^16` ordinary two-stage completion;
2. `rho^(8+8k)` conservative frozen-curl helper;
3. `rho^(11−9e)` force-eleven constants and fourth channels;
4. `rho^(12−10e)` same-endpoint force-twelve direct boundary.

Every other entry is one of `rho^(8+z)`, `rho^(4+z+k−3e)`, `rho^(4+z+2k−3e)`, `rho^(5+z+k−4e)` and its higher-force successors, or `rho^(8+z+k−6e)`. The seed has the separate marked term `rho^(4+z)`. Their proofs explicitly retain one centered original VALUE norm, not a dimension-dependent norm of the completed Gaussian output. The replacement word and all-marked curl packets preserve this property through root-private, covariance-removal, old-root and nonlinear comparisons. See LOW30 `t30:thm:marked`, especially lines 1117–1148 and 1207–1224; the rectangular VALUE-tree mechanism is `vtree:thm:marked-tree`, lines 1971–2267.

Replacing that one supplied norm by its physical `sqrt(d)` norm proves the intrinsic bounds for every energy-bearing class with precisely the old powers and types. There is no additional small-energy inference about a derivative. It remains to check the four ordinary classes rather than merely calling them numerical.

### 3.1 Two-stage seed and ordinary completion

In LOW30 `b27:calc:ordered-action` the uncompleted action `T` is a linear combination of original `D` VALUES at rotated fresh inputs. Thus `T`, its polarization sources `U,V`, and its repeated two-stage raw sources all take values in S. Their input dimension can be arbitrary.

Run `b27:calc:selector-difference`, `b27:calc:pure-fourth`, and `b27:calc:reserve-sixteen` on these rectangular S-valued sources, using `P` to identify S with `R^d`, and add an exact independent perpendicular Gaussian when an ambient output is required. This is the same literal constructor: matrix-vector actions are still actual first/adjoint actions of the rectangular VALUE program, not matrix columns or an averaged operator.

For a selector, the HS estimate is `C sqrt(d) r^2`, while the matrix-test scale remains `Cr^2`. The mixture proof retains the product of that scale and its one HS energy and therefore gives `C sqrt(d)r^4`. In the pure-fourth comparison the three bounded factors multiply one centered source norm. Iterating only on the uncompleted action of radius `r'=L r^2` gives `L sqrt(d)(r^4 delta+r^16)`. Conditional variance reduction and the tower property are exactly those of the native two-stage proof. Thus both the unmarked `rho^16` and marked seed are physical-dimensional.

The rotation `A0^T H+(I−A0^T A0)^(1/2)Z` may use an ambient input Gaussian. It has exactly the source's conditional Gaussian marginal; its ambient norm is never substituted for the energy of the subsequent S-valued source.

### 3.2 Frozen-curl helper

Keep LOW30 `b27:marked:frozen-curl` literally. Its `F=G_C+H` contains an ambient linear Gaussian `G_C`; do not assign that Gaussian a physical-dimension energy bound. The nonlinear VALUE part H and its independent helper H′ are S-valued. The ordinary sixteenth-power helper is applied to `(H′,0)`, not to `G_C`. By §3.1 its ordinary price is

    C sqrt(d) L^16 = C sqrt(d)(r kappa)^8.

This proves the physical prefactor for the conservative `rho^(8+8k)` row. The two other substantive curl allowances are already

    C r kappa^2 E_D,   C r^2 kappa E_D.

The first comparison roots the unexpanded Riesz chain at H, while `G_C` has a bounded first operator radius. The second keeps the actual selector HS energy `C kappa E_D`; its matrix-test scale is `Cr^2`. Neither comparison spends `||G_C||_2`. The identity Gaussian part cancels exactly in the compensated covariance. The possible curl HS factor also obeys the rank-`2d` estimate in §2. Numerical resolvent/root calibration can conservatively keep `sqrt(n)` and be tightened separately.

### 3.3 Force-eleven and force-twelve forest terms

This is the substantive extra argument missing from a bare appeal to the theorem-as-stated.

The native ideal forward vertex has

    Y_v=a_v M_v I_v+(I−a_v^2 M_v M_v^T)^(1/2)Z_v.

Every `M_v` is the selected first mean of a scalar-affine original D VALUE source. Hence `QM_v=M_v`. Consequently

    Y_v=Z_v + S-valued perturbation,
    range((M_v M_v^T)^j) subset S,  j>=1.

A readout adds fixed scalar Gaussian publics and independent Gaussian fills. Its zero-amplitude linear part is the baseline `LP`. After that baseline is removed, every positive-amplitude coefficient remains S-valued. This remains true after expanding the numerical-gap roots and substituting child polynomials: the outer output slot of each nonconstant forward term stays in S. A reverse representation is used only to identify coefficients; it does not change these same forward endpoint coefficients. A transposed slot is not a new isotropic source. More robustly, at the projected COMPLETE endpoint every nonconstant whole-polynomial coefficient is `Q B_(N,s)`, regardless of the ranges of intermediate reverse helpers. The fixed projector preserves its proper output cut, while the fresh perpendicular Gaussian contributes only to the exact zero-amplitude baseline. Thus the output-cut argument below needs no claim that reverse paths remain in S.

In LOW30 `b27:nc:polynomial` the native proof already supplies, conditionally on the physical old query roots, the proper-cut bounds for each whole Gaussian-polynomial coefficient. The last step there says “isolating the output gives sqrt(n) HS norm.” For this variant its output lies in S, so precisely that step gives `sqrt(d)`, by §2. This applies to the **whole** root-square coefficient; it does not assume individual cyclic Wick contractions are executable trees. The Gaussian degree and the operator moment bounds are unchanged. Thus

    ||B_(N,s)||_p <= L sqrt(d),
    ||D_P B_(N,s)||_(Lp,op) <= L,
    ||D_U B_(N,s)||_(Lp,op) <= L/s.

The root Taylor tail obeys the same estimate directly: it is `M M^T` times a bounded scalar function of `M M^T`, and its action is on an independent fresh Gaussian. Its HS norm is at most `sqrt(d)` times the existing operator remainder. Even without using a helper range invariant, the projected tail satisfies `||Q R(M)||_HS<=sqrt(d)||R(M)||op`; independent fresh-Gaussian isometry supplies its vector bound. The numerical source-field calibration uses its original one-Hilbert factor. Spectral caps can preserve range by singular-value capping, or their whole discrepancy can remain numerical. They do not create an ordinary ambient row.

The same-endpoint Riesz proof (`b27:compiler:forest`) then starts from these physical energy bounds and always multiplies one Riesz Hilbert factor by an endpoint operator derivative. Its centered/main step, finite polynomial base steps, and positive-grade stopping argument never introduce an identity Hilbert norm. It gives physical-dimensional force-eleven currents at their actual ranks and the physical-dimensional force-twelve direct boundary. The tensor-rank cancellation and force grades are unchanged.

One must not instead use full-input Poincare on the original-root derivative `L/s`: that would lose a floor and does not prove this assertion. The conditional whole-polynomial coefficient argument above is the native argument with its output rank localized.

## 4. Complete profiles and induction

Induct simultaneously over complete mean and pair calls, with the complete-return projection built into their definitions. Each new mean uses only older complete projected calls. Each pair then uses the newly constructed mean and older pairs. Their ideal references, variance allocations, signs and replay identities are unchanged.

For every current, channel and direct entry retain its own positive sum of conditional norms of centered original VALUE replays or common-tape differences. Energy-bearing entries keep their supplied physical mark. Ordinary entries use the fixed-S norm in §§2–3, before averaging the prescribed old roots. Conditional Minkowski and the tower property yield `sqrt(d)`, with exactly one energy factor per summand. No smallness conditional on every earlier public is inferred from a Gaussian average.

A complete distinguished spine is idealized before a side field changes. A distinguished input change is multiplied by its actual enclosing affine contraction. For a side change,

    ||(M(Y)−M(Y*))Z||_2 = ||M(Y)−M(Y*)||_(L2,HS)

uses an independent distinguished Gaussian. The existing input-to-matrix-HS Lipschitz estimate is dimension-free and multiplies the incoming `sqrt(d)` error; it adds no Hilbert dimension factor. The gapped-root Sylvester estimate has the same property. A complete return projection is a contraction and cannot enlarge either estimate.

For the actual incoming profile, the original query rows and source Lipschitz radius give the same ordinary slope `L rho`:

    ||e_tau(X)||_2 <= ||e_tau(X*)||_2 + L rho ||X−X*||_2.

Here `X*` is standard only after integrating the same earlier spine publics as in the original proof. Those publics remain identical in the coupling. The output projector never differentiates, conditions on, or replaces the retained incoming label. Numerical profiles keep their absolute path budgets. This is a full current/channel/direct and incoming-profile induction, not merely contraction of a final W2 error.

Every exponent calculation in LOW31 and Exact32 is now unchanged, with `sqrt(d)` factored out throughout. In particular the cap-lift proof for `B=21/2` implies the working `B=52/5`; its unmarked helper check is valid on all descendants `k>1/7`. The 10849-level schedule and all direct shifts remain intact.

## 5. Numerical precision and original-query cost

An ambient root implementation may still have numerical action error `sqrt(n) epsilon`, including a scalar error in its zero-amplitude identity term. Choose `epsilon` below the assigned physical numerical budget divided by `sqrt(n)`, after all weights and inverse floors are enumerated. Since `n<=C d A^(−c)` at fixed source order, this changes `log(1/epsilon)` only by `O(log(d+1)+c log(1/A))`. It does not add a power of `A^(−1)` to query work.

Where a polynomial root on a supported matrix is to preserve an exact scalar baseline, replace `p(x)` by `p(x)−p(0)+sqrt(nu)`. This has the same degree and at most twice the scalar approximation error. Its nonconstant action stays source-bearing. Either this implementation or explicit numerical payment is valid.

Freeze one projected version of each complete marked call before the final precision census, and use that same new named completed program in every matched occurrence and comparison. The wrapper owns a fresh Q-complement Gaussian block after the incoming variable and captured caller have been exposed. Preserve and reuse this block exactly whenever the same complete tape is aliased; redraw it only for an independent complete copy. The block is independent of the retained incoming variable and all captured labels, but is not resampled inside a same-tape difference.

For a captured caller theta, the new program's actual all-private-zero value is `QY(theta;0)`, since its owned complement block is then zero. Compute and charge that value when requested. It need not equal the old ambient zero; its physical projection `P QY(theta;0)=P Y(theta;0)` is preserved. Use the original full-caller/global-origin or portable-anchor centering convention and restore that new finite zero consistently in every matched occurrence. Retain all actual caller paths in numerical budgets; no caller-dependent offset is silently treated as a global constant.

The original raw F, retained source, proxies, same-tape raw subtractions, recorded coisometry, full-caller ports, and their original zero-tape restorations are unchanged. The raw source hierarchy is unaffected because these wrappers are added only after completion of marked priors. The returned marked kernels remain terminal: no first bound for their completed outputs is inferred. On any already admitted actual first-action port v, linear postcomposition has the formal executable rule `D_v W_Q=Q D_vY` and its adjoint first applies Q; the newly owned complement port has first/adjoint Qperp. This only extends an existing legal port and never creates D_xY or D_omegaY for a terminal kernel, nor differentiates a saved HVP through its captured location. Their existing linear norm-growth estimate is preserved by `|QY+(I−Q)Z|<=|Y|+|Z|`, uniformly in rough retained labels.

The extra complete-return maps add fixed known linear operations and one Gaussian block per complete call. They introduce no original gradient/HVP query and no HVP-valued nonlinear source. Existing complete calls and all requested first/adjoint directions still have their original bill. At fixed family depth the additional block count is a polynomial-logarithmic multiplier. Thus the hidden source heat-query exponent is preserved.

The growing ambient dimension still matters for unrelated full-gradient child comparisons and fixed-moment child replacement; their larger-order correction remains necessary. This intrinsic lemma fixes the substantive capped-parent rows instead of trying to absorb them by that unrelated order.

## 6. Exact32 consequence and scope

At the final mixed call `u=3/4`, the parent cap is `rho^B=A^(39/5)` with `B=52/5`. Without this extension the ambient prefactor gives only `sqrt(d) A^(39/5−c/2)`, and its rank-four interior grade is `83/5−c=16.6−c`. This is a real proof loss; increasing a gradient-child order does not change it.

With the intrinsic complete family, the same row is `sqrt(d) A^(39/5)` and has the original interior grade `83/5`. The self-prior cap rows are restored the same way. No algebraic Exact32 row is changed, so its stated covariance minima `(1607/100,141/10,22603/1250)` are preserved, subject to all the other admitted source/proxy/child/numerical contracts. This is not by itself an outer sampler theorem or an all-order chronological theorem.

The conclusion is certified for the explicitly described fixed-output-subspace source and its native affine VALUE descendants. It should be stated as this new intrinsic/projection variant when used; LOW30/LOW31/Exact32 themselves print ambient prefactors. Low pointwise rank alone, a low-dimensional final readout alone, a marginal W2 approximation, or a generic source with unverified proper-cut/native-topology contracts does not justify this extension.

## Source record

- LOW30: `LOW30 (external foundational source; not bundled)`, SHA256 `7767ac389d6b677eb9517b2585f7a24cb9f85ede47d099a653c635b7be750fb8`.
- LOW31: `LOW31 (external foundational source; not bundled)`, SHA256 `3de62d349338304cce43af4e49c8862e8c260b1f5816234d63d6abfbc1d004ba`.
- Exact32: `research/exact-slack/exact32-proof-v2.tex`, SHA256 `c7487ca24ffa44ea73b52e42e4b6f715ad964171fa5bd0a61cf9323eb646cbff`.

The companion diagnostic tests finite-dimensional algebra and the exact cap grade. It is not a numerical implementation of the 10849-level family and does not replace the proof.
