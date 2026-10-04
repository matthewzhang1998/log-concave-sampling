# Independent audit: finite VALUE reverse-OU curvature reserve

Date: 2026-10-04.

## Verdict and admitted scope

**PASS, conditional on the explicitly retained finite LOW30 radius, clock, filter, and numerical guards.** The final source is an actual finite VALUE construction of a Gaussian curvature reserve, not just a covariance identity.

Audited source: `../THIRD-ORDER-VALUE-CURVATURE-RESERVE.md`.

Final SHA256: `56c81676cecd5a76ae073178bc84338e2e9b0cb643788ceff58a554b980adc3d`.

The initially supplied pin was `692296739c0849b0e6a23312660a49a98bdb562864b934281c935a27c00750ee`. Its algebra and rates were already correct. The final version makes the signed-action gap constructive, tightens the raw numerical weight of the balanced lift, and supplies the independent-buffer interpolation proof explicitly. Its newly added diagnostics do not replace any analytical law argument.

The admitted conditional result, for the actual finite raw posterior and its restored state floor eta_Q, is

    W2(Law(R | a,r,t), N(0,(1-Delta)I-Delta H))
       <= Lambda sqrt(D) Delta A^3 s^4
          + C Delta (A/s) eta_Q + absolute finite-compiler floors,

where H is the exact conditional curvature, s²=1-r², 0<Delta=t²-r²<=min(s²,1/4), and 0<A<=1/2. All compiler guards are imposed at the actual scaled sources. The condition A<=1/2 alone is not permission to call a small-radius compiler.

This audit admits neither an exact non-Gaussian reverse transition nor an all-order endpoint recurrence. Higher conditional cumulants remain outside the result. The imported finite LOW30 programs are used as proved ports; this audit does not re-prove or numerically instantiate that entire compiler.

## 1. Literal imports and dimensions

The original source `30_low_acc.tex`, SHA256 `7767ac389d6b677eb9517b2585f7a24cb9f85ede47d099a653c635b7be750fb8`, was inspected directly at lines 180–420. The relevant ports are:

- `t30:eq:seed-data`: actual C1 square VALUE source, complete first ell, square curl ell*a, centered Gaussian energy, and separate caller first.
- `t30:eq:first-response` / `t30:eq:first-contract`: finite shared-root signed responses with Chebyshev filters, not an exact averaged-Jacobian oracle.
- `t30:lem:actions`: the forward covariance action, with actual nested VALUE response formula, positive unexpanded covariance clocks, mean calibration Lambda ell e mu, action energy Lambda ell e, complete first Lambda ell²/sqrt(mu), and caller Lambda ell L_Q/sqrt(mu).
- `t30:lem:value-mean`, especially its conditional-private-action interpolation calculation: error and actual residual/caller firsts with an independent Gaussian buffer.
- `t30:thm:value-pair`: the incoming p is retained in the conditional comparison to N(Mp,I-MM^T), with M=s0 E Dh. Its source-zero carrier is independent of incoming p; it is not the Gaussian chosen by a later law coupling.

The derivative source h is square on R^(2D); p and Y are both 2D-dimensional. The full gradients J_+ and J_- are also on R^(2D), and their covariance actions retain the entire 2D source. P is D by 2D. No D-dimensional fragment of a full-gradient source is improperly passed to the full-gradient action.

No high-order full-gradient pair or full-gradient mean compiler is needed by this reserve. Full-gradient status is used to make the analytical target of the LOW30 covariance action exactly Cov J. The near-gradient pair supplies the different retained-input derivative block.

The raw posterior input is the already audited finite conditional-velocity source, SHA256 `44294f14c4713882e656fbb835e03f90db61dfed224260a28d9501e39cf529c5`. Its finite-mode residual and original linear-tilt restoration remain in eta_Q. The proof does not assume an exact conditional mode.

## 2. Exact target and dimension-safe restoration

For the true conditional density proportional to exp(-|X-a|²/(2s²)-U(X)), integration by parts gives

    H=E Dg-Cov(g)=s^(-2) Cov(X,g(X)).

This matrix is symmetric. The same formula remains valid with an additional constant linear tilt in the potential; its covariance with g is zero. The coupling must nevertheless be to that tilted conditional law, which is why its state restoration stays in eta_Q.

For centered b and arbitrary a,

    ||E[a b^T]||HS <= ||a||2 sqrt(||Cov b||op).

One proof uses the operator from L2 to the coordinates of b: its operator norm is the square root of the covariance operator norm. Apply it separately to the rows of a and sum squares. Consequently, decomposing a covariance difference into Cov(A1-A2,B1)+Cov(A2,B1-B2) gives the source's inequality (1) exactly. There is only one Euclidean energy; a second sqrt(D) is not introduced.

Gaussian Poincare gives Cov(Q)<=C s²I and Cov(g(Q))<=C A²s²I. Conditional Poincare gives Cov(X)<=s²I. Couple Q and X and use the actual pointwise Lipschitz bound on g. Dividing by s² yields

    ||H_Q-H||HS <= C(A/s)W2(Q,X)
                 <= C A³s⁴ sqrt(D)+C(A/s)eta_Q.

This proves true-covariance restoration rather than differentiation of a law error. It uses no regularity rate for the Hessian.

Under the source's covariance convention, Cov(u,f)=E(D_u f)^T. Thus the un-symmetrized covariance has this transpose, while the source's symmetrized identity is correct:

    H_Q=Sym E D_u f-Sym Cov(K,f).

The omitted term is exactly Sym Cov(K,E/s), and is paid using Cov(K)<=alpha²I and ||E/s||2<=C A²s²sqrt(D). Its price is C A³s⁴sqrt(D). It is neither zero nor a derivative-small remainder.

## 3. Full-gradient lift and its complete source graph

The potential certificates in Section 3 differentiate to the displayed G_H and G_B, including all s and r_i factors. U is only an analytical certificate: the producer evaluates g, never U.

The exact equalities P G_H=K and P G_B=b give

    P[Cov(sG_B+G_H/s)-Cov(sG_B-G_H/s)]P^*
       =4 Sym Cov(K,b).

All K nodes and baseline values within one occurrence use the same W. Making these components independent would erase the desired mixed covariance. Independence is introduced only between complete action banks after their full raw sources have been assembled.

The positive dyadic rule has no r=0 node. Its smallest r is at least 1/(4m²), from Markov's polynomial derivative bound, so sum w_i/r_i<=4m². This is a public-log coefficient amplification and finite-precision obligation, not inverse-alpha replication. Balancing gives the actual full radius and energy

    ell_J<=Lambda As,    e_J<=Lambda As sqrt(D).

Without the s and 1/s balancing, the subsequent action allowance would not retain the required s powers.

The full square curl of P^*f includes both the upper-left commutator and the off-diagonal D_v f block. With D_u f=Dg(Q)(I-D_u K), D_v f=-Dg(Q)D_v K, its operator norm is O(A alpha). This is valid for noncommuting Hessians. No conclusion that g(Q) itself is a full gradient on W is used.

## 4. Retained-input pair and heat exponents

For h=Delta P^*f/(s0 v_D), the imported pair supplies the joint comparison with incoming p retained. Its reference has Cov(p)=Cov(Y)=I and Cov(p,Y)=M^T. Therefore

    Cov[sqrt(v_D/2)P(p-Y)]
       =v_D I-Delta Sym E D_u f.

This computation uses the entire joint comparison. Replacing Y by a fresh independent Gaussian would delete the required linear covariance. No forward square is treated as a transposed square.

At the actual scaling, ell_h<=Lambda Delta A, relative curl <=Lambda alpha, e_h<=Lambda Delta A sqrt(D), and mu=alpha=As². The imported pair error is therefore bounded by

    Lambda sqrt(D)[Delta² A³s²
        +Delta^4 A^4+Delta^4 A^(7/2)/s].

Each term is <=Lambda Delta A³s⁴sqrt(D): the respective ratios are Delta/s², Delta³ A/s⁴, and Delta³ sqrt(A)/s⁵. The last two are <=A s² and sqrt(A)s. The inequalities use Delta<=s², rather than an unjustified cancellation of 1/s.

The pair's actual residual first is Lambda[Delta A+Delta² A²/sqrt(alpha)]<=Lambda Delta A. Its normalized caller-a first is Lambda Delta A/s. These are imported graph firsts combined with literal chain rules, not differentiated Wasserstein estimates.

## 5. Signed covariance actions: positive law proof, not mixture matching

For a full gradient J, every smoothed Jacobian A_t is symmetric, so the analytical forward target B_J=int E A_t² equals Cov J. This licenses the exact target matrix, not an exact executed covariance action. The action still pays calibration, private randomness, finite clocks, and numerical error.

Write c=sigma Delta/(8 eta). Conditional on p, center the private action X=C_cov(J;p)-E_private C_cov(J;p). Its private first is at most Lambda ell²/sqrt(mu), and its conditional energy squares integrate to at most the squared full action energy.

The final source's standalone interpolation proof is valid. A centered Gaussian-source image X has a Stein matrix tau satisfying ||tau||L2(HS)<=Lip(X)||X||2. Along zeta Z+t c X, the continuity velocity has bound

    (t c²/zeta) Lip(X)||X||2.

Integrating and then averaging conditional energies proves the Delta² ell³e/sqrt(mu) law allowance. Crucially, it has two powers of c and one energy, and the independent Z buffer stays present. This is a direct smoothed-law comparison; no arbitrary Gaussian mixture is replaced by its mean covariance.

Mean calibration costs Delta ell e mu. After calibration the deterministic linear Gaussian has covariance

    v I+sigma (Delta/4)Cov J+c²(Cov J)².

The positive quadratic term is paid, not canceled. Since ||Cov J||op<=ell² and ||Cov J||HS<=ell e, its HS size is bounded by C Delta² ell³e. Fixed-gap Gaussian square-root stability gives the remaining law term.

The constructive guard Delta ell_J²<=2v_sigma ensures the minus-reference gap using only a declared deterministic source bound; no covariance query is needed. The near-gradient pair and all intermediate finite response/variance guards must additionally hold at their actual normalized arguments.

With ell_J,e_J above and mu=As², the three signed-action errors are

    Lambda sqrt(D)[Delta A³s⁴
       +Delta² A^4 s⁴+Delta² A^(7/2)s³].

Relative to Delta A³s⁴sqrt(D), the last two ratios are Delta A and Delta sqrt(A)/s, bounded by numerical constants using Delta<=s². The complete action residual first is <=Lambda Delta A^(3/2)s, and its normalized caller first <=Lambda Delta A^(3/2).

## 6. Independent joining, gaps, and exact centering

The three completed banks are independent conditional on the captured caller. Their Gaussian comparison covariances add to

    (1-Delta)I-Delta H_star.

Its gap is uniform: the displayed bound ||H_star||op<=A(1+2alpha)<=1 gives a lower bound at least 1/2. The exact target has lower bound at least 5/8 from 0<=H<=AI. Therefore restoration from H_star to H costs the HS matrix discrepancy times C Delta. The total error is exactly the admitted result stated above.

Optional antisymmetrization using two independent complete copies has mean exactly zero and preserves the Gaussian target. For independent identically distributed coupling errors e1,e2,

    E||(e1-e2)/sqrt(2)||²=E||e1||²-||Ee1||².

Thus the coupling allowance does not increase. This is legitimate exact centering without querying or estimating an expectation.

Multiplication by sqrt(v0) multiplies every reserve law and numerical floor by that factor, including C Delta(A/s)eta_Q. Independent joining with the already completed mean service can match the reverse transition's Gaussian mean/covariance reference. That statement supplies no control of its higher conditional cumulants.

## 7. Finite original query, numerical, caller, and zero ledger

A raw f occurrence evaluates the n_alpha K leaves and one terminal g(Q). A raw J occurrence evaluates the n_alpha K leaves and one baseline g(x+su). All reuse the captured g(x) only under the identical recorded key. Thus the source count

    M+1+(n_alpha+1)(N_pair+N_cov,+ +N_cov,-)

is a safe complete count when each N counts every substituted raw-source occurrence inside the finite responses. It cannot be replaced by three abstract calls. New raw arguments require full new K graphs. The pair's h_p substitutions, both signs of each response, all covariance clocks, and discarded-primal replays stay charged.

The private tape likewise includes the complete pair/action/root/clock/fill banks, not just the original two D-roots. Finite filter/clock parameters and tolerances are frozen before differentiation. The polylog claim is restricted to the stated admitted clock regime; log(1/s) remains in the literal bill for tiny positive widths.

Numerically, f has an absolute 1/s normalized g-value weight, while G_H/s has the tighter absolute weight L_Q. These are precision factors. No floor is divided by a possibly zero source energy. Intermediate inverse response widths and the complete outer readouts must be propagated before leaf tolerances are assigned.

Every nonlinear raw consumer source f, K, E, J_+, J_- has a literal zero using same-site recorded-anchor reuse. The uncentered helper values F and B at W=0 are g(x), not zero; the source's zero claim is correctly interpreted for its anchored raw consumer sources. The pair's known source-zero carrier, including retained incoming p in the final readout, remains an executed carrier. It is never identified with an analytical coupling Gaussian.

Actual mode derivatives obey ||D_a x||<=2. Literal differentiation gives D_a K=O(As), the anchored f caller O(A/s), and the balanced J caller O(Lambda A). The imported caller port then produces the bounds recorded in Sections 4–5. With the original physical source caller O(sqrt(A)), the complete reserve caller is <=Lambda Delta sqrt(A)/s<=Lambda sqrt(A)s, plus its separately retained finite-mode and numerical profiles. Only first/adjoint sweeps at original VALUE sites are used; no emitted HVP is differentiated.

## 8. Independent diagnostics and remaining boundaries

The independently written `check_reserve_algebra_independent.py` imports neither the author's checker nor prior diagnostic implementations. It passed **830 assertions**, seed 620041004. The saved JSON records:

- Exact dyadic-rule mass and first moment, r_min and lift-coefficient guards.
- 48 anisotropic rotated-matrix quadratic cases in dimensions 2,5,13,31, including widths down to 1e-4; full balanced polarization and retained-pair readout identities.
- Exact one-energy covariance inequalities on arbitrary Gaussian couplings, through dimension 64.
- Every displayed A,s,Delta exponent inequality down to s=1e-8.
- A different genuine C2/non-C3 fixture with scalar Hessian min(sqrt(|t|),1), three nonorthogonal ridges, and 45 raw tape tests. Directional finite differences check the VALUE graph's firsts; full lift Jacobians are symmetric, while source Hessians demonstrably do not commute. The maximum sampled Hessian commutator norm was about 0.00396.
- A negative control showing that a random Gaussian covariance mixture has a different fourth moment from the Gaussian with its mean covariance.

The maximum raw first finite-difference discrepancy was about 2.31e-8. These are diagnostics, not substitutes for the uniform C2 proofs. They do not instantiate the huge imported finite LOW30 action/pair programs.

The author's separate diagnostics and final Section 8 were inspected. Their exact-target integration uses the actual reweighted conditional law, and their covariance orientation and shared-root products are consistent. Their finite integration ratios are explicitly diagnostic, and no Hessian quadrature accuracy is silently used to prove the law bound.

No substantive construction gap remains within the stated guarded finite-order Gaussian-reserve scope. A full positive reverse-transition sampler still needs control of the genuine higher conditional cumulants and its own complete composition/cost theorem.
