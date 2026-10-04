# Independent audit: conditional OU velocity VALUE source and buffered law service

Date: 2026-10-04.

## Verdict

**PASS for the repaired source, strictly within its stated conditional, fixed-buffer, fixed-order scope and the explicit LOW30 consumer guards.**

Audited source SHA256:

`44294f14c4713882e656fbb835e03f90db61dfed224260a28d9501e39cf529c5`

The original submitted pin `a83ee91c695f5182ff1cb487bed9bd0e600a184306b9f2b9dba2b3c0ee230ed3` omitted a literal imported small-radius gate and suppressed the imported public-log factors. Its raw VALUE source was valid, but its unqualified completed-consumer claim was not admitted. The repaired bytes explicitly add the gradient-compiler threshold, a near-gradient covariance gap, the remaining finite-clock/filter guards, and the public-log factor. They also make restoration of the original global finite-mode linear tilt explicit.

For the repaired source the admitted bound is

    W2(Law(M_r | z), N(m_r(z), I))
       <= C A^3 s^5 sqrt(D) + Lambda_k A^4 s^3 sqrt(D)
          + A |R_M| + absolute numerical floors,
    s=sqrt(1-r^2).

It implies the source's slightly looser displayed bound with Lambda_k multiplying both intrinsic terms. The true original finite-mode posterior drift additionally uses the known shift ell and costs `A s^2 |ell|`. No exact expectation, conditional posterior sample, density, potential value, covariance matrix, or Hessian action is part of the VALUE producer.

This PASS does not supply a small-noise mean oracle, a strongly accurate random estimate, a joint law retaining private source roots, a native retained/protected/proxy certificate, a same-endpoint posterior sampler, or a complete sublinear order-to-cost recurrence.

## 1. Source pins and literal imported ports

All of these bytes were inspected and their hashes verified:

- Conditional service: `../THIRD-ORDER-CONDITIONAL-OU-VELOCITY-LAW-SERVICE.md`, hash above.
- Positive endpoint theorem: `../../POSITIVE-SECOND-ORDER-LAW-AND-QUADRATIC-AMPLIFIER.md`, `a085dfdef3f66f69208612034c17e4e18045fc00a21be3aed64ccfcb5e9b7085`.
- Its independent audit: `../../independent-audit/INDEPENDENT-POSITIVE-LAW-AUDIT.md`, `df32f63c88ec075180b231bd997f13524792296aaa161b1879c987bfcdd28512`.
- Complete-mean recipe: `no-copy-rank-20261004/host-observer-audit/COMPLETE-MEAN-FIRST-ORDERING-AND-REENTRY-RECIPE.md`, relative to the research root, `a7cde258f1c86defbfdc681f473076293cfe65253479493c6ef5baa8b08a88a0`.
- Its independent audit: `INDEPENDENT-COMPLETE-MEAN-REENTRY-AUDIT.md` in the same directory, `d9b4d45c276e72adc4b2590c1c17140bed2fe94dd2a6de5d8f7cc1b3a0a7d21a`.
- Original LOW30: `LOW30 (external prerequisite)`, `7767ac389d6b677eb9517b2585f7a24cb9f85ede47d099a653c635b7be750fb8`.

The decisive original LOW30 ports are:

1. `t30:eq:seed-data`, approximately lines 202–210: actual square C1 VALUE source, complete first ell, square curl ell*a, centered Gaussian energy, and separate caller first.
2. `t30:lem:value-mean`, lines 315–369: positive nine-replay VALUE mean, error `Lambda E {ell(a+mu)+ell^3(1+mu^(-1/2))}`, actual residual first `Lambda(ell+ell^2 mu^(-1/2))`, caller `Lambda(1+ell mu^(-1/2)) L_Q`, and known source-zero carrier.
3. `b27:compiler:mean`, lines 5704–5735: the exact full-gradient VALUE mean, its `Lambda_k sqrt(n) rho^k` error, frozen original labels, and explicit threshold `rho<=r_*(k,Delta,Lambda_k)`.
4. Lines 6458 onward: the threshold is imposed after every variance/readout normalization; no source outside it is admitted.
5. `b27:compiler:paths` and `b27:compiler:serial`: actual complete-graph first/caller/path bounds and full recursive query expansion.

The complete-mean recipe is therefore not being used to invent a new mean oracle. Its actual interfaces match those printed in LOW30. This audit imports those established compiler theorems; it does not re-prove the entire high-order LOW30 compiler.

## 2. Conditional normalization and finite mode

Fix z,r and all finite versions before sampling any private roots. For `0<=r<1`, set `a=rz`, `alpha=As^2`. The map `x -> a-s^2 g(x)` is an alpha-contraction. Its frozen M-step iteration and the recorded final anchor give

    R_M=x_M-a+s^2 g(x_M),
    |R_M|<=s^2 alpha^M |g(a)|,
    ||D_z x_M||<=r/(1-alpha),
    |x_M|<=r|z|/(1-alpha).

These follow from the actual iteration and need no derivative of its convergence error. In particular the residual remains caller-dependent, since `|g(a)|<=Ar|z|`.

The anchored conditional field

    gbar(xi)=s[g(x_M+s xi)-g(x_M)]

is an exact gradient with `gbar(0)=0` and Hessian between zero and `alpha I`. An analytical potential is

    Ubar(xi)=U(x_M+s xi)-U(x_M)-s g(x_M).xi.

After the change of variables `q=x_M+s xi`, the true conditional density has potential

    |xi|^2/2 + Ubar(xi) + (R_M/s).xi + constant.

Thus the sign and size of the finite-mode tilt in the source are correct. The positive endpoint theorem, with quadrature error at most alpha, gives standardized state error `C alpha^2 sqrt(D)`. Multiplication by s gives `C A^2 s^5 sqrt(D)`. Strong convexity bounds the standardized linear-tilt restoration by `|R_M|/s`; after physical conditional scaling it is exactly at most `|R_M|`. This division by s is analytical only.

Consequently the actual endpoint `Q_fin=x_M+s(u-H_alpha(u,v))` satisfies

    W2(Law(Q_fin | z), nu_(r,z))
       <= C A^2 s^5 sqrt(D)+|R_M|+e_state.

Every node of H_alpha retains the same u,v. It is not a mixture of independent node-root laws or a sampled conditional expectation.

## 3. True mean target and covariance signs

The executed identity is

    F=g(Q_fin)=B+E,
    B=g(x_M+s u),
    E=g(Q_fin)-g(x_M+s u).

Transporting the preceding state coupling through the A-Lipschitz g proves

    |E F-m_r(z)|
       <= C A^3 s^5 sqrt(D)+A|R_M|+e_force.

The right comparison is the true conditional mean under the tilted posterior density. It is not the mean of an arbitrarily fitted Gaussian or Gaussian mixture. The exact posterior and its expectations occur only in the proof. This argument does not require a pointwise Taylor remainder for g or a continuity rate for Dg.

The mode-translation correction to `P_r g(z)` is exactly the mean of

    g(x_M+s u)-g(a+s u)+E(u,v).

The translation chord is a genuine gradient in u, because both of its u-Jacobians are symmetric. Its finite caller dependence is retained rather than absorbed into a uniform sqrt(D) quantity.

The analytical amplitude-two identity also has the correct sign and factors:

    -Cov_gamma(U(a+s u),g(a+s u))
      =-s^2 integral_0^1 E[Dg(a+s u)
                         g(a+s(tu+sqrt(1-t^2)v))] dt.

Apply scalar/vector Gaussian covariance integration by parts to U(a+s u) and each component of g(a+s u). Each Gaussian derivative supplies one s; the Hessian factor is symmetric, so the displayed ordering is valid. There is no extra r and no factor t. This identifies a coefficient, not an executed HVP-valued source. The full mean estimate above is stronger and avoids expanding that coefficient.

As a separate sign check, differentiation of the exact conditional mean gives

    D_z m_r=r[E_nu Dg-Cov_nu(g)].

This is consistent with the upstream OU continuity velocity `-m_r`. It is not a producer-side derivative oracle.

## 4. Actual private first, square curl and physical energy

Write `B_1=Dg(Q_fin)`, `B_0=Dg(x_M+s u)`, `H_u=D_u H_alpha`, `H_v=D_v H_alpha`. The literal chain rule is

    D_u E=s[(B_1-B_0)-B_1 H_u],
    D_v E=-s B_1 H_v.

The two B matrices are symmetric and lie between zero and AI, so `||B_1-B_0||<=A`. Also H_u is symmetric, `0<=H_u<=alpha I/2`, and the complete H first has norm at most alpha. Therefore

    ||D_(u,v) E||<=As(1+alpha).

For the coisometry `P_u=(I,0)`, the square lift is `P_u^* E=(E,0)` on R^(2D). Its derivative skew has upper-left block

    s(H_u B_1-B_1 H_u)

and upper-right block `-s B_1 H_v`, with the lower-left negative transpose. Bounding the commutator without commuting its factors gives, for example,

    ||D(P_u^*E)-D(P_u^*E)^*||<=2 A^2 s^3.

The leading difference `B_1-B_0` is symmetric and contributes no curl. No smallness of that Hessian difference is assumed. This is exactly why the private first can be order As while the square curl has order A^2 s^3.

The pointwise force-chord estimate is

    |E|<=As|H_alpha|.

Every inner argument `r_i u+s_i v` is marginally standard Gaussian, irrespective of cross-node dependence. Positivity and total weight one give `||H_alpha||_p<=C_p alpha sqrt(D)`, hence

    ||E||_p<=C_p A^2 s^3 sqrt(D).

There is one physical Gaussian energy, not two factors sqrt(D). For u=v=0 the inner packet vanishes, and both force arguments equal x_M. Reusing its recorded original gradient gives E(0,0)=0 exactly, including in the stated coherent numerical convention. Recomputing the same site with separately perturbed values would instead create a numerical zero debt.

B is an exact gradient in its own D-dimensional u variable: `B=grad_u[U(x_M+s u)/s]` for s>0. Subtracting its recorded B(0) preserves exact gradient status. Its complete first is As and its anchored energy is at most `C_p As sqrt(D)`. The finite conditional mode is a captured constant for this gradient statement, so no approximate-implicit-gradient restoration is required here. The complete F or E on (u,v) is not being declared a full gradient.

## 5. Caller bounds and finite original derivatives

With z exposed before private roots,

    D_z gbar(xi)=s[Dg(x_M+s xi)-Dg(x_M)]D_z x_M,
    ||D_z H_alpha||<=2sAr.

The actual B and E caller firsts are therefore O(Ar). The source does not claim a derivative gain from E's order-A^2 VALUE energy.

For a physical caller y, the global finite mode record has `||D_y b_M||<=1/(1-A)`. At fixed normalized q,

    D_y g(q)=sqrt(A)[D^2 V(b_M+sqrt(A)q)-D^2 V(b_M)]D_y b_M.

This costs O(sqrt(A)), by the Hessian sandwich alone. The conditional recurrence then gives `D_y x_M=O(s^2 sqrt(A))`, and the actual H, B, E chain rules yield the source's O(sqrt(A)) physical-caller bounds. All explicitly captured additional labels must be included separately, as stated in the source.

These derivatives are original HVP/adjoint sweeps at already recorded original-gradient VALUE sites. Dense Jacobians are analytical notation, not one oracle query. No HVP output is differentiated. Finite VALUE restoration never proves convergence of Hessians: numerical first/adjoint actions retain the original-first accuracy contract and the actual finite graph's own radius bounds.

## 6. Positive complete-output consumer join

The B branch meets LOW30's full-gradient port on its D-dimensional tape. Its normalized radius is `rho_B=sqrt(2)As`. The repaired source explicitly requires `rho_B<=r_*(k,Delta_B,Lambda_k)`; A<=1/2 alone is insufficient. Its error after the fixed variance rescaling is

    Lambda_k sqrt(D)(As)^k+e_B.

For E one can make the rectangular adapter explicit: execute the square-source mean on `(E,0)/sqrt(v_E)` on R^(2D), then project by P_u and multiply by sqrt(v_E). Projection of its unit Gaussian target has exactly covariance I_D before the variance multiplier. The source energy is unchanged by the isometric lift. The dimension is 2D at this raw port, not an unpriced infinitely large tape.

The normalized complete first is at most `sqrt(2)As(1+alpha)`. The repaired numerical guard at 1/4 supplies a covariance gap: the nine independent main replays have `sum_i w_i^2=1/3`, so Gaussian Poincare gives

    Cov(S)<=ell_E^2 I/3,
    I-t^2 Cov(S)>=(47/48) I.

The finite clock/filter/numerical guards remain those of the imported circuit, at the actual normalized source and output share. Applying the near-gradient theorem with a declared first O(A), relative curl O(A), mu=A and energy `e_E<=C A^2 s^3 sqrt(D)` yields

    Lambda_k e_E [A(A+A)+A^3(1+A^(-1/2))]
       <= Lambda_k A^4 s^3 sqrt(D).

No actual source energy is queried to choose the program. Its energy parameter is an analytical bound, and absolute numerical floors are never divided by a possibly zero e_E.

Each branch is completed on an independent complete bank, while every occurrence of E internally preserves its same u,v, packet ancestors and terminal/base values. Conditional independent convolution adds their target means to EF and their positive variances to one. For k>=4, `(As)^k<=A^4 s^3`, proving equation (8) of the source and then its result (R).

The two Gaussian comparison targets are analytical couplings. They are not Gaussian draws substituted into the executed programs, and they are not identified with any actual source-zero carrier. Only the completed output and captured labels can be passed through a subsequent law-only comparison. Native retained/protected/proxy reentry is a different graph obligation.

## 7. Endpoints, global finite-mode tilt, domain and complete bill

At r=0, x_0=0 and every conditional mode iterate is exactly zero since g(0)=0. The conditional mean can nevertheless be nonzero; the force-chord service does not erase it. The z-caller derivatives vanish, consistently with the displayed r factors.

At r=1 use the explicit endpoint branch `g(z)+N`. It has exact anchored target mean and one original VALUE after the global anchor is available. Do not evaluate log(1/alpha), construct the positive quadrature at alpha=0, or divide by s there. As r approaches one through positive widths, the law powers and caller bounds remain valid; the node count retains its actual dependence on log(1/s).

If the global finite mode leaves linear tilt ell, the anchored target differs from the original normalized target. At fixed r,z the conditional physical-q Hessian is at least `s^(-2) I`. Adding ell.q moves its law by at most `s^2 |ell|`, so its g mean moves by at most `A s^2 |ell|`. Adding the known deterministic ell to M_r restores the constant part of the true drift. This is an explicit caller-dependent floor and requires no posterior expectation query.

With the already priced global mode/anchor supplied, the incremental service count is

    M+1 + N_B +(n_alpha+2)N_E + actual numerical/known-row work.

Starting instead from V,y requires adding the global finite-mode and anchor bill. Each g VALUE after that anchor is one original gradient query. A raw E occurrence uses all n_alpha original inner values plus its terminal and baseline values; its captured conditional anchor can be shared. Different raw arguments created by the compiler require their complete ancestor replays. Semantic same-site coincidences permit only exact-key reuse.

The literal private tape is proportional to D times the complete source occurrence count, plus all consumer clock/fill blocks. It is not 2D for the completed service. Every requested original HVP direction and every discarded-primal replay is separately charged.

The inner quadrature count is `O(log^2(1/(As^2)))`. Zero inverse-A exponent follows only under the supplied fixed-order public-log outer-clock hypothesis. Exponentially tiny widths are not free. The original-gradient query domain remains all Euclidean space with the recorded finite caller/displacement envelope; no algebraic D-versus-A guard or bounded-Gaussian-support fiction is introduced. Every finite parameter is frozen before differentiation.

## 8. Independent checks

`check_conditional_velocity_independent.py` imports neither the author's tests nor upstream diagnostic programs. It passes **1,685 assertions** with seed 620041027. The manifest pins it, the output JSON, this audit and all source inputs.

Coverage includes:

- 72 actual finite physical-source cases in dimensions 2, 5 and 13, A in {.5,.16,.035}, r in {0,.4,.93,.999}, and conditional mode depths 0 and 4.
- A genuinely C2, non-C3 fixture with continuous non-Lipschitz Hessian: nonorthogonal ridges have scalar Hessian `min(sqrt(|t|),1)`. Its global Hessian sandwich follows from the known ridge Gram norm, and sampled Hessian commutators are nonzero.
- Exact private first, full square-lift curl, mode residual/caller, physical caller, query-census and same-site zero checks, with independent directional finite differences in u,v,z,y.
- 40 scalar asymmetric exact-target quadrature cases, at two integration resolutions, checking the full conditional mean against the actual shared-root force graph. The maximum observed ratio `|EF-m_r|/(A^3 s^5)` is approximately 0.05262, with residual and resolution floors recorded separately.
- Direct covariance-integration checks for equation (3), exact conditional-mean derivative sign/r-scaling checks, and a negative control showing that independent node roots change the shared-root covariance.

The maximum observed square-curl ratio to `A^2 s^3` is approximately 0.485. These finite values are diagnostics, not universal constants or proofs. The analytical argument supplies the C2-uniform and dimension-safe statements. The checker does not execute LOW30's high-order mean compiler and does not claim to enumerate its complete graph.

## Final admitted boundary

The repaired construction supplies an actual reusable second-order shared-root force chord, and the imported guarded positive consumers turn it into a third-order conditional unit-buffer Gaussian-law service. The force mean is tied to the exact posterior conditional target by a state-law coupling, with both conditional and global finite-mode debts restored explicitly.

The next same-endpoint composition, smaller output variances, retained private observers, and all-order general-potential cost claims remain unproved. They cannot be obtained by subtracting a coupling Gaussian or treating this buffered output as a strongly accurate conditional mean.
