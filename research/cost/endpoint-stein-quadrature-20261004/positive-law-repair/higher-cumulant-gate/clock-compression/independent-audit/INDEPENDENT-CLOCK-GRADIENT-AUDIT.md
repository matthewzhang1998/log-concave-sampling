# Independent audit: protected OU-clock gradient and conditional covariance quadrature

2026-10-04. Independent source and theorem audit. No author files were changed.

## Verdict

**PASS for the finite protected-endpoint full-gradient lift, its actual first/origin/one-energy bounds, explicit PSD calibration, and the stated direct-predictor cosine separator.** The mathematical source is admissible for a frozen-endpoint full-gradient mean/pair port after its actual variance/readout normalization and the imported small-radius guard. This does not by itself construct a positive covariance-compensated path law, preserve private observers through a law-only comparison, or remove the printed consumer's ambient dimension factor.

A further independently verified analytical identity simplifies the proposed conditional covariance quadrature:

    v(z) = integral_0^1 P_r g(z) dr,
    B(z) = Dv(z) = integral_0^1 r P_r[Dg](z) dr,
    C_cont(z) = Cov(integral_0^1 g(X_r) dr | X_1=z)
              = 2 integral_0^1 q P_q[B^2](z) dq.

This gives a positive finite analytical covariance approximation with error at most `(3/2) delta A^2 sqrt(D)` in `L2(gamma;HS)`, using the existing uniform-Hermite positive rule. No dimension-dependent accuracy reduction is needed for that result. Expectations and Hessians in this identity are analysis, not newly authorized producer queries.

## 1. Inspected sources and pins

- `../POSITIVE-CLOCK-CALIBRATION-AND-DIAGONAL-GATE.md`: SHA256 `003476b645b0d83f49ece98d8a7d4e1bbd09cede992859e4d7714f950e46376d`.
- Original LOW30, `external:LOW30`: SHA256 `7767ac389d6b677eb9517b2585f7a24cb9f85ede47d099a653c635b7be750fb8`.
- The endpoint positive dyadic quadrature theorem and the existing complete-mean/caller recipe were read as supporting context. The checker imports none of their code or the author's checker.

The literal LOW30 interfaces inspected are `b27:compiler:mean` and `b27:compiler:pair` around lines 5696–5740, the variance/readout admission paragraph around 6458 onward, the complete source-path bounds `b27:compiler:paths`, and the serial occurrence expansion `b27:compiler:serial`. The covariance-only construction and `b27:covariance:module` were also inspected. These remain imported theorems, not re-proved or executed here.

## 2. Brownian normalization and the exact readout

Order distinct clocks `1>r_1>...>r_N>0`, and let

    q_i=(r_i^(-2)-1)/2, q_0=0, Delta q_i=q_i-q_(i-1),
    L_ij=sqrt(2) r_i sqrt(Delta q_j) 1_{j<=i}.

For `i<=k`,

    (LL^T)_ik=2 r_i r_k q_i=r_k/r_i-r_i r_k.

Thus `r_i Z+L_i G` has the claimed OU covariance, standard marginal, and endpoint covariance `r_i I`. The factor `sqrt(2)` and the clock's factor `1/2` are consistent; dropping either without changing the other is wrong.

The stable equivalent readout formula is

    t_j=sqrt(Delta q_j)/(sqrt(q_j)+sqrt(q_(j-1))).

It equals the displayed difference-quotient formula. Telescoping gives

    L_i t=sqrt(2) r_i sqrt(q_i)=sqrt(1-r_i^2)=c_i.

The first squared component is one. Every later squared component is

    t_j^2 = (sqrt(q_j)-sqrt(q_(j-1)))/(sqrt(q_j)+sqrt(q_(j-1)))
          = tanh((log q_j-log q_(j-1))/4).

Consequently `||t||^2 <= 1+(1/4)log(q_N/q_1)`. For the stated dyadic rule the terminal midpoint is the largest r clock, so `q_1` is comparable to the terminal panel length. The smallest r is the first-panel Gauss endpoint node, yielding the stated polynomial dependence of `q_N` on Gauss order. Alternatively the actual finite logarithm can be computed directly and used as the public guard.

This is an exact finite covariance/readout identity, not an inverse-heat formula or a comparison to an unexecuted path draw.

## 3. Executed full gradient and numerical normalization

At frozen endpoint z,

    Psi_z(G)=sum_i (w_i/c_i) U(r_i z+L_i G),
    J_z(G)=sum_i (w_i/c_i) L_i^T g(r_i z+L_i G).

The private variable and output both have dimension `N D`, or `M D` for a concatenated calibrated factor with M scalar private roots. This is a genuine full private gradient. An executed source occurrence makes N original gradient VALUES; it makes no U-value or Hessian-valued producer query.

A numerically cleaner equivalent expression is

    u_i=L_i/c_i, ||u_i||=1,
    J_z(G)=sum_i w_i u_i^T g(r_i z+c_i u_i G).

For the OU factor, `u_ij=sqrt(Delta q_j/q_i)1_{j<=i}`. Thus the apparently large `1/c_i` in the analytical potential need not appear as an uncompensated VALUE multiplier. Absolute source-VALUE path mass is one. Floating-point errors in clocks, rows and coefficients still need their ordinary absolute budget.

The readout is exactly

    (t^T tensor I_D)J_z(G)=sum_i w_i g(X_i)=H_Q.

This source is a gradient in G conditional on z. It is not automatically a full joint gradient after adjoining arbitrary caller coordinates while leaving its output unchanged.

## 4. Actual firsts, origins, and one physical energy

Write `A_i=Dg(X_i)`, so `0<=A_i<=A I`. Literal differentiation at recorded original VALUE sites gives

    D_G J=sum_i (w_i/c_i)(L_i^T tensor I) A_i (L_i tensor I),
    D_z J=sum_i (w_i r_i/c_i)(L_i^T tensor I) A_i.

Therefore

    0<=D_GJ<=A(sum_i w_i c_i)I<=A I,
    ||D_zJ||<=A sum_i w_i r_i=A/2.

There is also the useful raw Hilbert first bound

    ||D_GJ||_HS <= A sqrt(D) sum_i w_i c_i.

This follows term by term from rank at most D, not from the ambient output dimension. It is additional source information, not a replacement for an imported theorem's stated error factor.

The caller-only origin is

    J_z(0)=sum_i w_i u_i^T g(r_i z),
    ||J_z(0)||<=A|z|/2.

It needs N caller-only original VALUES, which may be cached once for the same complete frozen caller and finite version. It is generally nonzero. At `z=G=0` it is exactly zero under the host's coherent anchored same-site convention. Independent perturbed evaluations would instead incur a numerical zero debt.

Centering `Jc_z(G)=J_z(G)-J_z(0)` preserves the private-gradient property, with potential `Psi_z(G)-J_z(0) dot G`. Its private Hessian is unchanged. Its actual caller derivative is

    D_zJc=sum_i w_i r_i u_i^T [Dg(X_i)-Dg(r_i z)].

Because both Hessians lie between zero and A I, their difference has operator norm at most A. Hence the sharper centered bound is still `||D_zJc||<=A/2`. The origin paths must be retained to obtain this formula; they cannot be silently dropped because an origin was subtracted.

For every fixed p, conditionally on every z,

    ||Jc_z(G)||_Lp <= sum_i w_i ||g(X_i)-g(r_i z)||_Lp
                   <= C_p A sqrt(D) sum_i w_i c_i.

Every `L_i G` is one D-dimensional Gaussian with covariance `c_i^2 I`. The triangle inequality is applied to the literal weighted vector terms. There is no `sqrt(ND)` charge in this source energy. At A=0 all these source quantities and anchors vanish without dividing by A or by the actual source energy.

Additional physical callers inherited through the finite mode must be differentiated through their own actual original-VALUE/anchor paths. The displayed z-first is not a certificate for unspecified physical labels.

## 5. Calibrated concatenated factor

With the explicit independent or common-root auxiliary factor, let

    Lcal=(sqrt(1-eta)L_OU, sqrt(eta)L_aux),
    tcal=(t/sqrt(1-eta),0).

The calibration bracket is strict, so `0<=eta<1`. Its rows still have norm c_i, and

    Lcal tcal=c,
    ||tcal||^2=||t||^2/(1-eta).

Thus all of Sections 3–4, including source energy and exact gradient status, remain unchanged. Only the readout norm and literal private-root count change. Under the extra public guard `eta<=1/2`, the norm-squared bound worsens by at most two.

That guard is substantive. For `r=(.53,.51,.49,.47)` and equal weights, the calibration bracket holds but `eta=0.8782834769264231`. One cannot deduce the factor-two claim from the bracket alone. The source explicitly states the guard, and all fine rules checked here satisfy it.

Both auxiliary branches were independently checked. The simple clocks `(.99,.5,.01)` with weights `(.4,.2,.4)` have `S=0.4472404040404041<1/2` and use the one-common-root branch with `eta=0.2648961106697252`. The usual dyadic rules use the independent auxiliary branch. No diagonal ridge was added to the audited explicit factor.

## 6. Conservative join to the printed full-gradient ports

Let `T=t^T tensor I_D` (or its calibrated version), `beta=||t||`, `P=T/beta`, so `PP^T=I_D`. Put `H0=TJ_z(0)`. For a declared positive variance share v, the conservative mean recipe is

    f_z=(beta/sqrt(v)) Jc_z,
    Y_z=H0+sqrt(v) P N_B(f_z).

Its analytical target is `N(E[H_Q|z],v I_D)`. With

    n_G=M D,
    rho <= beta A(sum_i w_i c_i)/sqrt(v),

the literal printed theorem gives the safe error

    sqrt(v) Lambda_B sqrt(M D) rho^B + restored absolute floors.

Require its fixed-order threshold `rho<=r_*(B,Delta,Lambda_B)` after all positive reserve shares, readouts, filters and finite queues are enumerated. The raw source's `sqrt(D)` energy does not alone authorize changing `sqrt(M D)` in this imported statement. Since M and beta are public logarithmic quantities for the admitted rule, the conservative extra factor remains logarithmic, but it must be written and counted.

The mean target integrates the source-private G roots while keeping z and the original captured callers fixed. One may subsequently apply only a complete-output/captured-caller consumer justified by its real Lipschitz and law ports. Appending the old private roots, an old H, or a different output that still reads them is not licensed by this marginal law comparison.

Likewise the pair port targets its declared selected mean Jacobian at the normalized gradient source. It is not a covariance-matrix oracle. A covariance reserve needs the original finite clock/gradient-only square construction, a positive known variance gap, its complete output comparison, and its actual replay/first/zero ledger. The complicated specialized outer covariance theorem cannot be imported with its numerical grades merely by renaming J as its F.

There is a potentially cheaper separate route: LOW30 also states an admitted source class `f:R^(k n)->R^n` whose n-dimensional partial Jacobians are symmetric. Directly,

    D_{G_j}H_Q=sum_i w_i L_ij Dg(X_i)

is symmetric, and this remains true for scalar concatenated calibrated factors. However k=M grows with the clock rule, while the source text allows constants depending on fixed block counts. A direct D-output consumer avoiding the square-lift ambient factor needs that block-count dependence priced explicitly. This audit does not silently supply that uniform extension.

## 7. Positive calibration and cosine separator

The explicit mixture preserves endpoint correlations, unit outer marginals and PSD conditional covariance, while imposing `w^T Rcal w=1/2`. Appending independent local Markov innovations produces each named endpoint/outer/inner Gaussian triple exactly. This deliberately does not preserve the original whole path law.

For `g(x)=Bx`, the complete literal Gaussian covariance has

    Cov(Y)=I-B+B^2-2 Cov(H0,L0) B^3+Var(L0) B^4.

The coefficient of B squared is exactly one for arbitrary symmetric B; no diagonalization or omitted shared-root cross term is needed.

The smooth fixture `f_k(x)=a x+(b/k)(cos(kx)-1)` has derivative in `[1/4,3/4]` at `a=1/2,b=1/4`. Its centered cosine covariance is nonnegative at every correlation. Thus every standard-marginal PSD calibration retains the positive self-diagonal, independently of off-diagonal signs. The debt and current calculations in Sections 4 and 7 of the source are correct.

To remove typographical ambiguity, every isolated `e^-u/2` in its current/target formulas (15), (16), (19) means **exp(-u/2)**, not exp(-u)/2. The code uses the correct expression. Independent Gaussian integration verifies

    E[Z f_k'(X)f_k(Y)] = r[a^2 tau+b^2 F_u(tau)],
    F_u(tau)=exp(-u/2)-exp(-u)cosh(u tau)
                          +tau exp(-u)sinh(u tau),

and the centered target variance coefficient

    a^2+b^2[exp(-u/2)+(exp(-2u)-exp(-u))/u].

Subtracting the squared first-order mean is essential for the target coefficient. With `u=4/s2` and the specified inner accuracy, the executed nested predictor has coefficient debt at least `s2^2/256>=1/(256N^2)`. The family-specific remainder and the condition `A N^(5/2)->0` justify the polylogarithmic-N finite-A separator. Nothing here proves an unrestricted lower bound for arbitrary positive law compilers or arbitrary power-growing N.

## 8. Audit of the proposed L4 pair-quadrature argument

On a dyadic distance-to-one interval, the Bernstein ellipse of parameter two lies in the complex unit disk in the r plane. Hermite spectral calculus makes `r -> P_r` a bounded holomorphic operator on `L2(gamma;K)` for any Hilbert K, with norm at most one on that ellipse.

Let I_d be real-node Chebyshev interpolation and `E_r=P_r-I_dP_r`. Uniformly over the real panel,

    ||E_r||_(2->2) <= C 2^(-d),
    ||E_r||_(infinity->infinity) <= 1+Lambda_d,
    Lambda_d=O(log(d+1)).

The second line is a real Mehler contraction estimate. Bochner-space Riesz–Thorin gives

    ||E_r||_(4->4) <= C sqrt(1+Lambda_d) 2^(-d/2)=epsilon_d.

For vector-valued g,h, the degree-2d product interpolant obeys

    ||(P_rg) tensor (P_rh) - (I_dP_rg) tensor (I_dP_rh)||_L2(HS)
        <= (2 epsilon_d+epsilon_d^2)||g||_4||h||_4.

Multiplication by r raises the degree to 2d+1, exactly integrated by Gauss order `m>=d+1`. Positivity of the quadrature gives the same error times at most twice panel mass, and the terminal panel is bounded directly by its mass. Hence the proposed `epsilon A^2 D` relative-energy estimate is valid. Taking `epsilon=delta/sqrt(D)` only adds log D to the quadrature count.

For a bounded Hessian matrix field a stronger route is available. Its L2(HS) interpolation error is `O(2^(-d) A sqrt(D))`, while its real-node interpolant has L-infinity operator norm at most `Lambda_d A`. Multiplying an HS error by the bounded other factor gives `O((1+Lambda_d)2^(-d)A^2 sqrt(D))`, without L4 interpolation. If directly discretizing the proposed s-cubed-weighted product of two degree-d factors, use `m_s>=d+2`; the previous `d+1` degree count applies only to a linear weight.

## 9. Conditional covariance identity and a simpler positive square quadrature

Let `H=Dg`, and retain the actual continuous OU path conditional on `X_1=z`. On the ordered pair `X_r,X_(r tau)`, its conditional cross-covariance is

    P_r[g(P_tau g)^T](z) - (P_rg)(z)(P_(r tau)g)(z)^T.

The Mehler carré-du-champ identity, with `D(P_tau g)=tau P_tau H`, expresses this as

    2 tau integral_r^1 s ds
        P_(r/s)[(P_sH)(P_(s tau)H)^T](z).

Multiply by the ordered-pair measure `2r dr d tau`, symmetrize, exchange the triangle `0<r<s<1`, and set `q=r/s`. Since `r dr=s^2 q dq`, this is exactly

    C_cont(z)=4 integral_0^1 s^3 ds integral_0^1 q dq integral_0^1 tau d tau
                   P_q[Sym((P_sH)(P_(s tau)H)^T)](z).

Only bounded continuous H is required. One may first use smooth approximation to justify the semigroup product rule, then pass by bounded Hessian domination; no modulus of continuity or third derivative enters the result.

The two inner Hessian factors are symmetric but need not commute. Their symmetrized product is not termwise PSD. Nevertheless, integrating the ordered triangle factors it into an exact square. Set

    B(z)=integral_0^1 s P_sH(z) ds=Dv(z).

Then

    B^2=2 integral_0^1 s^3 ds integral_0^1 tau d tau
                 Sym((P_sH)(P_(s tau)H)^T),
    C_cont=2 integral_0^1 q P_q[B^2] dq.

Pointwise `0<=B<=A I/2` and `0<=C_cont<=A^2 I/4`.

Now take independent positive clock rules `(r_i,w_i)` and `(q_j,v_j)`, each with exact mass one and first moment one-half, and each with uniform Hermite moment error at most delta. Define only analytically

    B_Q=sum_i w_i r_i P_(r_i)H,
    C_Q=2 sum_j v_j q_j P_(q_j)[B_Q^2].

Both are pointwise PSD, with the same upper bounds as B and C_cont. The shifted-degree uniform moment estimate proves

    ||B_Q-B||_L2(HS) <= delta ||H||_L2(HS)
                     <= delta A sqrt(D).

For noncommuting matrices use the exact identity

    B_Q^2-B^2=B_Q(B_Q-B)+(B_Q-B)B.

It gives `||B_Q^2-B^2||_L2(HS)<=delta A^2 sqrt(D)`. The positive outer operator has total mass `2 sum_j v_j q_j=1`, so it is an L2(HS) contraction. The outer shifted-degree quadrature error is at most `2 delta ||B^2||_L2(HS)<=delta A^2 sqrt(D)/2`. Therefore

    ||C_Q-C_cont||_L2(gamma;HS) <= (3/2) delta A^2 sqrt(D).

This is a dimension-safe, one-energy, positive *analytical covariance approximation*. Its inner count is `O(log^2(1/delta))`; the outer count has the same order. A literal expansion has two independent inner clocks and one outer clock. Implementing their mean-Jacobian square action with legal original-VALUE producers remains a separate compiler step; no P, H, B_Q or C_Q evaluation is made available as an oracle by this theorem.

## 10. Independent diagnostics and stopping boundary

`check_clock_gradient_independent.py` passes **1,060 assertions**, seed 620041004. It checks exact root factors and both calibration branches; logarithmic readout identities; genuine C2, generally non-C3 ridge potentials; independent potential/first directional differences; centered and noncentered caller paths; zero and origin bounds; exact quadratic L2 energy without sqrt(N); full matrix Gaussian covariance; positive cosine diagonal debt; and independently integrated nonlinear currents and target centering.

Maximum observed potential-gradient directional error is `2.62e-11`; private/caller directional errors are below `2.46e-11`. At 401 outer nodes the cosine coefficient debt is `2.5564035911301176e-6`, with lower bound `1.5703380954884698e-6`. Its calibrated readout norm squared is `7.469403701210473`.

`check_conditional_covariance_square.py` passes **93 assertions**. It compares direct conditional Gaussian pair covariance to the positive square identity for nine bounded sinusoidal fixtures at two numerical resolutions. The maximum identity discrepancy is `5.56e-17`; the 145-node positive square-rule diagnostic differs by at most `6.49e-12`. It also checks the exact degree-two polynomial normalization. Its explicit expectation/Hessian calculations are diagnostic evaluation of analytical formulas, not source producer instructions.

The checker does not execute LOW30's high-order compiler, enumerate its covariance grammar, establish native retained/proxy reentry, or prove a new order-to-complete-cost recurrence. The independent proof and tests establish the raw gradient and analytical quadrature ports stated above, with all further consumer gates left explicit.
