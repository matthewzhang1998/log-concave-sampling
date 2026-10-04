# Independent audit: finite third-order full reverse-OU law

Date: 2026-10-04.

## Verdict and bounded scope

**PASS for the finite positive full-law endpoint construction, under the explicitly imported fixed-order mean and square-action guards.** The standardized endpoint error is

    Lambda(A) A^3 sqrt(D) + original mode/numerical floors,

and its original physical normalization has error `Lambda(A) A^(7/2) sqrt(D)` plus the correspondingly scaled original floors. Both complete original VALUE work and Gaussian-root count are fixed public-log polynomials at this grade. The public-log factor is retained; this audit does not turn the covariance service's `Lambda A^3` error into a log-free bound.

The verdict is a bounded composition of the admitted original-VALUE services. It does not re-prove the large imported mean and LOW30 response programs, admit an arbitrary higher-order recurrence, provide a strong conditional-mean oracle, or establish `c(P)/P -> 0`.

Inspected integrated source: `../THIRD-ORDER-FULL-REVERSE-OU-LAW-AND-ENDPOINT.md`.

Final source SHA256: e83d737da13d05054ff0dda129a05e72d350c6e3e33b019d548cc517b816e702.

All component and supporting pins appear in `MANIFEST.json`. The component mean pin is `b056e213ac28c60869727589d3ee6d1b2c53d6666908cdbe1961fa3121161ced`; its separate independent audit is PASS. The covariance service pin is `108af2124e3da0e2a04af92e9f0403bfce2f82f7e2912ba6b7e1975ef1c8b18b`; its separate source-to-VALUE-action audit is PASS. The Gaussianization source pin is `f347c88a3987e50542dbb97e1e9d2cb1d8eb684cb6327509e89154d75efbb9aa`, with an independent proof of the full-law lemma.

The independent checker passes **400 assertions**. It tests affine/variance ledgers, finite conditional-OU Gram and one-energy first mechanisms with noncommuting PSD matrices, contraction-weighted schedule and terminal/floor sums through A=10^-40, and a scalar quadratic full-schedule reference through A=10^-12. The largest tested quadratic `W2/A^3` is 0.034126. These diagnostics do not instantiate the imported finite compilers and do not replace the nonlinear proof below.

## 1. The continuous source can legitimately enter Gaussianization

For fixed fresh standard carrier Z, approximate the analytical OU path only in the proof by finite positive quadratures of total mass one. Use its exact finite conditional Gaussian Gram, rather than independent point noises or a common-root substitute. Every conditional row has operator norm at most one. The finite inner and outer force averages therefore satisfy, uniformly in their Gaussian input dimension,

    Lip(H_N) <= A,
    Lip(F_N) <= A(1+A).

For a D-output L-Lipschitz function of an arbitrary finite standard Gaussian bank, the output derivative has Hilbert-Schmidt norm at most `sqrt(D) L`. Gaussian Poincare hence gives centered energy at most `sqrt(D) L`, with no factor from the size of the private bank. Thus

    ||F_N - E(F_N|Z)||_2 <= A(1+A) sqrt(D).

The uncentered feedback chord satisfies

    ||F_N-H_N||_(2|Z) <= C A^2 (|Z|+sqrt(D)),

by the A-Lipschitz outer force and Minkowski for the inner force. Only one D-dimensional energy is used. These bounds are uniform over the finite approximations and over the external source parameters after anchoring.

On intervals away from zero the OU path is mean-square continuous. Lipschitz force and positive quadratures therefore converge in L2; the omitted interval near zero costs its length times the displayed uniform moment profile. A Gaussian point at clock zero is unnecessary. Nested approximations can be chosen to converge in L2 to the exact second-substitution force. Apply the finite-dimensional Gaussianization theorem before taking this limit. Means and covariance matrices converge, and the positive-buffer Gaussian reference laws converge in W2. The displayed uniform bounds pass to the limit. This argument never executes an infinite path or charges its approximation as the finite algorithm.

The buffered Gaussianization is applied to `-h(F_cont-E[F_cont|Z])` with buffer variance `1-h^2 >= 3/4`. It therefore prices every higher private-path cumulant by `C h^3 A^3 sqrt(D)/(1-h^2)`. No assumption about the covariance service can substitute for this full-law step.

## 2. The covariance substitution has one physical energy

Let `E=F_cont-H_cont`, centered when necessary. Exactly,

    Cov(F_cont|Z)-Cov(H_cont|Z)
      = Cov(E,F_cont|Z)+Cov(H_cont,E|Z).

For any centered vectors U,V, not necessarily Gaussian,

    ||E[UV^T]||_HS <= ||U||_2 ||Cov(V)||_op^(1/2).

To see this, apply the operator from scalar L2 to R^D, `a -> E[a V]`, to each coordinate of U; its squared operator norm is `||Cov(V)||_op`. Gaussian Poincare gives covariance operator norm `O(A^2)` for H and F. Hence the difference above is bounded by

    C A^3 (|Z|+sqrt(D))

in Hilbert-Schmidt norm. This is dimension-safe. Multiplying two D-sized vector energies would have been too weak.

The completed reserve restores the continuous covariance C(Z), while the mean service restores the exact continuous conditional mean m_cont(Z). Their restoration targets agree on the same Z and the same source. The square-clock rule is never used as though it reproduced the full conditional force variable.

## 3. Positive reserve, varying covariance, and complete-bank join

For `0<h<=1/2`, put `d=1-2h^2 >= 1/2` and `eta=zeta=sqrt(d/2)`. The executed reserve is

    R = eta p + h^2 D_Q(p)/(2 eta) + zeta z.

The full action D_Q retains its common p across the square-clock banks. Its other roots, including the owned coarse X clocks, remain inside its private bank. This is a VALUE response program, not evaluation of a covariance matrix. The action is odd in p, so the reserve is exactly centered under exact symmetry (with numerical symmetry errors assigned their own floor).

Conditional on Z, the reserve targets `N(0,dI+h^2 C(Z))`. Its error is

    Lambda [h^2 A^3 + h^4 A^(7/2) + h^4 A^4] sqrt(D)
       + C h^2 delta A^2 sqrt(D) + floors.

The h^4 terms include the private fluctuation and the positive quadratic covariance term of the linearized reserve; neither is silently erased. Multiplication by h^2 is applied after the completed action and does not change its normalized source radius or introduce inverse-h replay counts.

The mean and reserve programs are independent COMPLETE banks after Z and every external caller are captured. Their Gaussian reference covariance is exactly

    h^2 I + d I + h^2 C(Z) = (1-h^2)I+h^2 C(Z).

Its dependence on Z presents no positivity problem: C(Z) is positive semidefinite, and all comparisons have a fixed positive gap. The Gaussian-root/Sylvester bound works without matrix commutation. Conditional product couplings can therefore be formed at the same Z and then integrated over the fresh standard carrier.

Only whole completed outputs and retained callers are read. In particular no private clock, old path sample, source-zero Gaussian, or reference coupling noise is reattached after its integration. Consequently the compiled output reaches the buffered continuous P2 law, and then `h mu_U + sqrt(1-h^2) gamma`, with error `Lambda h A^3 sqrt(D)` plus floors.

## 4. Conditional scaling and the finite-mode tilt

At an arbitrary entering caller z, write `s^2=1-r^2`, `Delta=t^2-r^2`, `alpha=A s^2`, and use the recorded finite conditional mode x_M. The executed source

    f(y)=s[g(x_M+s y)-g(x_M)]

is anchored and has Hessian between zero and alpha I for every such caller. Its exact normalized target differs from the true conditional posterior only by the linear tilt `R_M/s`. Strong convexity bounds the normalized W2 shift by `|R_M|/s`.

The affine readout obeys

    B_rt s = sqrt(v0) h,
    v = v0(1-h^2),
    sqrt(v0) h alpha^3 = (Delta/t) A^3 s^5.

Thus the full local law error is exactly of the stated scale, plus `B_rt |R_M|`. This is the proper mode-restoration factor; the residual is not folded into a uniform fixed-caller sqrt(D) bound.

The mean and covariance quadrature guarantees are L2 in their standard carrier, rather than pointwise bounds in a deterministic carrier. This poses no mismatch here: every call draws its own standard Z after exposing the current external caller. None of the local clock theorems is being applied at the entering non-Gaussian state z.

## 5. Exact kernel contraction and the endpoint sum

The exact conditional density has potential

    |x-rz|^2/(2s^2)+U(x).

Its strong-convexity constant is at least 1/s^2. Changing z to z' adds constant drift of magnitude `r|z-z'|/s^2`. The synchronous strongly convex Langevin comparison yields

    W2(nu_(r,z),nu_(r,z')) <= r |z-z'|.

Using that coupling and the same independent affine Gaussian noise, the exact reverse kernel has W2 Lipschitz constant

    A_rt+B_rt r = r/t.

This extends directly to distributions of entering callers. At r=0 the reference kernel is independent of its input. At the terminal bridge its contraction is r_J. This is a law-kernel fact, not the derivative of the actual compiled stage.

For `rho=3/4`, `s_j^2=rho^j`, the positive-buffer stages have `Delta_j=s_(j-1)^2/4` and h=1/2. The final contraction weight of a local error produced at r_j is r_j. Consequently

    sum_j r_j [(Delta_j/r_j) A^3 s_(j-1)^5]
      = (A^3/4) sum_(k=0)^(J-1) rho^(7k/2)
      <= A^3 / [4(1-rho^(7/2))].

The constant on the right is approximately 0.393921. There is no J or inverse-heat loss in the intrinsic error sum, apart from the already explicit common public-log majorant.

Choose J minimally with `s_J^5<=A`. The separately admitted old positive conditional packet then costs at most `C A^2 s_J^5 sqrt(D) <= C A^3 sqrt(D)`. It is valid without the Gaussian reserve used on prior stages. Thus the zero-buffer terminal step is handled by an actual separately priced law program, rather than extrapolating Gaussianization to h=1.

## 6. Actual caller moments, finite modes, and firsts

For M>=3, the local mode floor is

    B_rt |R_M| <= (Delta/t) A^(M+1) s^(2M) r |z|.

This is evaluated at the actual entering law. Put `e_j=W2(actual_j,rho_(r_j))` and `u_j=r_j e_j`. Since exact OU marginals have second moment at most D, the entering moment is at most `sqrt(D)+e_(j-1)`. For M=3 the weighted recursion has the form

    u_j <= [1+C A^4 Delta_j s_(j-1)^6] u_(j-1)
             + C Lambda A^3 Delta_j s_(j-1)^5 sqrt(D)
             + C A^4 Delta_j s_(j-1)^6 sqrt(D)
             + weighted numerical floors.

The first r=0 mode residual vanishes. The heat sums of both coefficients are bounded. Discrete Gronwall gives `u_j<=C Lambda A^3 sqrt(D)` plus floors. Every j>=1 has r_j>=1/2, so actual entering moments are `O(sqrt(D))` under the smallness guards. The terminal mode floor is also of smaller order. This closes the caller-moment loop without assuming the approximate chain already has exact OU marginals.

For the actual first port, let the external stage input be z. The finite mode obeys

    ||D_z x_M-rI|| <= r alpha/(1-alpha),
    Lip_z f <= C A s r.

The imported mean caller port gives `C Lambda A s r`; the square-action parameter port gives `C Lambda sqrt(alpha) A s r` before its h^2 readout. Therefore the actual compiled stage satisfies the sufficient bound

    ||D_z K_new|| <= (r/t)
       [1+C Lambda A Delta+C Lambda A^(3/2) Delta^(3/2)].

At r=0 its input is absent from the source graph. The terminal packet has corresponding incoming bound `r_J[1+C A s_J^2]`. Since `sum Delta_j<=1`, the products of actual incoming first bounds are stable under the existing small-A/public-log guards. These inequalities use literal first chain rules, recorded origins, and captured parameters; they do not differentiate a W2 error. Original physical caller derivatives are separately propagated through the source and mode records with their original normalization.

The source-zero stage is an actual known Gaussian row with incoming coefficient r/t and fresh covariance `v0 I=(1-(r/t)^2)I`. Induction, followed by the terminal row, gives a unit coisometric endpoint carrier. Direct fresh-root residual derivatives contribute at most `C Lambda A sum_j sqrt(Delta_j) s_(j-1)^2`, which is a fixed geometric sum. One must also include the source-dependent incoming-Jacobian corrections above; their final sum is bounded by `C Lambda A sum_j Delta_j` times the already bounded entering first. The corresponding derivative recursion and discrete Gronwall therefore give a literal carrier-subtracted endpoint first `O(Lambda A)`, and full private first `1+O(Lambda A)`. The terminal's direct and incoming corrections have the same bound. This argument concerns the executed Gaussian tape and never replaces it by a law-coupling noise.

Every requested JVP/VJP uses only original HVPs at stored original VALUE sites, including the complete nested ancestors. A discarded primal pays its full replay; no saved HVP is differentiated. No stronger retained/protected/proxy interface is inferred from this law theorem.

## 7. Complete cost, precision, and source-zero audit

Minimal stopping gives

    rho A^(2/5) < s_J^2 <= A^(2/5),
    rho A^(7/5) < alpha_J <= A^(7/5).

Thus all local `log(1/alpha)` values are O(log(1/A)), and `J=O(log(1/A))`. The two raw mean sources cost `n_out` and `n_out(n_in+2)` original calls per complete occurrence. The square-clock service costs `N_r N_q(2K_f+5)` before explicitly requested restoration/replay. Each conditional mode, saved anchor, private origin and changed source argument is paid under its own semantic key. The terminal packet's separate mode and n_alpha source values are included. Summing fixed-order public-log occurrence counts over J remains public-log.

Independence is required only between complete banks at the specified cuts; the common p inside a covariance aggregate and the common inner/outer roots within a mean-source occurrence remain literal. The full Gaussian dimension counts all replay/filter/action banks, not merely the schematic displayed carrier roots. No root count is absorbed into a dimension-free first bound.

All quadratures, response widths, padding values, finite-mode counts, coefficient encodings and numerical versions are frozen before differentiation. Precision floors are absolute and propagated through their actual readouts and reference-kernel weights. The standard existing precision model can allocate them below the claimed grade using public-log bit precision; nothing divides an error by a realized source energy.

Raw anchored zeros are literal under recorded anchor reuse, while nonzero fixed-carrier origins are captured and restored. Completed mean and reserve programs retain their own known Gaussian source-zero carriers. Exact target modes and coupling Gaussians are not substituted for executed records. The original physical global-mode linear tilt and its caller profile remain separate restoration terms.

## Conclusion

The integrated proof closes the bounded next-law-order gate: a finite positive original-VALUE nonlinear endpoint program reaches standardized public-log order three and physical order seven-halves. The analytical continuum is only a reference; the finite mean and covariance branches plus Gaussianization supply an actual full-law bridge. The zero-buffer terminal step, actual entering-caller moments, exact contraction, complete replay costs, and source-zero records are all accounted for. No all-order or eventual-sublinear recurrence is certified.
