# Independent exact-text audit of the parallel radial-shell proof

2026-10-05. Reviewer: single-history-current analyst.

Reviewed source: `../RADIAL-SHELL-REFUTES-THE-UNSHIFTED-MEAN-COVARIANCE-REDUCTION.md`, SHA256 `9f8037bd618c82e188acbbc807935ee2b606c011d405603e31c27649d0dd0060`.

Verdict: PASS for the stated genuine same-source counterexample to the UNSHIFTED single-history covariance reduction. This review does not mark the reviewer's separate canonical author proof independently audited; that source has SHA256 `2807f7307646cc16c79049319b574259874b09e731a13f686696bdd60a05b035` and awaits a fresh third reviewer.

1. The radial bump yields a globally smooth anchored gradient. The stated fixed c,d choice makes both radial and tangential Hessian eigenvalues nonnegative and below A for large D.
2. The exact linear conditional OU integral has mean x/2 and covariance I/4. The actual common-root packet has linear private coefficient beta_Q. Both nonlinear K terms are bounded in norm by one. Their means differ by at most delta in L2 by the admitted one-variable multiplier guarantee. No joint-law or strong path approximation enters.
3. The finite-rule covariance gap follows directly from sqrt(1-tau^2)>=1-tau^2 and the n=2 multiplier error. It does not rely on convergence to a continuous common-root packet.
4. The strong radius expansion in source equation (7) is valid. After retaining the radial linear and tangential quadratic terms, the remainder is bounded on the good event by C[|n.zeta||zeta|^2/R^2+|zeta|^4/R^3]. Its Lp norm is O(A^3), because the radial Gaussian component is O(A), the displacement norm is O(1), and R=1/A. The bad-event tails are negligible. The expected radial displacement is (c^2 sigma^2/2)A+O_L2(A^2), with the fixed c/2 mean shift retained in the bump argument.
5. The direction increment is O_Lp(A); its conditional mean attenuation is O(A^2). Thus the scalar radial Taylor calculation yields the stated vector expectation with an O_L2(A^2) remainder.
6. The nonlinear K contribution is not silently discarded. Bounded D2f and K give a valid Taylor remainder. Df(y_a)-Df(y0)=O_Lp(A), so replacing that derivative while retaining K costs O(A^2) before multiplying by the outer A. The conditional K means then cancel up to delta. Correlation between G_a and K_a causes no omitted term.
7. Conditional first-chaos Bessel gives ||Cov(G_a,K_a)||HS<=||K_a-EK_a||2<=1 even without independence. Together with ||Cov(K_a)||HS<=1, this establishes the actual covariance remainder O(A^2), uniformly in x and D.
8. The displayed radial D2f contraction is correct. Its HS-to-vector operator norm is bounded independently of D, including outside the transition shell. Multiplication by A makes the covariance remainder current O(A^3). The isotropic part gives the stated unshifted psi'(S_D) coefficient.
9. The two leading terminal arguments differ by the fixed c/2 radial shift. The vector residual is k_D A^2 h(S_D)n+O_L2(A^3+A^2 delta), with k_D uniformly bounded away from zero.
10. The linear test Ax has norm one, and self-adjoint R1 sends it to Ax/2. The shell CLT and strict Gaussian convolution maximum of the symmetric bump make the limiting inner product nonzero. This proves the claimed A^2 lower bound after R1, beyond every fixed polylogarithmic A^3 allowance.

Scope is correctly bounded: the proof rejects the unshifted covariance-at-x formula using the SAME original gradient in all branches. It does not refute a centered/resummed current or close the native same-carrier m3 consumer. No ordinary Markov path-grid algorithm is executed.

Independent scalar diagnostics for the more explicit canonical choice c=d=1/4 are in `../radial_shell_separator_constants.json`: normalized bump maximum 0.41428441993455334, limiting shell E[h(S)]=-0.0018295069327453506, and a positive uniform asymptotic witness lower coefficient 6.451711731144085e-7. Those numbers are sanity checks, not a substitute for the analytical strict-sign proof or a finite-dimensional convergence-rate theorem.
