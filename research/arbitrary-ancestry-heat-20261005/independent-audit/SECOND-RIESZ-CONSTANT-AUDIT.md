# Additional independent check: Section 7 second-Riesz constant

2026-10-05. Scope: only Section 7 of `../FINITE-ANCESTRY-COMPILER-AND-BOUNDARY.md`.

## Verdict

**PASS.** The final Section 7 now explicitly includes the finite-dimensional independent standard Gaussian banks and the base condition `b(P) in L2`; the initially identified minor hypothesis omission is resolved. Assume finite-dimensional independent standard Gaussian Q and P, F centered in Q for almost every P, physical polynomial degree at most k, the stated first Sobolev regularity, and

    ||F||_2 <= E,
    ess sup_Q ||D_Q F||_(L4(P);operator) <= L.

For the ordinary Wasserstein space P2, also assume `b(P) in L2`. No derivative of b or C in P is needed. Under these assumptions the stated constant is valid:

    W2(Law(W1),Law(W0))
       <= [2 sqrt(2)/(3 kappa)] 3^(3k/2) E L².

There is no second Q derivative of F in the proof. The polynomial identity tests are correctly labeled as identity tests, not witnesses for the uniform derivative hypothesis.

## Norm estimates

All unmarked tensor norms below are Hilbert/HS norms. Let `R=D_Q N_Q^-1` componentwise. For a centered Hilbert-valued f, its Q-chaos expansion gives

    ||R f||_2² = sum_(n>=1) ||f_n||_2²/n <= ||f||_2².

This operation leaves physical degree unchanged. A Hilbert-valued P-polynomial of degree at most j obeys

    ||H(Q,.)||_L4(P) <= 3^(j/2) ||H(Q,.)||_L2(P)

at almost every Q. Apply this before integrating Q. Writing `B=D_Q F` and `A=R F`, matrix multiplication gives pointwise

    ||A B^T||_HS <= ||A||_HS ||B||_operator.

Conditional Holder in P, the uniform-in-Q derivative bound, and hypercontractivity therefore imply

    ||tau||_2 <= L 3^(k/2) ||R F||_2
               <= 3^(k/2) E L.

Gaussian integration by parts yields `C(P)=E_Q tau=Cov_Q F`. Conditional centering is an orthogonal projection, so with `S=tau-C`,

    ||S||_2 <= ||tau||_2.

S has physical degree at most 2k. Treat `R S` as a matrix with combined physical indices (i,j) on one side and the Q derivative index on the other. For

    U_ijl = sum_a (R S_ij)_a B_la,

the same HS/operator multiplication bound, now with physical hypercontractivity factor `3^k`, gives

    ||U||_2 <= L 3^k ||R S||_2
            <= 3^(3k/2) E L².

These estimates never exchange a Q supremum with a P expectation. The stated `sup_Q L4(P)` hypothesis is sufficient precisely because it is paired at each Q with the other factor's physical L4 norm, followed by Q integration. Only one Hilbert energy is used.

## Endpoint identity and transport factor

Set `s=kappa/2`. The analytical Gaussian path is

    W_t=b(P)+tF(Q,P)+sqrt(s)G0
          +[s I+(1-t²)C(P)]^(1/2)G1.

Differentiation in t and one Q integration by parts give

    d E phi(W_t)/dt = t E[S_ij partial_ij phi(W_t)].

S is Q-centered. The weak Riesz identity applies directly to S and the endpoint test, even when S itself has no Q derivative. Since C(P), b(P), and the Gaussian covariance root are Q-independent,

    d E phi(W_t)/dt
       = t² E[U_ijl partial_ijl phi(W_t)].

In particular the derivative lands on the endpoint only. This is why no `D_Q tau` or `D_Q²F` appears.

U is independent of G0. Transfer the j and l derivatives through this fixed Gaussian:

    d E phi(W_t)/dt
       = (t²/s) E[V_i partial_i phi(W_t)],
    V_i = sum_(j,l) U_ijl (G0_j G0_l - delta_jl).

Conditional on U, Gaussian Hermite isometry gives

    E_G0 |V|²
       = sum_(i,j,l) U_ijl² + sum_(i,j,l) U_ijl U_ilj
       <= 2 ||U||_HS².

No symmetry of U is required. The conditional velocity

    v_t(z) = (t²/s) E[V | W_t=z]

satisfies the continuity equation and has

    ||v_t||_(L2(Law(W_t))) <= (sqrt(2)/s) t² ||U||_2.

Integrating its metric-speed bound from zero to one yields the factor

    (sqrt(2)/s) integral_0^1 t² dt
       = 2 sqrt(2)/(3 kappa).

For first-Sobolev coefficients the identities can be obtained by Q smoothing/finite chaos approximation and truncation, followed by the displayed uniform integrability bounds. The covariance root remains uniformly positive because of sI. Adding `b(P) in L2` places the path in the standard P2 setting needed by this last transport argument.
