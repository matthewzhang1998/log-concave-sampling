# Rank-aware common auxiliary heat restoration for the actual K3 source

New source-qualified theorem, 2026-10-04. It uses P's exact K3 derivative and rank-aware energy observation, independently checked below. It supplements the held b03b0c2f physical source screen. Only the explicit K3 source is covered here; no arbitrary finite-K derivative extension or stationary nonlinear-constraint theorem is asserted.

## 1. Source and complete physical heat

Use the same complete source as c8b20bf4:

    x=cS+sU, Delta=g_t(S)-g_t(S+epsilon Z),
    x1=x+aM*Delta,
    t0=S+aM g_0(x), t1=S+aM g_0(x1),
    E=r[g_t(t0)-g_t(t1)].

The original S,U,Z are independent standard d-vectors, c²+s²=1, and M is a supplied known contraction. Assume BOTH original gradient Hessians lie in [.4I,.6I]. They may be the same primitive, and no derivative or continuity modulus of their Hessians is used. Write

    Q=MM*,  kappa0=ra,
    e=||E-EE||_2.

The physical first of E is at most C r for bounded a. Add independent standard V and execute the SAME displacement in both terminal calls:

    E_sigma=r[g_t(t0+a sigma MV)-g_t(t1+a sigma MV)].              (1)

All original feedback, ancestors and repeated query records remain intact.

## 2. The exact actual-energy lower bound

Differentiation of the literal K3 source gives

    D_Z E=r a² epsilon H_out M H_0 M* H_plus,                   (2)

where every H is the actual original Hessian at its recorded query. The feedback derivative has the stated positive sign. Each H=.5I+R with ||R||op<=.1.

The middle-factor Hilbert inequality is

    ||M R M*||HS <=||R||op ||MM*||HS.                            (3)

For example diagonalizing M's singular values gives the square sum sum_ij s_i²s_j² R_ij². The inequality s_i²s_j²<=(s_i^4+s_j^4)/2 and the operator row/column bounds prove(3). This is also valid for nonsymmetric R, although the actual R is symmetric.

Expanding the seven nonconstant product words in(2), using(3) for the middle occurrence and ordinary operator bounds for the two outside occurrences, gives

    ||H_out M H_0 M* H_plus-.125Q||HS
                    <=[(.5+.1)^3-.5^3]||Q||HS=.091||Q||HS.       (4)

Thus, for Q nonzero, the FIXED test Q/||Q||HS has pairing at least .034||Q||HS with every derivative product in(2). This is pointwise even though the product itself need not be symmetric.

At each fixed S,U, Gaussian Z integration by parts and first-chaos Bessel yield

    ||E_Z[D_Z E]||HS <=||E-E_Z E||_(L2Z),
    ||E_Z[D_Z E]||HS >=.034 r a² epsilon ||Q||HS.

Integrating the conditional variance proves the actual CENTERED energy bound

    e>=.034 r a² epsilon ||MM*||HS.                              (5)

No invertibility or lower singular-value assumption is needed. This is a same-version source statement, rather than a lower bound on a nominal envelope. If M=0 structurally then E=E_sigma=0 and no ratio below is formed.

## 3. Finite VALUE restoration and coefficient restoration

The two actual terminal increments in(1) give, simply by the original Lipschitz bound,

    ||E_sigma-E||_2 <=2r a sigma||M||HS.                         (6)

The same inequality holds for m_sigma=E_V E_sigma. Define the known geometry factor

    R_M=||M||HS/||MM*||HS.

Combining(5)-(6), with delta_sigma denoting either of these L2 differences, gives

    delta_sigma <=C [sigma R_M/(a epsilon)] e.                   (7)

This estimate uses neither a terminal third derivative nor a Gaussian displacement norm unrelated to the actual mark. The apparent dimension factor in(6) is paired with the source-specific lower bound(5).

Here is the exact bounded-first continuity estimate needed for the physical orientation. Let F,G be physical vector fields on the SAME complete standard tape, with actual complete firsts L_F,L_G. Let

    O(F)=Cov(F)-Sym integral_0^1 E[J_(F,S)(u)^2]du,

with the original physical S row, padded by zero on any added variables. Then

    ||O(F)-O(G)||HS <=2(L_F+L_G)||F-G||_2.                       (8)

For the covariance part, split into the two centered cross products and use Gaussian Poincare to give the ordinary row frames L_F,L_G, while the difference supplies the Hilbert mark. For the forward part, split J_F²-J_G² into J_(F-G)J_F+J_G J_(F-G). The remaining factors have operator bounds L_F,L_G, and Gaussian heat-gradient isometry gives

    integral E||J_(F-G,S)||HS du <=||F-G||_2.

This proves(8), including nonsymmetric physical forward coefficients. It does not execute a Jacobian or assume derivative convergence of the VALUE difference.

Both E and E_sigma have complete first at most C r when sigma<=1; the extra V first is at most2r a sigma||M||. Jensen gives the same old-root first for m_sigma. Hence(7)-(8) imply

    ||O(E_sigma)-O(E)||HS,
    ||O(m_sigma)-O(E)||HS
        <=C [sigma R_M/(a² epsilon)] kappa0 e.                  (9)

Therefore, for a desired relative coefficient tolerance b>0, the literal known choice

    sigma <= c b a² epsilon / R_M                              (10)

gives error at most b kappa0 e, with a sufficiently small numerical c.

For an orthogonal projection M of any rank, R_M=1. Thus at r=a=A and epsilon=A^.9, sigma=c b A^2.9 suffices UNIFORMLY in dimension and rank. For an arbitrary known M, its actual geometry factor remains in the width. It cannot be omitted or replaced by a dimension-independent constant without a supplied spectral premise. The target is the physical complete-source orientation; no mixed original-mark curl assertion beyond(9) is required for this theorem.

## 4. Actual implementation and precision

The auxiliary program still has the original six K3 gradient VALUES, one new Gaussian block and a known M action. Width(10) introduces no inverse-sigma replica or extra callback count. All old first/caller paths remain ordinary, and the actual pointwise uncentered mark from b03 is preserved by strong terminal monotonicity.

The width is small in arithmetic terms. Computing R_M or a certified conservative bound, multiplying a sigma M, generating the new Gaussian row, original-oracle accuracy and source encodings must use their actual precision/cost. If only an upper bound Rbar>=R_M is supplied, substitute Rbar in(10). Do not test an undecidable exact-real zero or divide by an uncertified tiny coefficient; use the structural-zero branch or a supplied certified positive geometry bound. No cheap dense-matrix arithmetic is assumed.

There is no inverse width in evaluating the physical source or the unnormalized auxiliary gradient h_W=M*E_sigma. An imported mean or pair compiler must retain its own complete call/variation bill and any physical inverse-M readout. The present count is for the literal modified K3 source, not an unproved complete stationary producer.

## 5. What this closes

This gives a constructive finite-width, one-actual-energy restoration of the constant orientation for a broad explicit native matrix class, even with merely bounded Hessians. It supplements the sharper source-specific drift test, and explains why its aligned auxiliary heat cannot be rejected using the unrelated isotropic primitive-heat counterexample.

It does not remove the original-root orientation from m_sigma: on the held coherent-drift family that mean retains its order-kappa0 e orientation. Nor does it make the nonlinear ambient constraint vector Gaussian or license a marginal comparison while another block observes the same auxiliary roots. Those are separate current/source-return obligations.
