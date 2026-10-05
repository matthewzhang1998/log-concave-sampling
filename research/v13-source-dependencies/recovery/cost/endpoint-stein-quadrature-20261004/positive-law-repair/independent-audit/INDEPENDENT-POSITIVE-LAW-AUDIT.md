# Independent audit: positive endpoint Stein law and quadratic amplifier

Date: 2026-10-04.

## Verdict and pins

**PASS, strictly for the stated bounded component.** The actual finite two-root endpoint packet satisfies the dimension-safe general-C2 marginal law bound

    W2(Law(Z-H_Q(Z,G)), mu_U) <= 10 (A^2 + delta A) sqrt(D),

under `g=grad U`, `g(0)=0`, `0<=Dg<=A I`, `0<A<=1/2`, and the audited positive quadrature's conditional-mean certificate. This yields physical order 5/2 when `delta<=A`, after separately restoring finite-mode and numerical errors. Its literal producer is unchanged and uses only original gradient VALUES. No conditional expectation, flow, entropy, density, or Hessian action is part of that producer.

The matrix-quadratic specialization also passes: the displayed recurrence uses exactly three additional anchored gradient VALUES per layer, preserves the actual common-G covariance, and attains standardized error `O(A^(2m+2) sqrt(D))`. A universal law constant `4/3` is available for that exact specialization. The corresponding physical power is `2m+5/2`.

No mathematical blocker was found within these claims. This verdict does **not** establish an all-order general-C2 law source, a joint law comparison retaining G, arbitrary later observers, a native retained/proxy port, a strong fine/coarse pair, or a sublinear complete-cost recurrence. The quadratic recurrence is not a nonlinear higher-current compiler.

Inspected source:

- `../POSITIVE-SECOND-ORDER-LAW-AND-QUADRATIC-AMPLIFIER.md`
- Final inspected SHA256: `a085dfdef3f66f69208612034c17e4e18045fc00a21be3aed64ccfcb5e9b7085`
- Original Sections 1–6 pin: `6073cd34d5bd01d2a17a632c540d3fab087fcb2d05b2a808a7114627c0f7b0fe`
- The final source's prefix before Section 7, with one terminal newline, hashes exactly to the original pin. The only added material is diagnostics/sharpness.
- Underlying endpoint packet SHA256: `4f0c182c36bf9dfbacc0f189899f689de52a2bb9bc37839e307eb76c4bda519f`.

The independent program `check_positive_law_independent.py` does not import the author's diagnostics or the earlier audit program. It passes **793 assertions**. Its detailed output is `positive_law_independent_checks.json`; the manifest pins the source, program, output, and this audit. Finite diagnostics supplement the proof below and are not a numerical proof of an infinite-dimensional or all-order statement.

## 1. Full joint entropy: cancellation and exact hypotheses

Let `Y_t=Z-tH(Z,G)`, with independent standard Gaussian roots `Z in R^d` and `G in R^n`. Fix G. If `||D_ZH||<=a<1`, the equation `z=y+tH(z,G)` is a contraction. It has a unique solution for every y. The inverse-function theorem, `||(I-tD_ZH)^(-1)||<=1/(1-ta)`, and the contraction construction give a global C1 diffeomorphism. Its determinant is positive: it starts at one when t=0 and cannot cross zero. Symmetry of `D_ZH` is unnecessary.

The joint transformation `(Z,G)->(Y_t,G)` is triangular and has determinant `det(I-tD_ZH)`. Comparing its density to the product Gaussian density gives exactly

    KL(Law(Y_t,G) || gamma_d x gamma_n)
      = E[-t Z.H + t^2 |H|^2/2 - log det(I-tD_ZH)].

This is the full private-bank joint law, not a nodewise computation. All common-root correlations remain present.

Gaussian integration by parts yields `E[Z.H]=E tr(D_ZH)`. For any real square B with `||B||<1`, the real analytic logarithm gives

    -log det(I-B)-tr B = sum_(k>=2) tr(B^k)/k.

For nonnormal B as well,

    |tr(B^k)| <= ||B||^(k-2) ||B||_HS^2.

One proof is to write `tr(B^k)=tr(B B^(k-1))`, apply Frobenius Cauchy–Schwarz, and use `||B^(k-1)||_HS<=||B||^(k-2)||B||_HS`. Bounding `1/k<=1/2` gives

    -log det(I-tB)-t tr B <= t^2 ||B||_HS^2/[2(1-ta)].

Thus the advertised bound is correct:

    KL_joint <= (t^2/2) [E_H^2 + J_H^2/(1-ta)].

No `D^2` factor is introduced by the trace. The derivative energy is essential; amplitude energy alone would not certify the entropy remainder.

In the packet, H is C1, all its firsts are bounded, and H has linear growth. These facts justify the change of variables, integrability, and integration by parts directly. The proof need not rest on the source's more general informal Sobolev extension. This audit admits the displayed C1 lemma and its C2-potential application, not every unspecified possible weak-regularity extension.

## 2. Conditional Gaussian comparison: the lost independence is paid

The marginal of G is still exactly gamma_n, though G and Y_t are generally dependent. The entropy decomposition is

    KL_joint = KL(nu_t || gamma_d) + I(Y_t;G),
    I(Y_t;G) = E_(Y_t) KL(Law(G|Y_t) || gamma_n).

All terms are finite by Section 1. Applying the Gaussian quadratic transport inequality separately to the conditional distributions gives

    E_(Y_t) W2^2(Law(G|Y_t), gamma_n) <= 2 I(Y_t;G) <= 2 KL_joint.

The external mathematical input is the Gaussian `W2^2<=2 KL` theorem, originally proved in [Talagrand, Transportation Cost for Gaussian and Other Product Measures](https://geodesic.mathdoc.fr/item/GFA_1996__6_1_58238/). Conditioning does not change its constant because the reference measure is the same standard Gaussian for every y.

Now let `h_t(y)=E[H(Z,G)|Y_t=y]`, and `v(y)=E_G H(y,G)` with the latter G independent. The exact splitting is

    h_t(Y_t)-v(Y_t)
      = E[H(Z,G)-H(Y_t,G) | Y_t]
        + E[H(Y_t,G)|Y_t]-E_G H(Y_t,G).

The first term has L2 norm at most `ta E_H`: `Z-Y_t=tH(Z,G)`. For the second, couple `G|Y_t=y` to a fresh standard Gaussian at fixed y and use the uniform G-Lipschitz constant b. Jensen and the preceding inequality then give

    ||h_t-v||_L2(nu_t)
      <= t [a E_H + b sqrt(E_H^2 + J_H^2/(1-ta))].

This is a valid dimension-free vector estimate. It does not sum coordinatewise Wasserstein bounds. It does not assume conditional independence. It does not infer strong closeness of H to v.

The motion term cannot be dropped. Already for `H=aZ` with no G dependence, mutual information is zero while `h_t(y)=a y/(1-ta)` differs from `v(y)=ay`. The independent exact-Gaussian diagnostics include this negative control and general nonsymmetric, affine, rectangular-root examples.

## 3. Continuity equations and the flow comparison

For a compactly supported smooth test, differentiating `E phi(Z-tH)` shows that `nu_t` solves the continuity equation with velocity `-h_t`. Conditional Jensen bounds its action by `E_H^2`. No smoothness or global Lipschitz assumption on h_t is needed.

The field v is C1 and globally a-Lipschitz. Let its remaining flow from t to 1 be `Phi_(t,1)`. The pushed-forward curve `(Phi_(t,1))#nu_t` has velocity

    D Phi_(t,1)(-h_t+v),

with derivative norm at most `exp(a(1-t))`. The standard finite-action Wasserstein length bound therefore gives source equation (5). Importantly, this is a **marginal continuity-equation argument**. Directly coupling the raw random paths would retain `H-v`, which can be first order, and would not prove the result.

For the deterministic flow itself,

    ||Phi_t(Z)-Z||_2 <= t exp(at) ||v(Z)||_2 <= t exp(a) E_H.

Integrating `|v(Phi_t(Z))-v(Z)|<=a|Phi_t(Z)-Z|` proves the Euler defect `a exp(a) E_H/2`. Thus equations (4)–(6) are closed, with no unimplemented producer hidden in the comparison.

## 4. Apply the lemma to the common-G packet and check the constant

Each `q_i=r_i Z+s_i G` is marginally standard Gaussian. Their joint independence is neither true nor needed. Positivity, total mass one, the exact first moment, and the Hessian sandwich give

    0 <= D_ZH_Q <= A sum_i w_i r_i I = (A/2) I,
    Lip_G(H_Q) <= A beta_Q <= A,
    ||H_Q||_2 <= sum_i w_i ||g(q_i)||_2 <= A sqrt(D),
    ||D_ZH_Q||_L2(HS) <= (A/2) sqrt(D).

The HS bound is pointwise before integration and uses the same derivative already permitted in analysis. It does not insert an HVP in the producer.

For `A<=1/2`, `1-tA/2>=3/4`, so

    KL_joint <= (2/3) t^2 A^2 D,
    ||h_t-v_Q||_2 <= (1/2+2/sqrt(3)) t A^2 sqrt(D).

Using `exp(a)<=exp(1/4)` and integrating `t` gives a convenient bound

    W2(Law(Z-H_Q), Law(Z-v_Q(Z)))
      <= exp(1/4) (1/2+1/sqrt(3)) A^2 sqrt(D)
      < 1.4 A^2 sqrt(D).

In particular the source's stated constant 3 is safe. The audited same-Z mean quadrature costs exactly `delta A sqrt(D)`. Neither quadrature accuracy nor marginal decoupling permits retaining G in the comparison: for quadratic g, the actual conditional fluctuation remains `beta_Q B G`, with L2 size of order `A sqrt(D)`.

## 5. Exact-target OU interpolation: signs, moments, and C2 regularity

Subtract the harmless constant U(0), so `0<=U(x)<=A|x|^2/2`. Then the target is normalizable, has Gaussian tails, and `|g(x)|<=A|x|`. Integration by parts against the target gives

    D = E_mu[|X|^2 + X.g(X)] >= E_mu |X|^2.

The Mehler symmetry identity proves that `rho_r = [P_r(e^-U)/E_gamma e^-U] gamma` is the law of `rX+sqrt(1-r^2)N`. Hence `E_(rho_r)|x|^2<=D` for every r.

For `0<r<1`, the OU parameter convention is crucial:

    partial_r f_r = -(1/r) L f_r,
    grad log f_r = -r m_r.

Because `L f gamma = div(gamma grad f)`, these identities yield

    partial_r rho_r = div(rho_r m_r),

so the continuity-equation velocity is indeed **-m_r**, with the sign in the source. The r=0 and r=1 limits are continuous. The drift at r=0 need not vanish: `m_0=E_mu g`, consistent with `E_mu X=-E_mu g`.

For fixed r,x, the auxiliary G-tilt has potential

    W(q)=|q|^2/2 + U(rx+s q),  Hess W >= I.

Its gradient is `q+s g(rx+s q)`. Synchronously couple its Langevin diffusion to a stationary standard OU coordinate N_t. Write the drift difference by evaluating the perturbation on N_t and use monotonicity of grad W. If R_t is the L2 separation, then

    R_t' <= -R_t + s ||g(rx+sN)||_L2(gamma).

Start with finite second moments and let t tend to infinity. This proves

    W2(nu_(r,x),gamma) <= s ||g(rx+sG)||_2.

The perturbation norm is under **gamma**, so no translated-mode moment estimate is being assumed. Since `q -> g(rx+s q)` has Lipschitz constant sA,

    |m_r(x)-P_r g(x)|
      <= s^2 A^2 sqrt(r^2|x|^2+s^2D).

After integrating against rho_r, the bound is `s^2 A^2 sqrt(D)`.

Differentiation in x uses only Dg and gives exactly

    D_x m_r = r [E_nu Dg(rx+sG) - Cov_nu(g(rx+sG))].

There is no additional r or s missing in the covariance term. The strongly log-concave tilt's Poincare inequality bounds, for every unit u,

    Var_nu(u.g(rx+sG))
      <= E_nu |s Dg(rx+sG) u|^2 <= s^2 A^2.

Thus `Lip(m_r)<=r(A+s^2A^2)` is safe. This use of the inverse-Hessian covariance inequality is consistent with [Carlen, Cordero-Erausquin and Lieb, Asymmetric Covariance Estimates of Brascamp-Lieb Type](https://arxiv.org/abs/1106.0709). No third derivative or Lipschitz modulus of Dg enters. Compact-x domination, Gaussian moments, and the positive denominator justify differentiation and the endpoint limits under the exact C2 assumptions.

The comparison flow `T_r'=-P_r g(T_r)` has `Lip(P_r g)<=rA`. Stability with the exact rho_r curve therefore costs at most `exp(A/2)` times the integrated drift error. The source's equation (11) follows (one may retain the extra integral factor 2/3). Gaussian invariance yields `||P_r g(Z)||_2<=A sqrt(D)`, and Gronwall gives

    ||T_r-Z||_2 <= r exp(A/2) A sqrt(D).

Its frozen-input Euler defect is consequently at most `exp(A/2) A^2 sqrt(D)/3`. Combining all terms leaves substantial room below the advertised constant 10. In particular the source's looser constants give

    W2 <= [3+(4/3)exp(1/4)] A^2 sqrt(D) + delta A sqrt(D),

which implies the claimed bound directly.

The complete argument is uniform over continuously varying noncommuting Hessians. It never interchanges products of local Hessians.

## 6. Physical scaling, finite caller, exact zero, and numerical scope

The ideal standardized theorem applies at every finite anchor b_M because

    g_M(z)=sqrt(A)[grad V(b_M+sqrt(A)z)-grad V(b_M)]

still has g_M(0)=0 and the same Hessian sandwich. The actual posterior at that anchor includes the linear tilt

    ell = [b_M-y+A grad V(b_M)]/sqrt(A).

With the finite fixed-count contraction iteration,

    |b_M-y+A grad V(b_M)| <= A^(M+1)|grad V(y)|,
    |ell| <= A^(M+1/2)|grad V(y)|.

Strong monotonicity bounds the standardized tilt restoration by `|ell|`; physical scaling costs `sqrt(A)|ell|`. This is a caller-dependent residual, not an intrinsic `sqrt(D)` error. Numerical mode errors remain additional. At delta<=A, the ideal law error scales by sqrt(A) to physical power 5/2.

The finite graph and original query sites do not change. Its complete VALUE count is `M+n+2` including the anchor and requested terminal canonical force; Gaussian dimension is exactly 2D. A directional first/adjoint sweep uses original HVPs at those recorded VALUE sites. An HVP is not itself differentiated, and a dense Jacobian is not being counted as one query.

Writing `B_M=D_y b_M`, the mode recurrence gives `||B_M||<=1/(1-A)` and `||B_M-I||<=A/(1-A)`. The literal caller derivative is

    D_y g_M(q)=sqrt(A)[Hess V(b_M+sqrt(A)q)-Hess V(b_M)] B_M.

Since both Hessians lie between zero and I, their difference has operator norm at most one. Thus `D_y X_Q=I+O(A)` with no modulus of Hessian continuity. The private first is `sqrt(A)[(I,0)+O(A)]`. Composing with the terminal original gradient gives canonical-force private first O(A) and caller first O(sqrt(A)). These match the earlier audited ledger.

Exact zero is genuinely conditional on the stated same-site rule. At zero roots the graph must reuse the recorded anchor, so every anchored gradient subtraction is bitwise identical and X_Q=b_M. Our diagnostic initially recomputed one site through a differently shaped floating-point matrix operation and observed a tiny zero mismatch; that was a **diagnostic implementation error**, not a contradiction of the source. The final diagnostic implements literal recorded-anchor reuse. Separate perturbed evaluations cannot be called exactly zero.

Rounded force/nodes/weights need not preserve an exact Hessian sandwich or exact quadrature moments. The theorem applies to the exact reference VALUE graph, and the actual numerical graph is coupled to it using the stated finite original-oracle allowances. No new numerical regularity assumption is inferred here. In particular, C2 alone does not give a quantitative Hessian-error modulus under perturbed query locations; requested first approximations retain the earlier original-first numerical contract.

A useful query envelope is unchanged:

    |b_j-y| <= A|grad V(y)|/(1-A),
    |q_i| <= sqrt(|Z|^2+|G|^2),
    |H_Q| <= A sqrt(|Z|^2+|G|^2).

Thus no maximum over independent node banks, hidden D-power heat guard, or inverse-A original query count appears. Random callers must be exposed before Z,G; then the same conditional theorem is integrated with the caller's actual residual/numerical profile. Dependence of a caller on this private bank is not licensed.

## 7. Matrix quadratic amplifier: complete VALUE graph and sharper law bound

For `g(x)=Bx`, `0<=B<=A I`, the shared-G packet is exactly

    Y=(I-B/2)Z-beta_Q B G,
    C(B)=Cov(Y)=I-B+cB^2, c=1/4+beta_Q^2.

The cross-node products are exactly those in `beta_Q^2`. Replacing them by `sum_i w_i^2 s_i^2` changes the law. Since `s(r)>=1-r` and s is concave, the exact first moment gives `1/2<=beta_Q<=sqrt(3)/2`, hence `1/2<=c<=1`.

The identity

    (I+B)C(B)=I+D(B), D(B)=(c-1)B^2+cB^3

is exact. On every eigenvalue `0<=b<=A<=1/2`,

    |D(b)| <= A^2 <= 1/4,  0<C(b)<=1.

The binomial coefficients `a_k=binom(-1/2,k)` satisfy `|a_k|<=1`. Therefore, uniformly in m,

    |sum_(k=0)^m a_k D(b)^k - (1+D(b))^(-1/2)|
      <= |D(b)|^(m+1)/(1-|D(b)|)
      <= (4/3) A^(2m+2).

Couple `Y_m=q_m(B)Y` to `Y_*= (I+D(B))^(-1/2)Y` using the same actual Y. All matrices are polynomials or spectral functions of the **single fixed symmetric B**, so

    Cov(Y_*) = (I+B)^(-1),
    E|Y_m-Y_*|^2 <= (16/9) A^(4m+4) tr C(B)
                         <= (16/9) A^(4m+4) D.

This proves equation (14) with constant 4/3. The diagonalization is legitimate for arbitrary dense B; no scalar-only covariance trace argument is used.

The execution is also literal. If `T_(k-1)=D(B)^(k-1)Y`, then three gradient VALUES return `BT_(k-1)`, `B^2T_(k-1)`, and `B^3T_(k-1)`, whose prescribed combination is T_k. No action of B is provided as a producer oracle. The complete count including finite mode, all quadrature nodes, the anchor, and a terminal canonical force **at the amplified endpoint** is

    M+n+3m+2.

The root count remains exactly 2D. If the terminal canonical force is not requested, subtract one. The terminal force at the unamplified endpoint need not be evaluated and is not being silently reused at a different endpoint.

For completeness, the omitted elementary first/zero/numerical ledger is safe, even for the executable nonlinear graph (without asserting its law):

- `|g(x)|<=A|x|` implies `|T_k|<=A^2|T_(k-1)|`, since `(1-c)A^2+cA^3<=A^2`.
- Its private first has the same factor `A^2` by the product rule and `||Dg||<=A`; no commuting local Hessians are assumed for this **bound**.
- The coefficients satisfy `|a_k|<=1`, so all fixed-m first and moment constants are finite. The output correction to the base private first is O_m(A^2); the full physical private first remains `sqrt(A)[(I,0)+O_m(A)]`.
- The direct caller derivative of every anchored g is at most `sqrt(A)/(1-A)`. Applying the finite chain rule at the three calls per layer gives an O_m(sqrt(A)) standardized caller correction and hence an O_m(A) physical correction. In a globally quadratic original potential the anchored g has no direct caller dependence at all.
- At zero roots, every T and every intermediate gradient argument is zero; recorded-anchor reuse makes all differences exactly zero. The physical output is b_M.
- Every added query has a displacement no larger than a fixed-m multiple of the earlier envelope. There is no new Gaussian block or D-dependent maximum.
- If each anchored VALUE has an added norm error at most eta, the layer recurrence has error bound `e_k<=A^2 e_(k-1)+(1+A+A^2)eta`. Since A<=1/2, this is at most `A^2 e_(k-1)+(7/4)eta`. Consequently `e_k<=A^(2k)e_0+(7/3)eta`, and the final finite sum is bounded by `(1-A^2)^(-1/2)e_0+(7m/3)eta`. Mode, node, and original-query errors must still be inserted with their actual scaling and correlations. There is no fictitious independent-error cancellation.

These graph facts do not transfer the quadratic law identity to nonlinear g. In general `g(g(x))` is not `Dg(x)g(x)`, and Hessians at successive query sites need not commute. The source states that limitation correctly.

## 8. Diagnostic coverage and interpretive limits

The independent program uses seed 73482109 and covers:

1. Exact nonsymmetric affine Gaussian entropy/conditional-velocity checks for `(d,n)=(1,2),(3,2),(5,7)`, four derivative scales, and four interpolation times. It checks joint entropy, mutual information, the motion-plus-information velocity bound, and deterministic Euler/flow error.
2. Twelve literal physical C2 source cases in dimensions 2, 5, and 11, two heats, and two finite-mode depths. It checks positive Z-first, full private first, caller and canonical firsts, derivative energy, mode residual, query census, exact recorded-anchor zero, and independent directional finite differences.
3. Genuine C2 but non-C3/non-Lipschitz-Hessian ridge potentials. The scalar Hessian is `sqrt(|t|)/(1+sqrt(|t|))`, and nonorthogonal ridge directions yield explicitly nonzero Hessian commutators. The normalized Hessian is globally bounded by .85 I by construction, rather than merely on tested points. Its Hölder cusp rules out a hidden Lipschitz-Hessian premise.
4. Thirty two-dimensional exact-posterior-drift quadrature cases, with r ranging from 0 to 1. These check the m_r derivative identity, covariance bound, and drift comparison on translated inputs. The quadrature is diagnostic; the analytical proof supplies the infinite-domain bounds.
5. Seventy-two dense matrix-quadratic heat/dimension/order cases through m=5 and D=47. The program evaluates the actual three-VALUE recurrence and matrix covariance, checks a deliberately different independent-node-root covariance, and verifies the spectral error bound. When the theoretical error drops below floating precision, the numerical tolerance is a floor, not evidence of asymptotic order; the binomial proof establishes those cases.
6. Direct one-dimensional CDF inversion for the actual two-root nonlinear positive law against its exact target density, independent of the entropy/velocity derivation. A monotone conditional Z map is inverted for each G quadrature node. Two resolutions are compared. At A=.5,.25,.125,.0625, the higher-resolution W2/A^2 values are approximately .02682,.03079,.03579,.03886. Corresponding resolution differences are about 3.90e-5, 1.85e-5, 8.21e-6, 3.85e-6 in W2. These are finite numerical observations, not claimed exact values or an asymptotic proof.

Section 7's quadratic sharpness statement is also algebraically correct: for fixed c<1 the unamplified standard-deviation error is `(1-c)A^2/2+O(A^3)`. The increasingly accurate dyadic quadrature has `c -> 1/4+pi^2/16<1`, so a real order-two law debt persists without amplification.

## 9. Final admitted boundary

The useful new admission is a positive, finite, original-VALUE, dimension-safe marginal W2 source at one fixed general-C2 physical law grade, together with a reusable joint-entropy/conditional-decoupling lemma and an exact matrix-quadratic amplifier.

The next nonlinear order still requires an actual finite grouped-current realization and a C2-uniform remainder. The source's formal smooth second-order field is not executed here; its inverse OU operator is understood on centered functions, and the constant does not affect the displayed gradient. Its inclusion supplies no higher-order producer authorization.

An arbitrary retained G observer, an inherited fine/coarse source pair, or a complete order-to-cost recurrence would require additional theorems. None is inferred from this PASS.
