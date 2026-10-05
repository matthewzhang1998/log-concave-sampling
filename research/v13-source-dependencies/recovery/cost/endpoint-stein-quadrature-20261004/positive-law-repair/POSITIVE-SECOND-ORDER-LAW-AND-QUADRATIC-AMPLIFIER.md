# Positive second-order endpoint law and a finite quadratic amplifier

2026-10-04. New result, pending independent review. This is not an all-order general-C2 posterior source, a paired-source theorem, or a sublinear complete-cost recurrence.

## Result

The already executed endpoint Stein packet has a stronger dimension-safe law certificate than its original first-order weak identity:

    W2(Law(Z-H_Q(Z,G)), mu_U) <= C (A^2 + delta A) sqrt(D).

Here `g=grad U`, `g(0)=0`, `0<=Dg<=A I`, `0<A<=1/2`, the exact target is `mu_U proportional to exp(-|z|^2/2-U(z))`, and the same two independent D-dimensional roots Z,G and the same positive dyadic quadrature Q are used. Its conditional mean `v_Q` has the previously proved Gaussian L2 error `||v_Q-v||_2<=delta A sqrt(D)`. No conditional expectation, derivative action, density, or entropy is evaluated by the finite program.

Thus choosing delta<=A gives physical posterior law error `C sqrt(D) A^(5/2)`, plus the explicitly restored finite-mode and numerical allowances, with only polylogarithmic original-gradient VALUE work. All actual first/caller/zero/query facts in the audited packet remain unchanged. This upgrades a weak first-order calculation to a positive W2 theorem, but is still only one fixed general-C2 law grade.

For quadratic potentials there is additionally an explicit finite VALUE-only amplifier that reaches arbitrary fixed order with O(J) extra original gradient values. It corrects the actual shared-G covariance, including every cross-node product. The nonlinear continuation of this amplifier is not proved.

## 1. Reusable Gaussian conditional-decoupling lemma

Let Z~N(0,I_d), G~N(0,I_n) be independent. Let H:R^(d+n)->R^d be C1, with

    ||D_Z H||_op <= a < 1,
    Lip_G H <= b,
    E_H = (E|H|^2)^(1/2) < infinity,
    J_H = (E||D_Z H||_HS^2)^(1/2) < infinity.

Put `v(z)=E_G H(z,G)`, `Y_t=Z-tH(Z,G)`, and let nu_t=Law(Y_t), for 0<=t<=1. Every fixed-G map `z -> z-tH(z,G)` is a globally invertible C1 map, because its departure from the identity has Lipschitz constant at most ta<1. Its determinant is positive, even if D_ZH is not symmetric.

The exact Gaussian change-of-variables identity is

    KL(Law(Y_t,G) || gamma_d tensor gamma_n)
      = E[-t Z dot H + t^2 |H|^2/2 - log det(I-t D_ZH)].       (1)

Gaussian integration by parts in Z gives `E[Z dot H]=E tr(D_ZH)`. Expanding logdet and using

    |tr B^k| <= ||B||_op^(k-2) ||B||_HS^2, k>=2,

therefore yields

    KL <= (t^2/2) [E_H^2 + J_H^2/(1-ta)].                     (2)

This is a full joint calculation. It does not replace common roots by independent node roots. By the entropy chain rule,

    I(Y_t;G) <= KL(Law(Y_t,G) || gamma_d tensor gamma_n).

Conditional Gaussian Talagrand transport gives

    E_(Y_t) W2^2(Law(G|Y_t),gamma_n) <= 2 I(Y_t;G).           (3)

These inequalities remain valid by approximation whenever the displayed C1/Sobolev and integrability assumptions hold. In the application below all derivatives are globally bounded and H has linear growth, so the change of variables and Gaussian integrations have no tail or invertibility gap.

Define the ACTUAL marginal velocity coefficient

    h_t(y)=E[H(Z,G)|Y_t=y].

The curve nu_t satisfies the continuity equation with velocity -h_t. Its difference from the mean field v is bounded by splitting

    H(Z,G)-H(Y_t,G),
    E[H(Y_t,G)|Y_t]-E_G H(Y_t,G).

The first has L2 norm at most `ta E_H`; the second is bounded using the b-Lipschitz G dependence and (2)-(3). Consequently

    ||h_t-v||_L2(nu_t)
      <= t [a E_H + b sqrt(E_H^2+J_H^2/(1-ta))].             (4)

Let Phi_t solve `d Phi_t(z)/dt=-v(Phi_t(z))`, Phi_0(z)=z. Since Lip(v)<=a, the usual continuity-equation/Wasserstein stability estimate gives

    W2(nu_1, Law(Phi_1(Z)))
      <= integral_0^1 exp(a(1-t)) ||h_t-v||_L2(nu_t) dt.     (5)

For clarity, this is only an analytical comparison: pull nu_t back through the deterministic remaining flow Phi_(t,1); its velocity is the remaining-flow derivative applied to `-h_t+v`, whose operator norm costs at most exp(a(1-t)). The dynamic W2 bound then gives (5). It is not necessary to execute h_t, conditional sampling, or the flow.

Also, because `||v(Z)||_2<=E_H`, the deterministic flow differs from its Euler endpoint by

    ||Phi_1(Z)-(Z-v(Z))||_2 <= (a exp(a)/2) E_H.              (6)

Equations (4)-(6) are a closed reusable positive-law lemma. In particular if `a,b<=C epsilon` and `E_H,J_H<=C epsilon sqrt(d)`, marginalizing the full private G bank costs `O(epsilon^2 sqrt(d))` in W2. The derivative energy J_H is a real hypothesis. Small VALUE energy alone does not make it small, and this lemma does not grant arbitrary joint retention of G after the marginal comparison.

## 2. Apply to the actual endpoint packet

The packet is

    H_Q(Z,G)=sum_i w_i g(r_i Z+s_i G), s_i=sqrt(1-r_i^2),
    w_i>0, sum_i w_i=1, sum_i w_i r_i=1/2.

All q_i use the SAME G. Its actual bounds are

    0 <= D_Z H_Q <= (A/2) I,
    Lip_G H_Q <= A beta_Q <= A,
    beta_Q=sum_i w_i s_i,
    ||H_Q||_2 <= A sqrt(D),
    ||D_Z H_Q||_L2(HS) <= (A/2) sqrt(D).

The last bound uses the actual full Hessian matrices only in analysis. No matrix is formed by the VALUE producer. Thus (2) gives, for A<=1/2,

    KL(Law(Y_t,G)||gamma tensor gamma) <= (2/3) t^2 A^2 D,
    ||h_t-v_Q||_L2(nu_t) <= (1/2+2/sqrt(3)) t A^2 sqrt(D).

Equations (5)-(6) prove

    W2(Law(Z-H_Q),Law(Z-v_Q(Z))) <= C A^2 sqrt(D).           (7)

One can use C=3 in (7); retaining an unspecified universal C is sufficient here. This pays the shared-root fluctuation through its actual conditional information, rather than asserting that H_Q is strongly close to v_Q. In the quadratic test that strong claim is false: `||H_Q-v_Q||_2=beta_Q ||B||_HS=O(A sqrt(D))`.

The audited quadrature then gives the literal same-Z coupling

    ||(Z-v_Q(Z))-(Z-v(Z))||_2 <= delta A sqrt(D),
    v(z)=integral_0^1 P_r g(z) dr.                          (8)

## 3. Compare the Gaussian Stein Euler map to the exact C2 target

Here is a complete dimension-safe first-order transport comparison; no high derivative of U is required. Let

    f_r(x)=P_r(e^-U)(x)/E_gamma e^-U,
    rho_r(dx)=f_r(x) gamma(dx), 0<=r<=1.

Then rho_0=gamma, rho_1=mu_U. Equivalently rho_r is the law of `rX+sqrt(1-r^2)N` for independent X~mu_U and standard N. Since `x dot g(x)>=0`, target integration by parts gives

    E_mu |X|^2 <= D,
    E_(rho_r)|x|^2 <= D.                                  (9)

Let s=sqrt(1-r^2), and define, analytically,

    m_r(x)= E[g(rx+sG) exp(-U(rx+sG))]
                /E[exp(-U(rx+sG))].

The OU differential identity shows that rho_r satisfies the continuity equation with velocity `-m_r`: indeed `partial_r f_r=-(1/r)L f_r` and `grad log f_r=-r m_r`. The r=0 identity is understood by its continuous limit. At r=1, m_1=g.

For fixed r,x, let nu_(r,x) be the tilted law of G with density proportional to exp(-U(rx+sG)) relative to standard Gaussian. Its potential has Hessian `I+s^2 Dg>=I`. Synchronous Langevin coupling to standard Gaussian, placing the perturbing force on the Gaussian coordinate and using monotonicity of the tilted drift, gives

    W2(nu_(r,x),gamma) <= s ||g(rx+sG)||_L2(gamma).

The same bound follows by the standard strong-convexity drift comparison; importantly it uses the perturbation under gamma, so there is no unproved moment bound at a translated mode. Since g is A-Lipschitz,

    |m_r(x)-P_r g(x)|
      <= s A W2(nu_(r,x),gamma)
      <= s^2 A^2 sqrt(r^2 |x|^2+s^2 D).                    (10)

After (9), this is at most `s^2 A^2 sqrt(D)` in L2(rho_r).

The bounded firsts needed for flow comparison also follow from the literal tilted expectations. Differentiation once gives

    D_x m_r = r [E_nu Dg(rx+sG)-Cov_nu(g(rx+sG))].

The strongly log-concave nu has Poincare constant at most one, so its covariance term has operator norm at most s^2 A^2. Hence `Lip(m_r)<=r(A+s^2A^2)`. Only Dg, which exists and is bounded under C2, is used. Likewise `Lip(P_r g)<=rA`.

Let T_r solve `dT_r/dr=-(P_r g)(T_r)`, T_0=Z. Stability against rho_r and (10) gives

    W2(Law(T_1),mu_U) <= exp(A/2) A^2 sqrt(D).              (11)

Finally compare this time-dependent flow with its frozen-input Euler endpoint `Z-v(Z)`. Minkowski, the Gaussian invariance of P_r, and `Lip(P_r g)<=rA` give

    ||T_r-Z||_2 <= r exp(A/2) A sqrt(D),
    ||T_1-(Z-v(Z))||_2
       <= integral_0^1 r A ||T_r-Z||_2 dr
       <= [exp(A/2)/3] A^2 sqrt(D).                       (12)

Combining (7),(8),(11),(12) proves the result. For example the safe absolute constant C=10 covers all terms for A<=1/2. The proof is uniform in all genuine C2 Hessian-sandwich potentials, including continuously varying, noncommuting Hessians without any Lipschitz modulus.

## 4. Exact finite-mode, original-query and source scope

Use exactly the finite-mode anchoring and numerical version from the audited endpoint packet. At b_M the anchored target has `g_M(0)=0` and the same theorem applies uniformly. The true standardized posterior has the residual linear tilt ell, so its physical restoration cost remains `sqrt(A)|ell|`, with

    |ell| <= A^(M+1/2)|grad V(y)|

before the separately propagated original numerical allowances. The exact packet zero is b_M by same-site anchor reuse. It is not the exact mode.

No program changes occur in this general-C2 theorem. The complete VALUE count is still `M+n+2` when the canonical terminal force is requested, and the Gaussian dimension is still exactly 2D. A requested first/adjoint still uses only original HVPs at the recorded VALUE sites, with the same complete count. The physical endpoint has caller first I+O(A), private first sqrt(A)[(I,0)+O(A)], and the canonical force has private first O(A), caller first O(sqrt(A)). The positive finite output remains the same literal map. Numerical floors are coupled to it, not to an unimplemented conditional mean.

No derivative of a saved HVP occurs. No protected late Gaussian, native retained/proxy premise, fine/coarse higher-order pair, or arbitrary later observer of G is supplied by this marginal law theorem. In particular (7) is not a joint comparison retaining G. If a caller is random, expose it before Z,G, use this conditional theorem, and integrate its separately restored numerical profile at that actual caller.

Taking quadrature delta<=A is enough for this new physical grade 5/2. Choosing delta=A^(J+2) for a future common fixed-target packet remains polylogarithmic, but does not promote the O(A^2) nonlinear law term to order J.

## 5. A VALUE-only all-order amplifier for the exact matrix quadratic test

Let `g(x)=Bx`, with B symmetric, `0<=B<=A I`. The already audited packet has

    Y=(I-B/2)Z-beta_Q B G,
    C(B)=Cov(Y)=I-B+c B^2,
    c=1/4+beta_Q^2.

There is no independence substitution inside beta_Q^2. Since the concave function sqrt(1-r^2) lies above 1-r and Jensen bounds its weighted mean from above,

    1/2 <= beta_Q <= sqrt(3)/2,
    1/2 <= c <= 1.

Define the scalar matrix polynomial

    D(B)=(I+B) C(B)-I=(c-1)B^2+cB^3,
    q_m(B)=sum_(k=0)^m binom(-1/2,k) D(B)^k.

For A<=1/2, `||D(B)||<=A^2<=1/4` (use `|(c-1)+cb|<=1` for 0<=b<=A). The binomial series is uniformly convergent there, and its tail obeys

    ||q_m(B)-(I+D(B))^(-1/2)|| <= C_m A^(2m+2).             (13)

Execute `Y_m=q_m(B)Y` by the following literal VALUE graph:

    T_0=Y;
    for k=1,...,m:
        u_1=g(T_(k-1)); u_2=g(u_1); u_3=g(u_2);
        T_k=(c-1)u_2+c u_3;
    Y_m=sum_(k=0)^m binom(-1/2,k) T_k.

Every displayed g call is the original anchored gradient VALUE query, not B as a known oracle. On this quadratic specialization the graph equals the polynomial exactly, with three extra original values per layer. The anchor is reused and each root/ancestor remains in the graph. No new Gaussian block is introduced. Query displacements, finite precision and the same caller/mode error are separately restored; the coefficients depend only on the already public quadrature.

The ACTUAL covariance is

    Cov(Y_m)=q_m(B)^2 C(B).

Since every matrix here is a polynomial in the same symmetric B, simultaneous spectral diagonalization is legitimate. Coupling the two Gaussian laws along common eigenvectors yields

    W2(Law(Y_m),N(0,(I+B)^(-1))) <= C_m A^(2m+2) sqrt(D).   (14)

The corresponding physical law order is `2m+5/2`. This is an exact matrix test in every dimension, not a scalar-only trace match. Original symmetric B may be arbitrary; all root covariance descendants are present. The polynomial firsts and same-site zero have the usual bounded actual VALUE-graph certificates for fixed m. This specialization does not make B a producer-side HVP.

For a general nonlinear C2 g the same displayed graph is executable, but g(g(x)) is not the required local Hessian action Dg(x)g(x), and different local Hessians need not commute. Neither (13) nor (14) applies there. A quadratic covariance repair is not a general second-order-current compiler.

## 6. What a further order gain would need

The law comparison in Section 1 is reusable: a higher current with BOTH amplitude/derivative scale epsilon and a conditional-mean field of size epsilon can be marginalized at W2 cost O(epsilon^2 sqrt(D)). It preserves its actual common-root covariance through conditional information.

It does not furnish that current. For the exact tilt interpolation a formal smooth second-order transport correction contains

    (1/2)[Dv(z) v(z) + grad(-L)^(-1)(g dot v)(z)],
    v=grad(-L)^(-1)(U-EU).

This formula identifies where local Hessian-vector and shared-root products arise; it is not an authorized HVP-valued producer and is not used in the theorem. A successful finite continuation must realize its entire grouped current, its later descendants, and a C2-uniform remainder with actual VALUES and priced finite firsts. Replacing this by g(g(z)), estimating a conditional mean strongly, or keeping only a covariance trace does not do so.

The proved general-C2 grade is fixed. Therefore the present result establishes no recurrence `c(P)/P -> 0`. It does establish a complete polylogarithmic positive known-center law source at physical order 5/2, and an exact all-order quadratic falsification standard for candidate continuations.

## 7. Finite diagnostics and sharpness

`check_positive_stein_law.py` writes `positive_stein_law_checks.json`. All 468 author checks pass. They cover positive quadrature identities; literal shared-root, non-diagonal matrix covariance; 120 quadratic dimension/heat/order cases through m=5; the exact three-query polynomial recurrence; covariance and Gaussian W2; a genuine C2 non-C3 two-dimensional ridge potential with continuously varying noncommuting Hessians; actual Z-first and one-energy bounds; the pointwise matrix-logdet remainder ledger; and scalar conditional-CDF W2 diagnostics for the actual positive map.

The scalar C2 fixture uses `h'(t)=sqrt(|t|)/(1+sqrt(|t|))`; at A=1/4,1/8,1/16 its computed W2/A^2 values are approximately 0.00864, 0.01662, 0.02100. These numerical CDF calculations are not a proof of asymptotic order. The tensor-Gaussian quadrature's nonzero integration-by-parts residual at the C2 cusp is explicitly reported (largest about 0.000913), rather than silently set to zero; the exact cancellation in the theorem is analytical. Firsts used to invert the scalar conditional map are diagnostic calculations, not extra producer actions.

For isotropic quadratic B=A I, the unamplified packet has exact W2

    sqrt(D) |sqrt(1-A+cA^2)-(1+A)^(-1/2)|
      = ((1-c)/2) A^2 sqrt(D) + O(A^3 sqrt(D)).

With an increasingly accurate quadrature `c -> 1/4+pi^2/16 < 1`, so the new general bound has the correct worst-case power already on quadratics. Finer endpoint quadrature alone cannot raise it. The polynomial amplifier removes that specific debt, while the nonlinear extension remains open.
