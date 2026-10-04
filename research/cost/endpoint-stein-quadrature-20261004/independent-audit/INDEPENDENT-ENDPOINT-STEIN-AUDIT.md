# Independent audit of the endpoint Stein packet

Date: 2026-10-04.

## Verdict and exact scope

**PASS for the bounded reusable component at the source pin below.** The positive dyadic time quadrature has a uniform moment error with logarithmic panel count and logarithmic Gauss order. Its Hermite multiplier bound is genuinely dimension-safe. The literal original-gradient packet has the stated private/caller firsts, numerical-zero qualification, finite-mode residual, complete query count, and two full Gaussian blocks. The first-order weak identity and uncancelled quadratic covariance ledger are correct.

This audit does **not** admit an all-order posterior source, a strong mean oracle, a strong fine/coarse pair, a retained/proxy theorem, or a sublinear complete-cost recurrence. Those are explicitly outside the inspected result. No counterexample or remaining mathematical blocker was found within its stated bounded scope.

Inspected source:

- `../ENDPOINT-STEIN-PACKET-AND-POLYLOG-QUADRATURE.md`
- SHA256: `4f0c182c36bf9dfbacc0f189899f689de52a2bb9bc37839e307eb76c4bda519f`

The audit assumes the normalized original Hessian sandwich `0 <= Hess V <= I`, `V` is C2, and `0 < A <= 1/2`. It also uses the stated original-gradient/directional-HVP numerical model. The baseline can retain its stricter heat guard. Earlier wording needed an explicit small-heat range and an exact numerical-zero convention; both are present in the audited pin.

The independent companion program is `check_endpoint_stein_independent.py`. It does not import, execute, or reuse the author's diagnostics. It passes **171,265 assertions**. Results are in `endpoint_stein_independent_checks.json`. These floating-point diagnostics supplement the analytical proof and are not an exhaustive numerical proof over infinitely many degrees, callers, dimensions, or potentials.

## 1 Uniform time quadrature with explicit constants

Write `h = 2^(-K)` with integers `K,m >= 1`. Every nonterminal panel is `[a,2a]`, with `0 < a <= 1/2`. On the Bernstein ellipse of parameter 2, write

    t = 3a/2 + (5a/8) cos(theta) + i (3a/8) sin(theta).

Its distance from its center is at most `5a/8`. Since `3a/2 <= 3/4`, the elementary triangle inequality already gives

    |1-t| <= |1-3a/2| + 5a/8 = 1-7a/8 <= 1.

This also verifies the source's equivalent real-endpoint argument. The same ellipse parameter works for all panels, including the largest one. Consequently `f_k(t)=(1-t)^k` is holomorphic and bounded by one throughout that ellipse for every integer `k >= 0`. No derivative of the force is being bounded here.

After affine transport of the panel to `[-1,1]`, its Chebyshev coefficients satisfy `|c_j| <= 2*2^(-j)`. Truncation through degree `2m-1` therefore has uniform error at most

    sum_(j=2m)^infinity 2*2^(-j) = 4*4^(-m).

The true integral and the positive Gauss rule each have mass `a` and agree on that truncating polynomial. Their difference is thus at most `8 a 4^(-m)`. Summing panel lengths gives at most `8(1-h)4^(-m)`, and hence at most `8*4^(-m)`.

The terminal interval and its positive midpoint approximation each have absolute value at most `h`, so their difference is at most `2h`, as stated. Since the real integrand is nonnegative, a sharper `h` bound is also available but unnecessary. Thus the exact advertised sufficient bound is

    sup_(k>=0) |sum_i w_i r_i^k - 1/(k+1)|
      <= 8*4^(-m) + 2*2^(-K).

One completely explicit schedule for `0 < delta <= 1` is

    m = ceil(log_4(16/delta)),
    K = ceil(log_2(4/delta)),
    n = Km+1.

Each of the two error terms is at most `delta/2`. Constants and linear polynomials are integrated exactly because even one-point Gauss and the terminal midpoint are exact on them. Therefore the weight and first-moment identities are correct.

For fixed maximum order J and `delta=A^(J+2)` times a prescribed precision/log factor, the count is logarithmic squared in the relevant precision parameter. It carries no inverse-A query exponent. This statement counts original queries, not constant-time arithmetic, storage, root computation, or Gaussian generation.

An additional check, not needed by the proof, is that exact Gauss-Legendre integration on these panels underestimates every monomial `(1-t)^k` once its degree exceeds exactness; the terminal midpoint does likewise for `k>=2`. The diagnostics verify this sign at their tested integer degrees. The infinite-degree claim rests on the ellipse proof above, not those samples.

## 2 Hermite vector norm and its precise probability space

On `L2(gamma;R^D)`, Mehler is the scalar Hermite multiplier `r^k` on each total-degree chaos, componentwise. Both the integral and finite positive sum define bounded operators, so their difference T has

    T g_k = e_k g_k,
    e_k = sum_i w_i r_i^k - 1/(k+1),
    ||T||_(L2 -> L2) = sup_k |e_k|.

Orthogonality of vector chaoses means that one sums the squared Euclidean coefficient norms once. There is no coordinatewise union bound or extra D multiplier. Therefore

    ||v_Q-v||_L2 <= delta ||g_M||_L2 <= delta A sqrt(D).

The anchored inequality `|g_M(z)| <= A|z|` holds for every finite center, so the estimate is uniform in the *exposed physical caller* through that anchored norm. It is not a pointwise or uniform-in-Z approximation statement.

**Essential interpretation:** `v_Q(Z)=E_G[H_Q(Z,G) | Z]` is an exact analytical conditional expectation. The displayed norm integrates the remaining standard Gaussian Z. Neither that conditional expectation nor v is queried by the finite program. Q2 is not a guarantee that the sampled vector H_Q is strongly close to v. In the quadratic example,

    H_Q = (B/2) Z + beta_Q B G,
    v = (B/2) Z,
    ||H_Q-v||_L2 = beta_Q ||B||_HS.

The last quantity stays of order `A sqrt(D)` for isotropic `B`, even as quadrature precision tends to infinity. This is a concrete obstruction to interpreting Q2 as a strong sampled-mean port. Sampling many independent G values to approximate that mean would require its own separately paid mean-estimation cost.

## 3 Stein identity and regularity actually used

Anchoring gives `U_M >= 0`, `U_M(z) <= A|z|^2/2`, and `|g_M(z)| <= A|z|`. These bounds suffice for Gaussian integrability. The centered inverse Ornstein-Uhlenbeck potential, with the sign convention `L=Delta-z dot grad`, has gradient

    v = integral_0^infinity exp(-t) P_(exp(-t)) g_M dt
      = integral_0^1 P_r g_M dr.

Gaussian integration by parts gives

    E_gamma[v dot grad(phi)] = Cov_gamma(U_M,phi)

for the asserted sufficiently regular integrable tests, with extension by usual approximation where appropriate. Differentiating only once in z gives

    Dv = integral_0^1 r P_r(Dg_M) dr,
    0 <= Dv <= (A/2) I.

Gaussian invariance and Minkowski give the stated `C_p A sqrt(D)` moment bound for each fixed finite p. No smoothness of V beyond the original C2 Hessian bound is used, and no analyticity of the original force at the endpoint is implied.

When the implemented center is b_M, both the Stein identity and W must use the finite-center `U_M,g_M`. They refer initially to the anchored target with its linear tilt removed. The true posterior is restored as in the next section; it is not identical to that intermediate target.

## 4 Finite mode and the missing linear tilt

Let `T(x)=y-A grad V(x)`, `b_(j+1)=T(b_j)`, and `b_0=y`. Its Lipschitz constant is at most A. The exact finite residual is

    R_M = b_M-y+A grad V(b_M) = b_M-b_(M+1),
    |R_M| <= A^(M+1) |grad V(y)|.

Thus the source's constant can be taken as one for exact original-gradient evaluations:

    |ell| = |R_M|/sqrt(A) <= A^(M+1/2)|grad V(y)|.

Original finite-query errors add their separately propagated numerical allowances.

Expanding the original posterior potential at b_M gives, after deleting constants,

    |z|^2/2 + U_M(z) + ell dot z.

This verifies the sign and scaling of the residual linear term. Set `F(z)=|z|^2/2+U_M(z)`, so `Hess F >= I`. Synchronously couple overdamped diffusions with drifts `-grad F` and `-grad F-ell`. Their separation satisfies `d|Delta|/dt <= -|Delta|+|ell|`. Passing to stationary laws yields a coupling with separation at most `|ell|`, in particular

    W2(anchored standardized target, true standardized target) <= |ell|.

After physical scaling the price is at most `sqrt(A)|ell|=|R_M|`. The same strong monotonicity of `x-y+A grad V(x)` gives `|b_M-b|<=|R_M|`. The packet's deterministic zero is b_M, not the exact mode. This residual must retain its caller-dependent `|grad V(y)|` factor; the source does so rather than turning it into an intrinsic `sqrt(D)` law error.

## 5 Actual finite private and caller firsts

Let `B_M = D_y b_M` and `H_x=Hess V(x)`. The finite recurrence gives

    B_(j+1) = I-A H_(b_j) B_j,
    ||B_M|| <= 1/(1-A),
    ||B_M-I|| <= A/(1-A).

Write `R_i=(r_i I,s_i I)`, with `s_i=sqrt(1-r_i^2)`. Each row R_i has operator norm one. Since

    D_(Z,G) H_Q = A sum_i w_i H_(b_M+sqrt(A)q_i) R_i,

positivity and unit total mass imply `||D_(Z,G)H_Q|| <= A`. This proves the full private first bound without treating different nodes as independent Gaussian inputs. The exact joint moment estimate follows from Minkowski and `q_i~gamma` separately; node independence is unnecessary.

For a physical caller direction,

    D_y g_M(q_i) = sqrt(A) [H_(b_M+sqrt(A)q_i)-H_(b_M)] B_M.

Both Hessians lie in `[0,I]`, so their difference has operator norm at most one, regardless of its continuity modulus. Therefore

    ||D_y X_Q-I|| <= 2A/(1-A) <= 4A,
    ||D_(Z,G) X_Q-sqrt(A)(I,0)|| <= A sqrt(A).

For `F_Q=sqrt(A) grad V(X_Q)`, explicit sufficient constants are

    ||D_(Z,G) F_Q|| <= A(1+A),
    ||D_y F_Q|| <= sqrt(A) [1+2A/(1-A)].

Each requested forward or adjoint action differentiates the actual finite graph once and calls original HVPs at stored gradient sites. Nothing requires a third derivative or differentiation of a saved HVP. Bounds on actual firsts should not be confused with arbitrarily accurate derivative *approximation* when query locations are numerically perturbed; the original first/numerical model supplies those separate allowances.

## 6 Zero query count weights and numerical profiles

With M mode updates, there are M mode VALUE calls, one b_M anchor, n node calls, and one terminal canonical-force call: `M+n+2` in the non-deduplicated complete ledger. The output X_Q without its canonical force omits the last call. An original directional forward/adjoint sweep has the same list of possible HVP sites, subject to legitimate same-site reuse. An adjoint may first accumulate multiple outgoing contributions before acting once at a shared stored site. No dense Jacobian is being counted as one HVP.

The implementation uses exactly two D-dimensional independent Gaussian blocks after exposing y. All q_i are deterministic reads of those blocks. At zero input, every q_i is zero. The revised source explicitly freezes a deterministic numerical-gradient version and reuses the same recorded anchor at the identical site. Then every difference is literally zero, and the finite numerical packet returns b_M. Separately biased evaluations at the same site do not cancel exactly; the revised source correctly prices that alternative as a numerical zero residual.

Positive mass one propagates common force biases with their actual absolute weights. It gives no `1/sqrt(n)` noise reduction. A useful explicit pathwise query envelope is

    |b_j-y| <= [A/(1-A)] |grad V(y)|,
    max_i |q_i| <= sqrt(|Z|^2+|G|^2),
    |H_Q| <= A sqrt(|Z|^2+|G|^2).

Thus all node, mode, and terminal query displacements are bounded by a constant times

    A|grad V(y)| + sqrt(A) sqrt(|Z|^2+|G|^2).

Lipschitzness of the original gradient gives the asserted caller profile after Gaussian moments. There is no dimension-dependent maximum over n independent blocks and no algebraic D-versus-A domain guard.

Encoded quadrature identities have to be priced as the source states. For example, total weight error eta_w and encoded nodes in `[0,1]` satisfying `|rhat_i-r_i| <= eta_r(1-r_i)`, `eta_r<1`, change every moment by at most

    eta_w + eta_r/(1-eta_r),

using `k u^(k-1) <= 1/(1-u)`. Since the smallest `1-r_i` is at least `2^(-K-1)`, enough absolute node precision requires only `O(K+log(1/eta_r))` bits. Endpoint square roots can likewise be assigned stricter finite precision using their elementary half-Hölder bound. This does not create an inverse-heat original-query exponent.

## 7 Weighted retained-coordinate factorization

For `d_i=sqrt(w_i)`, stack the rows `d_i R_i`. For every `(z,g)`,

    sum_i w_i |r_i z+s_i g|^2 <= |z|^2+|g|^2.

So the stack norm is at most one. The terminal row with blocks `d_i I` has norm one. In weighted node coordinate `u_i=d_i q_i`, the node output `d_i g_M(u_i/d_i)` is the gradient of `d_i^2 U_M(u_i/d_i)` and has Hessian bounded by A. Its original physical gradient evaluation has the stated inverse-width factor `1/d_i` relative to the standardized node.

These are genuine elementary graph facts. They do not turn the rectangular field `H_Q:R^(2D)->R^D` into a gradient on the entire Gaussian record and do not prove any native retained, proxy, protection, or later-observer theorem. The revised source explicitly respects this limitation.

For completeness, an elementary conservative width estimate follows from `|P_m|<=1` and the Markov derivative bound `|P'_m|<=m^2` on `[-1,1]`: the mapped Gauss weights satisfy `w_i>=2^(-K)/m^4`. Hence `1/d_i<=m^2 2^(K/2)`. Such widths can be inverse fixed-order heat powers, but their logarithms are a precision overhead. Pricing a future use of these widths is a separate obligation of that future graph.

## 8 Weak identity and quadratic remainder

Taylor's exact integral formula for `psi(Z-H_Q)` yields W once one replaces `E_G H_Q` by `v+(v_Q-v)` and applies the Stein identity. The sign of the covariance and of the quadrature-defect term is correct. The rank-two remainder is a same-record current containing all shared-G cross-node products.

For `g(z)=Bz`, the first moment identity gives exactly

    H_Q = (B/2) Z + beta_Q B G,
    Cov(Z-H_Q) = I-B+(1/4+beta_Q^2) B^2.

The reference covariance is `(I+B)^(-1)`. As K and m increase, `beta_Q` converges to the elementary integral `integral_0^1 sqrt(1-r^2)dr=pi/4`, so the limiting second-order coefficient is

    1/4+pi^2/16 = 0.8668502750680849,

rather than one. Refining time quadrature cannot remove that defect. Independent G at the different nodes would produce a different covariance involving `sum_i w_i^2(1-r_i^2)`, so that substitution would alter the implemented law.

The warning about a single uniform random r is also justified. For `B=aI`, the conditional scale is `s(r)=1-2ar+a^2`. As D tends to infinity, the normalized output radius tends to `sqrt(s(r))`, whereas the target normalized radius tends to `(1+a)^(-1/2)`. Their limiting squared radial discrepancy is

    integral_0^1 [sqrt(1-2ar+a^2)-(1+a)^(-1/2)]^2 dr
      = a^2/12 + O(a^3).

Thus one random time retains an order `a sqrt(D)` high-dimensional obstruction. The deterministic Q removes this random first-order scale, while preserving the displayed nonzero deterministic second-order covariance defect.

The crude weak remainder bound uses `E|H_Q|^2=O(A^2 D)`. Nothing in the admitted quadrature certificate alone improves this to the one-energy or strong-law estimates required by an all-order compiler.

## 9 Independent finite diagnostics

The companion checks cover:

- 42 combinations of K and m, up to 337 nodes, and sampled integer Hermite degrees through `10^12`, including endpoint-scaled degrees and exact low-moment identities.
- The common parameter-2 ellipse at four panel scales, including `2^(-50)`.
- Vector chaos multiplier identities in dimensions 1, 3, 17, and 257, including coefficient arrays attaining the sampled spectral norm.
- Literal nonlinear finite packets with noncommuting bounded Hessians in dimensions 1, 3, and 11; four heats; four mode counts; query census; private/caller analytic Jacobians; and independent directional finite differences.
- A genuinely C2 potential whose Hessian is continuous but not Lipschitz at zero.
- Exact quadratic mode-residual restoration in dimensions 1, 3, and 23.
- Shared-root quadratic covariance and deliberately different independent-node-root covariance in dimensions 1, 5, and 41.
- Weighted graph norms, finite inverse-width bounds, nonvanishing sampled bridge error, common numerical bias, and the failure of separately perturbed anchor values to cancel.
- A scalar nonlinear Stein identity and the full weak Taylor identity, evaluated by independent Gaussian quadrature. In the two fixtures, weak-identity residuals are 0 and `-2.22e-16`; the retained second-order currents are approximately `-0.0023651` and `-0.0022930`.
- Finite one-dimensional integrals displaying the random-time radial coefficient tending to `1/sqrt(12)`.

These checks make no empirical claim of posterior source order, general retained-port closure, or complete algorithmic cost. The analytical and implementation qualifications above are the audit boundary.
