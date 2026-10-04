# Independent audit: covariance-matched Gaussianization with a Gaussian buffer

**Verdict: VALID under the stated Gaussian-Lipschitz hypotheses.** The constant is `sqrt(2)/3`, and a scalar buffer of standard deviation `sigma > 0` gives `sqrt(2)/(3 sigma^2)`. No symmetry of the Stein matrix, Hessian of the map, higher-moment oracle, or smallness assumption is needed.

Audit date: 2026-10-04. This is an independent derivation, not a numerical certification of the theorem. The accompanying checker tests orientations, tensor constants, and examples.

Reviewed source: `CUBIC-GAUSSIANIZATION-AND-POSITIVE-BUFFERED-LAW-RETURN.md`, SHA-256 `f347c88a3987e50542dbb97e1e9d2cb1d8eb684cb6327509e89154d75efbb9aa`, including its anisotropic-buffer remark.

## 1. Precise statement and conventions

Let `V ~ N(0,I_n)`, let `f:R^n -> R^d` be globally `C^1` and `L`-Lipschitz for Euclidean norms, and set

\[
X=f(V)-\mathbb Ef(V),\qquad e^2=\mathbb E|X|^2,\qquad
\Sigma=\mathbb E[XX^\top].
\]

Let `Z ~ N(0,I_d)` be independent of `V`, and let `sigma > 0`. Then

\[
\boxed{\quad W_2\bigl(\mathcal L(X+\sigma Z),
  N(0,\sigma^2I_d+\Sigma)\bigr)
  \le {\sqrt2\over3\sigma^2}\,L^2e.\quad}
\]

The Lipschitz constant here is an **operator-norm** derivative bound: `||Df(v)||op <= L`. A coordinatewise derivative bound is not an interchangeable hypothesis without its dimension factor.

Write `(Dh)_{ij}=partial_j h_i` and `A:B=sum_ij A_ij B_ij`. Thus the Stein identity is

\[
\mathbb E[X\cdot h(X)]=\mathbb E[\tau:Dh(X)].
\]

All matrix `L^2` norms below use the Hilbert-Schmidt norm. A Stein matrix may depend on `V`, rather than measurably on `X` alone. This is sufficient throughout.

## 2. Stein matrix with the required dimension-free bound

Let `mathcal L=Delta-v dot grad` be the standard Gaussian Ornstein-Uhlenbeck generator. For each centered coordinate `X_i`, take

\[
u_i=\nabla(-\mathcal L)^{-1}X_i.
\]

The weak Gaussian Riesz identity gives

\[
\mathbb E[X_i g(V)]=\mathbb E[u_i\cdot\nabla g(V)],
\qquad
\mathbb E|u_i|^2\le\mathbb E X_i^2.
\]

For completeness, expand `X_i` in orthonormal Gaussian Hermites: a chaos coefficient of order `k >= 1` contributes `1/k` times its squared size to `E|u_i|^2`. This proves the estimate and constructs the field in `L^2`, without requiring second derivatives of `f`.

Let `U` be the `d by n` matrix with rows `u_i^top`, and define

\[
\tau=U(Df)^\top,\qquad
\tau_{ij}=u_i\cdot\nabla f_j.
\]

The Sobolev chain rule applied to `h_i(f(V)-Ef(V))` proves the stated Stein identity. In particular, testing a linear function (or using Sobolev approximation) yields `E tau_ij = Sigma_ij`. Moreover,

\[
\|\tau\|_{L^2(\mathrm{HS})}
\le L\|U\|_{L^2(\mathrm{HS})}\le Le,
\qquad
\|\tau-\Sigma\|_{L^2(\mathrm{HS})}\le Le.
\]

In fact, centering gives the exact identity

\[
\|\tau-\Sigma\|_{L^2(\mathrm{HS})}^2
=\|\tau\|_{L^2(\mathrm{HS})}^2-\|\Sigma\|_{\mathrm{HS}}^2.
\]

There is no need to symmetrize `tau`. In general `tau^top` does **not** satisfy the same vector-valued Stein identity under the displayed convention.

## 3. Centered-matrix smoothing lemma

Let `B(V)` be any centered `d by d` matrix in `L^2`, let `S` be symmetric positive definite, and put `Y=tX+SZ` for `t >= 0`. Then

\[
\boxed{\quad
\left\|\mathbb E[B S^{-1}Z\mid Y]\right\|_2
\le \sqrt2\,tL\|S^{-1}\|_{\mathrm{op}}^2
       \|B\|_{L^2(\mathrm{HS})}.
\quad}
\]

Take first `h in C_c^infinity(R^d;R^d)`. Gaussian integration by parts in `Z`, conditional on `V`, gives

\[
\mathbb E[(BS^{-1}Z)\cdot h(Y)]
=\mathbb E[B:F(V)],\qquad
F(V)=\mathbb E_Z Dh(tX+SZ).
\]

Because `E B=0`, Cauchy-Schwarz and scalar Gaussian Poincare, summed over the entries of `F`, imply

\[
|\mathbb E[B:F]|
\le \|B\|_2\|F-\mathbb EF\|_2
\le \|B\|_2\|D_VF\|_2.
\]

The ordinary first-derivative chain rule, not a Hessian of `f`, gives

\[
(D_VF)_{ija}
=t\sum_k\mathbb E_Z[\partial_k\partial_jh_i(Y)]
                 \partial_a f_k(V).
\]

Flatten the pair `(i,j)` into the row index. Right multiplication by `Df` contracts its last index, so

\[
\|D_VF\|_2^2
\le t^2L^2\,\mathbb E\|\mathbb E_ZD^2h(Y)\|_{\mathrm{HS}}^2.
\]

For each fixed `V` and coordinate `i`, twice integrating by parts gives

\[
\mathbb E_ZD^2h_i(tX+SZ)
=S^{-1}A_i(V)S^{-1},\qquad
A_i(V)=\mathbb E_Z[h_i(tX+SZ)(ZZ^\top-I)].
\]

The normalized second-Hermite basis consists of `(Z_j^2-1)/sqrt(2)` and `Z_jZ_k`, `j<k`. The Hilbert-Schmidt norm counts each off-diagonal coefficient twice. Consequently,

\[
\|A_i(V)\|_{\mathrm{HS}}^2
=2\|P_2[h_i(tX+SZ)]\|_{L^2(Z)}^2
\le 2\mathbb E_Z|h_i(tX+SZ)|^2.
\]

It follows that

\[
\mathbb E\|\mathbb E_ZD^2h(Y)\|_{\mathrm{HS}}^2
\le2\|S^{-1}\|_{\mathrm{op}}^4\mathbb E|h(Y)|^2.
\]

Combining these estimates proves the dual bound against smooth compactly supported `h`. The conditional expectation already belongs to `L^2`, by independence and conditional Jensen:

\[
\mathbb E|BS^{-1}Z|^2
\le \|S^{-1}\|_{\mathrm{op}}^2\mathbb E\|B\|_{\mathrm{HS}}^2.
\]

Smooth compactly supported vector fields are dense in `L^2(Law(Y))` (the law has a positive smooth Gaussian-convolution density), so duality proves the asserted `L^2` bound. The dependence of `B` on `V` can be arbitrary; `B` is never differentiated.

## 4. Covariance-preserving interpolation and its velocity

Define, for `0 <= t <= 1`,

\[
C_t=\sigma^2I+(1-t^2)\Sigma,\qquad S_t=C_t^{1/2},
\qquad Y_t=tX+S_tZ,\qquad\mu_t=\mathcal L(Y_t).
\]

The endpoints are the two laws in the theorem, in reverse order. Every law on the path has covariance `sigma^2 I + Sigma`.

Here `S_t` and `Sigma` commute, since `S_t` is a spectral function of `Sigma`. Thus

\[
\dot S_t=-t\Sigma S_t^{-1}.
\]

For any smooth compactly supported scalar test function `phi`, differentiate using the common `(V,Z)` coupling:

\[
\begin{aligned}
{d\over dt}\mathbb E\phi(Y_t)
&=\mathbb E[(X-t\Sigma S_t^{-1}Z)\cdot\nabla\phi(Y_t)]\\
&=t\mathbb E[(\tau-\Sigma):D^2\phi(Y_t)]\\
&=\mathbb E[v_t(Y_t)\cdot\nabla\phi(Y_t)],
\end{aligned}
\]

where

\[
v_t(y)=t\mathbb E[(\tau-\Sigma)S_t^{-1}Z\mid Y_t=y].
\]

The first part of the middle equality uses the Stein identity applied to `x -> partial_i phi(tx+S_tz)`; the factor `t` is essential. The last equality uses Gaussian integration by parts in `Z`. Nonsymmetry of `tau` creates no problem: the full Stein identity is correctly oriented, and the scalar Hessian only sees its symmetric part.

**Terminology caveat:** this is a velocity satisfying the continuity equation. It need not equal the direct conditional derivative `E[dot Y_t | Y_t=y]`. Equality of their pairings with gradients is exactly what the transport argument needs.

The centered-matrix lemma with `B=tau-Sigma`, and `S_t >= sigma I`, gives

\[
\|v_t\|_{L^2(\mu_t)}
\le {\sqrt2\,t^2L\over\sigma^2}\|\tau-\Sigma\|_2
\le {\sqrt2\,t^2L^2e\over\sigma^2}.
\]

## 5. Wasserstein passage and constant

The common coupling shows `t -> mu_t` is continuous in `W_2`; `X` is square-integrable and the matrices `S_t` are continuously differentiable even at the endpoints, since `sigma > 0`. The preceding gradient identity is the distributional continuity equation, and its velocity has finite time-integrated `L^2` norm. The dynamic Wasserstein length inequality therefore gives

\[
W_2(\mu_0,\mu_1)
\le\int_0^1\|v_t\|_{L^2(\mu_t)}dt
\le {\sqrt2 L\over3\sigma^2}\|\tau-\Sigma\|_2
\le {\sqrt2\over3\sigma^2}L^2e.
\]

One may first apply the inequality on interior subintervals and pass to the endpoints using `W_2` continuity. No inverse of `Sigma` is used, so singular covariance and zero coordinates are harmless. If `e=0`, both endpoint laws are identical.

## 6. Applicability and boundaries

- The assumptions can be weakened to a globally Lipschitz map with its almost-everywhere weak derivative; the stated `C^1` assumption is more than sufficient.
- The theorem is conditional pointwise when, after the retained variables are frozen, the remaining innovation is a standard Gaussian `V`, `Z` remains independent, and the conditional map has the specified Lipschitz constant. Freezing an arbitrary variable correlated with `V` does not automatically preserve this hypothesis.
- The theorem concerns the **actual law** `f(V)`, not an approximate mean/covariance oracle. Approximation errors in the selected mean or covariance must be added separately.
- With a nonzero mean `m`, translation gives the same bound between `f(V)+sigma Z` and `N(m,sigma^2I+Sigma)`.
- A componentwise Lipschitz estimate, an unbuffered target (`sigma=0`), or a non-Gaussian innovation requires additional work; this proof must not be cited for those cases as written.
- The scalar-buffer interpolation formula for `dot S_t` uses commutation. Do not replace `sigma^2I` by an arbitrary covariance matrix in that formula without changing the argument.

**Bottom line:** the proposed lemma does provide a dimension-free cubic-amplitude Gaussianization error, `O(L^2 e / sigma^2)`, for a genuine Gaussian-Lipschitz map with a positive Gaussian buffer and exact mean/covariance matching.

## 7. Valid noncommuting-buffer extension

For any fixed symmetric `Q >= q I`, `q > 0`, the same proof yields

\[
W_2(\mathcal L(X+Q^{1/2}Z),N(0,Q+\Sigma))
\le {\sqrt2\over3q}L^2e.
\]

Set `C_t=Q+(1-t^2)Sigma` and `S_t=C_t^(1/2)`. The displayed scalar-buffer formula for `dot S_t` need not hold, but it is unnecessary. Gaussian covariance differentiation gives

\[
{d\over dt}\mathbb E\phi(tX+C_t^{1/2}Z)
=\mathbb E[X\cdot\nabla\phi(Y_t)]
 -t\mathbb E[\Sigma:D^2\phi(Y_t)].
\]

Equivalently, differentiate the square-root coupling and use

\[
(\dot S_t S_t):D^2\phi
={1\over2}(\dot S_t S_t+S_t\dot S_t):D^2\phi
={1\over2}\dot C_t:D^2\phi=-t\Sigma:D^2\phi.
\]

The Hessian is symmetric, which justifies the symmetrization even when `Q` and `Sigma` do not commute. The same velocity and centered-matrix estimate apply, with `||S_t^(-1)||op^2 <= 1/q`. There is no condition-number factor. A retained-caller-dependent `Q` is allowed after that caller is frozen.
