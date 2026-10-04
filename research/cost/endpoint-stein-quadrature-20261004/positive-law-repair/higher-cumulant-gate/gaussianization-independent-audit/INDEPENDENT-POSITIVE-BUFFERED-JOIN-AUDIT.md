# Independent audit: positive buffered raw-source join

**Verdict: VALID as a source-qualified, positive, full-law return for the actual raw two-stage source.** It does not prove a posterior-law order-three result, a terminal-step result, or a complete improved-cost recurrence.

Audit date: 2026-10-04.

Reviewed source SHA-256: `f347c88a3987e50542dbb97e1e9d2cb1d8eb684cb6327509e89154d75efbb9aa`. The verdict applies to that frozen version.

## 1. Sources and scope

The audited composition is `CUBIC-GAUSSIANIZATION-AND-POSITIVE-BUFFERED-LAW-RETURN.md`, Sections 4–7, against:

- `TWO-STAGE-VALUE-REPAIR-OF-THE-THIRD-CUMULANT-WORD.md`, Sections 3–4: literal conditional-`v` source, mean-law service, caller and occurrence bounds.
- `reverse-ou-curvature/THIRD-ORDER-VALUE-CURVATURE-RESERVE.md`, Sections 5–6: full-gradient covariance action and its independent-buffer fluctuation/calibration estimate.
- The standalone Gaussianization lemma, independently proved in the companion audit.

The imported finite mean/action programs are treated as their established ports with all original guards. This audit verifies that this new composition meets the ports and propagates their errors; it does not re-prove or numerically execute those large imported programs. The manifest records hashes of the inspected versions.

Fix the exposed carrier `u`, external callers, finite versions, and all deterministic coefficients. Write `R_u=|u|+sqrt(D)`. The source facts used are

\[
\operatorname{Lip}_v K,\operatorname{Lip}_v K_*\le C\alpha,
\quad \|K_*-\mathbb E_vK_*\|_2\le C\alpha\sqrt D,
\quad \|E\|_2\le C\alpha^2 R_u,\quad E=K_*-K.
\]

The Lipschitz estimate is global in the fresh Gaussian innovation and uniform in the frozen `u`, as required by the new lemma. The source of the mean and covariance action is the executed finite VALUE graph, not an inferred function sharing its moments.

## 2. Covariance substitution is dimension-safe

For any vector-valued Gaussian-Lipschitz map `G`, scalar Gaussian Poincare applied to every projection gives `||Cov G||op <= Lip(G)^2`. The covariance difference decomposes exactly as

\[
\operatorname{Cov}(K_*)-\operatorname{Cov}(K)
=\operatorname{Cov}(E,K_*)+\operatorname{Cov}(K,E).
\]

For centered vectors `A,B`,

\[
\|\mathbb E[AB^\top]\|_{\mathrm{HS}}
\le \|A\|_2\|\operatorname{Cov}B\|_{\mathrm{op}}^{1/2}.
\]

Apply this to the two terms (transposing the second). Consequently

\[
\|\Sigma_*-\Sigma_K\|_{\mathrm{HS}}
\le C\alpha\|E-\mathbb EE\|_2
\le C\alpha^3R_u.
\]

Only one Euclidean marked energy occurs; there is no extra `sqrt(D)` loss.

## 3. Exact variance budget and reserve-action error

Let `0<h<=1/2`, `d=1-2h^2>=1/2`, and `eta=zeta=sqrt(d/2)`. Use the genuine-gradient anchored source

\[
J(v)=K(u,v)-K(u,0).
\]

For fixed `u`, subtracting a constant vector preserves the gradient property. Its origin is a literal caller-only VALUE record, and covariance is unchanged. Its radius and both its anchored and centered energies obey

\[
\ell\le C\alpha,\qquad e_J\le C\alpha\sqrt D.
\]

Set `mu=alpha`, use a fresh complete action tape and independent standard roots `p,z`, and execute

\[
R=\eta p+c\,C_{\rm cov}(J;p)+\zeta z,
\qquad c={h^2\over2\eta}.
\]

The imported action has target `Sigma_K p`, calibration error `Lambda ell e_J mu`, action energy `Lambda ell e_J`, and complete private first `Lambda ell^2/sqrt(mu)`.

1. Condition on `p`. Center only the private action fluctuation, keeping its private conditional mean. The imported independent Gaussian buffer estimate, with independent `zeta z`, costs
   `C c^2 ell^3 e_J / sqrt(mu)` after squaring and integrating over `p`.
2. Replace the actual conditional mean by `Sigma_K p`, at cost `C |c| ell e_J mu`.
3. The resulting deterministic linear Gaussian has covariance

\[
(\eta I+c\Sigma_K)^2+\zeta^2I
=dI+h^2\Sigma_K+c^2\Sigma_K^2.
\]

The last positive term is not dropped for free. The fixed-gap Gaussian-root bound prices it by `C c^2 ||Sigma_K^2||HS <= C c^2 ell^3 e_J`.

Because `eta,zeta` are bounded below, this proves

\[
W_2(\mathcal L(R\mid u),N(0,dI+h^2\Sigma_K))
\le\Lambda\sqrt D\{h^2\alpha^3+h^4\alpha^4
                         +h^4\alpha^{7/2}\}+\text{floors}.
\]

For `alpha<=1` and `h<=1/2`, both higher displayed terms are bounded by a constant times `h^2 alpha^3 sqrt(D)`. This is the claimed reserve estimate. Positive `Sigma_K` makes this reserve's linearized target automatically gapped by `d`; the stricter compiler-radius and finite-clock guards remain necessary.

## 4. Join with the mean law, then restore the actual raw law

The imported mean service is a whole output with comparison

\[
W_2(\mathcal L(M\mid u),N(\mathbb E_vK_*,I))
\le\Lambda\alpha^4R_u+\text{floors}.
\]

Execute the complete `M` and `R` banks independently after exposing the same `u`. Then

\[
T_{\rm comp}=hu-hM+R
\]

has Gaussian reference with mean `h(u-EK_*)` and covariance

\[
h^2I+dI+h^2\Sigma_K=(1-h^2)I+h^2\Sigma_K.
\]

The error to this reference is at most `h epsilon_M + epsilon_R`. The minus sign on `M` changes its mean as needed, without changing its covariance.

For positive definite matrices `A,B`, the same-root Gaussian coupling and Sylvester estimate give

\[
W_2(N(m,A),N(m,B))
\le\|A^{1/2}-B^{1/2}\|_{\mathrm{HS}}
\le {\|A-B\|_{\mathrm{HS}}\over
       \sqrt{\lambda_{\min}(A)}+\sqrt{\lambda_{\min}(B)}}.
\]

Thus replacing `Sigma_K` by `Sigma_*` costs at most `C h^2 alpha^3 R_u`, since the full covariance gap is `1-h^2>=3/4`. This bound does not require the covariance matrices to commute.

Finally use the independently verified Gaussianization lemma on the actual map

\[
X=-h(K_*-\mathbb E_vK_*),\qquad \sigma^2=1-h^2.
\]

Its Lipschitz constant and centered energy are `C h alpha` and `C h alpha sqrt(D)`. The error between its Gaussian reference and

\[
T_{\rm raw}=h(u-K_*(u,v))+\sqrt{1-h^2}N
\]

is at most `C h^3 alpha^3 sqrt(D)/(1-h^2)`. Since `h/(1-h^2)<=2/3` over the admitted interval, this is absorbed into `C h^2 alpha^3 R_u`.

All together,

\[
\boxed{W_2(\mathcal L(T_{\rm comp}\mid u),
                    \mathcal L(T_{\rm raw}\mid u))
\le\Lambda(h\alpha^4+h^2\alpha^3)R_u
       +\text{fully propagated absolute floors}.}
\]

Integrating the conditional comparisons while retaining the same Gaussian carrier `u` changes `R_u` to its `L^2` norm, at most `2 sqrt(D)`. It does not license re-attaching the private `v`, action tapes, or reference coupling noises.

## 5. Scaling, costs, and stopping boundary

For `s^2=1-r^2`, `Delta=t^2-r^2`, `alpha=As^2`, `v0=Delta/t^2`, and `h=sqrt(Delta)/s`, the admitted condition is exactly `Delta<=s^2/4`. Then

\[
h^2\alpha^3=\Delta A^3s^4,
\qquad
\text{physical error}\le
\Lambda\sqrt{v_0}(h\alpha^4+h^2\alpha^3)\sqrt D
+\text{scaled floors}.
\]

Both terms must be retained: `h alpha^4` cannot be absorbed into `h^2 alpha^3` uniformly for arbitrarily small increments. No inverse-`h` replication has been introduced by the join.

The conservative VALUE count

\[
Q_{\rm captured}+nN_K+(2n+1)N_E+nN_{\rm cov}
+Q_{\rm known/numerical}
\]

is consistent with the source ports. The anchored covariance calls reuse only the identical caller-only `K(u,0)` record; a changed private input requires a whole new `K` graph. Captured `f(u)` may lower the conservative `E` count but is not needed to justify it. Actual first/caller bounds must continue to be imported from the execution graph, never differentiated from this law comparison. The reserve residual first scales as `h^2 ell^2/sqrt(mu)=O(h^2 alpha^(3/2))`, as stated.

At `t=1`, `h=1` and the buffer/positive split used here disappears. That endpoint is excluded. The explicit raw rank-five posterior defect is also unaffected: Gaussianization compares the raw law to its moment-matched buffered Gaussian, and does not turn the raw law into the posterior. These restrictions are correctly retained in the source note.
