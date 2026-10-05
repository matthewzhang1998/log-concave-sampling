# Independent exact-text audit: canonical same-g radial-shell separator

Date: 2026-10-05 (UTC).

## Verdict and exact scope

**PASS.** The canonical proof establishes a genuine same-original-gradient counterexample to the **unshifted single-history formula**

\[
R_1(j_H-j_Q)=\tfrac12R_1[(C_H-C_Q):D^2g(x)]
                 +O(\Lambda A^4\sqrt D).
\]

For its smooth convex-gradient family, with \(A=D^{-1/2}\), the discrepancy is bounded below by \(c_*A^2\) for all sufficiently large \(D\), uniformly over the finite positive rules satisfying the stated assumptions. The proposed nonlogarithmic allowance is \(A^3\). A finite positive rule with polynomially small tolerance and polylogarithmic node count exists, so fixed public logarithmic factors do not rescue the formula.

There is no mathematical defect requiring correction to the reviewed canonical text. In particular, its bounded nonlinear history terms are retained and controlled; neither their correlation with the linear Gaussian history nor their contribution to the actual covariance is omitted.

This is **not** a refutation of the full \(m_3\) comparison, the literal finite \(F_{3,Q}\) source, every positive compiler, a centered covariance evaluation, or an exact shifted/resummed current. It does **not** prove a uniform centered fourth-order theorem or an original-VALUE consumer. No ordinary Markov path-grid algorithm is executed or revived.

## Exact reviewed text and independence

The primary source read and audited directly is:

- `../RADIAL-SHELL-REFUTES-UNSHIFTED-SINGLE-HISTORY-COVARIANCE.md`
- SHA256 `2807f7307646cc16c79049319b574259874b09e731a13f686696bdd60a05b035`

The result above applies to those exact bytes, not to a summary or to the parallel author's text. The canonical identities, remainder scales, Hessian contraction norm, and witness were checked before reading the parallel proof and its cross-review for comparison.

Comparison-only sources subsequently read:

- `../RADIAL-SHELL-REFUTES-THE-UNSHIFTED-MEAN-COVARIANCE-REDUCTION.md`, SHA256 `9f8037bd618c82e188acbbc807935ee2b606c011d405603e31c27649d0dd0060`.
- `../independent-mean-separator-audit/INDEPENDENT-RADIAL-PARALLEL-PROOF-AUDIT.md`, SHA256 `6d6d6885ef698a71758ef498424b4ac735b29ccf23b49cabf18ab2f3df622e3e`.

The positive quadrature existence premise was also checked against Section 3 of `../../../../ENDPOINT-STEIN-PACKET-AND-POLYLOG-QUADRATURE.md`. The full-\(m_3\) scope was checked against `../P3-MEAN-RAW-SOURCE-AND-SECOND-DECOUPLING-GATE.md`. All input hashes are recorded in `MANIFEST.json`.

## 1. Global smoothness, anchoring, and Hessian sandwich

Write \(R=\sqrt D\ge4\), \(A=R^{-1}\), \(c=d=1/4\), \(b=\psi'\), and

\[
f(x)=\psi(|x|-R)\frac{x}{|x|},\qquad g(x)=cAx+dAf(x).
\]

Because \(\psi=0\) on \((-\infty,-2]\), \(f\) is identically zero in the ball of radius \(R-2\ge2\). Thus the formal radial singularity at zero is absent, including every derivative. Flatness of the bump at its endpoints gives global \(C^\infty\) regularity. The displayed radial potential differentiates to \(g\), so the field is genuinely a gradient and \(g(0)=0\).

The normalization estimate is valid: \(\eta\ge e^{-4/3}\) on \([-1,1]\), whence

\[
\|b\|_\infty\le e^{-1}/(2e^{-4/3})=\tfrac12e^{1/3}<1.
\]

For \(n=x/r\), \(P=I-nn^\top\), differentiation gives

\[
Dg=cAI+dA[b(r-R)nn^\top+\psi(r-R)P/r].
\]

Both eigenvalues are nonnegative. The radial eigenvalue is at most \((c+d)A=A/2\). Wherever the tangential perturbation is nonzero, \(r\ge R-2\ge2\), so its eigenvalue is at most \((c+d/2)A=3A/8\). Inside the inactive ball it is just \(cA\). Thus \(0\preceq Dg\preceq AI\) globally, with ample margin. Any stricter small-\(A\) guard is met by restricting to sufficiently large \(D\).

All later smooth Taylor expansions are on this particular admissible family. The argument does not impose bounds on higher derivatives for arbitrary members of the class being refuted.

## 2. Exact source decompositions and finite clock feasibility

For the conditional OU process starting at \(x\), its linear integral has mean \(x/2\). Its conditional covariance is exactly \(I/4\), since

\[
\int_0^\infty\!\int_0^\infty e^{-t-s}
  (e^{-|t-s|}-e^{-t-s})\,dt\,ds
 =2\int_0^\infty t e^{-2t}dt-\tfrac14=\tfrac14.
\]

Consequently the canonical decompositions are exact:

\[
I_a=cA(x/2+\sigma_aG_a)+dAK_a,
\quad \sigma_H=1/2,\quad \sigma_Q=\beta_Q.
\]

For \(H\), the positive measure \(e^{-t}dt\) has total mass one, so the genuine nonlinear OU history has \(|K_H|\le1\). For \(Q\), positive weights of mass one give \(|K_Q|\le1\). The same original \(g\), including its nonlinear radial part, is used throughout. There is no requirement that \(G_a\) be independent of \(K_a\); only its conditional standard Gaussian marginal and independence of the endpoint Gaussian \(X\) are used.

The admitted one-variable Hermite operator certificate applies componentwise and sums without an extra dimension factor:

\[
\|\mathbb EK_H-\mathbb EK_Q\|_{L^2(\gamma_D)}
 \le\delta\|f\|_2\le\delta.
\]

The linear clock means coincide exactly, so

\[
\|\mathbb EI_H-\mathbb EI_Q\|_2\le dA\delta.
\]

This is a conditional-mean estimate for a bounded value field; it is not a strong path approximation, a mean oracle, or a claim that the two genealogies agree.

The rule assumptions are nonvacuous. On dyadic panels for \(t=1-\tau\), positive Gauss–Legendre weights plus a terminal midpoint preserve mass one and first moment one half exactly. The upstream Bernstein-ellipse argument gives the all-Hermite-degree bound

\[
8\,4^{-m}+2^{1-K}.
\]

Taking \(K=\lceil\log_2(4/\delta)\rceil\), \(m=\lceil\log_4(16/\delta)\rceil\) makes this at most \(\delta\), with \(Km+1=O(\log^2(1/\delta))\) interior positive nodes. This uses only analytical quadrature, not a discretized Markov path. The checker constructs explicit finite instances at \(\delta=1/D=A^2\).

The Hermite degree-two certificate gives

\[
\sum_jv_j\tau_j^2\le1/3+\delta,
\quad\beta_Q=\sum_jv_j\sqrt{1-\tau_j^2}
 \ge2/3-\delta\ge7/12,
\]

using \(\delta\le1/16<1/12\). Therefore

\[
\beta_Q^2-1/4\ge13/144,
\quad |k_Q|\ge\kappa:=13/18432>0,
\quad k_Q=\tfrac{dc^2}{2}(1/4-\beta_Q^2)<0.
\]

The loose \(7/12\) bound is correct; no convergence of \(\beta_Q\) to \(\pi/4\) is required.

## 3. Radial \(L^p\) geometry and the canonical remainder sizes

Let \(S_D=|X|-R\), \(\lambda=1-cA/2\), \(y_0=\lambda X\), \(\rho=|y_0|\), and \(w=-cA\sigma G\), with \(0\le\sigma\le1\). Uniform Gaussian norm concentration gives \(\|S_D\|_p\le C_p\). On the good event \(|X|\ge R/2\), \(|w|\le\rho/2\), one has \(\rho\asymp R\) and the standard expansions

\[
|y_0+w|-\rho=n\cdot w+\frac{|P w|^2}{2\rho}
                  +O(|w|^3/\rho^2),
\]
\[
n(y_0+w)=n+Pw/\rho+O(|w|^2/\rho^2).
\]

Here \(\|n\cdot w\|_p=O(A)\), \(\|w\|_p=O(1)\), and inverse powers of \(\rho\) on the good event contribute corresponding powers of \(A\). Thus both the radius increment and the direction increment have \(L^p\) size \(O(A)\).

The endpoint-small-radius complement has probability \(e^{-c_0D}\). On the endpoint-good event, \(|w|>\rho/2\) forces a Gaussian norm of order \(D\), rather than its typical order \(\sqrt D\), and has an even smaller tail. Global boundedness of \(f,Df,D^2f\), together with Gaussian polynomial moments for the geometric expressions, makes the discarded terms smaller than every needed power of \(A\). No expansion near the origin is used.

The canonical cubic radial remainder is only \(O_{L^p}(A^2)\), and that is entirely sufficient. This audit does not borrow the stronger \(O(A^3)\) refined radial remainder appearing in the parallel proof.

Using the exact Gaussian identities \(\mathbb E(n\cdot w)=0\) and \(\mathbb E|Pw|^2=c^2A^2\sigma^2(D-1)\), including exponentially small truncation corrections, gives

\[
\mathbb E_G(|y_0+w|-\rho)
  =\frac{c^2A^2\sigma^2(D-1)}{2\rho}+O_{L^2(X)}(A^2),
\]
\[
\|\mathbb E_G(|y_0+w|-\rho)^2\|_{L^2(X)}=O(A^2).
\]

The first directional increment has conditional expectation zero; its remainder is \(O_{L^2(X)}(A^2)\). The radial/directional cross term is also \(O(A^2)\), by joint \(L^4\) bounds. Scalar Taylor expansion of the fixed smooth \(\psi\) therefore yields

\[
\mathbb E_G f(y_0+w)
 =f(y_0)+\frac{c^2\sigma^2}{2}A
       b(S_D-c/2)n(X)+O_{L^2(X)}(A^2).
\]

The two necessary replacements are justified by

\[
\rho-R=S_D-c/2-(cA/2)S_D,
\quad A^2(D-1)/\rho=A+O_{L^2(X)}(A^2).
\]

The radial eigenvalue of \(Df\) is \(b(r-R)\), and its tangential eigenvalue is \(\psi(r-R)/r\). Their variation and the variation of the projections, using the same radius/direction estimates, give

\[
\|Df(y_0+w)-Df(y_0)\|_{L^p(X,G);\mathrm{op}}\le C_pA.
\]

This establishes the canonical equations (8)–(9) with constants uniform in \(\sigma,D,Q\). The Euclidean displacement itself has order-one norm; no false small-Euclidean-displacement premise enters.

## 4. The actual nonlinear history and correlation cancellation

Global boundedness of the bilinear \(D^2f\) norm and \(|K_a|\le1\) give the pathwise expansion

\[
f(y_{\sigma_a}-dAK_a)
 =f(y_{\sigma_a})-dA\,Df(y_{\sigma_a})K_a+O(A^2).
\]

Replacing \(Df(y_{\sigma_a})\) by \(Df(y_0)\) costs \(O_{L^2}(A^2)\) before multiplying by the outer \(dA\). Indeed,

\[
\|\mathbb E[(Df(y_{\sigma_a})-Df(y_0))K_a\mid X]\|_2
 \le\|Df(y_{\sigma_a})-Df(y_0)\|_{L^2(X,G);\mathrm{op}}
 \le CA.
\]

This bound remains valid for arbitrary genuine correlation between \(K_a\) and \(G_a\), and even when \(K_H\) uses Brownian randomness beyond \(G_H\). It makes no decoupling assumption.

After expectation the retained term is exactly \(-d^2A^2Df(y_0)\mathbb EK_a\). Its difference between \(H\) and \(Q\) is \(O(A^2\delta)\). The outer linear force contributes \(-cA(\mathbb EI_H-\mathbb EI_Q)\), likewise \(O(A^2\delta)\). The common \(f(y_0)\) term cancels exactly. Hence

\[
j_H-j_Q=k_QA^2b(S_D-c/2)n(X)
                   +O_{L^2}(A^3+A^2\delta).
\]

In particular, the conditional-mean quadrature floor is \(A^2\delta\le A^4\), below both the \(A^3\) remainder and the \(A^2\) obstruction. An unpriced dimension-sized linear mean error is absent because the first clock moment is exact.

## 5. Actual covariance: conditional Gaussian Bessel, not independence

At fixed \(x\), put \(M_a=\operatorname{Cov}(G_a,K_a\mid x)\) and \(V_a=\operatorname{Cov}(K_a\mid x)\). The coordinate functions \((G_a)_i\) form an orthonormal set in conditional \(L^2\), so applying Bessel to each component of \(K_a-\mathbb EK_a\) gives

\[
\|M_a\|_{\mathrm{HS}}^2
 =\sum_j\sum_i|\mathbb E[(G_a)_i(K_{a,j}-\mathbb EK_{a,j})]|^2
 \le\mathbb E|K_a-\mathbb EK_a|^2\le1.
\]

This is valid whether or not \(K_a\) is a function of \(G_a\) alone. Also \(V_a\succeq0\) and \(\|V_a\|_{\mathrm{HS}}\le\operatorname{tr}V_a\le1\). The exact covariance expansion is

\[
C_a=c^2A^2\sigma_a^2I+
 A^2[cd\sigma_a(M_a+M_a^\top)+d^2V_a].
\]

Consequently its remainder has HS norm at most \((2cd+d^2)A^2=3A^2/16\), uniformly in \(x,D,Q\). The source's \(O(A^2)\) claim is therefore valid pointwise in the endpoint. These are covariances of the actual nonlinear sources, not covariances substituted from a different quadratic example.

## 6. Exact Hessian contraction and the dimension-safe norm

For any matrix \(B\), including nonsymmetric \(B\), direct differentiation gives

\[
D^2f(x):B=b'(r-R)n(n^\top Bn)
 +a_r\{n\operatorname{tr}(PB)+P(B+B^\top)n\},
\quad a_r=b(r-R)/r-\psi(r-R)/r^2.
\]

This was independently derived and checked symbolically against all entries of a general \(3\times3\) matrix. Rotational covariance gives the invariant expression in all dimensions.

A sharper norm check is available. Rotate \(n=e_1\). The radial output uses \(B_{11},B_{22},\ldots,B_{DD}\), with coefficients \(b',a_r,\ldots,a_r\). Tangential output \(i>1\) uses the disjoint pair \(B_{i1},B_{1i}\), each with coefficient \(a_r\). These input subspaces are orthogonal in HS norm. For \(D\ge2\), the exact operator norm is therefore

\[
\|D^2f(x)\|_{\mathrm{HS}\to\mathbb R^D}
 =\max\!\left\{\sqrt{b'(r-R)^2+(D-1)a_r^2},\sqrt2|a_r|\right\}.
\]

In particular, the canonical upper bound \(|b'|+(\sqrt D+2)|a_r|\) is correct. On the transition shell \(r\ge R-2\ge R/2\), \(\sqrt D/r\le2\) and \(\sqrt D/r^2\le4/R\). Beyond the shell \(b=b'=0\) and \(a_r=-1/r^2\), which is smaller. Inside the inactive ball the tensor is zero. Thus the norm is globally \(O(1)\), and the corresponding norm of \(D^2g=dAD^2f\) is globally \(O(A)\).

This verifies two uses at once: the bilinear Taylor bound in Section 4, and

\[
\|(C_H-c^2A^2\sigma_H^2I-C_Q+c^2A^2\sigma_Q^2I):D^2g\|_2
 =O(A^3).
\]

There is no missing factor \(\sqrt D\) in this contraction.

Setting \(B=I\) gives

\[
\Delta g=dA[b'(S_D)+(D-1)(b(S_D)/r-\psi(S_D)/r^2)]n.
\]

The \(b'\) and \(\psi/r^2\) contributions are \(O(A)\). On the support of \(b\), \(r=R+O(1)\), so \(A(D-1)/r=1+O(A)\). Hence

\[
\Delta g=d\,b(S_D)n+O_{L^2}(A),
\]

and in fact the error is uniformly bounded pointwise by \(CA\). Contracting the isotropic covariance difference and adding the actual covariance remainder yields

\[
\tfrac12(C_H-C_Q):D^2g(x)=k_QA^2b(S_D)n(X)+O_{L^2}(A^3).
\]

## 7. Resolvent witness and the fixed-log separation

The pre-resolvent discrepancy is

\[
k_QA^2h(S_D)n(X)+O_{L^2}(A^3+A^2\delta),
\quad h(s)=b(s-c/2)-b(s).
\]

Choose \(T_D(X)=AX=X/R\). Its vector \(L^2\) norm is exactly one. The componentwise Gaussian Mehler resolvent is self-adjoint, contractive, and acts by \(1/2\) on first chaos. Therefore

\[
\langle R_1[h(S_D)n],T_D\rangle
 =\tfrac12\mathbb E[h(S_D)|X|/R].
\]

The ordinary chi-radius central limit theorem gives \(S_D\Rightarrow S\sim N(0,1/2)\). Because \(h\) is bounded and continuous and \(|X|/R\to1\) in \(L^2\), the right side tends to \(\frac12\mathbb Eh(S)\). No joint shell/history limit or conditional asymptotic independence is needed.

The limit is strictly negative by an explicit one-dimensional argument. Every nontrivial superlevel set of the symmetric strictly unimodal bump \(b\) is an interval \([-a,a]\), \(a>0\). If \(\varphi(s)=\pi^{-1/2}e^{-s^2}\), its translated mass satisfies, for \(t>0\),

\[
\frac{d}{dt}\int_{-a}^a\varphi(u+t)du
 =\varphi(a+t)-\varphi(a-t)<0.
\]

The inequality follows from \(|a+t|>|a-t|\), evenness, and strict radial decrease of \(\varphi\). Integrating positive layer-cake levels proves

\[
\mathbb E b(S-c/2)<\mathbb E b(S).
\]

Put \(m=-\mathbb Eh(S)>0\). For sufficiently large \(D\), the absolute leading witness is at least \(m/4\). Since \(|k_Q|\ge\kappa=13/18432\) and the remainder is uniformly \(O(A^3+A^2\delta)\), taking \(D\) still larger makes that remainder at most \(\kappa mA^2/8\). Thus a possible eventual constant is

\[
c_* =\kappa m/8>0.
\]

This yields the claimed norm lower bound by testing against the unit vector \(T_D\). The threshold is not quantified by the canonical text or by this audit, and none is needed for the uniform asymptotic refutation.

At \(A=D^{-1/2}\), \(A^4\sqrt D=A^3=D^{-3/2}\), while the obstruction has size \(A^2=D^{-1}\). Their ratio is \(\sqrt D\). For the concrete positive rules at \(\delta=A^2\), all their accuracy and node-count logarithms are \(O(\log D)\). Therefore the ratio exceeds every fixed polynomial in the public logarithms along this admissible sequence. No claim about an arbitrarily chosen superpolynomial accuracy parameter is needed.

## 8. Centering survives; the full target is not decided here

For the same family,

\[
\mu_H(x)=cAx/2+dA\mathbb EK_H(x),\qquad |dA\mathbb EK_H|\le dA.
\]

Thus \(z=x-\mu_H(x)\) has radius coordinate \(S_D-c/2+O_{L^p}(A)\). Its direction agrees with \(n(X)\) to the needed order; on the usual endpoint-good event the bounded nonlinear mean perturbation changes direction by \(O(A^2)\). The global Hessian contraction bound remains valid at \(z\), and the uniform Laplacian expansion from Section 6 gives

\[
\tfrac12(C_H-C_Q):D^2g(x-\mu_H(x))
 =k_QA^2b(S_D-c/2)n(X)+O_{L^2}(A^3).
\]

This matches the actual single-history force bias to \(O(A^3+A^2\delta)\) for this family. It demonstrates precisely why evaluating at \(x\) is wrong and why this example does not contradict coherent centering or exact resummation.

The reference P3 note has additional nested sources and mean/innovation currents. A failure of the specialized single-\(I\) unshifted formula cannot, without another comparison, be relabeled a proof that its complete full-\(m_3\) equation is false. The canonical source explicitly limits that inference, and this audit preserves it.

## 9. Reproducible checker and its limits

Run:

`python check_independent_radial_mean_audit.py`

It checks the exact canonical hash, symbolic OU mean/covariance constants, finite positive quadrature examples, the general-matrix radial derivative identity, the HS-to-vector norm formula, the scalar Gaussian convolution witness, and the algebraic \(A^2/A^3\) separation. Results are in `independent_radial_mean_checks.json`.

Independent numerical diagnostics give approximately:

- \(\|\psi'\|_\infty=0.4142844199345533\).
- \(\mathbb Eh(S)=-0.00182950693274534\).
- Limiting linear-resolvent witness \(=-0.00091475346637267\).
- Uniform absolute leading witness coefficient floor \(\kappa m/2\approx6.45171173114406\times10^{-7}\).
- One possible eventual lower-bound constant \(\kappa m/8\approx1.61292793278602\times10^{-7}\).

The strict sign is proved analytically above; these floating-point integrals are not a formal interval certificate. Finite-\(D\) scalar values concern only the leading-field witness. The checker does not simulate the true OU history, certify the full finite-\(D\) remainder threshold, prove the analytical moment estimates by computation, or validate a full-\(m_3\) consumer.

`MANIFEST.json` records exact input and artifact hashes. `SHA256SUMS` seals the audit, checker, results, and manifest. No source file outside this independent audit directory was modified.
