# Independent audit: two-replica conditional-Gram law

Date: 2026-10-05. Scope: the finite positive pushforward construction, its marginal transport bound, cap bias, and source dependencies. This does not certify a separate full-history coupling or an unprovided response-factorization argument.

## Verdict

The construction is valid. The proposed noncommuting Gaussian overlap identity is correct, and the stated dimension-dependent transport estimate follows. A rank-two observation gives the stronger bound

\[
 W_2\bigl(\mathcal L(Z),N(0,vI+C_R)\bigr)
 \le 5\sqrt v\,h^2,\qquad h=R^2/v\le 1/8,
\]

with no separate condition on \(nh^2\) and no explicit factor of \(n\). Dimension still enters through the cap radius \(R\).

For the original uncapped Gram target, the complete bound is

\[
\boxed{\quad
 W_2\bigl(\mathcal L(Z),N(0,vI+C)\bigr)
 \le {5\ell^4\over v^{3/2}}(\sqrt n+\sqrt{2L})^4
 + {\ell^2\over\sqrt v}\left(1+\sqrt{n/(2L)}\right)e^{-L},
 \quad L>0,\quad R^2/v\le1/8.
\quad}
\]

Thus optimizing the cap for the uncapped target gives a quartic rate with a logarithmic-squared loss. The displayed argument alone does not give a literal uniform \(O(\ell^4/v^{3/2})\) uncapped-target rate. For fixed \(L\), that literal rate holds for the capped target.

## 1. Roots, conditioning, and positivity

Write \(m(O)=\mathbb E_E\Psi(O,E)\) and \(C=\mathbb E_O[m(O)m(O)^T]\). Draw \(E_1,E_2\) independently conditional on the same \(O\), with the original conditional nuisance law. In the intended root model they are fresh independent copies of the nuisance roots, independent of the shared outer roots.

The supplied intended shared readset consists of the outer midpoint pair \((\xi,\xi')\), inner midpoint pair \((\zeta,\zeta')\), local endpoint data, and shared discrete geometry. The Gaussian centering subset \(U=(\zeta,\zeta')\) stays shared. Each nuisance copy has its own Gaussian bridge residuals, smoothing roots, and local integration times. A final standard Gaussian \(G\in\mathbb R^n\) is fresh and independent of all those roots.

Any conditioning used to prove Gaussian centering or concentration must freeze the remaining independent primitive roots, rather than inadvertently conditioning on derived data that constrain \(U\). Antisymmetry under \(\zeta\leftrightarrow\zeta'\), as specified by the parent, supplies conditional centering. Radial capping preserves this antisymmetry. If instead \(U\) were integrated out as part of \(E\), the stated conditional centering would force \(m=0\) and make the original Gram target trivial.

Let \(a=\operatorname{cap}_R\Psi(O,E_1)\), \(b=\operatorname{cap}_R\Psi(O,E_2)\), with \(\operatorname{cap}_R x=x\min(1,R/|x|)\), defined as zero at zero. Set

\[
 H=(ab^T+ba^T)/2,\qquad
 Z=\sqrt v\,(I+H/(2v))G.
\]

This is one genuine positive probability law, using two nuisance evaluations and one Gaussian output root. It is a finite-source pushforward, not necessarily a finite-support measure. Although \(H\) need not be positive semidefinite, the conditional output covariance is always positive semidefinite. For \(h<2\), the output scale matrix is itself positive definite.

Writing \(m_R(O)=\mathbb E_E\operatorname{cap}_R\Psi(O,E)\), conditional independence gives

\[
\mathbb E[H\mid O]=m_R(O)m_R(O)^T,\qquad
C_R=\mathbb EH\succeq0,
\]
\[
\mathbb EZ=0,\qquad
\operatorname{Cov}(Z)=vI+C_R+\mathbb E H^2/(4v).
\]

The last fourth-order covariance term is real and cannot be silently omitted.

## 2. Rank-two control

Direct calculation gives

\[
 \|H\|_F^2
 ={ |a|^2|b|^2+(a\cdot b)^2\over2}\le R^4,
 \qquad \|H\|_{\rm op}\le R^2.
\]

Consequently \(\|C_R\|_F\le R^2\) and \(\operatorname{tr}\mathbb EH^2\le R^4\). Define

\[
 \Sigma=vI+C_R,\qquad
 K=\Sigma^{-1/2}\left(H-C_R+H^2/(4v)\right)\Sigma^{-1/2}.
\]

Then the normalized law \(\Sigma^{-1/2}Z\), conditional on the roots, is Gaussian with covariance \(I+K\), and

\[
 \|K\|_F\le b_h:=2h+h^2/4,\qquad
 \|K\|_{\rm op}\le b_h,
\]
\[
 \mathbb EK=\Sigma^{-1/2}\mathbb EH^2\Sigma^{-1/2}/(4v)\succeq0,
 \qquad \|\mathbb EK\|_F\le\operatorname{tr}\mathbb EK\le h^2/4.
\]

These are dimension-free Frobenius bounds. An operator-norm-only argument discards this useful low-rank information and unnecessarily introduces powers of \(n\).

## 3. Noncommuting overlap and transport

Let \(r_K=dN(0,I+K)/d\gamma\), with \(\gamma=N(0,I)\), and let \(J\) be an independent copy of \(K\). For symmetric \(K,J\) with operator norms less than one, the Gaussian integral is finite. Indeed,

\[
 B=(I+K)^{-1}+(I+J)^{-1}-I\succ0,
 \qquad (I+K)B(I+J)=I-KJ.
\]

Taking determinants in this identity yields, without any commutation assumption,

\[
 \int r_Kr_J\,d\gamma=\det(I-KJ)^{-1/2}.
\]

The determinant is positive; the square root is its positive scalar square root. The product \(KJ\) need not be symmetric. Since \(\|KJ\|_{\rm op}<1\), the real trace-log series is nevertheless valid:

\[
 A(K,J):=-\tfrac12\log\det(I-KJ)
 =\tfrac12\operatorname{tr}(KJ)
 +\tfrac12\sum_{q\ge2}{\operatorname{tr}((KJ)^q)\over q}.
\]

For \(b=b_h\le1/3\), the nuclear/Frobenius product inequality gives

\[
 |\operatorname{tr}((KJ)^q)|\le b^{2q},\quad
 |A|\le {b^2\over2(1-b^2)},\quad
 |A-\tfrac12\operatorname{tr}(KJ)|\le {b^4\over4(1-b^2)}.
\]

Using \(|e^x-1-x|\le x^2e^{|x|}/2\), one obtains

\[
 \left|e^A-1-\tfrac12\operatorname{tr}(KJ)\right|\le b^4/2.
\]

Tonelli's theorem and independence of \(K,J\) now imply, for the normalized mixture \(\nu\),

\[
 \chi^2(\nu\mid\gamma)
 =\mathbb E e^{A(K,J)}-1
 \le\tfrac12\|\mathbb EK\|_F^2+\tfrac12 b_h^4
 \le9h^4\quad(h\le1/8).
\]

The leading expected trace is not zero: it equals \(\operatorname{tr}((\mathbb EK)^2)\). It is of fourth order, which is exactly enough.

Gaussian Talagrand transport, \(W_2^2\le2\operatorname{KL}\le2\chi^2\), and \(\|\Sigma\|_{\rm op}\le v(1+h)\) give

\[
 W_2\bigl(\mathcal L(Z),N(0,\Sigma)\bigr)
 \le\sqrt{18v(1+h)}h^2\le 5\sqrt v\,h^2.
\]

This concerns only the marginal law. It does not claim a small error under the original rootwise coupling \(Z\leftrightarrow\sqrt\Sigma G\).

## 4. Cap bias

Conditional Gaussian Poincare and concentration, from the stated rootwise centering and \(\ell\)-Lipschitz assumption, imply

\[
 \mathbb E|\Psi|^2\le n\ell^2,\qquad
 \Pr\{|\Psi|>\ell\sqrt n+t\}\le e^{-t^2/(2\ell^2)}.
\]

For an especially useful cap estimate, retain the two-replica representation of both Gram matrices. If \(X=|\Psi(O,E_1)|\), \(Y=|\Psi(O,E_2)|\), then

\[
 XY-\min(X,R)\min(Y,R)
 \le\tfrac12\left((X^2-R^2)_++(Y^2-R^2)_+\right).
\]

This is equivalent to the scalar cap being 1-Lipschitz. Applying the triangle inequality in nuclear norm gives

\[
 \|C-C_R\|_*\le\mathbb E(|\Psi|^2-R^2)_+.
\]

For \(R=\ell(\sqrt n+\sqrt{2L})\), integration of the Gaussian norm tail and the elementary Gaussian Mills bound give

\[
 \|C-C_R\|_*\le2\ell^2\left(1+\sqrt{n/(2L)}\right)e^{-L},\qquad L>0.
\]

No independence between the two replica norms is required for this estimate. There is no asserted Loewner ordering between \(C\) and \(C_R\).

Finally, for positive definite \(A,B\succeq vI\), the matrix-square-root Sylvester identity gives

\[
 W_2(N(0,A),N(0,B))
 \le\|A^{1/2}-B^{1/2}\|_F
 \le {\|A-B\|_F\over2\sqrt v}.
\]

Applying this with \(A=vI+C\), \(B=vI+C_R\) proves the cap term in the boxed bound.

## 5. Rate interpretation and remaining scope

With \(t=\ell^2/v\), taking \(L\) comparable to \(\log(1/t)\) makes the cap bias fourth order and yields

\[
 W_2\lesssim_n {\ell^4\over v^{3/2}}\bigl(1+\log(v/\ell^2)\bigr)^2
\]

in the small-\(t\) regime, subject to \(t(\sqrt n+\sqrt{2L})^2\le1/8\). Replacing this by a log-free rate requires an additional argument or a different construction.

The proof uses no full-history join and does not establish a native bounded-first theorem or a general-\(D\) improvement. It establishes the finite positive conditional-Gram law, with its exact mean/covariance, marginal transport estimate, cap correction, and explicit shared-versus-replicated root requirements.
