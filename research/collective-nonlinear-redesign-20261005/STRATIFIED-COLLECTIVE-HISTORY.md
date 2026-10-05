# A linear-cost, first-order strong approximation of the collective OU history

## Result and scope

Let \((X_t)_{t\geq0}\) be stationary standard OU in \(\mathbb R^D\), with
\(\mathbb E[X_sX_t^\top]=e^{-|t-s|}I_D\). Let \(g:\mathbb R^D\to\mathbb R^D\) be globally \(A\)-Lipschitz. The intended class \(g=\nabla U\), \(0\preceq Dg\preceq AI\), is a subclass. Anchoring \(g(0)=0\) suffices for all integrability statements below; the approximation bound itself only uses increments.

Define
\[
 H_t=\int_0^\infty e^{-h}g(X_{t+h})\,dh,\qquad
 B=\int_0^\infty e^{-t}g(X_t-H_t)\,dt,
\]
\[
 T=g(X_0-B),\qquad \psi(x)=\mathbb E[T\mid X_0=x].
\]

For every integer \(M\geq1\), the construction below uses exactly \(2M+1\) evaluations of \(g\), \(M\) independent scalar clock labels, and \(MD\) independent standard Gaussian coordinates. Conditional on all clock labels, it is a finite Gaussian-input graph. It satisfies
\[
 \left(\mathbb E_{U,X}\|T_U(X_0;Z)-T\|^2\right)^{1/2}
 \leq \delta_M:=24A^2(1+A)\frac{\sqrt D}{M}.
 \tag{1}
\]
In particular, if \(\mu_U(x)=\mathbb E_ZT_U(x;Z)\), then
\[
 \left(\mathbb E_U\|\mu_U-\psi\|_{L^2(\gamma_D)}^2\right)^{1/2}
 \leq\delta_M.
 \tag{2}
\]
The expectation defining \(\mu_U\) is a mathematical description of the source's own mean, not an executed expectation oracle. The executed object is the finite source graph below. The Gaussian-own-mean service, if used, must operate on that graph.

Every fixed-clock source has Gaussian-bank Lipschitz constant at most
\[
 \operatorname{Lip}_Z(T_U(x;\cdot))\leq A^2(1+A),
 \tag{3}
\]
uniformly in \(x,U,M,D\). No derivatives with respect to clocks, and no evaluations of \(Dg\) or \(D^2g\), occur in this construction or its error proof.

## Executable finite graph

Use the nonuniform deterministic bin edges
\[
 b_i=2\log\frac{M}{M-i},\quad 0\leq i<M,\qquad b_M=+\infty,
\]
and weights
\[
 w_i=e^{-b_i}-e^{-b_{i+1}}
 =\frac{2(M-i)-1}{M^2}.
\]
For each \(i\), independently draw \(U_i\in[b_i,b_{i+1})\) with density \(e^{-u}/w_i\). Equivalently, draw
\[
 q_i=e^{-U_i}\sim\operatorname{Unif}\left(\frac{(M-i-1)^2}{M^2},\frac{(M-i)^2}{M^2}\right).
\]
The list is automatically ordered: \(0<U_0<\cdots<U_{M-1}<\infty\), almost surely. Freeze this entire list as labels across the completed Gaussian bank and across all compared/replayed evaluations. These labels are not Gaussian coordinates and are never differentiated.

Given the input \(x\), set \(G_{-1}=x\), \(q_{-1}=1\). For independent \(Z_i\sim N(0,I_D)\), set
\[
 G_i=\frac{q_i}{q_{i-1}}G_{i-1}
       +\sqrt{1-\left(\frac{q_i}{q_{i-1}}\right)^2}\,Z_i.
 \tag{4}
\]
This finite chain is exactly the true conditional joint law of \((X_{U_i})_i\) given \(X_0=x\), not a time-stepping approximation.

Evaluate \(Y_i=g(G_i)\). Compute suffix vectors
\[
 S_i=\sum_{j>i}w_jY_j
\]
in one backward pass, and put
\[
 \widehat H_i
 =\left(1-\frac{e^{-b_{i+1}}}{q_i}\right)Y_i+\frac{S_i}{q_i},
 \qquad
 \widehat B=\sum_{i=0}^{M-1}w_i\,g(G_i-\widehat H_i),
 \tag{5}
\]
\[
 T_U(x;Z)=g(x-\widehat B).
 \tag{6}
\]
For the last bin, \(e^{-b_M}=0\), \(S_{M-1}=0\), so \(\widehat H_{M-1}=Y_{M-1}\); no division of a nonzero suffix by an arbitrarily small last clock is needed.

All weights inside each \(\widehat H_i\) are nonnegative and have total absolute weight exactly one:
\[
 1-\frac{e^{-b_{i+1}}}{q_i}+\sum_{j>i}\frac{w_j}{q_i}=1.
\]
Likewise \(\sum_i|w_i|=1\). Thus each inner history is a collective convex combination of the sampled nonlinear responses; it is not a Taylor expansion. Vector additions, scalar multiplications, suffix evaluation, and the exact OU chain cost \(O(MD)\). The graph uses \(M\) calls for \(Y_i\), \(M\) calls for the nonlinear outer responses, and one terminal call.

## Proof of the strong bound

All norms in this section are joint \(L^2\) norms over the stationary OU path and any displayed independent clock variables. An infinite OU path is only a proof coupling; the implementation uses (4).

### 1. Stationary increments and stratified variance

OU increments give
\[
 \|g(X_t)-g(X_s)\|_2
 \leq A\sqrt{2D(1-e^{-|t-s|})}.
 \tag{7}
\]
By Minkowski and time stationarity,
\[
 \|H_t-H_s\|_2
 \leq A\sqrt{2D(1-e^{-|t-s|})}.
\]
Writing \(V_t=g(X_t-H_t)\), it follows that
\[
 \|V_t-V_s\|_2
 \leq A(1+A)\sqrt{2D(1-e^{-|t-s|})}.
 \tag{8}
\]

For any path functional \(V_t\), conditional on its entire realized path, the independent stratified samples satisfy
\[
 \mathbb E_U\left\|\sum_iw_iV_{U_i}-\int_0^\infty e^{-t}V_t\,dt\right\|^2
 =\sum_i w_i^2\operatorname{Var}_{U_i}(V_{U_i}).
 \tag{9}
\]
For an independent duplicate \(U_i'\) in the same bin,
\(\operatorname{Var}_{U_i}(V_{U_i})=\tfrac12\mathbb E_{U_i,U_i'}\|V_{U_i}-V_{U_i'}\|^2\).
Combining with (8) yields
\[
 \left\|\sum_iw_iV_{U_i}-B\right\|_2^2
 \leq A^2(1+A)^2D\sum_iw_i^2\eta_i,
 \quad \eta_i:=\min\{b_{i+1}-b_i,1\},
 \tag{10}
\]
where the infinite last width has \(\eta_{M-1}=1\).

For convenient estimates, index bins in reverse by \(k=M-i\). Then
\[
 w^{(k)}=(2k-1)/M^2,\quad
 \Delta_k=2\log(k/(k-1))\quad(k\geq2),\quad \eta_1=1.
\]
For \(k\geq2\),
\[
 (w^{(k)})^2\eta_k
 \leq \frac{2(2k-1)^2}{M^4(k-1)}
 =\frac{8k+2/(k-1)}{M^4}
 \leq\frac{9k}{M^4}.
\]
The same last bound holds for \(k=1\). Consequently, for \(1\leq L\leq M\),
\[
 \sum_{k=1}^L(w^{(k)})^2\eta_k
 \leq\frac{9L(L+1)}{2M^4}\leq\frac{9L^2}{M^4}.
 \tag{11}
\]
In particular,
\[
 \left\|B-\sum_iw_iV_{U_i}\right\|_2
 \leq3A(1+A)\frac{\sqrt D}{M}.
 \tag{12}
\]

### 2. Inner-history error

Split \(H_{U_i}-\widehat H_i=L_i+Q_i\), where
\[
 L_i=\int_{U_i}^{b_{i+1}}e^{-(t-U_i)}[g(X_t)-g(X_{U_i})]\,dt
\]
is the same-bin remainder, and
\[
 Q_i=e^{U_i}\sum_{j>i}\left[\int_{b_j}^{b_{j+1}}e^{-t}g(X_t)\,dt-w_jg(X_{U_j})\right]
\]
is the later-bin quadrature error.

By (7), uniformly over \(U_i\) in its bin,
\[
 \|L_i\|_2
 \leq A\sqrt{2D}\int_0^{\Delta_i}e^{-h}\sqrt{1-e^{-h}}\,dh
 =\frac{2\sqrt2}{3}A\sqrt D\,(1-e^{-\Delta_i})^{3/2}.
 \tag{13}
\]
The equality follows from the substitution \(v=1-e^{-h}\), and also holds for the infinite last bin. In reverse index \(k\), \(1-e^{-\Delta_k}=(2k-1)/k^2\leq2/k\), so
\[
 \sum_iw_i\|L_i\|_2
 \leq\frac{8A\sqrt D}{3M^2}\sum_{k=1}^M\frac{2k-1}{k^{3/2}}
 \leq\frac{32}{3}A\frac{\sqrt D}{M^{3/2}}.
 \tag{14}
\]

For a bin in reverse position \(k\geq2\), its later bins have indices \(1,\ldots,k-1\). Conditional on the full path and \(U_i\), their sampling errors are independent and centered. By (7), (9), and (11), and by
\(e^{U_i}\leq e^{b_{i+1}}=M^2/(k-1)^2\),
\[
 \|Q_i\|_2
 \leq\frac{M^2}{(k-1)^2}\,A\sqrt D\,
       \frac{3(k-1)}{M^2}
 =\frac{3A\sqrt D}{k-1}.
 \tag{15}
\]
There is no later-bin error for \(k=1\). Thus, with \(H_0=0\),
\[
 \sum_iw_i\|Q_i\|_2
 \leq\frac{3A\sqrt D}{M^2}\,[2(M-1)+H_{M-1}]
 \leq9A\frac{\sqrt D}{M}.
 \tag{16}
\]

### 3. Propagation through the retained nonlinear responses

By Lipschitzness of \(g\), (12), (14), and (16),
\[
 \|B-\widehat B\|_2
 \leq3A(1+A)\frac{\sqrt D}{M}
 +A^2\sqrt D\left[\frac9M+\frac{32}{3M^{3/2}}\right].
 \tag{17}
\]
One more application of Lipschitzness gives
\[
 \|T-T_U\|_2
 \leq\frac{A^2\sqrt D}{M}
 \left[3+12A+\frac{32A}{3\sqrt M}\right]
 \leq24A^2(1+A)\frac{\sqrt D}{M},
\]
which proves (1). Jensen conditional on \((X_0,U)\) proves (2).

## Frozen-clock Gaussian-bank bound

For fixed clocks, expanding (4) gives
\[
 G_i=q_ix+\sum_{j\leq i}c_{ij}Z_j,\qquad
 \sum_{j\leq i}c_{ij}^2=1-q_i^2\leq1.
\]
Thus each map \(Z\mapsto G_i\) has Euclidean operator norm at most one. Consequently:

- \(Z\mapsto Y_i\) is \(A\)-Lipschitz;
- the absolute coefficient sum in \(\widehat H_i\) is one, so it is \(A\)-Lipschitz;
- \(Z\mapsto g(G_i-\widehat H_i)\) is \(A(1+A)\)-Lipschitz;
- the weights of \(\widehat B\) have absolute sum one, so it is \(A(1+A)\)-Lipschitz;
- the terminal source is \(A^2(1+A)\)-Lipschitz in the whole Gaussian bank.

These are bounds in the ordinary Euclidean norm of the full \(MD\)-coordinate bank. They do not hide a sum of \(M\) coordinatewise bounds, nor do they differentiate through the inverse-exponential clocks.

## Positive mixtures and what the mean estimate does not claim

Suppose a conditional service, for each frozen \(U\), provides a law \(\nu_{U,x}\) with
\[
 \left(\mathbb E_U\int W_2^2(\nu_{U,x},N(\mu_U(x),I_D))\,\gamma_D(dx)\right)^{1/2}\leq\varepsilon.
\]
Then the positive mixture \(\overline\nu_x=\mathbb E_U\nu_{U,x}\) satisfies
\[
 \left(\int W_2^2(\overline\nu_x,N(\psi(x),I_D))\,\gamma_D(dx)\right)^{1/2}
 \leq\varepsilon+\delta_M.
\]
This follows by coupling the component service to a Gaussian with its own mean, and then using the same standard Gaussian to shift \(\mu_U(x)\) to \(\psi(x)\). The clocks select a component of a positive mixture; no signed law and no extra expectation oracle are introduced.

Equation (2) is an average-over-clocks estimate. It is not a uniform guarantee for every possible frozen list. If one globally freezes a single random list, it provides the stated root-mean-square guarantee and the Markov bound \(\|\mu_U-\psi\|_{L^2}\leq\delta_M/\sqrt\rho\) with probability at least \(1-\rho\). A per-sample positive mixture uses (2) directly.

For the OU resolvent \(R_1=\int_0^\infty e^{-s}P_s\,ds\), its \(L^2(\gamma_D)\)-contraction implies
\[
 \left(\mathbb E_U\|R_1\mu_U-R_1\psi\|_{L^2}^2\right)^{1/2}\leq\delta_M.
\]
This statement presumes that the resolvent is actually integrated in the conditional mean or handled by a separately valid compiler. Freezing one additional outer clock \(S\) creates the mean \(P_S\mu_U\), not \(R_1\mu_U\). Its across-\(S\) spread is not bounded by \(\delta_M\); replacing the resolvent with that mixture needs a separate argument.

## Optional exact Gaussian linear-history control variate

This optional refinement is not needed for (1)–(3). It adds only a fixed-factor number of Gaussian coordinates.

Define exact Gaussian future integrals
\[
 J_t=\int_0^\infty e^{-h}X_{t+h}\,dh,\qquad
 K_t=\int_0^\infty h e^{-h}X_{t+h}\,dh.
\]
The stationary Gaussian process \((X,J,K)\), coordinatewise, has covariance
\[
 \Sigma=\begin{pmatrix}
 1&1/2&1/4\\
 1/2&1/2&3/8\\
 1/4&3/8&3/8
 \end{pmatrix}
\]
and forward Markov drift
\[
 C=\begin{pmatrix}-5&12&-8\\-1&1&0\\0&-1&1\end{pmatrix},
 \qquad Q=\operatorname{diag}(2,0,0).
\]
One direct derivation is to reverse time: then the process is the causal linear filter
\(dX=-X\,dt+\sqrt2\,dW\), \(dJ=(X-J)dt\), \(dK=(J-K)dt\). Its lower-triangular drift is \(B\); stationary Gaussian time reversal has drift \(\Sigma B^\top\Sigma^{-1}=C\). The displayed \(\Sigma\) solves its Lyapunov equation, as can also be checked by the defining integrals.

Given \(X_0=x\), initialize with independent \(Z_1,Z_2\sim N(0,I_D)\):
\[
 J_0=x/2+Z_1/2,\qquad K_0=x/4+Z_1/2+Z_2/4.
\]
Over any fixed time gap \(\Delta\), use the exact Gaussian transition with mean multiplier \(E=e^{C\Delta}\) and covariance \(\Sigma-E\Sigma E^\top\), tensorized with \(I_D\). This samples the true joint linear histories at all queried clocks in \(O(MD)\) work; it is a Gaussian sampler, not a nonlinear expectation oracle.

For any fixed scalar \(c\), write \(g=c\,\mathrm{id}+r\) and
\(R_t=\int_0^\infty e^{-h}r(X_{t+h})\,dh\). Then the following identity is exact:
\[
 B=cJ_0-c^2K_0+
 \int_0^\infty e^{-a}\left[r(X_a-cJ_a-R_a)-ca\,r(X_a)\right]da.
\]
It follows by \(H_a=cJ_a+R_a\) and stochastic Fubini;
\(\int e^{-a}J_a\,da=K_0\) and
\(\int e^{-a}R_a\,da=\int a e^{-a}r(X_a)\,da\).
Thus the linear response can be kept exact and only the residual terms stratified. For \(c=A/2\), the convex-gradient assumption gives \(\operatorname{Lip}(r)\leq A/2\). For \(g=c\,\mathrm{id}\), the residual vanishes and the canonical path functional is exact. No faster general-class weak rate is claimed from this refinement alone.
