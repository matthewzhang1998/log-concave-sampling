# Independent audit: same-carrier P3 contract and self-covariance counterexample

2026-10-04. **PASS for the exact analytical target contract, covariance replacement bound, and counterexample. The finite same-carrier fourth-order endpoint gate remains OPEN.**

Audited file: `../SAME-CARRIER-P3-MEAN-AND-FULL-COVARIANCE-GATE.md`.

Audited SHA256: `a63c38f14206ec278f7c12a56542968a7d1675af6f550610e4539e500a653eca`.

This review covers Sections 2–5 and the resulting bounded no-join conclusion. It does not certify later sibling proposals or supply the missing finite services. No source was modified.

## 1. The correct conditional mean

The analytical stationary-resolvent iteration has `P_(j+1)=X-integral g(P_j,past)` and contracts its stationary L2 error to the target by A at every substitution. The already pinned reference proves the contraction for every fixed j, so the P3 law error `A^4sqrt(D)` is justified as an analytical reference. It is not a query-count theorem for an executed continuous path.

For a fixed r, the Gaussian past before X_r is conditionally independent of the later endpoint Z given X_r. Rescaling its logarithmic time identifies its conditional law with the same original past-history law at endpoint x. Therefore

    chi2(x)=Law(x-F2,1 | X1=x),
    psi2(x)=E[g(x-F2,1)|X1=x],
    m3(Z)=int_0^1 P_r psi2(Z)dr.

The endpoint retained by chi2 is a Gaussian Markov endpoint. No posterior distribution occurs in this conditioning. A whole-output posterior endpoint coupling cannot be reused while retaining the integrated old Gaussian endpoint, so the previously audited posterior-force re-entry does not supply this target.

For g(x)=a x, direct conditional OU integration gives m2 coefficient `a/2-a^2/4` and m3 coefficient `a/2-a^2/4+a^3/8`. The covariance kernel for exponential linear path integrals is

    I(alpha,beta)=2/[(alpha+beta)(alpha+1)(beta+1)].

Its value and first/mixed derivatives at alpha=beta=1 give Var(H0)=1/4, Cov(H0,J0)=1/4 and Var(J0)=5/16. Hence the stated chi2 variance `a^2/4-a^3/2+5a^4/16` is exact. The fixed-heat posterior force is `a r z/[1+a(1-r^2)]`, already different at second order when r=1/2. Its conditional variance is order one at fixed heat, whereas chi2's is order a^2. The distinction is substantive, not notation.

## 2. The valid covariance replacement

On finite Gaussian path approximations, every conditional point row has norm at most one. Positive averaging and the actual chain rule give the stated private firsts `A(1+A+...+A^(j-1))`. Lipschitz chords give conditional energies `O(A^2(|z|+sqrt(D)))` for F2-H and `O(A^3(|z|+sqrt(D)))` for F3-F2. The uniform bounds pass to the analytical path limit.

For Gaussian-root U,V, dualizing their covariance against an HS-unit matrix M and applying Gaussian Riesz integration to the centered V gives

    |Cov(U,V):M|<=||DU* M||_(L2 HS) ||R(V-EV)||_(L2 HS)
                 <=Lip(U)||V-EV||2.

This proves the one-marked-energy covariance bound. Apply it to

    Cov(F3)-Cov(F2)=Cov(F3-F2,F3)+Cov(F2,F3-F2).

After transposing the first term when needed, the O(A) first belongs to F2 or F3 and the O(A^3) energy to their difference. Squaring and integrating over standard Z gives `O(A^4sqrt(D))` as claimed.

The exact F2 decomposition retains its self-covariance:

    Cov(F2)=Cov(H)+Cov(H,E2)+Cov(E2,H)+Cov(E2).

An O(A^2sqrt(D)) VALUE mark for E2 does not imply an O(A^2) first. The source's Hessian difference can still be O(A), so the one-marked bound for Cov(E2) is only O(A^3sqrt(D)). The counterexample below shows that this distinction cannot be removed uniformly by a better estimate of the same stated hypotheses.

## 3. Admissibility of the counterfamily

Let A=D^(-1/2), v=(1,...,1)/sqrt(D), c=1/4, epsilon=1/2, and

    g_D(x)=cA[b(x_1),...,b(x_D)]+cA v tanh(v·x),
    b(t)=t+epsilon log cosh(t).

This is a smooth anchored gradient. Its Hessian is the sum of a positive diagonal matrix and a positive rank-one matrix, bounded above by `c(2+epsilon)A I=(5/8)A I`. Thus the entire family satisfies the required convex Hessian sandwich with room to spare. No unbounded Hessian or hidden algebraic dimension-versus-A exception is being used.

The coordinate OU histories are independent. Put `K_r^i=int b(X_(r tau)^i)dtau` and `L_r=int tanh(U_(r tau))dtau`, with U_r=v·X_r. The vector inner shift has coordinates `I_r^i=cA(K_r^i+D^(-1/2)L_r)`. Its collective projection is

    v·I_r=c[D^(-1)sum_i K_r^i+D^(-1/2)L_r].

Since the K histories have uniformly finite second/fourth moments and stationary marginals, the ordinary coordinate L2 law of large numbers gives convergence to

    d=c epsilon E log cosh(N)>0

uniformly in r in L2, and therefore after integration in r. The common rank-one term is bounded and vanishes at rate D^(-1/2).

The exact projection formula for E2 in the source follows by substituting g_D. Taylor expanding only log cosh, whose second derivative is bounded by one, has summed L2 remainder

    D^(-1/2) sum_i O(||I_r^i||4^2)=O(D^(-1/2)).

The leading separable sum is `-epsilon c D^(-1)sum_i K_r^i tanh(X_r^i)` plus a vanishing rank-one term. It converges to a deterministic constant by the same coordinate law of large numbers. The tanh difference is Lipschitz in its shift, so replacing v·I_r by d costs another vanishing L2 error. Consequently the source's exact asymptotic form is correct:

    (v·E2)/A=c[C0+J_D]+o_L2(1),
    J_D=int_0^1 [tanh(U_r-d)-tanh(U_r)]dr.

The law of the one-dimensional OU path U, and therefore of J_D, is the same for every D. No coupling between different dimensions or conditional limit interchange is needed.

## 4. Conditional variance survives complete endpoint conditioning

Fix 0<a<1 and let V_D=U_a-aU_1. This is a Gaussian innovation independent of every coordinate of the complete endpoint Z_D, with variance 1-a^2. For `h(u)=tanh(u-d)-tanh(u)`, scalar Gaussian integration by parts gives

    E[J_D V_D]=E h'(N) int_0^1[min(r,a)/max(r,a)-ra]dr
              =E h'(N)(-a log a).

The covariance integral is exact. Moreover `E sech^2(N-d)<E sech^2(N)` for d>0: the even strictly decreasing function sech^2 is a positive layer-cake integral of centered interval indicators, and each shifted Gaussian interval mass is strictly smaller away from zero. Therefore E h'(N) is strictly negative, and this inner product is a fixed nonzero constant independent of D.

Conditional centering of v·E2 given Z_D does not change its inner product with V_D, because the latter is independent of Z_D and has mean zero. The unconditional L2 approximation above transfers that nonzero inner product. Cauchy–Schwarz therefore gives, for all sufficiently large D,

    E Var(v·E2 | Z_D)>=c_* A^2,

with c_*>0 independent of D. This argument directly avoids the potentially invalid step of passing an L2 limit through varying conditional distributions.

Since conditional covariance is PSD and v is unit,

    ||Cov(E2|Z_D)||_(L2;HS)>=E[v*Cov(E2|Z_D)v]>=c_* A^2.

But `A^4sqrt(D)=A^3`, while `A^3sqrt(D)=A^2`. Thus the proposed fourth-order self-covariance remainder fails by a factor of order 1/A, which dominates every fixed polynomial in the public dimension/heat logarithms. The existing one-marked allowance is sharp at this scaling. The lower-bound constant can be tiny; its fixed positive value is all the uniform asymptotic counterexample requires.

## 5. Consequence and diagnostics

The bounded target contract is correct: one still needs an actual same-standard-Z mean service for m3 and a full covariance service for Cov(F2|Z), including the E2 self-covariance, or an equivalent exact correction. The separately audited cubic-current packet does not supply either. The valid Cov(F3)->Cov(F2) replacement does not remove the terms inside Cov(F2).

The independent checker verifies the exact quadratic conditional moments, posterior target mismatch, bridge innovation covariance, Hessian coefficient and asymptotic scales. Gaussian quadrature only illustrates the nonzero constants; the lower bound is established analytically above, not inferred from a finite-dimensional simulation. Counts and hashes appear in the adjacent JSON and manifest.

**Conclusion:** the original pinned contract/counterexample passes. The dropped-self-covariance route fails under the stated uniform Hessian class, including fixed public-log losses. This is not an impossibility theorem for a positive fourth-order construction that actually realizes the missing covariance.
