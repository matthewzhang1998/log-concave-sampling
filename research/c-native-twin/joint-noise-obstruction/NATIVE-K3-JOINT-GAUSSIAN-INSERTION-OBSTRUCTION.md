# A source-qualified obstruction to joint Gaussian noise insertion

This NEW bounded obstruction uses the literal scalar K3 native twin c8b20bf4, with its complete retained terminal/twin and ancestor records. It rules out a high-order SINGLE-GAUSSIAN replacement of the joint base/shifted ancestor, and of that Gaussian insertion before the actual paired terminal observable. It does not rule out a weaker SAME-E coefficient/current estimate, a non-Gaussian retained kernel, or a Gaussian mixture carrying extra declared roots.

## 1. One original strongly convex primitive and actual weak scales

Fix m=1/2, d=1/10 and 0<epsilon<=1. Use the SAME original primitive at both nodes,

    g0(x)=gt(x)=g(x)=m x+d epsilon sin(x/epsilon).

It is the gradient of m x²/2-d epsilon² cos(x/epsilon), with

    lambda=m-d=.4 <=g'(x)<=m+d=.6=Lambda.

It is anchored at zero and has uniform unit first, regardless of epsilon. Take the exact original C,P rows and protected increment from c8b20bf4, scalar M=1 and 0<a<=1/16. In the independent coordinates(S,U,Z), define x=c0 S+s0 U and

    k_S(Z)=a[g(S)-g(S+epsilon Z)],
    E(S,U,Z)=r{g(S+a g(x))-g(S+a g(x+k_S(Z)))}.         (1)

This is the ACTUAL K3 remainder, not an independent surrogate. Write rho=a epsilon and theta=S/epsilon. Then exactly

    k_S(Z)/rho=-m Z-d[sin(theta+Z)-sin(theta)].         (2)

Its Z derivative is negative and lies between -Lambda rho and -lambda rho. Its full S derivative is bounded by2ad, and |k_S(Z)|<=Lambda rho|Z|. At the native weak scales r=a=A, epsilon=A^.9, rho=A^1.9.

## 2. Exact distance from every scalar Gaussian, uniformly in S

The map Z -> k_S(Z) is strictly decreasing. In one dimension the optimal transport to any Gaussian is therefore the matching anti-monotone affine quantile mu-sigma Z, sigma>=0. Consequently

    inf_(mu,sigma>=0) W2(Law(k_S|S),N(mu,sigma²))²
      =inf_(mu,sigma>=0) E|k_S-(mu-sigma Z)|².          (3)

The optimal mean and positive scale are

    mu/rho=d sin(theta)(1-exp(-1/2)),
    sigma/rho=m+d cos(theta)exp(-1/2)>0.

After subtracting the constant and first chaos, the exact residual is

 -d rho{sin(theta)[cos Z-exp(-1/2)]
          +cos(theta)[sin Z-exp(-1/2)Z]}.

The two brackets are even and odd and hence orthogonal. Their variances are

    v_even=(1+exp(-2))/2-exp(-1),
    v_odd =(1-exp(-2))/2-exp(-1)>0.

Since v_even-v_odd=exp(-2)>0, equation(3) is

    d² rho²[sin²(theta)v_even+cos²(theta)v_odd]
      >=d² rho² v_odd.                                (4)

Thus every phase is covered, including phases at which either sine or cosine vanishes. The UNIFORM lower constant is

    c_k=d sqrt(v_odd)=0.0253875790910144...,
    inf_Gaussian W2(Law(k_S|S),Gaussian)>=c_k rho.      (5)

No asymptotic phase equidistribution or exceptional-set argument is used.

## 3. Exact metrics on the jointly observed ancestor pair

Let U be standard independent of Z, conditional on S, and set the actual physical ancestor base x=c0 S+s0 U. Consider the LITERAL ancestor pair

    B_S=(x,x+k_S(Z)).

Give R² its ordinary Euclidean norm. For ANY bivariate Gaussian Gamma_S, including singular or nonstandard ones, its coordinate contrast is a scalar Gaussian. The map L(u,v)=v-u has norm sqrt(2), so

    inf_Gamma_S W2(Law(B_S|S),Gamma_S)
       >=c_k rho/sqrt(2)
       =0.0179517293331661... rho.                    (6)

If the base U (equivalently x at fixed S) is retained EXACTLY and the shifted coordinate is compared conditionally on it to any single Gaussian kernel, the stronger metric is

    inf_(mu(S,U),sigma(S,U)>=0)
      || W2(Law(x+k_S(Z)|S,U),
                    N(mu(S,U),sigma(S,U)²)) ||_(L2_U).

Translation by x and(5) show that this quantity is at least c_k rho. Allowing the Gaussian parameters to depend arbitrarily on the retained(S,U) does not remove the bound.

For the standardized pair (U,U+k_S(Z)/s0), the corresponding lower bounds are divided by s0. All norms above are therefore on the displayed physical pair, with no hidden change of metric.

This is the observation-level distinction: two jointly retained coordinates reveal their nonlinear contrast. A theorem for the buffered marginal U+k_S(Z) alone need not give either metric above. Conversely, retaining Z as an additional root and keeping the ACTUAL nonlinear mean k_S(Z) gives a conditional Gaussian mixture, which this obstruction does not prohibit; it simply has not eliminated the nonlinear insertion.

## 4. The same obstruction reaches the ACTUAL E mark

For every retained S,U define the scalar map

    T_(S,U)(t)=r{g(S+a g(x))-g(S+a g(x+t))}.

Its derivative lies between -r a Lambda² and -r a lambda². It is globally monotone and bi-Lipschitz, with lower slope r a lambda². Therefore for any scalar Gaussian G_(S,U), inserted before these SAME original terminal queries,

    W2(Law(T_(S,U)(k_S(Z))|S,U),
                    Law(T_(S,U)(G_(S,U))|S,U))
      >=r a lambda² W2(Law(k_S|S),Law(G_(S,U)))
      >=lambda² c_k r a rho.                          (7)

The first step follows either from one-dimensional quantile transport or by applying the Lipschitz inverse of T to any coupling. In particular this is not merely a lower bound for an unused hidden noise coordinate. Its left side is the actual paired E observable in(1), with its base and terminal roots retained.

The actual centered energy e=||E-EE||₂ is comparable to the native mark. Pointwise,

    |E|<=Lambda³ r a rho |Z|.

The exact positive Z derivative is the product of the three original derivatives, so D_ZE>=lambda³ r a rho. Conditional Gaussian integration by parts gives Cov_Z(E,Z)>=lambda³ r a rho at every S,U. Hence

    lambda³ r a rho <=e<=Lambda³ r a rho.              (8)

Equations(7) and(8) imply an integrated conditional error at least

    (lambda² c_k/Lambda³)e
        =0.0188056141414921... e.                      (9)

Thus no uniform gain A^eta e with eta>0 is possible for this SINGLE-GAUSSIAN insertion at the actual native weak scales. The source is genuinely the one-potential K3 twin with e comparable to A^3.9; no foreign large-energy marker is used.

## 5. Exact boundary of the negative result

The result excludes:

- a high-order joint Gaussian law for the complete base/shifted pair after Z is integrated, conditional on the retained terminal root;
- or replacing the shift by one conditional Gaussian before the SAME paired terminal observable, with S and the base U retained.

It does NOT exclude:

- the buffered marginal comparison after the base observer is removed;
- a reference keeping the actual nonlinear k_S(Z) or a properly declared Gaussian mixture over additional roots;
- a SAME-E weak coefficient/current comparison that preserves the covariance and every required observer rather than controlling the stronger metrics above;
- an executed adjoint/skew construction based on the paired source;
- the scalar direct-Gram construction, which retains the original paired VALUES and never makes this Gaussian substitution.

In particular the present stationary-pair marginal interface cannot be promoted into the excluded joint interface by identifying its source-zero carrier with an analytical innovation. Any weaker positive insertion port must state which joint observables it returns and account for the exact noise covariance, mean and spawned currents.

## 6. Independent diagnostics

The adjacent checker verifies the exact Gaussian projection for401 phases, including pure even/odd cases, and optimizes Gaussian insertions before the literal nonlinear E for several a,epsilon,S,U values. All946 checks pass. The numerical minimum is above the analytical uniform lower bound. The proof, rather than these optimizations, establishes uniformity.
