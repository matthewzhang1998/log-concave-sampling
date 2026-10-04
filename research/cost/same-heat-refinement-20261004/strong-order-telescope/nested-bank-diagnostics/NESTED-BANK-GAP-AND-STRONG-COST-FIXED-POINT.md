# Literal nested banks and the strong-allocation cost fixed point

2026-10-04. This diagnostic does not by itself admit the complete hidden source. It isolates the exact sampling and indexing issue, gives an original-C2 nonlinear fixture, and states the sufficient strong-allocation recurrence. The retained graph, full caller/anchor, and actual history topology are separate obligations in the continuation construction.

## 1. The exact nested-bank identity

Let E_1,...,E_N be independent complete copies, conditional on the same original caller theta, of a vector E with covariance V. All within-copy aliases are preserved. For M<=N, use the first M copies for the coarser mean. Then

    Cov(average_N E - average_M E | theta) = (1/M-1/N) V.

This is an exact identity: the first M coefficients are 1/N-1/M and the other N-M coefficients are 1/N. A common independent Gaussian output carrier cancels in this difference. It cannot reduce this sampling variance.

Consequently a law estimate with error proportional to 1/N does not establish a literal nested-bank strong estimate proportional to 1/N. Under nonzero variance, the strong error has scale 1/sqrt(M). A complete paired increment must pay this fact; a Wasserstein coupling selected after the fact is not the executed source coupling.

The identity holds without symmetry, Gaussianity, differentiability, or independence between coordinates. Independence is across complete copies only, not across records inside an occurrence. It also holds for covariance conditional on a captured caller.

## 2. A quadratic posterior fixture

Take V(x)=lambda x^2/2 and the exact posterior at center zero. Its canonical force is

    F(W)=sqrt(A) V'(sqrt(A/(1+A lambda)) W)
        =epsilon W,  epsilon=A lambda/sqrt(1+A lambda).

The positive finite empirical statistic Y_N=Z+average_N F has the exact law

    N(0,1+epsilon^2/N),

so its W2 error from N(0,1) is sqrt(1+epsilon^2/N)-1, asymptotic to epsilon^2/(2N). With the SAME Z and literal nested force samples, however,

    ||Y_N-Y_M||_2 = epsilon sqrt(1/M-1/N).

This proves the weak-versus-literal-strong gap even for a smooth quadratic source. It is not a lower bound against every legal coupling.

Indeed this special Gaussian source permits an explicit orthogonal carrier rotation. Put xi_N=epsilon/sqrt(N) and use independent U,V. Define

    Z_N=(U-xi_N V)/sqrt(1+xi_N^2),
    G_N=(xi_N U+V)/sqrt(1+xi_N^2).

Within each level, Z_N and G_N are independent standard Gaussians, while

    Z_N+xi_N G_N=sqrt(1+xi_N^2)U.

Completing the sample-mean direction G_N to an orthogonal iid sample bank preserves its full marginal. This achieves a cross-level strong gap of the law grade in this promised linear case. It depends on the known linear direction and covariance. It does not provide a rotation that Gaussianizes a nonlinear empirical mean, nor justify sharing a sampled hidden source as an external caller.

## 3. A nonlinear fixture requiring only original C2 VALUES

At each fixed A choose the globally smooth potential

    V_A(x)=x^2/4 - (A^6.8/4) cos(x/A^3.4).

Then

    V_A''(x)=1/2+(1/4)cos(x/A^3.4),

which lies in [1/4,3/4]. At the original Gaussian query x=sqrt(A)W its canonical force is

    sqrt(A)V_A'(sqrt(A)W)
      =(A/2)W + E_A(W),
    E_A(W)=(A^3.9/4)sin(W/A^2.9).

Thus E_A has first radius at most A/4, is odd and zero at the actual zero, and has exact variance

    Var E_A = A^7.8 (1-exp(-2 A^-5.8))/32.

The nested-bank identity still applies exactly. It cannot be improved merely because the residual is tiny and has a bounded original-Hessian implementation. Only original gradient VALUES and their ordinary HVP firsts are used; no third derivative is invoked. This is a local original-force Gaussian-query diagnostic, not a claim that sqrt(A)W is the exact posterior for V_A, and not a lower bound for the whole nonlinear hidden construction.

## 4. Strong allocations deliberately spend a factor two in order

Assume an ACTUAL finite hierarchy of complete paired canonical forces F_l, conditional on the original caller, with

    ||Delta_l-E Delta_l||_Lp <= Lambda sqrt(Dim) A^l, l>=1,
    ||F_0-E F_0||_Lp <= Lambda sqrt(Dim) A,
    Q_l <= Lambda A^-2(l-1)_+.

Every pair has the named fine and coarse finite marginals. No derivative gain of Delta_l is assumed: its complete private first may remain Lambda A and caller first Lambda sqrt(A).

For a target level j>=1 take independent complete banks across levels, with

    n_0(j) >= Lambda_(j,p) A^-2(j-1),
    n_l(j) >= Lambda_(j,p) A^-2(j-l), 1<=l<=j.

Round upward and include at least one sample at every level. The constants allocate the fixed number j+1 of terms and the required finite moment list. The vector-valued independent-sum moment inequality gives strong mean error

    ||hat m_j-E F_j||_Lp <= Lambda_(j,p) sqrt(Dim) A^j.

No oracle computes E F_j: this is the exact analytical target of

    hat m_j=average(F_0)+sum_(l=1)^j average(Delta_l).

The complete original query bill has, for every l>=1,

    n_l(j) Q_l <= Lambda A^-2(j-l) A^-2(l-1)
                 = Lambda A^-2(j-1).

The base term has the same exponent. Finite sums, actual zeros, all original first/adjoint sweeps and logarithmic finite restorations are included in Lambda only after their circuits are counted. The number of independent Gaussian records is genuinely inverse-heat and must be recorded separately; it is not polylogarithmic.

More generally, if a complete source hierarchy has exponent beta(l-1) and the same strong energy A^l, this allocation yields

    beta_next = max(2,beta).

In particular beta=2 is a fixed point. This is a sufficient upper bound, not an optimization or a lower bound. It is precisely why accepting slope two may close a linear complete family even if weak-law-optimal counts fail to maintain slope one.

## 5. The adjacent-level offset is restored by the physical CW7 rows

Freeze a final maximum J and one decoder length k_J for EVERY level j<=J. The positive statistic uses the same variance A/(k_J+1), with one independent Gaussian keep and the usual one shared observation pool. Couple the actual means at j and j-1 by literal nested banks for all common levels. Independence across different levels is preserved. At level j the added top increment has its own fresh complete record.

The coarser mean's strong error is A^(j-1), not A^j. Triangle inequality plus the exact finite mean difference therefore gives a safe adjacent normalized mean gap O(A^(j-1)); nesting is an explicit admissible coupling, not a claimed extra cancellation.

For the actual CW7 hidden center, write the original physical history rows as

    q=b+sum_s C_s E g_s,   sum_s ||C_s||<=Lambda A.

With canonical F_s=sqrt(A)(g_s-anchor_s), absorb the actual anchors into b and put B_s=C_s/A. Then

    q=b0+sqrt(A)sum_s B_s E F_s,   sum_s ||B_s||<=Lambda.

The normalized adjacent mean gap A^(j-1) consequently gives a physical statistic/state gap A^(j-1/2). The final canonical force readout is sqrt(A)grad V, so its adjacent force gap is A^j again. The two physical sqrt(A) factors restore the single offset exactly:

    (j-1)+1/2+1/2=j.

A source bias in physical state law is propagated through the same history row with factor A, because sqrt(A) times the canonical-force Lipschitz factor sqrt(A) is A. This is the necessary normalization, not a generic common-statistic contraction: the decoder response to a common T perturbation is order one.

For correlated force/history occurrences one must telescope and copy the COMPLETE history source. The independent-bank theorem does not authorize destroying its aliases. The construction needs its actual causal caller contract, bounded affine history resolvent and physical A predecessor. A depth-indexed conclusion requires checking these conditions for the selected consumer.

## 6. Anchors and known limits of this diagnostic

All-private-zero execution has an exact telescoping identity regardless of the sample counts:

    hat m_j(theta;0)=F_j(theta;0).

A centered implementation may therefore retain the actual finest zero plus differences with their actual offsets. Centered variance bounds do not require a small zero value. If an uncentered source energy or retained origin envelope is needed, prove its deterministic zero recursion separately; never evaluate an Lp Gaussian inequality at zero.

A complete-family claim still needs the actual retained graph and error, caller/anchor derivatives, protected late source, finite original numerical floors, enlarged independent Gaussian tape, and genuine-gradient reference/restoration ports. The empirical signal is kept in full. It is not declared a gradient merely because it approximates a mean.

The companion script checks 272 exact algebraic/numerical identities: nested variances, Gaussian-law and legal rotation fixtures, the nonlinear C2 variance/first envelope, all strong allocation exponents through level15, the physical adjacent offset, and the general cost-slope recurrence. It does not numerically prove the source-admission theorem.
