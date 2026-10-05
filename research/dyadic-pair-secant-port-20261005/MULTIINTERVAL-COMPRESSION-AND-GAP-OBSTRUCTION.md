# Fixed-cluster translation, and the relative-gap obstruction

2026-10-05. Analytical audit of the Gaussian compression step. This note supplies no original-VALUE source, native action, history grid, simulated path, or all-interaction compiler.

## Main findings

1. A fixed ordered cluster of any number of intervals, with all internal geometry frozen, has exactly the same dimension-free complex translation wedge as its convex-hull endpoint pair. Arbitrary interactions within the cluster survive this reduction. The reason is an exact Markov projection onto the two extreme endpoints, not deletion of pair or higher Hoeffding components.
2. Moving relative gaps changes the stationary law of the physical source coordinates. There is no corresponding dimension-free contraction theorem for arbitrary fixed multi-interval fields. This fails on the real axis in a fixed L2 source space and on every nontrivial complex gap neighborhood for bounded sources.
3. More strongly, a source-independent finite positive quadrature over a continuum of relative gaps cannot have a dimension-free node count uniformly over all bounded positive physical pair fields. An explicit empirical-correlation test gives a necessary order-sqrt(D) node count for fixed accuracy. This is an obstruction to a generic theorem, not a lower bound for the particular nonlinear history fields if they have additional supplied structure.

## 1. Exact ordered-cluster theorem

Let X be the stationary D-dimensional OU process with covariance

    E[X_s X_t^*] = exp(-|s-t|) I_D.

Retain only the original endpoints (X_0,X_T), where T>0. Fix ordered intervals

    0=u_1 < v_1 <= u_2 < v_2 <= ... <= u_r < v_r=S<T.

Coincident interval endpoints are one coordinate. All interval lengths, relative gaps, marked positions, and conditional-copy conventions are fixed. Let C be their union. A source may depend on finitely many sites in C, the complete path restricted to C, and additional independent roots with a fixed law. Let its fixed square-integrable Hilbert-valued field be Phi; matrices with Hilbert-Schmidt norm are allowed.

For 0<=t<=L:=T-S define the actual conditional averaging operator

    A_t Phi = E[Phi((X_(t+s))_(s in C), eta) | X_0,X_T].

The source space carries the fixed stationary cluster law. Let P_C be conditional expectation of the canonical source onto its two extreme endpoints (X_0,X_S), integrating the independent eta roots as well. Then

    A_t = T_t^S P_C.                                      (1)

Here T_t^S is precisely the fixed-span pair operator. Indeed, conditional on (X_t,X_(t+S)), the path inside [t,t+S] is independent of the path outside that interval. Its internal conditional law depends only on the frozen relative geometry. P_C is a contraction. Equation (1) preserves every interaction in Phi; it averages the complete interacting field before translating it.

Use stationary orthonormal coordinates

    retained: (X_0, (X_T-e^(-T)X_0)/sqrt(1-e^(-2T))),
    extreme pair: (X_t, (X_(t+S)-e^(-S)X_t)/sqrt(1-e^(-2S))).

Put q=e^(-T) and S_v=sqrt(1-e^(-2v)). Their cross-correlation is

    K_z = [ e^(-z)    2q sinh(z)/S_T          ]
          [ 0         e^(z-L) S_S/S_T       ].          (2)

Extend the two extreme-pair coordinates to an orthonormal Gaussian coordinate system for the cluster and its private roots. All remaining coordinates are internal bridge innovations or independent marks. Their cross-correlation with the retained pair is zero. Thus the full cross-correlation matrix consists of K_z and zero rows, and has the same operator norm as K_z, regardless of r, site count, or D. Equivalently, use the fixed contraction P_C in (1).

For z=x+iy,

    det(I-K_z K_z*)
       = [(1-e^(-2x))(1-e^(-2(L-x))) - 4q^2 sin^2(y)]
           /(1-q^2).                                    (3)

The lower diagonal entry of I-K_z K_z* is positive for 0<=x<=L. Since

    (1-e^(-2x))(1-e^(-2(L-x))) >= 4e^(-2L)x(L-x),
    e^(-2L)>=q^2,
    sin^2(y)<=y^2,

we obtain

    0<x<L,   y^2<=x(L-x)  =>  ||K_z||op<=1.              (4)

On Gaussian chaos n, the conditional averaging is the symmetric tensor power of this cross-correlation. Its second quantization therefore extends holomorphically in the open disk |z-L/2|<L/2 and satisfies

    ||A_z||_(L2(cluster;H) -> L2(retained;H)) <= 1.        (5)

The same holds on the smaller two-sided wedge |Im z|<=min(Re z,L-Re z). It is dimension-free and also independent of cluster cardinality. Infinite path restrictions follow by Gaussian-cylinder approximation, or directly by second quantization of the first-chaos Hilbert-space contraction. At the real boundary points t=0,L use the actual conditional expectation. If S=T there is only one translation and no clock to compress.

## 2. Positive translation quadrature

Let a positive translation measure be supported on a lattice t_j=j h, 0<=j<=M, with L=Mh. Keep endpoint atoms exactly. Divide the remaining indices into dyadic panels by distance from the nearer boundary. Every panel has width at most its distance from either endpoint. A fixed Bernstein ellipse of parameter 2 around such a panel lies in (4).

Apply m-node positive Gaussian quadrature to the measure on each panel, keeping smaller panels exactly. Analytic operator-valued polynomial approximation and degree-(2m-1) exactness give

    || integral A_t dmu(t) - sum_j omega_j A_(t_j) ||
        <= C 4^(-m) mu([0,L]).                           (6)

The number of nodes is O(m log(M+1)), independent of D and r. The weights are positive and their total mass is exact. Consequently a PSD source field remains PSD under the real quadrature. For a field of scale A^4, (6) costs only C 4^(-m) mu([0,L]) times its actual L2-HS norm; no interaction is dropped or demoted in A-grade.

This is a node-count theorem, not a generic moment-computation theorem. For geometric lattice weights exp(-beta t), panel moments are Taylor coefficients of the explicit geometric sum, as in the supplied dyadic pair note. For arbitrary inaccessible weights, one must separately supply their moments without enumeration.

For dyadic r-tuples, the differences among their indices are frozen in this theorem. Their common shift is compressed. Summing over every relative index pattern is still a separate task; this result does not turn the potentially many gap patterns into a public-logarithmic program.

## 3. Why a relative-gap clock is a different operator

Already take two physical sites with gap g>0. Their stationary law is mu_rho in R^(2D), where

    rho=e^(-g),    Sigma_rho = [ I   rho I ]
                                [ rho I  I ].

For a fixed physical field F(x,y), the gap-dependent scalar functional is

    L_rho F = integral F dmu_rho.                        (7)

It is a lower bound on any conditional-output operator: averaging that output over its retained endpoints gives (7), and averaging is a contraction.

Choose a fixed source reference mu_(rho_0), with 0<=rho_0<1, and put h=rho-rho_0. If |h|<1-rho_0, the exact likelihood-ratio calculation gives

    ||L_rho||_(L2(mu_(rho_0))->C)
       = [(1-h^2/(1+rho_0)^2)(1-h^2/(1-rho_0)^2)]^(-D/4).  (8)

This is strictly greater than one and exponential in D whenever rho!=rho_0. Outside the L2-integrability range the functional is unbounded. At rho_0=0 this simplifies to

    ||L_rho|| = (1-rho^2)^(-D/2).                       (9)

To prove (8), diagonalize into the normalized sum and difference coordinates. The two covariance ratios are

    lambda_+=1+h/(1+rho_0),
    lambda_-=1-h/(1-rho_0).

For a one-dimensional variance ratio lambda, the squared L2 norm of the Gaussian likelihood ratio is [lambda(2-lambda)]^(-1/2), finite exactly when 0<lambda<2. Multiply the two factors and then all D coordinate copies. Truncated nonnegative likelihood ratios show that the obstruction can be approached by bounded positive fields as well.

Whitening at each rho does not evade (8). It changes the physical field into

    F_rho(x,z) = F(x, rho x+sqrt(1-rho^2) z).

Its input Gaussian law is now fixed, but the source itself varies with the gap. Any theorem treating F_rho as one frozen field omits this dependency. A semigroup applied to y while keeping x captured ends with evaluation on a correlated diagonal; that diagonal evaluation is exactly where the claimed dimension-free L2 contraction fails.

## 4. Bounded fields also lack a dimension-free complex wedge

Let rho=a+ib satisfy |rho|<1 and b!=0. The complex Gaussian density defines the holomorphic extension of (7) for every bounded measurable fixed F. Its exact total variation, hence its complex L-infinity-to-scalar norm, is

    ||L_rho||_(Linf->C)
       = [ |1-rho^2| / (1-a^2) ]^(D/2).                (10)

For one coordinate pair, use the covariance eigenvalues 1+rho and 1-rho. The determinant of the real part of the inverse covariance is

    (1-a^2)/|1-rho^2|^2.

Integrating the absolute complex Gaussian density yields (10). Moreover,

    |1-rho^2|^2-(1-a^2)^2 = 2(1+a^2)b^2+b^4 > 0.

Thus every fixed nonreal gap point rho=e^(-z) gives exponential growth in D. The conclusion holds up to an absolute factor for real nonnegative fields bounded by one: rotate the complex measure, retain the positive half-plane of its real part, and average the rotation angle. Some such indicator has expectation modulus at least total variation divided by pi. Smooth bounded approximations give the same obstruction if smooth test fields are preferred. Multiplying the field by A^4 I_D makes it a PSD matrix field with real op bound A^4 and real HS bound A^4 sqrt(D).

This excludes a generic uniformly bounded complex wedge for gap variation on arbitrary bounded interacting sources. Positivity of the real source is insufficient.

## 5. Real positive gap quadrature needs dimension-dependent nodes in the unrestricted class

Let rho be mixed uniformly over a fixed interval [a,b] inside (0,1), of length ell=b-a. This is a legitimate positive gap measure after g=-log rho. Consider any m-node positive mass-one rule with correlations rho_1,...,rho_m in that interval. Define the physical statistic

    Z_D(x,y) = D^(-1) sum_(i=1)^D x_i y_i.

Under mu_rho,

    E Z_D=rho,      Var Z_D=(1+rho^2)/D <= 2/D.

Fix R>sqrt(2), let delta=R/sqrt(D), and take the bounded positive fixed source

    F_D(x,y) = indicator{ min_j |Z_D(x,y)-rho_j| <= delta }.

Every quadrature node satisfies E_(rho_j) F_D>=1-2/R^2, so the discrete mixture does too. For a single node, split the continuum integral into |rho-rho_j|<=2delta and its complement. The first part has length at most 4delta. On the complement Chebyshev bounds the probability by 2/[D(|rho-rho_j|-delta)^2]. Integrating both sides of the complement costs at most 4/(D delta). A union bound over the nodes therefore gives

    E_(continuous mixture) F_D
        <= (m/ell) [4delta + 4/(D delta)]
        = (4m/(ell sqrt(D))) (R+1/R).

Hence the quadrature error on this positive bounded field is at least

    1-2/R^2 - 4m(R+1/R)/(ell sqrt(D)).                 (11)

For example, take R=4. To make the uniform error <=1/4 requires m>=c ell sqrt(D), for an absolute c>0. With any fixed finite m the error tends to at least 7/8 as D increases. Indicators can again be replaced by smooth bounded fields with arbitrarily small finite-D loss. The PSD matrix field A^4 F_D I_D has the same obstruction at the A^4 sqrt(D) HS scale.

The test field depends on the proposed rule, as is appropriate when refuting a source-independent operator quadrature theorem. This is not a claim about a rule specifically tailored to one supplied history coefficient, nor a claim that all such indicators are realizable as true F2 Hoeffding covariances.

## 6. Exact readset and source boundary

The translation theorem requires all of the following.

- The ordered cluster geometry, all internal gaps and marked positions, source formula, conditioning/copy pattern, and independent private-root laws are fixed before the translation clock.
- The physical OU readset of the field is contained in the translated hull, except for genuinely independent owned roots. A field may use a deterministic captured parameter, with a uniform bound, but correlated external skeleton values cannot silently be treated as independent parameters.
- If the target conditions on additional coarse skeleton coordinates, those coordinates change the conditioning operator. The two-original-endpoint K_z above is then not its full cross-correlation matrix.
- If a purported cluster coefficient marginalizes a coherent history outside the hull, one must first prove that this produces one fixed endpoint/cluster field with the stated stationary law. Markov localization of the process does not localize an arbitrary nonlinear history functional.
- Conditional projection P_C is analytical. It does not supply an executing original-g-VALUE service. Native evaluation, original-site replay, ancestor and anchor replay, captured first bounds, common-input assembly, and occurrence/root counts must be supplied separately.
- Dimension-free contraction does not say that a literal source's derivative or replay cost is independent of r. Physical stacked endpoints, private roots, and actual pathwise dependencies retain their real row factors and counts.

The actual two-midpoint A^4 Hoeffding covariance from the inspected obstruction note is not removed by (1): if its full fixed field is supplied, it is retained inside Phi. The unresolved questions are whether the true coefficient has the required local readset, how its full environment is supplied, and how all relative-gap patterns are summed. The generic Gaussian theorem closes only common translation with frozen geometry.

## Source scope

Read in full:

- /workspace/shared/dyadic-bridge-square-20261005/DYADIC-POSITIVE-BRIDGE-COVARIANCE.md
- /workspace/shared/dyadic-marked-history-obstruction-20261005/DYADIC-MARKED-HISTORY-NONCLOSURE.md

No external literature, browser, external message, upload, or numerical experiment was used. The proofs above are independent Gaussian and Markov calculations based on the explicitly stated OU model.
