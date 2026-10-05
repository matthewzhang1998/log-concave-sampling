# Independent audit: genuine joint-bridge quadrature and the square route

2026-10-05. This is a new audit artifact. No previously sealed file is changed.

## Results

1. The positive, genuine jointly sampled bridge quadrature gives a rigorous finite RAW covariance port, with original VALUE leaves, the correct retained endpoints, private first at most `L=A w^(3/2)`, and retained-pair first at most `Aw`.
2. Its covariance error can be certified as `C L^2 sqrt(D)/M`, using stationarity and a spectral integration-by-parts argument. The simpler path-increment argument gives only `C L^2 sqrt(D)/sqrt(M)` and is valid but unnecessarily expensive.
3. A public-log deterministic-node quadrature does not uniformly buy any extra power of `w` by itself. For every such positive actual-bridge rule with a fixed positive fraction of its weight in the middle half, some smooth anchored monotone scalar source has conditional covariance defect at least `c A^2 w^3/J^2`. This is an obstruction to this specified construction, not an impossibility theorem for different finite sources, corrections, or LAW services.
4. The bridge martingale-square identity and its original-VALUE source bounds are correct. A public-log, dimension-uniform compression of its outer square integral remains an additional theorem. The exact clock change alone does not supply it.

All conditional covariance norms below are `L2` over the actual stationary Gaussian endpoint pair, with Hilbert-Schmidt norm inside. This qualification is important: no uniform pointwise-in-endpoints moment error is claimed.

## 1. Finite RAW port

Let `q=e^(-delta)>=1/2`, `w=1-q`, and

    P = integral_0^delta e^(-s) g(X_s) ds,

where the OU process is stationary, `X=X_0`, `Y=X_delta`, and `g(0)=0`, `0<=Dg<=A I`. Partition `[0,delta]` into panels of maximum length `Delta<=delta/M`. On each panel use positive quadrature weights with EXACT mass for the measure `e^(-s)ds`. Define

    Q = sum_j beta_j g(X_(s_j)),   sum_j beta_j=w,

where every node value belongs to the SAME jointly sampled true OU bridge conditional on `(X,Y)`. A scalar Markov-bridge factorization samples this finite Gaussian vector using finitely many independent D-roots. Its conditional rows have the exact bridge covariance; sharing one arbitrary root in place of this covariance is not this construction.

### Actual source firsts

The conditional bridge means are `a_s X+b_s Y`, with

    a_s=sinh(delta-s)/sinh(delta),
    b_s=sinh(s)/sinh(delta),
    a_s,b_s>=0, a_s+b_s<=1.

Every private Gaussian row has norm `sigma_s<=sqrt(w)`. The triangle inequality, positivity, and the derivative bound therefore give

    Lip_private Q <= A sum_j beta_j sigma_(s_j) <= A w^(3/2)=L,
    Lip_(X,Y) Q <= A sum_j beta_j sqrt(a_j^2+b_j^2) <= Aw.

The same private bound holds for the analytical full bridge integral. Also `D_X Q=sum beta_j a_j Dg` is symmetric PSD and at most `Aw I`, so the coherent shifted consumer keeps its contraction property when `Aw<=1`.

### Conditional mean accuracy and literal cost

Refine the existing positive analytic bridge-mean panels until their width is at most `delta/M`; retain the geometrically graded endpoint panels and positive degree-`2n-1` quadrature with `n=O(log(1/epsilon))`. Exact positive weighted Gaussian quadrature is one way to enforce each panel's scalar mass. The old ellipse proof remains valid on the refined panels and the panel masses sum to `w`, so

    ||E[P|X,Y]-E[Q|X,Y]||_2 <= epsilon A w sqrt(D).

There are `O(M+log(1/epsilon))` panels and

    J=O((M+log(1/epsilon)) log(1/epsilon))

original VALUE leaves. All scalar rows, clocks, and weights must be frozen and numerically charged as usual.

## 2. Improved strong and covariance bounds

The elementary stationary increment estimate yields

    ||P-Q||_2 <= A w sqrt(2D delta/M)
              <= 2L sqrt(D)/sqrt(M).

There is a stronger estimate requiring no second derivative of `g`.

Set `nu(ds)=e^(-s)ds-sum_j beta_j delta_(s_j)(ds)` and let `F(s)=nu([0,s])`. Exact panel masses make `F` zero at every panel boundary; positivity gives `|F(s)|<=Delta` everywhere. Thus `integral F^2<=delta Delta^2`. For each spectral frequency `lambda>=0`, Stieltjes integration by parts gives

    integral integral e^(-lambda|s-t|) nu(ds)nu(dt)
      = 2lambda integral F(s)^2 ds
        -lambda^2 integral integral e^(-lambda|s-t|)F(s)F(t)dsdt
      <= 2lambda delta Delta^2.                         (1)

The double integral subtracted in (1) is nonnegative because the exponential kernel is positive semidefinite. The formula is valid for the finite atomic measure: its cumulative function has bounded variation, the boundary terms vanish, and the diagonal Dirac term in the mixed distributional derivative is `2lambda delta_(s=t)`.

For the stationary reversible OU semigroup, apply its scalar spectral measure separately to each component of `g`. The total first spectral moment is the Gaussian Dirichlet energy,

    integral lambda dmu_g(lambda)
        = E ||Dg||_HS^2 <= A^2 D.

The constant component is killed by the exact total mass. Integrating (1) proves

    ||P-Q||_2 <= sqrt(2) A sqrt(D delta) Delta
              <= 4 L sqrt(D)/M,                       (2)

where `delta<=2w` was used in the last deliberately conservative bound. This is an integrated stationary-endpoint result, exactly the norm required by the port.

Conditional Gaussian Poincare gives `opCov(P|X,Y), opCov(Q|X,Y)<=L^2`. For any centered random vectors `U,V`,

    ||Cov(U,V)||_HS <= sqrt(opCov(U)) (E|V|^2)^(1/2).

Use `P-Q` in the second slot and expand the covariance difference. Integrating conditional squares gives

    ||Cov(P|X,Y)-Cov(Q|X,Y)||_(2;HS)
        <= 2L ||P-Q||_2
        <= 8 L^2 sqrt(D)/M.                            (3)

Hence the requested two-source conditional weak comparison has the certified bill

    epsilon A^2 w sqrt(D)
      + C A^3 w^3 sqrt(D)/(h M)
      + C A^4 w^(9/2) sqrt(D)/h^2.                     (4)

It pays both third Stein defects, compares at the original endpoint pair, and then pays the already justified coherent `A^3 w sqrt(D)` shift. The old common-root source cannot be substituted into (2)-(3).

### Cost-qualified grade check

At

    w=A^(4/5), h=A^(9/5), M=ceil(A^(-1/5)),
    eta=A^(8/5), variable-buffer mu=A^(1/2),

the covariance term in (4) has grade `19/5`; the third Stein term has grade `4`. The existing terms `A^2h`, `A^3w`, `A^5w^(-3/2)`, and `A^3sqrt(eta)` all have grade `19/5`. The mean floor has this grade when `epsilon<=A` and is smaller at the previously used `epsilon<=A^2`.

The smallest-buffer mean radius has exponent `1/5`; the native self-reserve radius has exponent `3/20`. Both have strict power margin. The pointwise physical residual term `A^2 mu^(-1/2)/sqrt(eta)` has exponent `19/20`, so the complete theorem must use the established weighted aggregate-first argument rather than claim every node has first `O(A)`. The corresponding aggregate term is `A^(7/4)/sqrt(w)=A^(27/20)`, which is smaller than `A`.

This checks the proposed RAW port and its scalar interface with the existing ledger. It does not replace a fresh audit of every imported LAW guard, all finite-precision floors, or the complete assembled graph. In particular `M` is explicitly an inverse-`A` power. Keeping the weaker increment bound instead uses `M=ceil(A^(-2/5))` for the same grade.

## 3. A quantitative obstruction for ordinary true-bridge quadrature

Fix any deterministic positive actual-bridge rule with `J` nodes, total weight `w`, and at least `c0 w` weight in `[delta/4,3delta/4]`, where `0<c0<=1`. The rule may have arbitrary joint Gaussian implementation, provided it is the genuine joint bridge law. There is an admissible scalar `g` for which

    ||Cov(Q_g|X,Y)-Cov(P_g|X,Y)||_2
        >= (c0^4/2048) A^2 w^3/J^2.                   (5)

### Exact scalar kernel

Let `R_st=e^(-|s-t|)`. Let `m_s=a_sX+b_sY` and `C_st=E[m_s m_t]`. Explicitly

    C_st=a_sa_t+b_sb_t+q(a_sb_t+b_sa_t).

Both `C_st` and the bridge covariance `K_st=R_st-C_st` are nonnegative. For `phi_ell(x)=cos(x/ell)`, direct Gaussian integration gives

    E Cov(phi_ell(X_s),phi_ell(X_t)|X,Y)
      = exp(-ell^(-2))
          [cosh(R_st/ell^2)-cosh(C_st/ell^2)] >=0.     (6)

For stable floating-point evaluation, write each `exp(-ell^(-2)) cosh(z/ell^2)` as

    (exp(-(1-z)/ell^2)+exp(-(1+z)/ell^2))/2.

For `s` in the middle half, `sigma_s^2>=3w/32`. Set

    theta=c0^2/32,   ell^2=theta w/J.

Then (6) on the diagonal is at least `1/4`. Cauchy-Schwarz over middle nodes gives

    E Var(sum_j beta_j phi_ell(X_(s_j)) | X,Y)
        >= c0^2 w^2/(4J).                             (7)

For the true integral, subtracting conditioning only decreases integrated variance, and

    Cov(phi_ell(X_s),phi_ell(X_t))
        <= (1/2) exp(-|s-t|/(2ell^2)).

Here `1-e^(-r)>=r/2` for `0<=r<=delta<=log 2`. Therefore

    E Var(integral e^(-s) phi_ell(X_s)ds | X,Y)
        <= 4w ell^2.                                  (8)

Now put `f_ell(x)=(Aell/2)(cos(x/ell)-1)`. Multiplying (7)-(8) by `A^2ell^2/4` shows that its expected conditional variance defect is at least

    [theta c0^2/16-theta^2] A^2w^3/J^2
       = (c0^4/1024) A^2w^3/J^2.                     (9)

The oscillatory function by itself need not be a monotone gradient, so use the THREE admissible smooth functions

    g0(x)=Ax/2,
    g+(x)=Ax/2+f_ell(x),
    g-(x)=Ax/2-f_ell(x).

They vanish at zero and each has derivative in `[0,A]`. If `d(g)` denotes the expected conditional variance defect, exact quadratic polarization says

    d(g+)+d(g-)-2d(g0)=2d(f_ell).

At least one of `|d(g0)|,|d(g+)|,|d(g-)|` is therefore at least half (9). The `L2` conditional covariance defect is at least its absolute expectation, proving (5).

In particular, for `J` bounded by a fixed polynomial in public logarithms, `J^(-2)` cannot be replaced uniformly by any positive power of `w`. This result excludes only the specified uncorrected jointly sampled positive quadrature. It does not address endpoint-only degeneracies, source-adaptive moment repairs, alternative Gaussian correlations, nonquadrature sources, or a separately proved positive square-LAW service. It also does not contradict excellent approximation for each fixed smooth `g`; the worst-case admissible family oscillates at scale `ell`.

## 4. What ordinary random clocks do and do not fix

Let `mu(ds)=e^(-s)ds/w`, take `N` independent clocks from `mu`, independent of the entire true bridge, and return `Q_clock=(w/N)sum_i g(X_(S_i))`. Conditional on the full path this is an unbiased estimator of `P`, but conditional total covariance gives the exact identity

    Cov(Q_clock|X,Y)=Cov(P|X,Y)
       +(w^2/N) E[Var_(S~mu)(g(X_S))|X,Y].

The second term is positive semidefinite and generally nonzero, including for a linear source. Unbiased clock integration therefore does not match the desired covariance. Independent stratification reduces its size but does not cancel it. A pair-clock U-statistic can estimate a moment without bias; such a matrix estimator is not automatically the covariance of an allowed finite positive vector source.

There are two further interface issues. A continuously random clock evaluates a nondifferentiable sample path at a random time, so a complete Gaussian-private-first bound does not follow from the fixed-time bridge row bound. Also a single random clock shared by every physical coordinate can generate a rank-one retained-endpoint covariance error whose HS norm scales like `D`; one cannot assume the required `sqrt(D)` ledger. These observations exclude no suitably designed, separately proved randomized-clock construction. They explain why ordinary sampling is not already the missing port.

## 5. Bridge martingale square: verified identity and VALUE bounds

Condition on `X,Y`. The pinned bridge satisfies

    dX_s=[-coth(delta-s)X_s+csch(delta-s)Y]ds+sqrt(2)dW_s.

Let

    V_s(x,Y)=integral_s^delta e^(-u)E[g(X_u)|X_s=x,Y]du,
    B_s=D_x V_s.

The martingale for the additive functional has stochastic differential `sqrt(2)B_s(X_s,Y)dW_s`. Since `Dg` is symmetric,

    Cov(P|X,Y)=2 integral_0^delta E[B_s(X_s,Y)^2|X,Y]ds. (10)

The identity follows first for smooth approximations and then by the bounded first and Dirichlet estimates; no modulus of continuity of `Dg` is needed.

For `r=delta-s`, `t=u-s`, the future bridge coefficients are

    a=sinh(r-t)/sinh r, b=sinh t/sinh r,
    sigma^2=2 sinh t sinh(r-t)/sinh r.

Thus

    B_s=integral_s^delta e^(-u)a E Dg(a x+bY+sigma N)du,
    0<=B_s<=A c_s I,
    c_s=e^(-s)[1/2-r/(exp(2r)-1)]<=e^(-s)r/2.          (11)

An analytical common-root primitive is

    F_s(x,Y,N)=integral_s^delta e^(-u)(a/sigma)
        [g(a x+bY+sigma N)-g(a x+bY)]du.              (12)

It is a genuine gradient in `N`, with `E D_NF_s=B_s`. Its source bounds are direct original-VALUE estimates:

    Lip_N F_s <= A c_s,
    |F_s| <= A c_s |N|,
    Lip_(x,Y) F_s <= (pi/sqrt(2)) A e^(-s)sqrt(r).     (13)

Indeed `a+b<=1`, and

    a/sigma <= sqrt((r-t)/(2rt)),
    integral_0^r a/sigma dt <= pi sqrt(r)/(2sqrt(2)).

Each of the two VALUE sites must be differentiated/replayed on its actual ancestry; their subtraction is part of the source. The bound in (13) uses two bounded `Dg` terms and does not differentiate `Dg`.

After substituting the true marginal bridge `x=X_s=a_sX+b_sY+sigma_sZ`, the additional private `Z` first is `O(A r)` because `sigma_s^2<=2r`. The retained endpoint first remains `O(A sqrt(r))`. These estimates support an integrable square-service candidate, but alone do not prove the finite positive outer quadrature, its LAW error, or the complete endpoint derivative of a proposed reserve implementation.

For the linear test `g(x)=Ax`, (10) agrees exactly with

    Var(P|X,Y)=A^2[(1-q^2)/4-q^2 delta^2/(1-q^2)]
             =2A^2 integral_0^delta c_s^2 ds.

### Exact clock change is insufficient to identify a fixed field

The pinned bridge also has the representation

    X_s=a_s[X+W_(tau_s)]+b_sY,
    tau_s=sigma_s^2/a_s^2
         =2sinh(delta)sinh(s)/sinh(delta-s).

The normalized centered bridge becomes stationary OU under clock `(1/2)log(tau_s)`. Nevertheless its integrand is

    g(a_sX+b_sY+sigma_s Z_clock),

which depends explicitly on the clock through three different scalar coefficients. This is not a fixed-field stationary covariance identity. A dimension-uniform analytic product/conditional-semigroup theorem, or another constructive outer-square compression, must be established rather than inferred from this reparametrization.

## 6. Numerical check

The companion `check_joint_quadrature_and_martingale.py` checks the exact cosine kernel and its predicted inverse-node covariance scale, the spectral cumulative-measure estimate, and the linear martingale identity. Its checks supplement the proofs above; they are not a substitute for them.

All 12 cosine-family cases, 36 spectral-bound cases, and 5 martingale-identity cases passed. The report is `joint-quadrature-and-martingale-checks.json`.

## 7. Reviewed inputs and audit scope

Reviewed source snapshots, SHA256:

- `../COVARIANCE-QUALIFIED-RAW-PREFIX.md`: `d28057b710d06140ce21b3d8cced4bd8822e2bd7da49f40576e23522e0286249`
- `../BRIDGE-MARTINGALE-SQUARE-GATE.md`: `2b82590fd3f103d50a2ed2b8ec6f45195a4f20604a2e25de805064ddffa38fed`
- `/workspace/shared/smoothed-bridge-prefix-20261005/SMOOTHED-BRIDGE-PREFIX.md`: `c26e44e932f0626cd7101d092d679f7b416d4d1e71adaf09e3cf4a48941e4258`
- `/workspace/shared/combined-mean-reentry-20261005/POST-SYNTHESIS-REENTRY-AND-NEXT-PORT.md`: `53748adb366902bc875e5d74871bcd1a9cda74ec22524bba699f38fc4c918f49`

The reviewed RAW note now explicitly chooses the endpoint cutoff at most `min(epsilon w/64,delta/M)`. This was the only required correction to its finite-prefix argument found in this audit; without the minimum, a generic choice `M>>1/epsilon` could leave an endpoint replacement interval wider than the claimed maximum mesh. It does not affect the displayed power-law parameters.

The RAW prefix proof, its two-source comparison assuming the imported Stein lemma, and the displayed scalar grade/cost calculations pass. The square note correctly labels the outer square-clock theorem and actual native caller paths as unresolved. Its martingale identity and literal primitive bounds pass. Imported full LAW-service theorems, their numerical guards, and a separate final own-mean compiler with the enlarged inverse-power tape have not been independently re-proved here.
