# Positive finite quadrature for the true conditional two-time ancestry

2026-10-05. New stage-two result; none of the sealed stage-one files is modified.

## Result

The required two-time quadrature exists under the declared bounded-block qualification. It does **not** require analyticity of the source. There is a source-independent positive rule, with O_b(log^4(1/delta)) ordered-pair nodes, approximating every conditional pair function with global linear growth |F_Y(u,v)|<=H sqrt(|Y|^2+|u|^2+|v|^2) in L2 of the retained Gaussian variables to error delta H sqrt(d), for d<=b. Applied blockwise, the actual ancestry correction has H=sqrt(11/8)A^3 and the actual Q correction has H=A^2/2. Thus delta_Q<=A^2 gives an order-four quadrature budget. A polynomial-growth extension also gives the optional sharper C_(b,B) delta A^3 sqrt(D) bound for Q under ||D^2g||<=BA. Constants are independent of A, D, and the particular source; b is a declared bound on block dimension.

Only two fresh D-dimensional Gaussian roots are needed in addition to x,N,M. Including the outer G gives **5D raw private coordinates**. Two different history times in the same node use the true pair covariance. Sharing these same two roots across different finite nodes makes no claim about cross-node history laws and is legitimate because only expectations of individual node outputs are used.

A fully rational certificate proves the uniform standardized-innovation covariance bound K >= (1/16)I. A completely explicit safe holomorphy radius is eta=2^-24. The constants are intentionally conservative; their size does not create an inverse-A power.

## 1. Ordered pair and its true conditional law

Write a>0 for the earlier time and h>0 for the gap. Let

    sigma_s = sqrt(1-exp(-2s)),
    f(s) = 2s/sqrt(exp(2s)-1),
    R(a,h) = [ f(a)                 f(a)(a-1)                 ]
             [ exp(-a) f(h)        exp(-a) f(h)(2a+h-1)      ],
    K(a,h) = I_2-R(a,h)R(a,h)^T.

The unconditioned Markov innovations

    Z1=(X_a-exp(-a)x)/sigma_a,
    Z2=(X_(a+h)-exp(-h)X_a)/sigma_h

are independent standard Gaussians and independent of x. Their cross-covariance with (N,M) is exactly R. Consequently, conditioned on Y=(x,N,M), generate

    (W1,W2)^T = R(a,h)(N,M)^T + K(a,h)^(1/2)(L1,L2)^T,
    U = exp(-a)x + sigma_a W1,
    V = exp(-h)U + sigma_h W2.

Then (Y,U,V) has exactly the law (Y,X_a,X_(a+h)). In particular its means are q(a)Y,q(a+h)Y and its covariance is

    Cov(U,V | Y)=exp(-h)-q(a).q(a+h),

with the corresponding exact marginal conditional variances. Each complete scalar row on (x,N,M,L1,L2) has Euclidean norm one.

A stable source-independent 2-by-2 square root is

    K^(1/2) = [K+sqrt(det K) I]/sqrt(tr K+2 sqrt(det K)).

The denominator and det K are uniformly separated from zero by the certificate below. No source Hessian, source-specific covariance root, or D-by-D matrix factorization is used. Near a=0 or h=0, sigma_s can be computed with expm1 and enters only as a scalar coefficient. All quadrature nodes below are strictly positive, so literal repeated times are unnecessary; the formula nevertheless has the correct continuous limiting law.

## 2. Certified uniform covariance margin

Put F(s)=f(s)^2=4s^2/(exp(2s)-1), with F(0)=0. There are no square-root singularities in

    R^T R = F(a) [1,a-1]^T[1,a-1]
            +exp(-2a) F(h) [1,2a+h-1]^T[1,2a+h-1].

The companion `certify_pair_covariance.py` uses only exact rational arithmetic. It subdivides [0,8]^2 into 2,792 dyadic boxes and certifies

    I-R^T R >= (1/10)I

on every box, including the axes. For exp(x), 0<=x<=16, it sums the rational Taylor polynomial and bounds the remaining positive tail by

    next_term / [1-x/(n+2)]

once n+2>x. The rational enclosure width is below 2^-72. For F(s), it uses

    F(s)=2s/E(2s), E(2s)=(exp(2s)-1)/(2s), E(0)=1,

and monotonicity of E to avoid division by zero at s=0. Interval matrix entries certify the positive leading diagonal and determinant of (9/10)I-R^T R. The supplied script is a finite check, rather than an unverified numerical eigenvalue assertion. Its result is recorded in `pair-covariance-certificate.json`.

Here are the analytic tails completing the global certificate. For s>=8,

    F(s)[1+(s-1)^2] <= 8s^4 exp(-2s)
                         <= 8*8^4*(2/5)^16 < 1/64.

We used exp(1)>5/2 and the fact that s^4 exp(-2s) decreases for s>=8. If a>=8, the first row has squared norm <1/64. Also F(h)<=2 and h^2 F(h)<=6 for every h>=0, by the positive exponential series. Hence the second row has squared norm at most

    exp(-2a)(16a^2+14) <= 1038*(2/5)^16 < 1/64.

If h>=8, the function exp(-2a)[1+(2a+h-1)^2] decreases in a>=0, since its derivative is -2 exp(-2a)(2a+h-2)^2. Thus the second row has squared norm at most F(h)[1+(h-1)^2]<1/64. The first row has squared norm <=9/10 when 0<=a<=8, by the compact certificate at h=0, and <1/64 when a>=8. Therefore outside [0,8]^2,

    ||R||^2 <= 9/10+1/64 < 15/16.

Since R R^T and R^T R have identical eigenvalues, this proves globally

    (1/16)I <= K(a,h) <= I.

## 3. Density holomorphy, with an explicit radius

This step handles the retained S(Y) and arbitrary nonanalytic source values directly.

For one block of dimension d<=b, let F_Y(u,v) be any measurable function with a uniform polynomial envelope

    |F_Y(u,v)| <= H (1+|Y|+|u|+|v|)^m.

Hilbert-valued functions are allowed. Define

    P_(a,h)F(Y)=E[F_Y(X_a,X_(a+h)) | Y].

For every real a,h>0, this L2(gamma_(3d))-valued function has a holomorphic extension on

    |z-a|<=eta a, |w-h|<=eta h, eta=2^-24,

bounded by C_(b,m) H. The function F_Y itself is never complexified. Only its conditional Gaussian density is complexified.

Here are explicit uniform estimates establishing that assertion. Let

    d(z)=sqrt(1-exp(-2z)),
    L(z,w)=[[d(z),0],[exp(-w)d(z),d(w)]],
    m(z,w)=(q(z)Y,q(z+w)Y),
    C(z,w)=L(z,w)K(z,w)L(z,w)^T.

Use the analytic square root positive on the positive real axis for d; it exists on Re z>0. In the larger relative polydisk of radius 1/4, the elementary estimates

    |f(z)|<=4, |z f(z)|<=8, |z exp(-z)|<=1

give ||R(z,w)||<=24. For example, 1-exp(-3a/2)>=3a/(2+3a), while |z|<=5a/4 and Re z>=3a/4; splitting a<=1 and a>=1 proves the displayed estimates. Cauchy's inequality on the radius-1/4 polydisk then gives, on radius eta<=1/8,

    ||R(z,w)-R(a,h)|| <=384 eta.

Writing L0=L(a,h), R0=R(a,h), K0=K(a,h), one has

    L0^(-1)L(z,w)
      =[[d(z)/d(a),0],
        [(exp(-w)-exp(-h))d(z)/d(h),d(w)/d(h)]].

The ratios d(z)/d(a) are bounded by 2 on radius 1/4, so their deviations from 1 are at most 16 eta by Cauchy's inequality. The lower-left entry is at most 4 eta; this follows from

    |exp(-w)-exp(-h)|<=eta h exp(-7h/8),
    h exp(-7h/8)/sqrt(1-exp(-2h))<=2.

Therefore

    ||L0^(-1)L(z,w)-I||<=24 eta,
    ||K(z,w)-K0||<=1152 eta,
    ||L0^(-1)C(z,w)L0^(-T)-K0||<=8192 eta.

The latter follows by inserting D=L0^(-1)L and expanding D K(z,w)D^T-K0. All these estimates hold for eta<=2^-24. They do not divide by a small conditional variance without first normalizing the Markov innovations.

The mean variation similarly satisfies

    |L0^(-1)(m(z,w)-m(a,h))|<=512 eta |Y|.

For the N,M part, use ||D R(z,w)-R0||<=432 eta. For the x part, the two normalized entries are

    [exp(-z)-exp(-a)]/d(a),
    exp(-z)[exp(-w)-exp(-h)]/d(h),

each bounded by 2 eta.

Now whiten using the real K0. The complex standardized covariance G and displacement beta obey

    ||G-I||<=16*8192 eta <=1/128,
    |beta|<=4*512 eta |Y| <=2^-13 |Y|.

In these real standardized coordinates r, the complex density has a canonical determinant branch near G=I and absolute value bounded by

    2^d (2pi)^(-d) exp(-|r|^2/8 + 8|beta|^2).

This is an integrable bound. The transformed real arguments u,v have magnitude at most C(|Y|+|r|), because ||L0||<=2, ||K0||<=1, and each real q row has norm at most one. Integrating the polynomial envelope gives

    |P_(z,w)F(Y)|
      <= C_(d,m) H (1+|Y|)^m exp(2^-23 |Y|^2).

Its square is integrable for Y~gamma_(3d), uniformly in a,h. Dominated integration proves scalar holomorphy, and the same domination, or Banach-valued Morera, proves L2-valued holomorphy. Local branches agree with the positive real Gaussian density, so patch over the relative-disk neighborhood. This proof requires neither analyticity nor a complex extension of g, S, or F_Y.

For the global linear envelope

    |F_Y(u,v)|<=H sqrt(|Y|^2+|u|^2+|v|^2),

one may use the fully explicit bound

    ||P_(z,w)F||2 <= C_b H sqrt(d),
    C_b = 12*8^b*(1-2^-21)^(-(3b+2)/4).

Indeed sqrt(|Y|^2+|u|^2+|v|^2)<=3|Y|+2|r|. Comparing exp(-|r|^2/8) to the N(0,4I_(2d)) density yields

    |P_(z,w)F(Y)|<=8^d H exp(2^-23|Y|^2)[3|Y|+6sqrt(d)].

The stated C_b follows by the exact Gaussian moment-generating function in dimension 3d and 3sqrt(3)+6<12. Consequently, use the quadrature below with internal accuracy delta/C_b to obtain output error at most delta H sqrt(d). The node count is O_b(log^4(1/delta)). For homogeneous degree-two envelopes, the analogous bound is C_(b,2) H sqrt(d), again because d<=b.

## 4. An actual positive finite tensor rule

Let nu_lambda(ds)=lambda exp(-lambda s)ds, for lambda in {1,2}. For target 0<delta<1, set epsilon=delta/32. Set ell=epsilon/8 and H=log(8/epsilon), so H>ell. Replace [0,ell] by an atom at ell/2 of its exact nu_lambda mass, and [H,infinity) by an atom at H+1 of its exact mass. For functions bounded by M on the positive real axis, these two replacements together have error at most (lambda+1)epsilon M/4.

Partition [ell,H] into dyadic panels [A,B], B<=2A. On every panel use the positive n-node Gaussian rule for the exact measure nu_lambda. It has positive interior nodes, degree-(2n-1) polynomial exactness, and exact panel mass. Choose

    rho=1+eta,
    n=ceil(log(C/(epsilon*(rho-1)))/(2 log rho)),

where C=128 suffices. The uncorrected marginal error is at most epsilon M: the head/tail error is at most 3epsilon M/4 and the middle error below is less than epsilon M/16. The Bernstein ellipse with parameter rho over a dyadic panel lies in the relative disks from Section 3. Banach-valued Chebyshev approximation and positivity therefore bound the middle error by

    4 rho/(rho-1) rho^(-2n) M.

This is the same positive weighted-Gaussian argument as the sealed one-time proof, now applied to the density-smoothed pair kernel. The number of marginal nodes is O(log^2(1/epsilon)).

Tensor the lambda=1 or lambda=2 rule in a with the lambda=1 rule in h. By a telescoping difference of the two probability quadratures, the L2 error is at most 2epsilon M before moment correction, and the number of pair nodes is

    J_pair=O(log^4(1/epsilon)).

All weights are positive and sum exactly to one. The source-independent rule approximates P_(a,h)F for every F in the stated polynomial-growth class, not just analytic source functions. A common list of nodes serves every block and every source with the declared b and envelope degree.

### Exact exponential moments

The marginal rule can additionally enforce

    sum p_i exp(-s_i)=lambda/(lambda+1)

with one positive atom. If its current moment m differs from the target m0, add by convex admixture a node with exponential value r=(m0+1)/2 when m<m0, or r=m0/2 when m>m0. Set alpha=(m0-m)/(r-m). This is positive, below one, and bounded by a constant times the moment error. The latter is at most epsilon by the scalar quadrature error for exp(-s). For lambda in {1,2}, alpha<=6epsilon. Since both old and added rules have exact mass one, the change in the F error is at most 2 alpha M. Thus each corrected marginal has error at most 13epsilon M, and the corrected tensor rule has error at most 26epsilon M<delta M for the specified epsilon=delta/32. The node-count order and positivity are unchanged.

For lambda=1 in both variables this gives exactly

    sum p_(i,j) exp(-(a_i+h_j))=1/4.

It preserves the desired exact linear-source cancellation of the added ancestry correction. The lambda=2 marginal can likewise enforce 2/3, although Q vanishes for a linear source and does not need this extra identity.

## 5. Application to the stage-two VALUE stencil

Keep exactly

    v=x/2+N/2, w=x/4+N/2+M/4,
    S=x-g(v)+g(g(w)),
    d_y=g(v)-g(y),
    e_(y,z)=g(y)-g(y-g(z))-g(g(w)).

Define C and Q exactly as in SECOND-CORRECTION-GAP.md. The two history-dependent conditional integrals are

    ancestry = int exp(-t-u) E[C(S,e_(X_t,X_(t+u))) | Y] dt du,

    quadratic = int exp(-t-u) E[Q(S,d_(X_t),d_(X_u)) | Y] dt du.

The first is directly the ordered-pair kernel under nu_1(da)nu_1(dh). For the second, Q(S,a,b)=Q(S,b,a), and the diagonal t=u has measure zero. Splitting into t<u and u<t gives exactly

    quadratic = int 2 exp(-2a-h)
                     E[Q(S,d_(X_a),d_(X_(a+h))) | Y] da dh,

which is nu_2(da)nu_1(dh). These are true pair laws, not a tensor product of incorrect one-time bridge laws.

The simplest application needs only g(0)=0 and ||Dg||<=A<=1/2 for the quadrature estimate. With

    R5=sqrt(|Y|^2+|y|^2+|z|^2),

one has the global bounds

    |C(S,e_(y,z))| <= sqrt(11/8) A^3 R5,
    |Q(S,d_y,d_z)| <= (A^2/2) R5.

The first uses |w|<=sqrt(3/8)|Y|. For the second, grouping opposite edges in the four-value difference gives |Q(S,a,b)|<=A min(|a|,|b|)/2. Then min(|v-y|,|v-z|)<=R5 because the linear map (Y,y,z)->(v-y,v-z) has squared norm 2 and averaging the two squared distances bounds their minimum. These are linear envelopes without B or higher source derivatives. The normalized linear-growth theorem therefore gives total quadrature error at most

    [sqrt(11/8) delta_anc A^3 + (delta_Q/2)A^2] sqrt(D).

Taking delta_anc<=A and delta_Q<=A^2 gives order A^4, with inverse-A exponent zero. Taking both at most A^2 is also valid. Quadrature nodes remain source-independent and their constants depend only on the declared b.

For an optional sharper A^3 Q envelope, impose on each declared block ||D^2g||<=B A. For arbitrary real y,z,

    |e_(y,z)| <= A^2(|z|+|w|),
    |C(S,e_(y,z))| <= A^3(|z|+|w|),
    |Q(S,d_y,d_z)| <= (B A^3/2)(|v|+|y|)(|v|+|z|).

The last inequality follows from the exact double integral of D^2g over the rectangle defining Q. These are polynomial envelopes of degrees one and two, with S(Y) retained. Section 3 therefore applies directly to the full finite VALUE stencil, rather than merely a derivative Taylor surrogate. The two finite positive rules yield per-block errors bounded by C_b delta A^3 sqrt(d) and C_b B delta A^3 sqrt(d), respectively. Squaring and summing over the orthogonal blocks gives

    total pair quadrature error <= C_(b,B) delta A^3 sqrt(D).

The C(S,d_t) term can retain the sealed one-time bridge rule and its independent remainder analysis; alternatively, this same pair theorem with an unused second history argument approximates it at C_b delta A^2 sqrt(D).

Using the optional degree-two envelope, choosing pair delta<=A gives an order-A^4 pair quadrature budget as well. In either version, J_pair=O_b(log^4(1/A)) and **inverse-A exponent zero**. This theorem does not itself prove the stage-two Taylor/stencil remainder, source first/curl guards, or final numerical floor budget; those are separate ledgers. It does close the true-two-time positive finite quadrature gate, including the retained Y and nonanalytic C3 sources.

## 6. Numerical implementation scope

The mathematical nodes and weights are real-arithmetic objects, as in the sealed theorem. Their actual finite-precision computation still needs an absolute tolerance ledger. All pair rows use source-independent scalar functions; the uniform K gap makes the displayed 2-by-2 square root stable. At each node one may compute scalar rows and renormalize with a declared absolute row error, then price the resulting source perturbation using the source Lipschitz bounds. No numerical clipping is part of the exact Gaussian theorem. Clipping a negative eigenvalue caused by rounding must be explicitly charged if used.

The same L1,L2 may be reused for all actual nodes and for both pair rules. Positivity here refers to time quadrature weights and to the final genuine Gaussian source law; signed arithmetic inside C and Q is not a signed probability law.
