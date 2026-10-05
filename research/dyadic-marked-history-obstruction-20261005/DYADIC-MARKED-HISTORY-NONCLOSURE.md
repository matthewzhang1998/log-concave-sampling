# A positive marked-history martingale and the exact failure of pair-local dyadic closure

2026-10-05. A new source-qualified obstruction and analytical replacement identity. No path grid, nested-history sampler, paused Picard program, cold source, external upload, or publication is implemented here.

## Result

The positive dyadic conditional-variance mechanism extends to the true coherent pair W=(F2,F3-F2) as a block martingale identity. It has a dimension-safe, one-Hilbert-Schmidt remainder bound, and logarithmically many analytical refinement levels suffice for any declared covariance accuracy. However, its coefficient sources do not close under the fixed-span, two-endpoint VALUE fields proved for an additive prefix.

This is a demonstrated obstruction to that specific transfer, not an impossibility theorem for another marked-history compiler. An actual scalar smooth convex gradient gives a nonzero TWO-INTERVAL conditional Hoeffding component of F2. Its covariance is positive and is absent from an additive sum of the exact singleton midpoint contributions. This issue precedes native admission, clock accuracy, and smoothing; it is neither the previously stated lack of a Delta first bound nor the already known missing cross block between F2 and Delta.

The full block target remains positive, including its history cross block. There is no finite public-logarithmic original-VALUE contract for it in this packet.

## 1. Sources and exact target

Use the original stationary OU process, g=grad U, g(0)=0, 0<=Dg<=AI, A<=1/2, and the true future histories

    F0,t=0,
    Fj,t=integral_0^infinity e^-s g(X_(t+s)-F_(j-1,t+s)) ds.

Let Y=X0 and W=(F2,0,F3,0-F2,0), in physical slot order in R^(2D). Its conditional block covariance B(Y)=Cov(W|Y) gives the true F3 covariance after the readout [I I]. All ancestors share one coherent future.

The inspected input manifests are pinned in INPUT-PINS.json. In particular the dyadic covariance packet has manifest 2365fc791230cd2d73be22964abf45c071825771ca9d6e860b483a7843a67a0a, and the full prefix-LAW join has manifest 130131317183c6f00ea3566c187296912db7850081904322215a3ec9d5eb37f2. Their positive prefix theorem is unchanged.

## 2. The exact positive block martingale

Fix T>0. Let S_k be the sigma field of Y and the OU samples at the dyadic endpoints jT/2^k, 0<=j<=2^k. Thus S_0 includes the terminal X_T. Put

    M_k=E[W|S_k],
    D_k=M_(k+1)-M_k.

For every finite K,

    B(Y)=Cov(M_0|Y)
          +sum_(k=0)^(K-1) E[D_k D_k^*|Y]
          +E[Cov(W|S_K)|Y].                             (1)

Every term is positive semidefinite. Martingale orthogonality is conditional on Y. Within each W, the F2 and Delta slots are coherent; no independent Delta bank appears. [I I] applied to the whole identity gives Cov(F3|Y), with every cross term.

Equation (1) is an analytical reference. In particular M_k is not an available conditional-mean oracle or a finite original-VALUE source.

### Full path sensitivity and one-HS remainder

For a deterministic path variation h with suitable compact support, differentiating the finite-depth history recursion and retaining matrix order gives

    |D Fj,0[h]| <= integral_0^infinity beta_j(s)|h(s)|ds,
    beta_j(s)=A e^-s sum_(r=0)^(j-1) (As)^r/r!.

The stacked pair W therefore has a valid sensitivity kernel

    beta(s)=2 beta_2(s)+beta_3(s)
           =A e^-s[3+3As+(As)^2/2].                    (2)

This is an O(A) complete first bound. The small L2 energy of Delta does not replace it by an O(A^3) first.

Conditioned on S_K, the bridges in separate mesh cells of length d=T/2^K and the OU future beyond T are independent Gaussian objects. A bridge point has private row norm at most C sqrt(d). Hence the derivative norm contributed by one cell I is at most C sqrt(d) integral_I beta. Gaussian Poincare, applied in an arbitrary physical unit direction of R^(2D), gives

    0<=Cov(W|S_K)
      <=C [d sum_I (integral_I beta)^2
                          +(integral_T^infinity beta)^2] I_(2D)
      <=C A^2[d^2+(1+T^4)e^(-2T)] I_(2D).              (3)

The first bound includes all within-cell ancestor effects through beta; it does not assume W additive over cells. The tail row norm is at most one at each time and its full contribution is paid by the displayed integral. Thus

    ||E[Cov(W|S_K)|Y]||HS
        <=C A^2[d^2+(1+T^4)e^(-2T)] sqrt(2D).           (4)

This is uniform in the retained skeleton values and Y, and remains valid after L2 integration. It is one operator bound followed by one HS conversion, with no product of two dimension-sized vector energies. It concerns the orthogonal conditional-mean approximation E[W|S_K]; it is not an error theorem for evaluating W on an interpolated skeleton path. For explicit scalar checks, integral beta^2<=8A^2 and (integral_T^infinity beta)^2<=80A^2(1+T^4)e^(-2T); the midpoint bridge row bound is sqrt(d/2).

Choose T=C log(1/epsilon) and d^2<=epsilon. Then the remainder is at most C epsilon A^2 sqrt(2D), with K=O(log(1/epsilon)+log log(1/epsilon)). The literal skeleton has 2^K cells. Logarithmic depth does not establish logarithmic node count. No such skeleton is executed in this packet.

## 3. What the level innovation actually contains

Conditioned on S_k, the 2^k new midpoint coordinates are independent Gaussians. Let xi=(xi_1,...,xi_N) denote their standardized roots and write

    H_k(xi)=E[W|S_k,xi],
    d_S(xi_S)=sum_(R subset S) (-1)^(|S|-|R|)
                                  E[H_k|S_k,xi_R]

for every nonempty S subset {1,...,N}. These are the ordinary conditional Hoeffding projections. Their physical slots still include the complete coherent W.

Orthogonality yields the further POSITIVE identity

    E[D_k D_k^*|S_k]
       =sum_(nonempty S) E[d_S d_S^*|S_k].              (5)

For the additive prefix integral in the sealed dyadic theorem, only singleton terms occur: changing a midpoint changes the conditional mean of its own bridge integral, and its source is a fixed function of that cell's two endpoints. That is exactly the locality permitting the 2-by-2 Gaussian operator Gamma(K_t) and fixed-span position compression.

For true F2 neither conclusion holds. A singleton conditional projection may already depend on the other coarse endpoints because F1 inside g sees the entire coherent future. In addition, the pair term d_{I,J} is genuinely nonzero in the following actual-history example.

## 4. An actual-history two-interval interaction witness

Take D=1 and, for 0<A<=1/100, the SAME original gradient

    g(x)=(A/2)x+(A/4) log cosh x.

It satisfies g(0)=0, is C-infinity, and

    l=A/4 <=g'(x)<=u=3A/4<A,
    g''(x)=(A/4) sech^2 x>0.                            (6)

Thus U'=g has the required convexity and Hessian bound. No enlarged function class, Hessian discontinuity, or unbounded higher derivative is involved.

Choose two ordered disjoint mesh cells I=[a,a+d] and J=[b,b+d], with a>=0 and b>=a+d. Let h and k be their OU midpoint regression hats, zero off their cells, positive in their interiors, and equal to one at their respective midpoints. For I, for example,

    h(s)=sinh(s-a)/sinh(d/2)           on [a,a+d/2],
    h(s)=sinh(a+d-s)/sinh(d/2)         on [a+d/2,a+d].

Use unstandardized midpoint values m,n; changing them shifts the conditional path by mh+nk. Freeze all coarse endpoints to zero, including Y, and condition the two selected midpoint values to zero. Other fine midpoint values and the remaining coherent future are integrated, not fixed. The resulting whole conditional path is centered Gaussian, with point variances at most one.

Write (Rf)(t)=integral_t^infinity e^(-(s-t))f(s)ds, so F1,t=R(g(X))(t). Let

    Z_t=X_t-F1,t,
    L_h(t)=D F1,t[h]=R(g'(X)h)(t),
    L_k(t)=D F1,t[k]=R(g'(X)k)(t).

Because hk=0 pointwise, D^2 F1[h,k]=R(g''(X)hk)=0. Ordered supports also give k L_h=0. Direct differentiation of the actual F2 gives EXACTLY

    D^2 F2,0[h,k]
      =integral_0^infinity e^-t g''(Z_t)
                          [-h(t)L_k(t)+L_h(t)L_k(t)]dt.             (7)

No expectation has been moved through g. All paths and all F1 ancestors are the actual ones. In particular the two terms in (7) have the displayed signs.

Under the zero conditioning, ||X_t||2<=1 and ||F1,t||2<=u, whence ||Z_t||2<=1+u. Markov's inequality gives

    E sech^2 Z_t >=(3/4) sech^2(2(1+u)) >=1/20.         (8)

The last scalar bound holds for A<=1/100. Also l Rh<=L_h<=u Rh and l Rk<=L_k<=u Rk pathwise. Put

    J0=R(h Rk)(0)>0,
    K0=R((Rh)(Rk))(0).

Tonelli and the ordered supports give

    J0=integral_I integral_J e^-t h(s)k(t) dt ds,
    K0=integral_I integral_J e^-t(1-e^-s)h(s)k(t)dt ds,
    0<=K0<=J0.                                         (9)

Let f(m,n)=E[F2,0|S_k=0, midpoint_I=m, midpoint_J=n], integrating every other path coordinate. Bounded derivatives and the exponential envelopes justify differentiation under this exact conditional expectation. Equations (6)-(9) imply

    partial_m partial_n f(0,0)
       <=(A/4)[-l J0/20+u^2 K0]
       <=-(11/6400)A^2 J0 <0.                         (10)

Indeed l/20=A/80 and u^2<=9A/1600. This is a quantitative, strictly nonzero mixed conditional derivative of a TRUE force history.

A square-integrable additive function f(m,n)=f_I(m)+f_J(n)+constant has zero mixed distributional derivative. Our f is smooth and has (10), so its two-variable Hoeffding component is nonzero. Therefore d_{I,J} has strictly positive F2-slot variance in (5). Continuity in the finitely many coarse endpoints extends this fact from the all-zero conditioning to an open set of positive Gaussian probability. The interaction is not a null-set conditioning artifact.

A convenient closed-form scalar check uses

    integral_I h(s)ds=2 tanh(d/4),
    integral_I e^-s h(s)ds=(d/2)e^(-(a+d/2)).

Thus J0=[2 tanh(d/4)](d/2)e^(-(b+d/2)) and
K0/J0=1-(d/2)e^(-(a+d/2))/[2 tanh(d/4)]. These verify (9) without constructing or simulating a path grid.

### 4.1 The omitted covariance is at the critical A^4 grade

Fix the horizon and a level with two ordered cells independently of A, and write g=A f with

    f(x)=x/2+(1/4)log cosh x.

Taylor expansion of the TRUE outer composition, without replacing an ancestor by its mean, gives

    F2,0=A Rf(X)(0)-A^2 B2(X)+O_Lp(A^3),
    B2(X)=R(f'(X) Rf(X))(0),                            (11)

for any fixed finite p under the actual stationary law. The remainder follows from bounded f'' and fixed Gaussian moments of Rf; conditional versions hold for fixed finite coarse values, and the integrated versions retain the actual Gaussian caller law. No derivative modulus beyond the explicit witness is being imported into the general class.

The pair Hoeffding projection of the first term is exactly zero, since Rf(X) is additive across disjoint bridge cells. For B2, direct differentiation along the two hats instead gives

    D^2 B2[h,k]=R(f''(X) h R(f'(X)k))(0)>0.              (12)

Every other product-rule term contains hk or k R(f'(X)h), both identically zero. Since f''>0 and f'>=1/4, (12) is strictly positive on every path. Its conditional expectation proves that the two-midpoint projection of B2 is nonadditive for every fixed coarse skeleton.

Let Psi be that nonzero pair Hoeffding projection of B2. Orthogonal projection is Lp-bounded at this fixed two-slot order, so

    d_(I,J)(F2)=-A^2 Psi+O_Lp(A^3),
    E[d_(I,J)(F2)^2|coarse]=A^4 E[Psi^2|coarse]+O(A^5), (13)

with strictly positive leading coefficient. Integrated L2 conditional-covariance versions follow using p=4. Thus dropping these terms loses covariance at order A^4, precisely the grade the missing true-history port must improve; this is not just an arbitrarily tiny nonzero defect.

After averaging the coarse endpoints while retaining Y, define q_pair(Y)=E[Psi^2|Y] and c=||q_pair||_(L2(Y))>0. Taking independent coordinatewise copies of the same separable gradient in D dimensions gives a diagonal omitted F2 covariance block. Its integrated L2(Y;HS) norm is c A^4 sqrt(D)+O(A^5 sqrt(D)), with c>0 independent of A,D. No dimension-sized product of vector energies is used. Even exact clocks and exact native actions on the singleton terms cannot make the resulting omitted-interaction covariance error o(A^4 sqrt(D)). This is a witness against deleting the interaction terms in the actual singleton decomposition, not a lower bound against an implementation that separately supplies or reorganizes them.

## 5. Consequence for the proposed fixed-span lift

The exact covariance of the level innovation is not the sum of its exact singleton midpoint covariances. The missing pair component in (10) contributes a nonzero PSD block in (5), already inside the F2 marginal. This is separate from, and survives proper inclusion of, Cov(F2,Delta).

Consequently the proof of the sealed local formula cannot be reused by merely replacing its local integrand g(X_u) by g(X_u-F1,u), or its local C_d by a nominal marked C_d. The latter integrand depends on future bridges outside the cell. Conditional independence of those bridges does not imply additive conditional means or vanishing cross-interval Hoeffding components after the nonlinear ancestor.

The fixed-span pair operator theorem remains valid for any genuinely fixed square-integrable field of the local endpoint pair. What is missing is precisely a finite original-VALUE supplier for such a field that includes the marginalized complete environment and all interaction terms. The 2-by-2 wedge does not, by itself, manufacture that supplier or compress dependence on the full skeleton. No claim is made that interactions could never be reorganized or allocated into a different family of effective fields.

Sequentially revealing one midpoint at a time removes explicit subset sums, but then each martingale increment conditions on all previously revealed coordinates. It is not the fixed two-endpoint source of the prefix theorem, and the source/readout dimension and replay ancestry grow with that reveal history. This is a different representation of the same unresolved coefficient obligation, not a public-log program.

A further exact warning against an averaging shortcut comes from this same strictly convex scalar g: for a nondegenerate conditional ancestral F1,

    E[g(x-F1)|labels] > g(x-E[F1|labels]).

Thus replacing the coherent ancestor by its conditional mean inside g changes the target even with smooth g. An average of finite shifted atoms cannot be identified with the nested true-history source by such a swap.

## 6. Native source, readset, caller, and finite-cost boundary

The complete analytical source at a level is H_k=E[W|S_k,xi]. Its readset contains every revealed skeleton coordinate and the integrated coherent future. In (5) the same physical F2/Delta history is used in every required projection; an independent complete native service would own a fresh entire such source tape, while copies INSIDE one projection must have precisely the conditional coupling specified by that projection.

Represent the finite retained Gaussian skeleton by orthonormal innovation coordinates, including Y. Every path point's regression row on these coordinates has operator norm at most one. Equation (2) then gives O(A) full firsts for M_k and H_k, and O(A) retained-coordinate firsts, uniformly in k. Conditional averaging can preserve these analytical first bounds, but it does not turn H_k into an executing original-VALUE graph. Physical endpoint coordinates with nonorthonormal rows must retain their actual row factors.

The mark still has only the available complete first O(A). Normalization of the pair at output buffer u has radius O(A/sqrt(u)), not O(A^3/sqrt(u)). The block physical dimension is 2D; [I I] has norm sqrt(2), and its buffers, errors, and caller derivatives receive that factor. The local square action from the prefix theorem consumes a genuine-gradient VALUE source with symmetric PSD Jacobian. The general rectangular derivative of H_k is not such a supplied source. A formal positive Jacobian Gram is not admission of that native leaf.

Any future construction must supply and charge: the complete conditional environment for all interaction coefficients; original-g VALUE evaluation of every nonlinear shifted ancestor; full caller/first and required marked cuts; an admitted rectangular/block Gram action or a valid reduction to the existing square action; independent complete service tapes; physical readouts and keeps; and literal occurrence/root counts. Changed coarse endpoints, midpoint coordinates, source mode, conditional-copy keys, or anchors replay every affected ancestor. A caller may retain unread labels but may not inspect a consumed skeleton bank after a whole-law comparison.

The current packet supplies no executing occurrences of H_k. Its Q and d_native are therefore UNSUPPLIED, rather than assigned unit cost or called polylogarithmic. Equations (1)-(5) and (10) are useful analytical information, not a new block-LAW compiler. There are no finite error floors concealing interaction covariance.

## 7. Smoothing and exact scope

The witness is already C-infinity with bounded g'' and every higher derivative bounded. Small Gaussian smoothing, followed by subtraction of the smoothed value at zero to retain g(0)=0, preserves strict convexity of g, the positive first-derivative interval, and the nonzero interaction by continuity. Thus smoothness alone does not restore additive pair-local closure. A smoothing-based compiler may still be possible, but it must pay smoothing bias and conditional ancestor replay, then explicitly represent the surviving interactions.

This result does not rule out a new multivariate fixed-span field construction, a joint current compiler, effective state augmentation, a controlled marked cluster expansion, or a source-qualified cancellation. It shows exactly why the current additive-prefix dyadic theorem plus the seven-VALUE local shifted atom is insufficient to claim the missing coherent (F2,Delta3) block port.

The full prefix-LAW join, canonical m3 theorem, and existing positive local marked current remain valid at their stated targets. No m4 accuracy increase, higher-rank cut theorem, whole-law same-carrier update, or new original-VALUE complexity theorem follows from this packet.
