# C2 stage-two transfer by analytical smoothing, with the original VALUE bill

2026-10-05. A finite bounded-order extension of the sealed nonlinear bridge repair. Smoothing is used only in the proof. The executing graph is the original stage-two graph on the original source.

## 0. Result and exact scope

Assume the original source is an anchored convex gradient

    g = grad U,  U in C2,  g(0)=0,  0 <= Dg <= A I,
    0 < A <= 1/36.

No modulus of continuity of Dg, bound on D2g, or third source derivative is assumed. The native compiler's other radius, caller, curl, active-dimension, clock, filter, mode, and precision guards remain mandatory.

If g acts in a fixed orthogonal decomposition into blocks of dimension at most b, with no requirement to supply that basis to the executor, the unchanged stage-two finite original-VALUE graph admits

    integrated W2 error to N(m3(g)(Z), I)
        <= {110 A^(8/3) b^(2/3)
             + delta1 A^2 + sqrt(11/8) deltaL A^3
             + (deltaQ/2) A^2
             + delta_out A C_F2 + Lambda_comp A^4} sqrt(D)
             + restored absolute floors.                         (0.1)

The constant 110 is a deliberately conservative explicit bound, derived below. All quantities are as declared in the sealed stage-two theorem, with the pair rule built using this actual b. Choosing the old delta1<=A^2, deltaL<=A, deltaQ<=A^2, delta_out<=A^3 keeps all quadrature and completion terms O(A^4 sqrt(D)).

This is uniform over the entire C2-potential/Hessian-interval class WITHIN the stated bounded-block structure. Neither the source shape nor its Hessian modulus must be fixed as A varies. The b-dependent target estimate holds at arbitrary ambient D. It is stronger than the earlier fixed-shape qualitative little-o statement, but is not order four.

For an unrestricted D-dimensional source, take the single block b=D. The same construction and proof give

    error <= 110 A^(8/3) D^(7/6)
             + quadrature/completion/floors.                      (0.2)

Thus a quantitative gain survives uniformly over sources at every fixed D, with its dimension loss exposed. There is NO claim of O(A^(8/3) sqrt(D)) for arbitrary high-dimensional sources. Against the existing O(A^2 sqrt(D)) coarse certificate, the ratio in this new sufficient bound is O((AD)^(2/3)); improvement is certified in the corresponding small-AD regime. An upper-bound comparison is not a lower bound or impossibility theorem outside that regime.

The private root dimension remains 5D. The complete original-VALUE bill is unchanged:

    Q_rule/certificate_setup + Q_captured + N_out N_B
       + (5+3J1+5JL+6JQ) N_out N_E
       + Q_known/numerical/replay.                                (0.3)

There is no exact smoothed-gradient leaf, finite smoothing grid, additional smoothing root, or multiplier of source calls in (0.3). All native repetitions, live caller captures, original HVP sites and complete replays in this bill are still present. The original resummed shift is retained.

## 1. The graph and the two nonlinear remainder maps

Use exactly the sealed graph, rules, original VALUE sites, and genuine conditional Gaussian rows. Its inner output is

    F_Q(g)=g(S)
          + sum_single p C_g(S,d_g(U))
          + sum_linear q C_g(S,e_g(U,V))
          + sum_quadratic r Q_g(S,d_g(U),d_g(V)),

where

    v=x/2+N/2, w=x/4+N/2+M/4,
    S=x-g(v)+g(g(w)),
    d_g(U)=g(v)-g(U),
    e_g(U,V)=g(U)-g(U-g(V))-g(g(w)),
    C_g(S,a)=[g(S+a)-g(S-a)]/2,
    Q_g(S,a,b)=[g(S+a+b)+g(S-a-b)
                    -g(S+a-b)-g(S-a+b)]/8.

All clocks have positive mass-one weights and the sealed source-independent true conditional one-time or ordered-pair rows. Each U and V has unit Gaussian row norm and its x coefficient has absolute value at most one. Let psi_Q,g(x)=E F_Q(g)(x) over the four inner Gaussian roots. Define, only as analytical functions,

    J_g = R1(psi_Q,g-g),
    R_g = m3(g)-R1 g,
    R1 f(z)=int_0^1 E f(tz+sqrt(1-t^2)G) dt.

J_g is the finite stencil's exact-outer residual mean. R_g is the genuine-history canonical nonlinear remainder. The target discrepancy satisfies the exact identity

    R1 psi_Q,g - m3(g) = J_g-R_g.                    (1.1)

This identity is why an unsuppressed baseline smoothing error need not be paid.

## 2. Gaussian weak transport, without higher derivatives

The detailed independent derivation is in BASELINE-SMOOTHING-TRANSFER.md. The common lemma is short enough to state here.

For X~N(m,c^2 I_d), c>0, let T(x)=x-F(x) with F in C1 and ||DF||<=l<1. T is a global orientation-preserving diffeomorphism, even when DF is nonsymmetric. The relative entropy identity, followed by Gaussian integration by parts, gives

    2 KL(T#Law(X) || Law(X))
       <= E|F(X)|^2/c^2 + d l^2/(1-l).

Indeed the linear displacement term cancels tr DF; the remaining log-determinant series begins at degree two. It requires no derivative of DF. Pinsker therefore gives, for any vector-valued f with sup|f|<=delta,

    |E[f(T(X))-f(X)]|
       <= delta sqrt(E|F(X)|^2/c^2+d l^2/(1-l)).       (2.1)

Condition on each original private Gaussian history/row root before applying this inequality. In the outer resolvent c=sqrt(1-t^2); the only singularity is integrable:

    int_0^1 dt/sqrt(1-t^2)=pi/2.

Consequently a source-uniform L2 displacement of order A sqrt(d) produces a weak bound of order A delta sqrt(d). This is a conditional mean comparison, not a pathwise perturbation theorem and not an improved numerical-error floor.

## 3. Two coherent-source stability bounds

Let g,h be anchored C1 gradients with the same Jacobian interval [0,A I], and let delta=sup|g-h|<infinity on one d-dimensional block.

### 3.1 Genuine-history remainder

Use the actual OU future X_s=e^-s x+xi_s. Its nested force F2,g has

    ||D_x F2,g|| <= L=A/2+A^2/4,
    ||F2,g||_2 <= A(1+A)sqrt(d),
    |F2,g-F2,h| <= (1+A)delta.

For f=g-h, the exact remainder difference is

    f(x-F2,g)-f(x)
       + h(x-F2,g)-h(x-F2,h).

Apply (2.1) to its first part and Lipschitzness to the second. This proves

    ||R_g-R_h||_2 <= A delta C_d(A),
    C_d(A)=(1+A)+sqrt(d)[(1+A)pi/2
                              +(1/2+A/4)/sqrt(1-L)]
               < 5sqrt(d)  for A<=1/2.              (3.1)

The entire genuine future path is used only in this proof, never as an executed producer.

### 3.2 Actual finite residual stencil

In F_Q(g), the signed coefficients of terminal g VALUES sum to one and their total absolute mass is 7/2. Subtracting its original g(x) baseline lets every coherent source difference f appear as f(T_g(x))-f(x). Each actual terminal map T_g is one of

    S, S+/-d_g(U), S+/-e_g(U,V),
    S+/-d_g(U)+/-d_g(V).

For every fixed inner root, these maps have

    ||DT_g-I|| <= L_*=(7/2)A+A^2/4<1  for A<=1/36.

Put k=sqrt(3/8), d0=1+1/sqrt(2), e0=1+k,

    H_sum=2+11/(2sqrt(2))+A[1+(9/2)k].

The absolute-weighted sum of their stationary L2 displacement norms is at most A H_sum sqrt(d). Coherently replacing g by h in the terminal arguments costs at most

    A delta [14+(11/2)A].

For example S_g-S_h is bounded by (2+A)delta, d_g-d_h by 2delta, and e_g-e_h by (3+2A)delta. These estimates retain all nested occurrences; they are not single-site replacements at frozen descendants. Applying (2.1) terminal by terminal yields

    ||J_g-J_h||_2 <= A delta H_d(A),
    H_d(A)=14+(11/2)A
       +sqrt(d)[(pi/2)H_sum
              +(7/2)(7/2+A/4)/sqrt(1-L_*)]
             < 37sqrt(d)  for A<=1/36.              (3.2)

The bound is uniform in every finite positive node list meeting the stated row laws, including lists whose lengths and times vary with A. No comparison of unretained multi-time histories across distinct quadrature nodes occurs.

## 4. Analytical mollification and the complete bias inequality

Define the anchored Gaussian mollification

    h_e(x)=E[g(x+eY)-g(eY)],   Y~N(0,I),  e>0.

It is used ONLY in the proof. On each original d-block,

    h_e(0)=0, 0<=Dh_e<=A I,
    sup|g-h_e|<=2 A e sqrt(d),
    Lip(Dh_e)<=A/[sqrt(2pi)e],
    Lip(D2h_e)<=A/[sqrt(2)e^2].                      (4.1)

The dimension-free induced tensor bounds follow by integrating against first and second Gaussian derivative kernels and centering Dg by (A/2)I. Isotropic convolution retains every original orthogonal block without revealing its basis.

Apply the smooth stage-two theorem to h_e, on the SAME finite clock rules used by g. Use the identity

    J_g-R_g=(J_g-J_he)+(J_he-R_he)+(R_he-R_g).         (4.2)

The first and last terms are bounded by (3.2) and (3.1). Summing squared block errors and using d_j<=b gives the following explicit bound for the unchanged original source:

    ||R1 psi_Q,g-m3(g)||_2
      <= {84 A^2 e sqrt(b)
          + A^4 [p_A sqrt(b+2)/e
                    + q_A sqrt((b+2)(b+4))/e^2]
          +delta1 A^2+sqrt(11/8)deltaL A^3
                    +(deltaQ/2)A^2}sqrt(D),          (4.3)

where

    p_A=[d0 e0+(A/2)e0^2+1/2]/sqrt(2pi),
    q_A=[(d0+A e0)^3+4d0^3+A^3 e0^3]/(6sqrt(2)).

The coefficient 84 is 2(37+5), with no smoothness constant hidden as A varies. For A<=1/36, p_A<1.32 and q_A<3. The displayed epsilon^-1 and epsilon^-2 growth is paid in full.

Choose the proof parameter

    e=A^(2/3)b^(1/6).

No restriction e<=1 is needed in (4.1)--(4.3). Since b>=1,

    sqrt(b+2)<=sqrt(3b),
    sqrt((b+2)(b+4))<=sqrt(15)b,
    A^(2/3)b^(-1/3)<=1.

The three smoothing terms in (4.3) are at most

    [84+1.32sqrt(3)+3sqrt(15)] A^(8/3)b^(2/3)
       < 98 A^(8/3)b^(2/3)
       < 110 A^(8/3)b^(2/3).                        (4.4)

This establishes the conservative main constant in (0.1). The smoother-source A^4 bound has not been relabeled as a general-C2 A^4 bound. The certified interpolation exponent is 8/3.

## 5. Finite outer rule, native completion, and full cost

All weak comparisons above are made under the exact outer R1. Only after that proof do we replace it by the original finite positive outer operator Q_out. Since the original finite inner graph still has

    ||psi_Q,g||_2 <= A C_F2 sqrt(D),
    C_F2=1+A[3/2+5/(2sqrt(2))+A(1+2sqrt(3/8))],

the finite outer difference contributes exactly the allowed

    delta_out A C_F2 sqrt(D).

There is no unproved uniform estimate at an outer quadrature node close to t=1. The sealed L2 Gaussian operator rule handles the outer replacement separately.

The executed original source requires only the Hessian interval for its ports: normalized private first<=9A, curl<25A^2, O(A^2)(|z|+sqrt(D)) residual energy, actual captured origins, baseline/residual complete-bank independence, coisometry, positive fills/buffers, and every original guard. Consequently the existing guarded own-mean service contributes Lambda_comp A^4 sqrt(D) plus its absolute floors. Combining means with that actual completed law gives (0.1).

Every original VALUE and HVP occurrence is exactly the one already listed by the sealed graph. In particular:

- Residual raw VALUES: (5+3J1+5JL+6JQ)N_out per occurrence.
- Private Gaussian dimension: 5D.
- Caller origins and all source/anchor/scale dependencies remain live.
- One original HVP is charged per recorded original VALUE site for each requested first/adjoint sweep, with complete replay of discarded primals.
- No HVP is differentiated. No mollified value or derivative is queried.
- Source-dependent or coefficient-version changes invalidate caches under their full original keys.
- Original absolute VALUE, row, weight, moment, mode, compiler, and arithmetic floors are unchanged. The weak coherent-source estimates above do not reduce arbitrary numerical floors.

The pair-rule constant must use the actual block parameter:

    C_b=12*8^b*(1-2^-21)^(-(3b+2)/4).

For unrestricted b=D its exponential size enters the tolerance logarithm of the explicit positive rule; it cannot be called a dimension-independent constant. The displayed sealed bound gives marginal node count

    O((D+log(1/delta))^2)

and pair count O((D+log(1/delta))^4), with its very large fixed holomorphy constants and lower-order logarithms bounded in this notation. Source-independent coefficient setup and numerical work remain separately charged. If b or D itself varies with A, its displayed bias factors and node/guard bills remain part of the A-dependent cost; it is not treated as a fixed constant. For fixed b and fixed grade these are polylogarithmic in 1/A, and the new transfer introduces zero inverse-A original-VALUE exponent. Vector arithmetic is still O(D) per original leaf, not dimension-independent work.

### 5.1 Conditional source normalization

For a live exterior anchor a and a scale s>=0, the actual normalized source

    f(y)=s[g(a+s y)-g(a)],   alpha=A s^2,

is anchored, is a gradient, obeys 0<=Df<=alpha I, and retains the same orthogonal block sizes. Therefore the present C2 transfer applies with A replaced by alpha whenever 0<alpha<=1/36 and every new native guard holds. When alpha=0, use the literal zero-source branch rather than the positive-radius or mollification formulas. No source modulus or smoothing derivative constant needs to be assumed fixed under this normalization.

This is a structural closure statement, not an endpoint induction. Every f VALUE must execute its original g(a+s y) and actual g(a) capture (with complete-key reuse only where valid); all a and s dependencies remain live in firsts/adjoints, including the scale prefactor. Changed anchors/scales and discarded primals pay their complete original-g replay. The physical output scaling and the outer program's actual caller/radius/precision bill must still be restored; the normalized alpha^(8/3) bound is not silently a physical-coordinate rate.

## 6. Why an exact smoothed oracle or finite smoothing grid is unnecessary

A naive implementation would replace each source VALUE by a Gaussian convolution. That is not a permitted finite producer. A deterministic positive grid can implement an approximate anchored convolution, but its source-node count and all anchor/HVP/replay costs must be charged; FINITE-SMOOTHING-IMPLEMENTATION.md supplies such an optional construction.

The present theorem avoids that cost honestly rather than omitting it. The single-point approximation QY=0 has projected W2 error sqrt(d), and its anchored macro source is exactly

    g(x+e*0)-g(e*0)=g(x).

Equivalently, use g and h_e directly in (4.2). Their uniform difference is enough because BOTH residual maps have the extra-A weak stability. No finite approximation is actually added to the graph. Thus smoothing is a comparison argument, while the implementation remains precisely the original finite graph whose own ports and bill were already sealed.

By contrast, a naive replacement of the whole m3 target pays O(Ae sqrt(D)) instead of the residual O(A^2e) transfer. Balancing Ae with A^4/e^2 only gives order A^2 at fixed dimension. That bottleneck is a limitation of that naive sufficient certificate, not a general impossibility result. Preserving and canceling the original baseline is the constructive escape.

## 7. Remaining limits and next precise gate

This is one improved canonical-m3 mean-law component. It neither completes the full endpoint join nor proves an all-order finite recurrence. It does not change the genuine target to a conditional-mean surrogate. No paused unsmoothed Picard/path-grid or cold-start route has been resumed.

For unrestricted high dimension the current proof has two explicit losses: Gaussian near-identity transport costs sqrt(d), and the smooth stage-two third-order remainder costs sqrt((d+2)(d+4)). Their balanced certificate gives the b^(2/3) factor. There is also a simple global cap, requiring no smoothness: the original residual energy and genuine-history Lipschitz estimate give

    ||J_g||_2 <= A^2[3/2+5/(2sqrt(2))+A(1+2k)]sqrt(D),
    ||R_g||_2 <= A^2(1+A)sqrt(D).

Their sum is below 5 A^2 sqrt(D) for A<=1/36. Thus the exact-outer target error may be bounded by the minimum of this cap and (4.3), before adding the actual finite-outer/completion/floor terms. An arbitrary source is not assigned a worse guarantee merely because D is large.

Removing the dimensional loss requires a new dimension-controlled conditional weak norm estimate or a sharper shifted remainder, with the actual retained row geometry and one-energy norm preserved. Merely freezing B_e or C_e as constants would be invalid.

Within the present bound, the exponent ceiling comes from balancing A^2e against A^4/e^2. A further finite accuracy step would need a new cancellation/stability estimate, not an invocation of an unproved all-order recurrence. The viable next construction is therefore to sharpen the SAME finite residual/current comparison under the true ancestry and existing native ports, rather than installing an unpriced smoothing oracle or replacing the Gaussian genealogy.

## 8. Provenance

Sealed stage-two entry:

    ../nonlinear-bridge-stage-two/STAGE-TWO-POSITIVE-VALUE-THEOREM.md
    SHA256 c7275de9261bce2eeaee1ba054351055b1e236590340ee7b1554ff064ffa2f13

Sealed stage-two manifest:

    SHA256 cd1111e28dc399643f27bcdd0f97dc735ada196f8a43f220211500bcb45be3ae

The separate proof, optional implementation, and independent audit are contained in this new packet. No sealed stage-one or stage-two file is modified.
