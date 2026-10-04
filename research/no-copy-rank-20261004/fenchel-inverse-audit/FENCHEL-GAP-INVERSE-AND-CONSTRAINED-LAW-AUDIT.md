# Fenchel gap: bounded first, finite inverse, and the actual constrained law

Independent bounded audit, 2026-10-04. Source: `p-native-twin/JOINT-GRADIENT-LIFTS-AND-CONSTRAINT-DESCENDANTS.md`, especially its (7)–(11). This is a test of the Fenchel-gap escape, not an impossibility theorem for all native producers.

## Result

The exact Fenchel gradient has a dimension-free bounded first and uses no third derivatives. However:

1. A finite inverse contraction gives an approximate inverse VALUE, generally **not an exact gradient field** in dimension two or higher.
2. A globally normalized positive Fenchel penalty changes the Gaussian roots, even for a quadratic primitive if there is an independent Gaussian auxiliary prior.
3. Conditional normalization preserves those roots but introduces a non-free conditional moment/normalizer, leaves the graph at every positive temperature, and recovers a singular graph law at zero temperature.
4. The exact graph law has a Hessian-dependent tangent current. Its current cannot be replaced by ordinary ambient Gaussian integration by parts.

The counterexamples and constants below are literal, and the adjacent Python checker verifies them. No HVP is used as an executed VALUE, and no HVP is differentiated.

## 1. Exact Fenchel gap and its actual first

Let g=grad V, with mI <= A(x)=Dg(x) <= LI, m=2/5, L=3/5. Assume g(0)=0; translation by g(0) gives the corresponding unanchored statements. Strong monotonicity makes g a bijection. Let h=g^{-1}=grad V* and

    Q(x,p)=V(x)+V*(p)-x·p.

Potential VALUES are only analytical here. Its gradient and Hessian are

    grad Q = (g(x)-p, h(p)-x),
    D²Q = [[A(x), -I], [-I, A(h(p))^{-1}]].

The exact inverse is C1, so Q is C2 under the stated C2-potential assumption. The two equivalent Bregman forms give

    (m/2)|x-h(p)|² <= Q <= (L/2)|x-h(p)|²,
    (1/(2L))|p-g(x)|² <= Q <= (1/(2m))|p-g(x)|².

Thus Q>=0 and Q=0 exactly on p=g(x). Both gradient blocks vanish on that graph.

The following sharp-for-the-block-envelopes bounds hold globally:

    -(sqrt(1261)-31)/30 I <= D²Q <= (31+sqrt(761))/20 I,
    lower endpoint = -0.15035206030431353,
    ||D²Q|| <= 2.9293114224133725.

For the lower bound compare to [[mI,-I],[-I,L^{-1}I]]; for the upper compare to [[LI,-I],[-I,m^{-1}I]]. Multiplying Q by κ multiplies this first bound by |κ|. The penalty field grad(Q/τ) therefore has first at most 2.929312/τ, not a τ-independent first.

On the graph, D²Q is PSD with kernel {(v,A(x)v)}. Its nonzero eigenvalues are λ+λ^{-1}, λ in spectrum A(x), and are at most 2.9.

### Concrete off-graph negative curvature

In one dimension set

    V(x)=x²/4 + (1-cos x)/10,
    g(x)=x/2 + sin(x)/10.

At (x,p)=(π,0), h(p)=0 and

    D²Q = [[.4,-1],[-1,5/3]],     det D²Q = -1/3.

Its negative eigenvalue is the lower endpoint above. Positive gap, strong separate convexity, and a minimum manifold do not imply joint convexity. The upper block envelope is also attained in this scalar example at (x,p)=(0,g(π)).

## 2. Finite original-VALUE inverse contraction: what is and is not obtained

Use the optimal constant step for the interval [m,L]:

    u_0=0,
    u_{n+1}=u_n+2(p-g(u_n)),
    h_N(p)=u_N.

The update contracts in u by q=1/5. It uses at most N original g VALUES; the anchored first query can be omitted if g(0)=0 is supplied. In particular,

    |h_N(p)-h(p)| <= (5/2) 5^{-N}|p|.

For |p|<=R and desired absolute error η, it suffices to take

    N=max(0, ceil(log((5R)/(2η))/log 5)).

The same formula with R=||P||_Lp gives an Lp VALUE error. For R=1, N=10 gives error <=2.56e-7; N=18 gives error <=6.5536e-13. The a posteriori test |g(u_N)-p|<=mη certifies inverse error <=η, with one extra VALUE if not already cached.

FIRST/adjoint verification may use the ordinary first sites of g. With J_n=Dh_n,

    J_0=0,
    J_{n+1}=(I-2A(u_n))J_n+2I,
    ||J_N|| <= (5/2)(1-5^{-N}).

Consequently the approximate field F_N=(g(x)-p,h_N(p)-x) has actual first at most 2.929312 by the corresponding two-by-two norm bound. These derivative products are FIRST verification, never executed VALUES. This is a valid small-first approximate vector field. It is not generally a genuine gradient.

### Exact smooth two-dimensional non-gradient fixture

Let ε=1/20, e=(1,0), v=(1,1), and

    V(z)=|z|²/4 - ε cos z_1 - (ε/2)cos(z_1+z_2).

Then g(0)=0 and

    A(z)=.5I+ε cos(z_1)eeᵀ+(ε/2)cos(z_1+z_2)vvᵀ.

The two perturbation norms sum to 2ε=.1, proving .4I<=A<=.6I everywhere. Take p=(π/4,π/4). The first two inverse iterates are

    u_1=(π/2,π/2),
    u_2=(π/2-2ε,π/2).

Writing A_i=A(u_i),

    A_1=.5I-(ε/2)vvᵀ,
    A_2=.5I+ε sin(2ε)eeᵀ-(ε/2)cos(2ε)vvᵀ,
    J_3-J_3ᵀ=8(A_2 A_1-A_1 A_2)
            =-4ε² sin(2ε) [[0,1],[-1,0]] != 0.

The off-diagonal skew is -0.000998334166468282. Thus replacing h by the finite inverse iterate destroys the claimed exact gradient, although VALUE accuracy and a bounded first remain valid. N=1 and N=2 happen to be integrable for the zero start; this does not extend to arbitrary N. Special commuting/affine cases are exceptions, not a general nonlinear construction.

### No uniform FIRST convergence count follows from bounded Hessians

The VALUE iteration count does not buy a corresponding FIRST approximation. For every prescribed N>=2, put δ=.1^(N+1)/2 and define

    f(t)=-t(1-t²)^3 for |t|<=1, and f(t)=0 otherwise,
    g_N(x)=.55x+.05δ f((x-1)/δ).

Here f is C2, -1<=f'<=32/49, g_N(0)=0, and .5<=g_N'<=.55+.05(32/49)<.6. At p=.55 the exact inverse is h(p)=1 and h'(p)=2. All first N queried iterates stay outside the bump; exactly

    u_j=1-(-.1)^j,
    Dh_N(.55)=(1-(-.1)^N)/.55,
    |Dh_N(.55)-Dh(.55)|>=.18  for all N>=2.

Proof: f'(t)=(1-t²)²(7t²-1) on [-1,1], whose minimum is -1 and maximum is 32/49; f, f', and f'' match zero at the endpoints. For each j<N the putative orbit has distance .1^j from 1, strictly greater than delta, so induction shows g_N(u_j)=.55u_j and g_N'(u_j)=.55 at every queried point. The exact linear recurrences therefore give both displayed formulas for u_j and Dh_N. At x=1, f(0)=0 and f'(0)=-1 give g_N(1)=.55 and g_N'(1)=.5. For even N the derivative error is at least 2-1/.55=2/11; for odd N>=3 it is at least 2-(1+.001)/.55=.18. The VALUE error is exactly .1^N at this point. The fixture may depend on N, which is precisely what disproves a uniform FIRST guarantee over the original bounded-Hessian class. Smooth compact variants can be made without changing this conclusion; the displayed C2 g already exceeds the required regularity for a C2 potential.

A modulus of continuity ω for A would instead yield the honest estimate

    ||J_N-Dh|| <= (5/2)q^N
      +5 sum_{j=0}^{N-1} q^{N-1-j} ω((5/2)q^j|p|).

No uniform such modulus was supplied. This paragraph is not needed for the exact-Q bound in §1; it explains why finite VALUE inversion does not preserve the literal derivative record automatically.

## 3. Positive penalties and normalization: exact root accounting

Without a confining root law, exp(-Q/τ) dx dp is not normalizable: integrating over each p fiber gives a positive lower bound independent of x. Define

    Zτ(x)=∫ exp(-Q(x,p)/τ) dp.

The Bregman inequalities give

    (2πmτ)^(d/2) <= Zτ(x) <= (2πLτ)^(d/2).

As τ decreases to zero, ordinary local Laplace expansion (C2 is enough for the leading term) gives

    Zτ(x)/(2πτ)^(d/2) -> sqrt(det A(x)).

The displayed uniform bounds justify dominated convergence against a Gaussian root density φ. Thus the globally normalized law

    Cτ^-1 φ(x) exp(-Q(x,p)/τ) dx dp

has limiting x marginal proportional to φ(x)sqrt(det A(x)), **not φ(x)** in general. For the scalar cosine fixture, fiber normalization at x=0 versus x=π tends to sqrt(.6/.4)=sqrt(3/2). This is already a nonlinear root-law change.

The same determinant follows from graph geometry: graph area is sqrt(det(I+A²))dx, while the determinant of the nonzero normal Hessian is det(A+A^{-1})=det(A)^{-1}det(I+A²). Their combination gives sqrt(det A)dx. A delta constraint in dp coordinates and surface area on the graph are different measures.

### Gaussian ambient-prior counterexample, even with affine g

Let d=1, g(x)=αx with α in [.4,.6]. Then

    Q=(p-αx)²/(2α).

Starting from independent standard Gaussian x,p and globally weighting by exp(-Q/τ) yields

    Varτ(x)=[1+α²/(1+ατ)]^{-1}.

For α=.5, this tends to 4/5 rather than 1. Merely making the penalty infinitely stiff therefore does not restore the original Gaussian root.

### Conditional normalization preserves the x root, but is an additional operation

The normalized conditional joint law is

    μτ(dx,dp)=φ(x) exp(-Q(x,p)/τ)/Zτ(x) dx dp.

Its x marginal is exactly φ. Let bτ(x)=Eμτ[p|x]. Exact differentiation gives

    grad_x log Zτ(x)=(bτ(x)-g(x))/τ,
    grad_x log μτ = -x+(p-bτ(x))/τ,
    grad_p log μτ = -(h(p)-x)/τ.

Equivalently, its negative log potential has x gradient x+(bτ(x)-p)/τ. Omitting grad log Zτ produces the wrong score. The conditional moment bτ or the normalizer is not a supplied finite original-gradient VALUE oracle. It must be constructed and its actual first and approximation errors paid for. The identity

    D² log Zτ(x)=-A(x)/τ + Cov(p|x)/τ²

is exact; no third derivative was used. At positive τ the conditional law has full-dimensional fiber noise, not a graph constraint. Elementary integration by parts and monotonicity of h give

    E[|p-g(x)|² | x] <= Lτd.

Hence μτ converges weakly to φ(x)dx δ_{g(x)}(dp); that limit is singular and is not an ordinary ambient Gaussian input law.

Even in the affine case where Zτ is constant and its score is free, the normalized negative-log joint Hessian contains

    [[1+α/τ,-1/τ],[-1/τ,1/(ατ)]],

whose largest eigenvalue grows as (α+α^{-1})/τ. A concentrating density does not come with uniformly small first at the fixed original source scale.

## 4. Literal K3 consequence: positive temperature changes the query record

The native source uses p_i=g_0(x_i), x_1=x_0+aM*Delta and the same Gaussian roots S,U,Z. Applying an unnormalized Fenchel penalty to p_i changes the root law by its fiber factors. Applying normalized conditional penalties retains those roots but randomizes p_i, and therefore changes the actual outer arguments S+aMp_i. At τ=0 the algebraic graph returns, together with its singular-law current.

A particularly sharp test uses g_0=g_t=x/2, M=1, and epsilon=0. The native record has x_1=x_0, p_1=p_0 and E=0 identically. Independent normalized Fenchel fibers instead have

    p_i | x_0 ~ N(x_0/2, τ/2), independently,
    Eτ = (ra/2)(p_0-p_1),
    E[Eτ² | roots]=r²a²τ/4 > 0.

Thus a positive-temperature product penalty creates actual energy where the literal K3 factor vanishes. Specially tying the auxiliary noises is an extra coupling construction, not a consequence of positivity or conditional normalization; at nonlinear arguments it still needs a same-record proof. This counterexample does not prohibit such a separately proved construction.

Adding a FULL extra exact κQ(x_i,p_i) to an ambient potential has zero gradient on its exact constraint graph, so that particular addition does not remove existing graph companions or make the graph Gaussian. This is distinct from the main artifact's useful DUAL COMPLETION: adding only the signed V0*(p0)-V0*(p1) terms to the prior source's already-existing V0(x_i)-x_i·p_i terms really does remove the -x0/+x1 p companions and the unmarked rM*d common-sum term. The resulting exact common sum is aM*E; the opposite baseline remains. Other K3 constraints, such as d=a[gt(S)-gt(S+epsilon Z)], also remain to be enforced. The independently reviewed final main artifact is pinned below.

## 5. Exact graph current: where the Hessian rows return

For the intended graph law μ(dx,dp)=φ(x)dx δ_{g(x)}(dp), a smooth test F satisfies the exact Gaussian-root integration-by-parts identity

    Eμ[∂_{x_i}F + sum_j A_{ji}(x) ∂_{p_j}F] = Eμ[x_i F].

This follows directly by applying Gaussian integration by parts to F(x,g(x)). The tangent current has frame (e_i,A(x)e_i). Its A(x) row is the actual retained Hessian correlation, not a known constant Gaussian row. Dropping it already fails for F=p_j. With g(x)=x/2 in dimension one, the falsely omitted-current identity would read 0=E[xp]=1/2.

The same statement for the entire feature map Gamma is E[D Gamma(Y)^* grad F(Gamma(Y))]=E[Y F(Gamma(Y))]. It contains exactly the graph derivative rows identified in source §5. Taking full Euclidean gradient pullbacks may reintroduce HVP-valued companions and higher-jet firsts; the weak current identity itself does not differentiate an HVP. A new retained-current theorem could be a useful descendant, but the stationary Gaussian theorem cannot be invoked without that operation.

## Scope and returned components

Returned: exact bounded-first Fenchel source; explicit inverse VALUE counts and actual first bounds; an exact smooth finite-iterate curl fixture; a finite-iterate FIRST-error fixture; exact conditional scores and normalization bias; affine Gaussian-root and zero-native-energy counterexamples; and the exact graph-current identity.

Not returned: a finite exact genuine-gradient inverse source for general g; a free conditional-normalization or exact conditional-sampling oracle; a Gaussian-compatible compiler for the singular nonlinear graph; or a retained one-energy bound for the extra channels. The audit excludes this naive Fenchel/positive-penalty escape as a literal native K3 realization, not all possible repairs.

## Final main-artifact review pin

The complementary dual-completed construction in `FENCHEL-INVERSE-LIFT-AND-LITERAL-K3-GRAPH-CURRENT.md`, SHA-256 `7357f4e314c2d8024baea2a6cb1fa8021762728a6456a951f1831e7c5b78cbe4`, was independently reviewed. Its exact marked common mode aM*E, opposite baseline, literal product-fiber W score, and common-noise displacement identities are correct. The unnormalized common-mode sqrt(2) readout factor and split finite inverse tolerance budget were corrected and verified in this pinned version. The sharper inverse-defect bound and its no-product counterexample appear in the accompanying `INDEPENDENT-REVIEW-OF-DUAL-COMPLETED-AMBIENT-LIFT.md`. Exact inverse reference plus coherent finite VALUE restoration is a legitimate component when the receiving theorem supports that restoration; the present audit does not deny this. It does not supply the still-needed constrained-law current or all-rank closure.
