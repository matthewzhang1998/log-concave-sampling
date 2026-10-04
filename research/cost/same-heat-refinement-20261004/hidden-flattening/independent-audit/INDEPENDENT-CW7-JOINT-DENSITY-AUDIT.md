# Independent CW7 joint-density and hidden-flattening audit

2026-10-04. Bounded mathematical audit of the actual named CW7 hidden provider and its ideal observation decoder. This is not an impossibility theorem for all hidden-mean methods.

## Result

The displayed CW7 construction does not presently furnish a known-gradient joint target to which the new known-center constructor can be applied once. There are three distinct objects:

1. The **actual finite CW7 graph** is a pushforward of its complete Gaussian tapes through original-gradient programs. Flattening those tapes preserves every source evaluation and ancestor; it does not create a cheap joint potential of original V evaluations.
2. The **ideal decoder reference** has an explicit conditional-posterior path density, but every conditional normalizer remains. Its gradient requires posterior means, while its observation mean q=E[H|theta] is itself unavailable.
3. A **Bayesian undirected joint** obtained by dropping those normalizers is a different distribution, generally tilting the ancestors. The available Gibbs identity cancels a child normalizer when the predecessor has already supplied the exact stationary marginal. CW7 starts with a finite mode, not that marginal.

There is a useful positive structural fact: with q known and a deterministic initial state, the ideal decoder's path potential is strongly convex at sufficiently small heat. Convexity is therefore not the principal obstruction for this affine-center path. Oracle closure is. Conversely, arbitrary nonlinear hidden-center edges need not preserve joint convexity.

A promised quadratic potential together with an explicitly affine hidden graph is a genuine positive special case. Its means and normalizers close by linear algebra. Dimension, history size, matrix actions, and actual-tape certificates must still be charged.

## 1. Sources and conventions

Read sources:

- `COMPONENTWISE-OLD-R7-HIDDEN-MEAN.md`, especially sections 1, 2, 4, 5.
- `exact-slack/cw7/01_retained_hidden_force_certificate.md`, sections 2--5.
- `no-copy-rank-20261004/host-observer-audit/ACTUAL-CW7-WEAK-E-HOST-SEPARATOR-AND-COMPARISON-ORDER.md`, especially sections 2, 3, 5--7.
- LOW30: `LOW30 / 30_low_acc.tex (external source; not bundled)`, sufficient-statistic decoder at lines 12676--12734 and finite modes at lines 13695 onward.
- The input signature of `cost/same-heat-refinement-20261004/protected-short-refresh/DETERMINISTIC-KNOWN-CENTER-SOURCE-WITH-PROTECTED-REFRESH.md`.

Condition throughout on the actual captured caller theta. This includes previously exposed physical labels and the once-refreshed phase momentum, but not freshly sampled histories. Write q=E[H|theta]. Let phi_b be the normalized d-dimensional N(0,bI) density. Define

    z_b(c) = integral exp(-V(x)) phi_b(x-c) dx,
    Q_(b,c)(dx) = exp(-V(x)) phi_b(x-c) dx / z_b(c).

If instead Z_b(c)=integral exp(-V(x)-|x-c|^2/(2b)) dx, then Z_b=(2 pi b)^(d/2) z_b. All center derivatives below are identical for Z and z. Normalized phi is used to avoid hiding constants in Bayes identities.

The original VALUE oracle is F=grad V, not a potential-value oracle. The supplied first/adjoint actions are HVPs under C2 regularity. Writing V in a mathematical density does not grant potential values, log normalizers, their gradients, arbitrary third derivatives, or a new joint-potential oracle.

## 2. Exact actual finite law before any posterior idealization

Put n=k+1, sigma^2=a/n and

    S(theta,Omega) = sigma [h0(theta) + sum_l sqrt(v_l) Y_l(theta,Omega_l)],
    T = S + sqrt(3/4) sigma Z,
    R_i = T + sqrt(a)(E_i-Ebar).

All Y_l are complete actual local-statistic programs with their original internal aliases. Their banks are independent conditional on theta. Z and the observation pool are new. The exact statistic density is

    p_T(t|theta) = E_Omega phi_(3 sigma^2/4)(t-S(theta,Omega)).             (A)

Let rbar=n^(-1) sum_i r_i. The exact observation density with respect to Lebesgue measure on R^(nd) is

    p_R(r|theta)
      = n^(-d/2) p_T(rbar|theta)
        (2 pi a)^(-k d/2)
        exp[-sum_i |r_i-rbar|^2/(2a)].                                  (B)

The Jacobian factor n^(-d/2) follows by using the orthonormal common coordinate sqrt(n) rbar. If, and only at this analytical comparison boundary, p_T is replaced by phi_(a/n)(t-q), (B) becomes product_i phi_a(r_i-q). The actual finite p_T is not being asserted Gaussian.

Keeping Omega explicit instead gives a perfectly known Gaussian conditional R|Omega: its mean is (S,...,S), and its block covariance is

    a I_n - (a/(4n)) 11^T,

tensored with I_d. Its common-mode eigenvalue is 3a/4 and its other eigenvalues are a. However, the mean S remains the entire nonlinear source program. No original query has disappeared.

Let M_J be the literal finite mode and K_i^fin the actual known-center child output kernel. The exact reduced actual law is

    p_R(r|theta) dr delta_(M_J(r0))(dx0)
      product_(i=1)^k K_i^fin((x_(i-1)+r_i)/2,dx_i).                     (C)

These finite children are not asserted to possess densities equal to Q. One can always write the complete program as independent standard-Gaussian tape measure with delta constraints at every deterministic node. That is a faithful tape representation, not a new cheap full-dimensional log-concave state density. Sampling its tapes still requires evaluating S, every history, the finite mode, and every child to get the output.

The sufficient-output separator justifies comparing (A) first and appending the host. It does not identify (A) with a Gaussian law, or permit replacing actual tapes in the independent retained/first/protected proof.

## 3. Exact ideal decoder density, including every partition factor

Make both idealizations explicit: take R_i iid N(q,aI), and replace each child by its exact Q law. With c_i=(x_(i-1)+r_i)/2 and b=a/2, the path law is

    [product_(i=0)^k phi_a(r_i-q) dr_i]
    delta_(M_J(r0))(dx0)
    product_(i=1)^k {
       exp[-V(x_i)] phi_b(x_i-c_i) dx_i / z_b(c_i)
    }.                                                               (D)

One may replace M_J by the exact prox for a stronger, simpler reference, paying the actual numerical error separately. This replacement does not remove any z_b(c_i).

The more general ideal ancestral posterior DAG has the same issue. For unique stochastic nodes v with parents pa(v), root density p_root and centers c_v,

    p(z) = p_root(z_root)
           product_v exp[-V(z_v)] phi_(b_v)(z_v-c_v(z_pa))
                     / z_(b_v)(c_v(z_pa)).                            (E)

Every occurrence carries its own parent-dependent denominator. Exact aliases are one shared variable, not independent factors. Removing the denominators multiplies the original joint density by product_v z_(b_v)(c_v), then renormalizes; this generally changes ancestral marginals. Equation (E) is an ideal posterior DAG, not a claim that every actual CW7 finite source already has this density.

## 4. The exact Bayes cancellation and its required predecessor law

Set y=z+r-q. The fresh observation r~N(q,aI) makes y|z~N(z,aI), and the decoder center is (q+y)/2. Gaussian multiplication gives

    phi_a(x-q) phi_a(y-x)
       = phi_(2a)(y-q) phi_(a/2)(x-(q+y)/2).                            (F)

If Z~Q_(a,q) and Y|Z~N(Z,aI), then

    p_Y(y) = phi_(2a)(y-q) z_(a/2)((q+y)/2) / z_a(q),                  (G)
    Law(Z|Y=y) = Q_(a/2,(q+y)/2).

Now (G) cancels the child denominator in (D), and one Gibbs update returns exactly Q_(a,q). This is the stationarity proof and yields detailed balance. It is not an identity replacing the whole finite decoder by one original-V joint without paying for stationarity.

For a nonstationary incoming distribution nu,

    p_Y^nu(y) = integral phi_a(y-z) nu(dz),

which is generally not (G). Appending its exact conditional Q child leaves the ratio p_Y^nu(y)/z_(a/2)((q+y)/2). Thus iid observations by themselves do not provide the cancellation. The CW7 finite-mode initialization is the concrete nonstationary boundary. A different construction might deliberately supply (G); that would require proving how to generate this V-dependent tilted marginal without already importing the target sampler/mean service.

One can also telescope likelihood ratios against a stationary path by retaining the initial density ratio. Its conditional expectation at the endpoint is precisely a nonstationarity correction; it is not identically one. No endpoint cancellation follows from merely renaming this ratio.

## 5. Scores, joint convexity, and oracle closure

Let m_b(c)=E_(Q_(b,c)) X and Sigma_b(c)=Cov_(Q_(b,c)) X. Gaussian differentiation and integration by parts give

    grad log z_b(c) = (m_b(c)-c)/b = -E_(Q_(b,c)) F(X),
    Hess log z_b(c) = Sigma_b(c)/b^2 - I/b.                            (H)

If 0<=Hess V<=L I, then

    -L/(1+bL) I <= Hess log z_b(c) <= 0.                               (I)

For the upper bound use the posterior covariance upper bound; the lower bound follows, for example, from Gaussian differentiation in the representation z_b(c)=E exp[-V(c+sqrt(b)G)]: Hess log z=-E Hess V+Cov(F), combined with the sharp posterior covariance lower bound. The weaker -L I lower bound, obtained directly from that representation, is already sufficient below.

The negative log conditional has center derivative

    grad_c[-log Q_(b,c)(x)] = (m_b(c)-x)/b.                            (J)

Consequently the ideal decoder path energy, with x0 fixed, is

    U = sum_i [ |r_i-q|^2/(2a) + V(x_i)
                + |x_i-(x_(i-1)+r_i)/2|^2/a
                + log z_(a/2)((x_(i-1)+r_i)/2) ].                     (K)

For example, its r_i gradient is

    (r_i-q)/a + (m_(a/2)(c_i)-x_i)/a.

An interior x_i gradient contains

    F(x_i) + (x_i-c_i)/b + (m_b(c_(i+1))-x_(i+1))/(2b).

Therefore evaluating grad U requires one posterior mean at each current child center, in addition to original gradients. Differentiating it requires posterior covariance actions. Formula (H) is an exact analytical identity, not an available original-F/HVP implementation. Applying the new known-center constructor to U while treating these terms as primitive would reintroduce the hidden service.

### Positive convexity fact for the affine decoder path

For fixed x0, write S for the backward shift on x_1:k, and u=x-(Sx+r)/2. The Gaussian part of (K) has quadratic second variation

    (|r|^2+2|u|^2)/a >= (|r|^2+|x|^2)/(4a).

Indeed, |x|<=2|u|+|r|. The center map c=(Sx+r)/2 satisfies |dc|^2<=(|dx|^2+|dr|^2)/2. Equations (I) and convexity of V imply

    Hess U >= (1/(4a)-L/2) I.                                         (L)

This bound is uniform in k. A known quadratic ridge can likewise be moved between the Gaussian part and the residual potential to make the latter convex with bounded Hessian, for small aL. Thus a formal Gaussian-plus-convex joint model exists if q and the posterior-score oracle are supplied. It has dimension 2kd before an initial random state. This is a structural positive result, not an executable original-gradient reduction.

### Why general nonlinear histories are different

For a scalar Gaussian ancestor h and one conditional Q_(b,c(h)), the h-h second derivative of its negative log joint, whenever c is twice differentiable, is

    1 + c'(h)^2 Sigma_b(c(h))/b^2
      + c''(h)(m_b(c(h))-x)/b.                                       (M)

If c''(h) is nonzero, the final term is unbounded below in one x direction. Smooth convex potentials with bounded Hessian already supply literal nonlinear force centers of this form: take V(t)=lambda t^2/2+epsilon(1-cos t), 0<epsilon<lambda and lambda+epsilon<=1, and c(h)=u-tau F(sqrt(a)h). At suitable h, c'' is nonzero. This one-edge density is not globally convex even with arbitrarily small positive tau. Any bounded-guard alternative needs its own proof; (M) is not a claim about an unexamined complete guarded CW7 joint.

Under only C2, a center built from F is generally only C1. Its first action can use HVPs, but differentiating that action again is not an available third-derivative operation. This matters if a proposed joint-density constructor asks for bounded Lipschitz gradients or Hessian actions of a new composed potential.

### Mode pushforward is not a free remedy

For the exact mode X0=prox_(aV)(R0), R0=X0+aF(X0), so its density is

    p_0(x) = phi_a(x+aF(x)-q) det(I+a Hess V(x)).                       (N)

The gradient of the logarithmic determinant requires third-derivative information when it exists; C2 alone does not make this density a C2 potential. For the finite mode, the corresponding Jacobian is built from the actual finite HVP chain and its derivative has the same regularity problem. Keeping R0 and the deterministic constraint avoids these derivatives but produces a singular constrained joint. Neither operation is automatically the known-center constructor's admitted interface.

Potential-value-based rejection, weighting, or normalizer-ratio methods would be new algorithms with new oracle and query analysis. The original VALUE notation does not license them for free. This does not rule out implementing potential differences from F with separately certified quadrature.

## 6. Scalar quadratic counterfixtures, valid at arbitrarily small heat

Let V(x)=lambda x^2/2, lambda>0, and put

    beta=1/(1+a lambda), alpha=1/(2+a lambda).

Then Q_(a,c)=N(beta c,a beta). For H~N(q,tau^2), three laws differ:

    Q_(a,E H)                 = N(beta q, a beta),
    E Q_(a,H)                 = N(beta q, a beta+beta^2 tau^2),
    dropped-normalizer joint  has X~Q_(a+tau^2,q).

The last equality comes from integrating exp(-V(x)) phi_a(x-H) against the original Gaussian H density before globally normalizing. It tilts H by z_a(H). With a=lambda=tau^2=1 and q=0 the three variances are respectively 1/2, 3/4, and 2/3. The discrepancy persists for every a>0, including the admitted small-heat regime; these values are merely readable numbers.

For the ideal decoder with exact-mode initialization and iid R_i~N(q,a),

    X0=beta R0,
    Xi=alpha X_(i-1)+alpha R_i+sqrt(a alpha) G_i.

Its mean is beta q at every index, but

    Var X_k = a beta - a^2 lambda beta^2 alpha^(2k).                    (O)

Thus no finite k gives stationarity for lambda>0. At a=lambda=1, k=1 gives 17/36 rather than 1/2. At a=1/10, lambda=1, k=1 gives 4751/53361 rather than 1/11. Finite CW7 modes and finite children add their separately paid errors; they do not justify replacing this residual by an exact identity.

A one-step normalizer-deletion fixture is even more direct. Fix x0=0, R~N(0,a), X|R~Q_(a/2,R/2). The true variance is a alpha+a alpha^2. Deleting z_(a/2)(R/2) and renormalizing gives

    Var R = 2a(2+a lambda)/(4+3a lambda),
    Var X = 3a/(4+3a lambda).

For a=lambda=1 these are 6/7 and 3/7; the original R variance is 1 and the original X variance is 4/9. This exhibits the ancestor tilt rather than merely an endpoint error.

## 7. Genuine positive quadratic/affine cases and their costs

Suppose it is promised that

    V(x)=x^T A x/2+h^T x+constant, 0<=A<=L I,

and that the relevant hidden graph is explicitly affine in its Gaussian tapes (or is an affine-center exact Gaussian conditional DAG whose means can be propagated). Put B_b=(I+bA)^(-1). Then

    Q_(b,c)=N(B_b(c-bh), b B_b),
    z_b(c)=det(I+bA)^(-1/2)
       exp[-c^T A B_b c/2-h^T B_b c+(b/2)h^T B_b h-constant].           (P)

All conditional partition factors are explicit quadratic functions. For an explicitly affine H=h0+M Omega, q=h0 is its zero-tape value. For an affine-center Gaussian DAG, means propagate once through its topological order. This licenses direct sampling from Q_(a,q), after computing q, without Gaussianizing H. It does not assert that every actual CW7 correction program is affine merely because V is quadratic; that reduction must be checked on the literal program.

For the ideal quadratic decoder define D=(2I+aA)^(-1), B=(I+aA)^(-1). Then

    Xi=D X_(i-1)+D R_i-aD h+sqrt(aD) G_i,
    E Xi=B(q-ah),
    Cov X_k=aB+D^k(aB^2-aB)D^k.

The finite endpoint itself can be generated as a d-dimensional Gaussian if only its law is required. The full ideal path is Gaussian after eliminating the affine mode, with dimension (2k+1)d. Its original triangular generation uses O(k) block operations; dense one-shot factorization can be more expensive. Changing to a law-equivalent sampler does not inherit the CW7 same-tape/protected certificate.

### Matrix-free fixed-order version: no dense reconstruction required

The d-HVP reconstruction below is optional. Under the explicit constant-Hessian promise, let ||bA||<=rho<1 and define

    P_m = sum_(j=0)^m (-bA)^j,
    S_m = sum_(j=0)^m (-1)^j [binom(2j,j)/4^j] (bA)^j.

These approximate (I+bA)^(-1) and (I+bA)^(-1/2). An executable conditional output is

    X_m = c - b P_m F(c) + sqrt(b) S_m G.

It uses one original F(c) call and at most m HVPs on each of the two vector power chains, at any fixed known reference since A is constant. No matrix is formed. Its errors satisfy

    |mean(X_m)-B_b(c-bh)| <= b (bL)^(m+1) |F(c)|,
    ||sqrt(b)(S_m-B_b^(1/2))G||_L2
       <= sqrt(bd) (bL)^(m+1)/(1-bL).                                 (Q)

The first bound follows from P_m-(I+bA)^(-1)=-(-bA)^(m+1)(I+bA)^(-1). The second uses binom(2j,j)/4^j<=1. The anchored mean expression avoids paying a spurious factor |c| when F(c) is small. It does not remove the genuine weighted caller-profile factor |F(c)|.

For b comparable to, or smaller than, the outer small heat a, a fixed target a-power needs degree O_R(1), together with the original public precision/logarithmic allocations and the declared caller-moment weights. More generally degree O(log(1/epsilon)) works uniformly for bL<=rho<1. This is a restricted zero-inverse-heat-query-exponent construction, not a general CW7 flattening theorem.

Apply these actions at each affine conditional node, and compute node means in topological order using the same finite polynomial versions. Independent centered Gaussian inputs are set to zero only after the whole relevant program has been certified affine; otherwise a zero tape is not an expectation. If a mean error at node v has root readout W_v, the root error budget is bounded by sum_v ||W_v|| delta_v. Use the actual finite coefficient/depth bounds to allocate delta_v, rather than assuming all histories contract.

For N distinct affine nodes and maximal chosen degree m, the elementary bound is O(N m) HVPs plus O(N) original gradient calls and the literal known affine operations. Straightforward retained storage is O(Nd), with each node's vector power sequence streamed using O(d) temporary storage; a specified reverse/caller schedule can have additional storage or paid replay. No dense d-by-d matrix reconstruction is necessary. Assembly/evaluation of unknown affine coefficients, complete mean/origin traversals, terminal force queries, and numerical source versions remain part of the bill. For independent calls, one cannot reuse a sampled ancestor unless exact equality of the needed deterministic affine summary has been established.

Accounting:

- If A is not explicitly given, a promised general quadratic is identified by one F(0) query and d HVPs on a basis, with O(d^2) storage and conventional O(d^3) dense factorization. This is dimension-dependent work, not a dimension-free polylogarithmic original-query theorem.
- Matrix-free inverse/square-root approximations can instead use polynomial approximation and paid HVP actions; bounded aL gives good conditioning, but the accuracy degree and every action must be counted.
- If A=lambda I is promised, one HVP identifies lambda and F(0) identifies h. In this special case covariance actions are scalar, and a direct output costs O(d) arithmetic/Gaussian generation after the hidden mean is known.
- An affine history DAG with N genuinely distinct nodes costs at least its actual required node/matrix evaluations to propagate a mean. Topological depth can be parallelized where allowed; it cannot silently erase N, coefficient assembly, or private zero-tape computations. A recursion that duplicates full ancestors must either pay those occurrences or prove exact reusable affine summaries.
- Caller derivatives of q and source-version changes, if consumed downstream, require their actual affine coefficient/caller computation. A law-only d-dimensional collapse gives no automatic retained-label or protected-row contract.

## 8. Admission conclusion

PASS: all conditional normalizers and the exact stationary Bayes cancellation are identified.

PASS: the actual statistic/observation law can be written explicitly as (A)--(C), preserving the shared T. A Gaussian-tape flattening is faithful but does not reduce source cost.

PASS: the fixed-initial-state ideal affine decoder is jointly strongly convex at small heat; a nonconvexity claim against that specific reference would be incorrect.

NOT SUPPLIED: an original-F/HVP implementation of the path score in (H)--(K), an accessible reference q, or a cancellation of every normalizer from the actual finite-mode/nonstationary chronology. These are concrete missing interfaces, not a universal lower bound.

POSITIVE SPECIAL CASE: promised affine hidden graphs and quadratic V admit explicit Gaussian closure, with the dimension/history/query costs in section 7.

The attached checker verifies Gaussian identities, covariance distinctions, normalizer-induced tilt, decoder residuals, matrix quadratic closure, and the path precision bound. These diagnostics do not certify an unprovided hidden sampler or its retained source graph.

Publication provenance: local absolute source location replaced by its named external-source identifier. Mathematical content and historical source hashes are unchanged.
