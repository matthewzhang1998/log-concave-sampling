# Exact CW7 hidden law and the joint-flattening gate

2026-10-04. Bounded test of flattening the selected hidden histories. No universal impossibility claim is made.

## Result

The selected CW7 hidden target is Q_(a,q), where q=E[H|theta] is an inaccessible deterministic mean conditional on the entering caller. It is not the mixture E[Q_(a,H)|theta]. Replacing that mean by an ordinary latent variable changes the target already for a quadratic potential.

The exact ideal-child decoder path does have a joint density, but it contains one conditional posterior partition factor for every child. The Gaussian sufficient-statistic observation identity does not cancel these factors. Even in an affine-initialized case where Gaussian whitening gives a well-conditioned strongly log-concave joint potential, its gradient contains posterior means through grad log Z. That is not an original-gradient/HVP VALUE oracle.

Thus this direct flattening does not establish c(R,depth R)=O(R). It either changes the target, retains the unavailable posterior-mean score in the joint gradient, or merely re-encodes the original finite Gaussian program and its full existing query count. A separately certified affine/quadratic case is a genuine zero-heat-power subproblem, described below.

## 1. The actual inaccessible center

At the captured caller theta, the selected componentwise program has

    H(theta;Omega)=h_base(theta)+sum_l C_l g_l(theta;Omega_l),
    q(theta)=E[H(theta;Omega)|theta].                   (1)

Each g_l is a complete original force/history occurrence. The conditional independence and all aliases are exactly those in `COMPONENTWISE-OLD-R7-HIDDEN-MEAN.md`; an ancestor does not become a caller cache merely because its evaluation is expensive.

The target of the hidden sampler is

    Q_(a,q)(dx)=exp(-V(x)) phi_a(x-q) dx / Z_a(q),
    Z_b(c)=integral exp(-V(x)) phi_b(x-c) dx.          (2)

Here phi_b is the normalized N(0,bI) density. Potential values in these formulas are analytical quantities; the original VALUE oracle returns grad V, not an evaluation of Z or an expectation.

Conditional on theta, q is fixed. A valid joint coupling may attach any auxiliary variables to the target, but its x marginal must still equal (2). The simplest product extension is p_Omega(Omega|theta) Q_(a,q)(x). It still contains the unknown q. The natural proposed replacement

    p_Omega(Omega|theta) Q_(a,H(Omega))(x)             (3)

is normalized but has a different x marginal. There is no Jensen identity equating (2) and (3).

### Exact quadratic test

Let V(x)=lambda x^2/2 and H~N(mu,tau^2). Then

    Q_(a,E H)=N(mu/(1+a lambda), a/(1+a lambda)),
    integral Q_(a,H) dLaw(H)
       =N(mu/(1+a lambda),
           a/(1+a lambda)+tau^2/(1+a lambda)^2).       (4)

Taking tau=a^(3/2) is compatible with the ordinary small center-fluctuation scale. The discrepancy in (4) is of order a^(5/2), not arbitrary order. Averaging n independent latent H samples divides that extra variance by n; to obtain order a^R by this particular route needs n of order a^[-(R-5/2)]. Such replicas are real full history evaluations.

This screens the latent-mean substitution, not all possible joint representations.

## 2. The actual sufficient-statistic observation identity

The finite statistic T is constructed so its conditional law approximates

    T_ref~N(q,aI/(k+1)).

With independent E_0,...,E_k, the observations are

    R_i=T_ref+sqrt(a)(E_i-E_bar).

The exact covariance identity makes the entire reference array independent:

    (R_0,...,R_k)~product_i N(q,aI).                  (5)

The one common T is essential. Replacing T_ref by H+the same Gaussian does not produce (5): every cross covariance acquires Cov(H|theta). The actual finite statistic is compared in complete-output law before this host; its private source variables do not appear as additional observations in (5).

Identity (5) removes the Gaussian correlation introduced by the common statistic. It does not replace q by a sampled H, and it does not evaluate the mean in (1).

## 3. Exact decoder-path density

Let M_a be the actual finite mode/prox map used at the first observation. The ideal-posterior-child reference sets

    x_0=M_a(r_0),
    b=a/2,  c_i=(x_(i-1)+r_i)/2,
    x_i | (x_(i-1),r_i) ~ Q_(b,c_i),  1<=i<=k.

Conditional on theta, its exact density with respect to dr_0...dr_k dx_1...dx_k is

    product_(i=0)^k phi_a(r_i-q)
       times product_(i=1)^k
           [exp(-V(x_i)) phi_b(x_i-c_i)/Z_b(c_i)].    (6)

If x_0 is included as a separate coordinate, add the deterministic constraint delta_(M_a(r_0))(dx_0). The actual finite child replaces its Q factor by its actual pushforward kernel; it is not assigned the analytic Q density for free. Its law comparison to (6) has the already charged child errors.

The negative log of (6), apart from constants and the initial constraint, is

    sum_i |r_i-q|^2/(2a)
      +sum_(i>=1) [V(x_i)+|x_i-c_i|^2/(2b)+log Z_b(c_i)]. (7)

All c_i use their actual parent values. No factor in the independent observation prior is a matching Z_b(c_i). Dropping the last terms multiplies the correct path law by product_i Z_b(c_i), followed by a new normalization. That is a tilt of the ancestor/observation law, not the same joint model.

### The exact Bayes cancellation and why it is unavailable here

Suppose, additionally, x_(i-1) already has the exact stationary law Q_(a,q). Put Y=x_(i-1)+r_i-q. Its density is

    phi_(2a)(Y-q) Z_(a/2)((q+Y)/2) / Z_a(q).          (8)

The Z factor in (8) cancels the corresponding conditional-child denominator, recovering the simple Gibbs joint law. This cancellation uses stationarity of the entering x_(i-1), not only the independent Gaussian law of r_i. The finite decoder begins from M_a(R_0) and approaches stationarity by its geometric contraction; (8) is not its initial input law. Assuming (8) at the outset assumes the hidden posterior already exists as a sampler.

Even where stationarity supplies (8), q is still the inaccessible mean from (1). The cancellation therefore does not turn the selected hidden service into a known-center original-gradient problem.

## 4. Quadratic check that dropping Z changes a real decoder

Take q=0, V(x)=lambda x^2/2, and a literal finite linear mode iteration

    M_a(r_0)=tau_m r_0,
    tau_m=sum_(j=0)^m (-a lambda)^j.

For one ideal child, c=(tau_m R_0+R_1)/2 is Gaussian with

    v_c=a(tau_m^2+1)/4,
    t_b=1/(1+b lambda), b=a/2.

The correctly normalized path has endpoint variance

    v_true=b t_b+t_b^2 v_c.

After deleting Z_b(c) and renormalizing the purported joint density, c is tilted by Z_b(c), so its variance becomes v_c/(1+lambda t_b v_c). The endpoint variance is instead

    v_drop=b t_b+t_b^2 v_c/(1+lambda t_b v_c).         (9)

These are unequal, even though the potential and every conditional are Gaussian. At small a their W2 discrepancy is of order a^(3/2). The factor is not an irrelevant scalar normalization.

For the simpler latent-mean model H~N(mu,tau^2), omitting Z_a(H) from (3) and integrating H gives Q_(a+tau^2,mu). This is a third law, generally different from both laws in (4).

## 5. Why a favorable whitened joint still lacks the VALUE oracle

Even grant that q is known and replace the initial mode by a fixed known point so every c_i is affine. The V=0 path in (6) is a known Gaussian linear system:

    R_i=q+sqrt(a)G_i,
    X_i=(X_(i-1)+R_i)/2+sqrt(a/2)Z_i.

Its triangular Gaussian map has uniformly bounded normalized condition numbers in k, because the inherited coefficient is 1/2. It has O(kd) coordinates. After whitening, the non-Gaussian potential is a sum of original V terms and log Z_b terms. At small a its Hessian can remain bounded above and below by fixed positive constants. A blanket nonconvexity objection would therefore be incorrect.

The missing oracle is explicit:

    grad log Z_b(c)=(E_(Q_b,c) X-c)/b
                  =-E_(Q_b,c) grad V(X),             (10)
    Hess log Z_b(c)=-I/b+Cov_(Q_b,c)(X)/b^2.

For one parent-to-child center map c=F(parent), its contribution to the joint score simplifies to

    DF(parent)* [E_(Q_b,c)X-x_child]/b.               (11)

The simplification does not remove the conditional posterior mean. A joint gradient evaluation would require these means at all actual centers, and a joint first/HVP would require their covariance/response actions. Treating (10) as an original gradient of a new known potential silently introduces the service being eliminated.

The favorable O(kd) latent dimension and stable Gaussian whitening do not fix this oracle problem. Calling finite posterior-mean programs inside every such gradient is a real recursive bill; it is not one application of the known-center constructor with only original V queries.

## 6. Actual nonlinear histories add a second source-grammar issue

The selected hidden centers are weighted original-force histories, rather than arbitrary known affine functions of latent posterior states. If a proposed joint factor uses

    |x_child-h_base-sum_l C_l grad V(x_l)|^2/(2b),

its gradient with respect to x_l already contains an original Hessian action. Its first/HVP generally differentiates that action and requires an original third derivative. C2 does not supply it. Adding force variables with exact constraints f_l=grad V(x_l) gives a singular measure instead of the smooth Gaussian-plus-potential model. A finite penalty relaxation has a separate bias, stiffness, and conditioning price; it is not an exact flattening.

The initial mode has the same distinction. For an exact prox and R_0~N(q,aI), its pushforward density is

    phi_a(x+a grad V(x)-q) det(I+a Hess V(x)).         (12)

Differentiating the determinant term generally needs information beyond original C2 VALUE/HVP access. Keeping R_0 and its deterministic mode constraint is exact, but does not make that constraint a known quadratic term. Actual finite mode maps have their own pushforward Jacobian, not the density exp(-V(x)) times a Gaussian.

These are additional difficulties. The normalization gate already remains in the more favorable affine-initialized case of Section 5.

## 7. Primitive-Gaussian flattening does not erase executed work

The actual finite CW7/weak-E program is already a deterministic VALUE graph on a finite collection of primitive Gaussians. One may write its joint law as the standard Gaussian density times deterministic constraints for all recorded queries and outputs. Sampling its Gaussian inputs and evaluating those constraints is exactly executing the original graph.

This representation preserves all partition effects implicitly, but it preserves all original queries too. An empirical bank with N complete private histories still has N such records. A changed source point still needs its changed original query. Repeated incoming theta labels and identical zero queries can be cached only according to the established exact-key rule; no hidden sample becomes a deterministic caller because it is placed in a joint-variable list.

If the auxiliary constraints are instead relaxed to a smooth density, the resulting model and its conditioning are new and need a complete law/source/cost proof. No such relaxation is supplied by the Gaussian sufficient-statistic identity.

## 8. A genuine restricted zero-power subproblem

Under an explicit **quadratic-potential and certified affine-history** promise, the normalization issue closes algebraically. Write

    V(x)=x* A x/2+g* x+constant,
    B_b=(I+bA)^(-1).

Then

    Q_(b,c)=N(B_b(c-bg), b B_b),
    log Z_b(c)=constant_b
      -c* A B_b c/2-g*B_b c+(b/2)g*B_b g.            (13)

Every conditional partition correction is an explicit quadratic. Affine history means close by a deterministic affine recursion; setting private centered Gaussians to zero computes a mean only when that whole supplied history map is actually certified affine. Quadratic V alone does not prove that every arbitrary literal old correction tape is affine.

A dense reconstruction uses d original HVP directions to recover A and then real matrix factorization work. A matrix-free alternative is often preferable: for b||A||<1, evaluate B_b and B_b^(1/2) on a requested vector by finite Neumann/binomial polynomials. Degree fixed by the desired source/mean/noise tolerances costs that many original HVPs per node, with no inverse-b power at fixed target. The original F(c)=A c+g form gives the mean c-b B_b F(c), making its actual caller envelope explicit. Every affine-node mean/origin traversal, Gaussian coordinate and matrix-vector action is charged.

For a decoder of k nodes, the latent dimension is O(kd) and the number of such matrix-function actions is O(k), times the fixed-target/precision degree. For a larger promised affine history DAG, replace k by its actual number of nodes and retain its known conditioning. If that node count is polynomial-logarithmic at fixed target/depth, this is a zero-heat-exponent subproblem. If the literal supplied DAG already has inverse-heat empirical banks, their count does not vanish without a separately proved affine simplification.

This special case does not remove the nonlocal factors (10) for general C2 potentials and nonlinear hidden histories.

## 9. The concrete identity a successful general flattening would need

A successful route must exhibit an executable joint density or transport such that:

- its terminal marginal is Q_(a,E[H|theta]) at the actual selected hidden target, rather than a mixture centered at H;
- every conditional Z factor is present or canceled by an explicitly supplied matching tilted marginal, without assuming the desired hidden stationarity first;
- its known-gradient and first/adjoint actions reduce to original grad V and original HVPs, with no posterior-mean, covariance, determinant-gradient, or original-third-derivative oracle;
- its Gaussian whitening, smallest variance, condition number, latent dimension and complete query count are explicit when hidden depth grows with R;
- its actual retained/caller/zero/protected source graph is returned independently of its marginal-law proof.

The existing CW7 observation identity proves none of those cancellations beyond (5). The precise live factor is log Z_b(c) in (7), whose gradient is (10), together with the inaccessible mean q in the observation prior itself. That identifies the route's missing source identity without ruling out a different positive construction.

## Sources

The read-set and formulas were checked against `COMPONENTWISE-OLD-R7-HIDDEN-MEAN.md`, `exact-slack/cw7/01_retained_hidden_force_certificate.md`, the current host-separator audit, and LOW30's sufficient-statistic decoder, hidden-center programs, and expanded serial accounting. The original-gradient/HVP oracle model is the one in `exact-slack/exact32-proof-v2.tex`.
