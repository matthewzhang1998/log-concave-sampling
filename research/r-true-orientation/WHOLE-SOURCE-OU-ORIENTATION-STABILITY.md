# One-energy stability of true orientation under whole-source OU heat

## 0. Bounded analytical return

Let f:R^n->R^n be the SAME finite square source, with centered Gaussian energy e, first at most A and pointwise curl C=Df-Df* bounded in operator norm by kappa. Define

    O_f=Cov(f)-Sym B_f,
    B_f=integral_0^1 E[J_u²]du,
    J_u(X)=E Df(sqrt(u)X+sqrt(1-u)Z).

Let P_t=exp(-t N) be the Gaussian OU semigroup on the ENTIRE input array, with every source alias preserved inside each original evaluation. Then

    ||O_f-O_(P_t f)||HS <= kappa e min(1,sqrt(2t)),         (1)
    ||O_f-O_(P_t f)||op <= A kappa min(1,sqrt(2t)).         (2)

More sharply replace min(1,sqrt(2t)) by

    eta(t)=sup_(k>=1) [1-exp(-2kt)]/sqrt(k).

This is an energy-relative target restoration. It does not require a small L2 VALUE difference between f and P_t f. It does not independently heat primitive Hessians, construct an averaged VALUE oracle, or certify an actual finite source/caller implementation. It is an analytical bridge for the common-input target only.

The operator convention is the same as `astra-coherent-gram/closed-curl-operator/CLOSED-CURL-SWAP-ROW-FRAME.md`, SHAba28738df7cc0724b2787a3cd59f4d75027ae10cd38db25a476cd64a60293911, sections1--2. In particular T_0=0, and Aop=I-E-T; the symbol Aop is distinct from the numerical first bound A.

## 1. Exact curl/divergence identity, with normalization

Write the degree-k Hermite coefficient of f as f_(i;a1,...,ak), symmetric in its k input indices. For k>=1,

    (T_k f)_(i;a1,...,ak)
      =(1/k) sum_l f_(a_l;i,a1,...,omit(a_l),...,ak).

The curl coefficient on degree k-1 is

    C_(ij;beta)=k[f_(i;j,beta)-f_(j;i,beta)].

Let delta be Gaussian divergence on the j index: delta_j=x_j-partial_j. It creates one Hermite input index. Symmetrizing that created index with beta yields EXACTLY

    delta(Curl f_k)=k Aop_k f_k.                        (3)

For k=1 this reads delta(F-F*)=(F-F*)x, agreeing with Aop_1 f=(F-F*)x. Thus the denominator in(3) is k, not k+1.

For a vector field U in degree k-1, its divergence satisfies

    ||delta U||2 <= sqrt(k) ||U||2.                      (4)

Indeed the squared input norm is(k-1)! times the coefficient norm; the output norm is k! times its symmetrized coefficient norm, and symmetrization is an orthogonal contraction. This is a fixed-chaos statement and is dimension-free.

Contract the physical output with a fixed u. Equations(3)--(4) give

    ||u* Aop_k f_k||2 <= k^(-1/2) ||u* C_(k-1)||2.       (5)

No derivative of the curl is evaluated or estimated.

## 2. The small-time row frame

Orthogonality of distinct Hermite degrees and(5) imply

    ||u* Aop(I-P_(2t))f||2²
      <= eta(t)² sum_(k>=1)||u* C_(k-1)||2²
      = eta(t)² E|u*C|²
      <= eta(t)² kappa² |u|².                           (6)

The elementary inequality1-exp(-x)<=min(1,sqrt(x)) gives

    eta(t)<=min(1,sqrt(2t)).

The same proof only needs the weaker Gaussian row-frame condition E[CC*]<=kappa_row² I. Therefore(6) also holds with kappa_row=sqrt(||E CC*||op) when the pointwise operator hypothesis is unavailable. The pointwise bound in section0 supplies that row frame immediately.

For smooth finite Hermite sums the calculation is algebraic. For a Lipschitz weakly differentiable source, its Gaussian Sobolev chaos expansion, orthogonal truncation and L2 closure justify(3)--(6). This does not differentiate a saved sampler Hessian or import a Hodge decomposition.

## 3. Orientation and the single Hilbert mark

The cited swap identity is

    O_f=Sym E[(f-Ef) outer Aop f].

Aop is a bounded self-adjoint degreewise operator and commutes with OU. Since OU is self-adjoint,

    O_(P_t f)=Sym E[(f-Ef) outer Aop P_(2t)f],
    O_f-O_(P_t f)=Sym E[(f-Ef) outer Aop(I-P_(2t))f].      (7)

View f-Ef as a Hilbert-Schmidt map from Gaussian L2 to the physical output; its Hilbert-Schmidt norm is e. The second factor in(7) has output-row operator norm at most eta(t)kappa by(6). The HS/operator product inequality gives(1).

Gaussian Poincare gives Cov(f)<=A²I from the actual first bound. Using that ordinary row frame for the first factor gives(2). Symmetrization does not increase either norm. Exactly one f occurrence supplies the energy mark.

The estimates apply conditionally on fixed exterior data whenever the stated conditional first/curl/energy premises hold. For a coisometric physical projection P, projecting the square orientation by P and P* cannot increase the norms. Integrating conditional energy profiles preserves the same one-energy bound.

## 4. Source and consumer qualifications

P_t f means the heat of the WHOLE original finite VALUE graph under one Gaussian input substitution. Every nonterminal node, repeated query, twin bank, origin and coefficient alias is evaluated on that same substituted array. The Gaussian expectation defining P_t f is analytical in this note.

For h=sqrt(2t), equation(1) supplies an O(h kappa e) restoration tolerance between two constant orientation targets. If both appear inside a fixed positive covariance gap, the usual covariance transport estimate can consume that HS tolerance. It requires no derivative of the restoration error. It does not yield a finite program for either target.

A product of independently heat-averaged primitive Hessians is generally different from the derivative/commutator of this whole-source heat. Equation(1) does not identify them. In particular, cancellations making E=F-H small must stay inside the same complete E before the OU substitution. A source-specific extraction identity, finite VALUE implementation, full retained-root/private/caller return, and the true-orientation adjoint/two-curl ordering are still needed for a positive native packet.
