# Independent audit of same-heat quarter-flow refinement

Date: 2026-10-04. Bounded review of the new same-heat alternative.

Audited report SHA256: `37efa7576cfac04fbe86141aa09963cd856f03cded90b74135d18929583a8ada`.

## Verdict

**Scoped PASS for the exact posterior refresh, contraction, legal finite VALUE/HVP DAG, and its additive query upper certificate.** No complete new source/hidden-host theorem or sublinear heat exponent is established.

One scope clarification is important: the blind-bump construction is a worst-case deterministic information lower bound for the **first-Picard integral at a fixed input**. It is not a lower bound for the exact nonlinear flow, Gaussian-input Lp strong approximation, randomized methods, or a posterior-law compiler. The report mostly states this limitation correctly; its summaries should preserve the fixed-input/first-Picard qualification.

The particular safe quadrature certificate is

    c_R ≤ max(c_P,R−3/2)=max(c_P,P−1), R=P+1/2.

This is a valid upper certificate for the displayed approximation, not a proof that any complete same-heat method must pay that exponent. The unproved protected/retained-host join prevents treating it as a new admitted family recurrence today.

## 1. Exact flow is a genuine posterior-preserving Markov kernel

For fixed y and 0<a<1, use

    x'=v,
    v'=−(x−y)−a grad V(x).

Its canonical Hamiltonian is a times

    H(x,v)=V(x)+|x−y|²/(2a)+|v|²/(2a).

Thus the displayed equations conserve H and preserve phase volume. The joint density proportional to exp(−H) is invariant. Its position marginal is Q_(a,y), and its independent velocity marginal is N(0,aI). This validates the refresh without a Metropolis decision. The scalar time normalization is internally consistent.

Bounded Hessian makes the vector field globally Lipschitz; V in C2 makes it C1. The finite-time flow and its first derivatives exist, with no third derivative assumption.

At T=π/2, variation of constants gives

    x(T)=y+sqrt(a)Z−a integral_0^T cos(s) grad V(x(s)) ds.

It contains one incumbent at heat a through its initial condition, with no reheated old call.

## 2. The O(a) coupling contraction is correct

Let d(t) compare two initial positions under the same fresh Z. With ||H_s||≤1,

    d(t)=cos(t)d(0)−a integral_0^t sin(t−s)H_s d(s) ds.

The kernel is nonnegative and its row mass is 1−cos(t)≤1. Therefore

    sup_t |d(t)|≤|d(0)|/(1−a),
    |d(T)|≤a|d(0)|/(1−a).

Couple the old approximation with a target initial position and choose the fresh Z independently of that entire coupling. Target invariance then gives exactly the stated W2 contraction. This is a contraction of posterior law error, including its mean and covariance defects, not merely of an old-source mean.

The same inequality applies to a same-record retained replacement of the initial position. Its old deletion error gains one factor a/(1−a). It does not certify the other deletions or observed-label contracts of a future host.

## 3. Actual first actions and leading carrier

Differentiating the integral equation uses only Hess V at actual path points. The endpoint starting-position derivative is O(a). The intermediate path derivative in the caller is bounded if the old actual D_y X0 is bounded. Because the direct initial-position term vanishes at T,

    D_y x(T)=I−a integral cos(s) Hess V(x(s))D_y x(s)ds
            =I+O(a).

This remains valid for a CW7 incumbent with actual I+O(sqrt(a)) center derivative. No law-error derivative is being taken.

If the complete old fresh first is O(sqrt(a)), all of its endpoint contribution is O(a sqrt(a)). The new momentum derivative is sqrt(a)I+O(a sqrt(a)). Consequently the complete leading row is the fresh momentum row, and the ordinary original-force leading symmetric block/curl conclusion follows in its stated normalized units.

At an intermediate time the harmonic row (cos(t)P_old,sin(t)I) is a coisometry, but the old actual remainder must still be included in the nonlinear record. At the endpoint its direct harmonic term vanishes. None of this creates an unread independent additive Gaussian keep: Z is used by the force path. The report correctly leaves protected chronology open.

The physical sqrt(a)-scale moment statement is understood around the posterior mode, with its dimension/moment factors. It is not a claim that |m−y| is uniformly O(sqrt(a)) for arbitrary entering y.

## 4. The finite DAG uses legal original queries

The original model in exact32-proof-v2.tex permits grad V(x) and directional Hess V(x)v at actual recorded points/directions. It does not supply potential values, a free dense Hessian, or derivatives of HVP outputs.

The finite layer

    q_i^[k+1]=q_i^[0]−a sum_(j<i)w_ij f(q_j^[k])

uses only original gradients. Store each actual point and gradient once. One directional forward or reverse sweep uses original HVPs there, with all incoming cotangents accumulated before the reverse action. There is no derivative of an HVP child and no hidden finite-difference Hessian modulus assumption.

One incumbent evaluation and its saved record feed all layers. This is lawful deterministic sharing of one input record, not aliasing nominally independent empirical banks. The new Z is sampled after the entering caller and independently of the incumbent's private records.

The query count Q_P+MN+Q_mode+Q_zero is correct for M layers and N grid intervals, up to fixed oracle conventions and separately priced numerical restoration. Mode and zero programs are not free. For the claimed exponent they require the admitted polynomial-logarithmic initialization/precision profiles, or their actual extra exponent must be added. Uniform caller domains and finite numerical encoding remain explicit return obligations, as the report states.

The weights are nonnegative and have row mass at most one. Finite encoded versions must preserve a bounded row mass and price their error; their exact positivity cannot be silently inferred from independently rounded cosine differences. This is a standard finite numerical task, not a newly identified heat-power obstruction.

### Optional arithmetic simplification of this exact DAG

The report honestly quotes O(MN²) vector arithmetic for dense weight multiplication. This is avoidable without changing the source or query count:

    w_ij = cos(t_i)[cos(t_(j+1))−cos(t_j)]
           +sin(t_i)[sin(t_(j+1))−sin(t_j)].

Maintain the two weighted prefix sums of f(q_j^[k]) at each layer. Each output row then takes two known scalar-vector products. This gives O(MN) vector arithmetic, in addition to the original-gradient evaluations and O(Dim) cost per vector operation. It does not improve the heat exponent of the original-query certificate.

## 5. Strong-error upper certificate

Centering at the exact posterior mode gives f(0)=0 and Lip(f)≤1. An admitted old sqrt(a)sqrt(Dim) moment profile and the bounded-Hessian ODE imply the same fixed-p moment scale for both q and its speed on [0,T].

Two errors must be separated:

1. Discrete Picard truncation: the finite-grid map is a contraction with factor a in the grid sup norm. Its initial displacement is O(a sqrt(a) sqrt(Dim)). After M layers, the tail is O(a^(M+3/2) sqrt(Dim)).
2. Product integration: replacing f(q(s)) on one interval by its left value loses at most Lip(f) times the path-speed envelope times T/N. The outer integral contributes a, and discrete/continuous stability contributes only 1/(1−a). The total is O(a^(3/2) sqrt(Dim)/N).

Thus the displayed estimate is valid, with mode, oracle and encoding errors separately added. Choose M fixed large enough for R and N at least a fixed-target multiple of a^(−(R−3/2)). The old W2 error gains order one, leaving room for the requested half-step. The original-query upper exponent follows.

This is an approximation to an exact posterior-invariant kernel. The finite DAG itself is not exactly posterior-invariant. Its output is a genuine positive probability law as a pushforward of Gaussian primitives; that fact alone does not give the native positive covariance reserves, protected increments or complete source grammar required for recursive use.

## 6. Exact scope of the deterministic bump obstruction

The report's first-Picard integral is

    I(h)=integral_0^1 h(u)du,
    x^[1](T)=sqrt(a)[1−a I(h)]

at the fixed normalized input (z0,Z)=(0,1).

Run any deterministic adaptive gradient/HVP integration algorithm on h0(u)=u/2. Collect ALL its queried physical points, normalized by sqrt(a), and use the points lying in [0,1], together with the endpoints, to partition the interval. If the total oracle budget is q, there are at most q+1 gaps. Queries outside the interval cause no difficulty because the bumps vanish there.

For a gap of length ell use psi(u)=ell r²(1−r)². At every queried point psi=psi'=0, so both gradient and every scalar HVP reply agree under h+=h0+psi and h−=h0−psi. Hence adaptive queries and stopping decisions reproduce the same baseline transcript. In dimension one all directional HVPs at a point reveal only a scalar multiple of h'(u), so additional directions do not defeat the construction.

The exact derivative maximum is

    max |2r(1−r)(1−2r)|=1/(3sqrt(3)),

so h' lies in [1/2−1/(3sqrt(3)),1/2+1/(3sqrt(3))], strictly inside [1/4,3/4]. The assembled h is C1. Integrating it gives a globally C2 potential, even though no higher derivative modulus is uniform. The physical scaling V_a(x)=a H(x/sqrt(a)) keeps this same Hessian sandwich.

The gap integrals give

    I(h+)−I(h−)=2 sum ell_i²/30 ≥1/[15(q+1)].

The two first-Picard endpoint values differ by at least a^(3/2)/[15(q+1)]. Any one common deterministic answer therefore has worst-case error at least half that separation. This is a valid information lower bound in the specified oracle model.

Important limits:

- The two potentials can depend on a and on the algorithm's completed baseline transcript. This is appropriate for a uniform worst-case Hessian-sandwich theorem, not a claim about one fixed potential's nonuniform asymptotic rate.
- The shown separation is for the first-Picard integral. A lower bound for the exact nonlinear flow would need to control its additional path-response terms.
- The selected Z=1 is one fixed input. No Gaussian-input Lp lower bound follows without an additional positive-measure/distributional argument.
- A deterministic adaptive transcript is covered. Choosing a different adversarial potential separately for every random seed does not prove a randomized lower bound.
- The original model has no potential-value query. If such queries were added, the indistinguishability construction would need modification.
- No posterior-law lower bound follows from a strong fixed-input integration obstruction.

Thus the fixture rejects a free uniform deterministic quadrature premise, but does not make the upper exponent P−1 necessary for the research program.

## 7. Quadratic and Stein identities

The exact quadratic map has coefficients c=cos(Tsqrt(1+a lambda)) and s=sin(Tsqrt(1+a lambda))/sqrt(1+a lambda). Its covariance identity is correct. Replacing the carrier's Z inside the additive correction by an independent Z' adds exactly 2a(1−s) to the variance, giving the claimed leading lambda a² variance error and (lambda/2)a^(3/2) W2 error. This is a useful explicit warning against independentizing a lawful shared endpoint record.

The finite-source Gaussian-IBP Stein identity is also correct with the stated full-owned-tape and Frobenius conventions:

    E[(S+b(S))·phi(S)−div phi(S)]
      =E[(r+b(S))·phi(S)+(P Dr*) : D phi(S)].

The scalar S=tau W test gives [(1+a lambda)tau²−1]E phi'(tau W), as claimed. Its rank-two coefficient is not a new differentiable VALUE oracle. Full matrix reconstruction and differentiating its HVP representation remain outside the free interface.

## 8. Adoption boundary

The exact posterior kernel and its incumbent contraction are solid. The literal finite approximation has the stated legal, additive old-call query structure. To admit it as R_R in the full induction still requires:

- Actual retained controlled-block enumeration and exact finite origins.
- The complete observed-label/caller and moving-center source-domain return.
- An unread/protected late Gaussian construction with its actual error, positive gap and cost.
- Whole finite tape, precision, replay and attached hidden-family verification.

Until that join exists, write the formula as a candidate law-level/query certificate. Do not call it a completed family improvement or a subpolynomial-dimension result.

The accompanying independent script passed 900 checks, including noncommuting finite matrix derivative fixtures, exact rational blind-gap areas, the prefix-sum identity, quadratic invariance, and verification of all three audited report output hashes. These diagnostics support the algebra and numerical fixtures; the analytic arguments above establish their scope.
