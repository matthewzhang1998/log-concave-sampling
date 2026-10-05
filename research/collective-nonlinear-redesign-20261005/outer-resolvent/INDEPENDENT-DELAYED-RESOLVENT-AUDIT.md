# Independent audit: delayed conditional-noise cancellation

2026-10-05. Audited construction: `/workspace/shared/delayed-resolvent-cancellation-20261005/DELAYED-CONDITIONAL-NOISE-THEOREM.md`.

## Verdict

The core conditional-noise lemma and its new original-VALUE construction are accepted. They give a dimension-uniform general-C2 fixed-grade bound

    ||R1 mu_Q-m3||_2 <= 5 A^(11/4) sqrt(D), 0<A<=1/12,

with positive inner rule tolerance A^(3/4), J=O(log^2(1/A)) original history VALUES, and two more original terminal VALUES. This is a new graph, with an explicitly changed genealogy. It is not a better estimate for the old graph ruled out by the radial counterexample.

The finite positive outer rule, existing native completion, all original numerical floors, actual radii, and all native occurrences must be added. The raw source is not itself a Gaussian at its own mean. This audit accepts the finite source/target theorem and its interface; it does not claim a fully expanded native compiler execution or an any-order induction.

## 1. The derivative and covariance identity have the correct live variables

Set h=R1 f for an arbitrary L2 vector test f. Hermite orthogonality gives

    ||h||_2<=||f||_2,
    ||Dh||_(L2,HS)<=||f||_2/2.

The derivative norm is summed over every output and input coordinate. There is no unexposed scalar-only restriction.

Condition on Y. In the intended application, X|Y is N(qY,sigma^2 I), and the complete future Gaussian noise bank Z is independent of X conditional on Y. For b(Y,Z), m=E_Z b, define e_Z=b-m. Differentiation in lambda of E F(X-m-lambda e_Z) and the Gaussian covariance formula yield

    lambda integral_0^infinity exp(-t)
      E sum_(j,k) C_jk partial_k partial_j F_i(X-m-lambda e_(Z_t)) dt,
    C_jk = sum_a partial_a b_j(Y,Z) partial_a b_k(Y,Z_t),
    Z_t=exp(-t)Z+sqrt(1-exp(-2t))Z'.

The coefficient is `Db(Z) Db(Z_t)^T`, and its operator norm is <=L^2. It has no derivative of Db. Its dependence on Y remains live. Its independence of X conditional on Y is the exact reason the next integration by parts is allowed without differentiating C.

Integrating in the k coordinate of X transfers the derivative to h and the conditional Gaussian density. The matrix contraction is `DF C`, because its (i,k) entry is sum_j partial_j F_i C_jk. The two required bounds are

    ||DF C||op <= K L^2,
    ||DF C||HS <= sqrt(D) K L^2.

Therefore the first term costs `sqrt(D) K L^2 ||Dh||_2`. The score term costs

    (K L^2/sigma) E|h(X)||eta|
       <= (K L^2 sqrt(D)/sigma)||h||_2,
    eta=(X-qY)/sigma.

This uses ordinary Cauchy-Schwarz on the original joint law. It does NOT pretend that h(X) and eta are independent. This is a legitimate one-energy sqrt(D) estimate rather than an induced higher-tensor norm multiplied by D.

The lambda integral contributes 1/2; the exp(-t) integral contributes one. Self-adjointness of R1 now proves exactly

    ||R1 e||_2 <= (K L^2 sqrt(D)/2)(1/2+1/sigma).

All covariance/matrix factors have their actual order. A nonsymmetric DF is allowed. Thus applying the lemma to F(u)=g(u-wg(u)) is valid, although F need not be a gradient.

For the true tail, finite-time finite-dimensional Gaussian path projections can first be used. Positive quadrature/integration bounds give the same bank operator first <=qA, uniformly in the projection. The resulting tails converge strongly in L2, and their conditional means converge by Jensen. Euclidean smoothing of F preserves its Lipschitz constant and converges with a vanishing uniform error. This supplies the infinite-bank/C1 passage without executing an infinite history or assuming a modulus of Dg.

## 2. Actual history and finite source have the required different roles

The true tail is

    V=q integral_0^infinity exp(-s) g(X_(delta+s))ds,
    Y=X_delta, q=exp(-delta).

By the OU Markov property it is independent of X_0 conditional on Y. Its complete future bank first is at most qA and its conditional mean is qR1g(Y).

The executed finite replacement is

    V_Q=q sum_j a_j g(r_jY+sqrt(1-r_j^2)N).

Given Y, it is independent of X_0, has bank first <=qA, and mean qQg(Y). Reusing N across nodes does not create a true OU future. None is needed: the proof compares the true and finite objects through their conditional means and the proved weak Jensen bound, rather than using a false multi-time identity.

The short prefix remains part of the proof, not a deleted contribution. Its replacement by wg(X_0) costs exactly at most

    (2sqrt(2)/3) A sqrt(D) w^(3/2)

before the final A-Lipschitz source. Moving that prefix into F(u)=g(u-wg(u)) costs at most qwA^3sqrt(D), with both source descendants changed coherently. The original nested F2-to-F1 replacement separately costs A^3sqrt(D). These three quantities cannot be omitted or charged to quadrature.

With w=sqrt(A), sigma=sqrt(w(2-w)), epsilon=A^(3/4), the complete exact-outer bias bound is

    sqrt(D)[A^3+(2sqrt(2)/3)A^2 w^(3/2)
               +qwA^3+qA^2 epsilon
               +q^2 A^3(1/2+1/sigma)].

Since sigma>=sqrt(w), this is <=5A^(11/4)sqrt(D) for A<=1/12. The numerical script tests the full coefficient over 1,000 logarithmically spaced values and observes a maximum below 2.905, but the analytical bound and the uniform range supply the theorem.

## 3. Positive finite quadrature

The author's normalized time-panel rule has a valid uniform Hermite-multiplier proof. Its final first-moment repair must use the diagonal bounds `||P_0-R1||<=1/2` and `||P_1-R1||<=1`, not a generic two-contraction bound of two. Starting at epsilon/3 is therefore sufficient for the repaired rule to have error <=epsilon.

An independently derived alternative is in Section 2 of `RESUMMED-COLLECTIVE-WEAK-GRADE-THREE.md`: dyadic correlation panels with Gauss-Legendre rules and one final midpoint have exact mass and first moment with no repair. Its certificate is

    2*2^(-K)+8*4^(-n) <= epsilon,
    K=ceil(log_2(4/epsilon)),
    n=ceil(log_4(16/epsilon)),
    J=Kn+1.

Its holomorphic operator proof is dimension-free, and all nodes are strictly in (0,1). Either certified rule may be used. The independent implementation uses this alternative and verifies positivity, mass, first moment and sampled high-chaos errors. Sampling up to chaos 10^11 supplements rather than replaces the uniform proof.

## 4. Native first/curl ports are uniform in the cutoff

Use an inner rule with exact first moment one-half, and use the same moment for the outer rule. Let w=1-q, sigma=sqrt(1-q^2), and denote the private delayed root by Z and the private inner quadrature root by N. Let b=V_Q, u=x-b, v=u-wg(u), H_v=Dg(v), H_u=Dg(u), H_0=Dg(x). Then

    b_x=q^2 sum a_j r_j Dg(U_j),
    ||b_x||<=q^2 A/2,
    ||b_(Z,N)||<=q A sqrt(1-q^2/4)<=A sqrt(3)/2.

The full-bank bound follows from the norm of each entire row and concavity of the square root, not coordinatewise summation. The actual residual derivative is

    E_x=(H_v-H_0)-w H_v H_u-H_v(I-wH_u)b_x,
    E_(Z,N)=-H_v(I-wH_u)b_(Z,N).

The leading difference is symmetric with norm <=A. The remaining x term has norm <=A^2 K, where

    K=w+q^2/2=(1+w^2)/2 <=13/24 for w=sqrt(A), A<=1/12.

The PSD source condition guarantees `||I-wH_u||<=1`. The product order above is unchanged.

Let beta_out=sum omega_i sqrt(1-t_i^2)<=sqrt(3)/2. The raw private first is at most

    A sqrt[beta_out^2(1+A K)^2
                      +A^2 q^2(1-q^2/4)].

The square lift into the outer G block has raw curl at most

    2 beta_out A^2 K+A^2 q sqrt(1-q^2/4).

The pinned half-variance adapter divides the source by sqrt(1/2), and therefore multiplies both first and curl by sqrt(2). This was independently verified against `nonlinear-bridge-repair/INDEPENDENT-REPAIRED-SOURCE-AUDIT.md`, Section 3, and the matrix-free backbone's source/adapter formulas. Thus normalized private first <=2A and curl <=3A^2 are valid uniformly at the actual cutoff. There is no hidden sigma^-1 in these source ports.

The raw caller first is <=A(1+A K)/2. Capturing and restoring the actual all-private-roots-zero origin doubles this bound at most; the normalized anchored bound is safely 2A. Additional exterior anchors/scales/source parameters require their own live firsts. If the cutoff varies with such an exterior parameter, its scalar row derivatives must also be retained and guarded; this proof does not silently freeze a live cutoff.

A safe conditional residual energy envelope is

    ||E_raw||_(Lp|caller=z)
      <= A^2 {[(1+w^2+wq^2 A)/4]|z|
                        +kappa_p(1+wqA)sqrt(D)}.

The actual caller origin obeys this caller term and is not discarded. Private dimension is exactly 3D after adding the finite outer G root. Literal zero source and all-zero anchored inputs produce literal zero.

## 5. Complete original-VALUE bill and precision

An inner terminal occurrence uses J original history VALUES and the two nested terminal VALUES. Its residual uses J+3. The full guarded baseline/residual positive own-mean bill is

    Q_rule/setup+Q_captured+N_out N_B
      +(J+3)N_out N_E+Q_known/numerical/replay.

`N_B,N_E` mean complete native occurrence counts at their new actual radii, tolerances and private dimension 3D. They include complete independent banks, native clocks, filters, marks, fills, modes, origins, replays and numerical versions. First/adjoint sweeps retain and charge every original HVP site in the literal source. No HVP is differentiated.

The independent COMPLETE baseline/residual banks, positive half-variance shares and pinned native coisometry are required before a completed law can be coupled. Declaring `ell_E=2A`, `a_seed=2A`, `padding mu=A` gives the existing bracket

    16A^2+8A^3+8A^(5/2)=O(A^2),

and hence O(A^4sqrt(D)) from the O(A^2sqrt(D)) residual energy, IF every literal native guard holds. The baseline keeps its separate gradient-order/radius/caller allowance. A<=1/12 does not certify all of those guards.

The source's original-VALUE absolute floor is at most `(1+A+wqA^2)nu` for the terminal and `(2+A+wqA^2)nu` for the residual, before separate caller-origin/replay charges. Scalar and HVP floors remain separate.

Use `sigma=sqrt(w(2-w))` numerically. No inverse sigma or time delta is executed in the graph. The singular factor sigma^-1 belongs to the analytical weak comparison only. Inner and outer square roots can be evaluated as `sqrt((1-r)(1+r))`; verify actual positive weights and moments at the used precision. q(A), w(A) and all row-coefficient errors retain absolute budgets. This establishes no hidden inverse-A VALUE replication; it does not claim zero known-scalar bit complexity or erase the costs of high precision and actual source dimensional arithmetic.

## 6. Explicit source-valid obstruction check

The old two-shell family is

    g_(A,D)(x)=A sqrt(D) phi_(A,A^2)(|x|/sqrt(D)) x/|x|,
    phi(q)=(1/2)H_(A^2)(q-1/2)
                     +(1/2)H_(A^2)(q-[1-A tau]),
    tau=(1-1/sqrt(2))/4.

It is zero near the origin, its potential is C2, and its radial and tangential Hessian eigenvalues lie in [0,A]. All the preceding proofs apply directly and uniformly even though this source depends on A and D. There is no frozen source shape or missing Hessian-modulus assumption.

For the old obstruction sequence `A=1/n`, `D=(3000n)^2`, the new exact-outer ratio is bounded by

    ||new_mean-m3||_2/[A^2 sqrt(D)] <=5 n^(-3/4) ->0.

This is an ACTUAL finite-D consequence of the uniform theorem. For instance n=10^6 gives an upper ratio below 0.000159 before separately added outer/completion/numerical floors, whereas the old graph had a lower ratio at least 0.002. It is therefore consistent only because the executed graph has changed.

The independent diagnostic also evaluates this exact smooth source formula on the deterministic Gaussian coefficient rows used in the old radial reduction. With A from 10^-2 through 10^-6, the outer first-chaos defect divided by A^2 decreases from approximately -4.04e-4 to -4.22e-8; divided by A^3 it remains near -0.04. This is a Gaussian-row-limit diagnostic, not a replacement for a finite-D path theorem. The actual finite-D assertion is the preceding uniform bound.

Noncommuting C-infinity convex-gradient 2D source checks differentiate every original site of the literal new graph. The sampled normalized curl/A^2 stays below 0.48 and first/A below 0.04, within the proved 3 and 2 bounds. Numerical examples are not the source qualification: the theorem applies to the full C2 interval class.

## 7. Quadratic accuracy, optional exact correction, and next type

The uncorrected source is not exactly quadratic-exact. For `g(x)=Kx`, its inner conditional mean is

    Kx-(w+q^2/2)K^2x+(wq^2/2)K^3x.

Relative to the genuine inner target, its defect is

    -(w^2/2)K^2x+(wq^2/2-1/4)K^3x.

At w^2=A this is O(A^3sqrt(D)), below the A^(11/4) target scale. The script verifies this exact polynomial identity.

If a later consumer requires exact quadratic means, a possible explicit correction is

    C(x)=(w^2/2)g(g(x))+(1/4-wq^2/2)g(g(g(x))).

Its two coefficients are nonnegative. It cancels the displayed defect on every anisotropic quadratic. At w^2=A it adds at most `(3/4)A^3sqrt(D)` to the general target bound. It uses two additional original VALUES when the residual baseline g(x) is already recorded. This is a separately named variant, not part of the accepted J+3 bill or unchanged first/curl declarations; those must be re-audited if the correction is installed.

The returned finite residual has the near-gradient input type of the existing own-mean compiler. It is not proved to be a convex gradient suitable for recursively applying this theorem. A normalized original anchored/scale source retains its original interval class, but any changing cutoff and physical scaling require their actual caller chains and guards. Consequently the accepted result is a cheap fixed-grade general-class advance with a genuine next consumer, not an any-order recurrence. Higher-grade cancellation and reusable closure remain open.

## 8. Independent check outputs

- `check_delayed_resolvent_independent.py`
- `delayed_resolvent_independent_checks.json`

The script passes positive-rule checks, high-chaos sampled multipliers, exact quadratic algebra, the actual smooth radial source row calculation, noncommuting-source first/curl calculations and a cutoff-constant sweep. None instantiates the full native completion or substitutes a diagnostic for the analytical proof above.
