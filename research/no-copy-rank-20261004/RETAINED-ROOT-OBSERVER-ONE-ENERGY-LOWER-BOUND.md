# A constant mixed correction cannot be reinserted under a retained-root observer

Exact finite scalar-output example, 2026-10-04. This is an obstruction to upgrading the established CONSTANT covariance service to a SAME-retained-root nonlinear-observer service. It is not an impossibility theorem for a constructor that explicitly retains and pays the original root/current. The example is an admissible small-full-gradient descendant, with an actual nonzero near-gradient curl; it does not claim to be the specific K3 lambda descendant.

## 1. An admissible marked source on a two-dimensional Gaussian tape

Fix 0 < g <= 1, 0 < A <= 1/2, and kappa = A^(1+g). Use X=(X1,X2) standard Gaussian, physical coisometry P=(1,0), and known descendant coisometry B=(1,0). Put

    E(X) = (A/2) sin X1 + (kappa/2) X2,
    G(X) = (kappa X1, 0),
    H(X) = B G(X) = kappa X1.

G is the exact genuine-gradient VALUE source of kappa X1^2/2. It has full first kappa, no changing caller, and uses the SAME complete tape as E. There are no higher-derivative or finite-to-exact oracles. The physical square lift is

    P* E = (E,0),
    D(P*E) - D(P*E)* = [[0,kappa/2],[-kappa/2,0]].

Thus the actual recorded curl is nonzero and <= kappa. Also Lip E <= sqrt(A^2+kappa^2)/2 <= A. Every original/root first is bounded explicitly. The source is centered, with actual L2 energy

    e = (1/2) sqrt(A^2 v_s + kappa^2),
    v_s = (1-exp(-2))/2.

In particular

    (sqrt(v_s)/2) A <= e <= (sqrt(v_s+1)/2) A.

All fixed higher Gaussian moments obey the analogous dimension-free bounds. The zero-source anchor is literal. These are a conservative first radius A and curl radius kappa exactly of the admitted mixed-return type. There is no reliance on an unmarked dimension factor.

If an original bounded-Hessian gradient VALUE representation is desired, define

    U(x) = |x|^2/4 - (A/2)cos x1 + (kappa/2)x1 x2.

Then E=P[grad U(x)-x/2], an exact one-gradient-VALUE construction with only a known linear subtraction. The Hessian perturbation of U from I/2 has operator norm at most (A+kappa)/2<=A. Hence for A<=1/10, its original Hessian lies in [.4 I,.6 I], and grad U(0)=0. Multiplying the two nonquadratic coefficients by eta gives the same representation for section 6. This strengthens the regularity/source certificate but still does not assert the six-VALUE K3 genealogy.

## 2. One finite clock node already computes the exact constant cross

Take one coarse standard Gaussian R=(R1,R2) and one independent fine standard Gaussian V. The complete query is X=(R+V)/sqrt(2). This has numerical positive root and innovation widths c=v=1/sqrt(2). The conditional Gaussian Jacobians are exactly

    J_E(R) = ((A/2) exp(-1/4) cos(R1/sqrt(2)), kappa/2),
    J_H(R) = (kappa,0).

Consequently the scalar mixed covariance field (scalar Sym is the identity) is

    S(R) = J_E(R) J_H(R)* = lambda q(R1),
    lambda = (A kappa/2) exp(-1/4),
    q(r) = cos(r/sqrt(2)).

Let m=E q(R1)=exp(-1/4). Then

    E S = lambda m = (A kappa/2) exp(-1/2)
        = Cov(E(X),H(X)).

The last equality follows directly from E[X1 sin X1]=exp(-1/2). Hence this single positive node is an EXACT finite quadrature for this particular mixed covariance, not an approximation to an omitted covariance tail. All widths, weights, dimensions and values are finite numerical choices. The point is a logical obstruction even at the ideal conditional Gaussian reference, before any pair-execution/numerical errors are added.

Define the two joint reference laws, retaining the same entire coarse root R,

    Y = sqrt(1-lambda q(R1)) Z,
    Ybar = sqrt(1-lambda m) Zbar,

where Z and Zbar are independent standard Gaussians independent of R. The square roots describe analytical Gaussian reference laws; they are NOT proposed executable random-covariance primitives. For lambda <= 1/8 all variances have a fixed positive gap. The constant reference has exactly the averaged covariance, and both Y and Ybar are centered.

## 3. Explicit bounded smooth observer and exact joint-law lower bound

Use the parameter-independent observer

    phi(r1,r2,y) = (cos(r1/sqrt(2))-exp(-1/4)) cos(y) / 4.

It is C-infinity, bounded in absolute value by 1/2, and globally 1-Lipschitz. Indeed its two nonzero first-derivative bounds are 1/(4sqrt(2)) and (1+m)/4, whose Euclidean square sum is <1. All higher derivatives are bounded by numerical constants independent of A,g. Its Hessian operator norm is also <1, as follows already from the entrywise bounds.

Let

    V_q = Var(q(R1)) = (1-exp(-1/2))^2/2 > 0.

Under the constant joint law E phi(R,Ybar)=0. Under the conditional law,

    Delta := E phi(R,Y) - E phi(R,Ybar)
           = exp(-1/2) Cov(q, exp(lambda q/2)) / 4.

For any independent copy q' and any lambda>0,

    Cov(q,exp(lambda q/2))
      = (1/2) E[(q-q')(exp(lambda q/2)-exp(lambda q'/2))]
      >= (lambda/2) exp(-lambda/2) V_q.

Here q,q' lie in [-1,1], so the displayed inequality is the ordinary mean-value theorem. Therefore the EXACT lower bound is

    Delta >= lambda V_q exp(-(1+lambda)/2) / 8.           (1)

In particular, for lambda <= 1/8,

    d_BL(Law(R,Y),Law(R,Ybar)) >= Delta,
    W1(Law(R,Y),Law(R,Ybar)) >= Delta
      >= [V_q exp(-9/16)/8] lambda
      = c A kappa >= c' kappa e.                        (2)

Consequently W2 also dominates Delta. These are joint-law lower bounds allowing arbitrary couplings of the root, so they also apply to the more restrictive SAME-root coupling requirement.

For completeness the exact same-root conditional W2 distance satisfies

    lambda sqrt(V_q)/(2sqrt(1+lambda))
      <= [E W2^2(Law(Y|R),Law(Ybar|R))]^(1/2)
      <= lambda sqrt(V_q)/(2sqrt(1-lambda)).              (3)

This follows by coupling the one-dimensional Gaussian noises monotonically and using

    sqrt(1-lambda q)-sqrt(1-lambda m)
      = -lambda(q-m)/(sqrt(1-lambda q)+sqrt(1-lambda m)).

Thus the retained-root gap is genuinely order A kappa, not merely an upper bound that a better coupling could remove.

## 4. The unobserved marginal has only a second-order gap

Dropping R really does remove the first-order discrepancy. For the 1-Lipschitz output observer cos y,

    E cos Y - E cos Ybar
      = exp(-1/2)[E exp(lambda q/2)-exp(lambda m/2)].

Taylor's theorem and the same fixed interval give

    lambda^2 V_q exp(-(1+lambda)/2)/8
      <= E cos Y-E cos Ybar
      <= lambda^2 V_q exp(-(1-lambda)/2)/8.               (4)

There is also a full marginal W2 upper bound of this second order. Write v=1-lambda m, delta=-lambda(q-m), and let p be the density of Y and gamma_v that of N(0,v). Gaussian integration gives, for independent delta,delta',

    1+chi^2(p|gamma_v)
      = E[1-delta delta'/v^2]^(-1/2).

The binomial expansion has positive coefficients <=1. Its n=1 term vanishes because E delta=0. Since |delta| <= 2lambda, eta=2lambda/v<1, and E delta^2=lambda^2 V_q,

    chi^2(p|gamma_v)
      <= lambda^4 V_q^2/[v^4(1-eta^2)].                 (5)

The elementary Gaussian transport bound W2(p,gamma_v) <= 2sqrt(v chi^2(p|gamma_v)) yields, when lambda <= 1/8,

    W2(Law(Y),Law(Ybar)) <= 3 lambda^2 V_q.             (6)

One way to obtain that transport bound without any special property of the mixture is to solve the mean-zero Gaussian weighted Poisson equation for p/gamma_v-1. Gaussian Poincare bounds its flux energy by v chi^2. Along the linear density interpolation (1-t)gamma_v+tp, the energy is at most v chi^2/(1-t); its square-root integral is 2sqrt(v chi^2). Approximation handles a merely square-integrable density ratio. Alternatively the usual Gaussian transport-entropy inequality gives a slightly better constant. Together (4) and (6) prove a genuine marginal order lambda^2 and a joint order lambda.

## 5. Exact powers and the missing current

For kappa=A^(1+g), the retained-root defect is

    Theta(A^(2+g)) = Theta(kappa e).

The unobserved marginal correction error is

    Theta(A^(4+2g)) = Theta(A kappa^2 e).

This exactly matches the one-energy grade of the admitted CONSTANT small-gradient mixed return, while the forbidden observer reinsertion loses the factor

    1/(A kappa) = A^(-2-g).

There are no logarithmic qualifications in this finite example. The observer, dimension, clock and Gaussian variance gap are all fixed. Even an arbitrarily accurate finite execution of the ideal conditional references cannot erase (1). For example if both executed joint laws are within o(kappa e) in W1 of the displayed reference laws, triangle inequality preserves a positive c kappa e lower bound.

The missing term can be written directly as the derivative at zero covariance amplitude:

    d/dlambda E phi(R,sqrt(1-lambda q(R1))Z)|_(lambda=0)
      = exp(-1/2) V_q / 8 > 0,

whereas replacing q(R1) by its mean before applying the observer makes this derivative zero. This is a root/observer covariance current. Integrating the root is legitimate for the constant reserve's marginal target; it does not license a later nonlinear observer of that same root.

Scope: the example rules out a generic upgrade based solely on the cited small-full-gradient/one-energy covariance frames. It neither proves that the K3 lambda defect itself attains this lower bound nor excludes a new SAME-root positive channel that explicitly carries the joint current. An actual K3-specific obstruction or construction requires its extra genealogy.

## 6. Arbitrarily small actual energy, including the weak A^3.9 profile

For any 0<eta<=1, replace E by E_eta=eta E and leave G,H unchanged. Keep the declared source envelopes A and kappa. The actual first and actual curl only become smaller, the complete sources/callers remain admissible, and

    e_eta = eta e,
    S_eta(R) = eta lambda q(R1),
    Cov(E_eta,H) = eta lambda m.

Apply the SAME parameter-independent observer phi from section 3. Every displayed bound repeats with lambda_eta=eta lambda. Consequently

    joint W1 and W2 gap = Theta(eta A kappa)=Theta(kappa e_eta),
    marginal W2 gap = Theta(eta^2 A^2 kappa^2)
                    = Theta(eta A kappa^2 e_eta)
                    <= O(A kappa^2 e_eta).

Here the joint W2 upper bound is supplied by the explicit same-root coupling in (3); the lower bound is joint-law and requires no forced coupling. The parent source radius is an envelope, not a demand that every example saturate it. In particular, choose

    eta=A^(29/10), so e_eta=Theta(A^(39/10))=Theta(A^3.9).

This gives the native weak nominal energy exponent in the fixed tape dimension two. Its exact comparison grades are

    retained-root gap: Theta(A^(49/10+g)),
    actual marginal gap: Theta(A^(49/5+2g)),
    admitted constant-service envelope: O(A^(69/10+2g)).

Thus the retained-root error is still larger than the admitted envelope by the factor A^(-2-g). For g=1 these three grades are A^5.9, A^11.8, and A^8.9, respectively. The eta=1 construction saturates the constant-service order; the eta<1 construction need not do so. Neither construction is asserted to have the specific K3 genealogy.

## 7. This is not an artifact of the single-node choice

For any finite positive clock with weights w_j summing to one and c_j^2+v_j^2=1, this particular source has the exact mixed-covariance identity at EVERY node. If the same coarse root is retained across nodes, its coefficient is

    S_clock(R) = (A kappa/2) Q(R1),
    Q(r) = sum_j w_j exp(-v_j^2/2) cos(c_j r),
    E Q = exp(-1/2),
    Var Q = exp(-1) sum_(i,j) w_i w_j[cosh(c_i c_j)-1].

Here |Q|<=1 and |Q'|<=1. The observer (Q-EQ)cos y/4 has the same uniform bounded-smooth profiles, and all of (1)-(6) repeat with lambda=A kappa/2 and V_q=Var Q. If a clock gives total weight at least theta to c_j>=c0>0, then

    Var Q >= exp(-1) theta^2[cosh(c0^2)-1] > 0.

Thus any such finite positive clock has the same dimension-free leading retained-root defect. If instead every node has its own retained independent root, the exact variance becomes exp(-1) sum_j w_j^2[cosh(c_j^2)-1]; the argument still applies with that explicit clock factor. The one-node instance avoids every clock logarithm and is exact for the displayed source. This analysis concerns the ideal reference after conditional comparisons, not a claim that a generic source clock can be replaced by one node.

## 8. Source connections

- r-native-twin/SMALL-FULL-GRADIENT-MIXED-DESCENDANT-RETURN.md: its ideal conditional coefficient and the A kappa^2 e owned-root marginal grade.
- p-native-twin/INDEPENDENT-SMALL-GRADIENT-MIXED-RETURN-AUDIT.md: the constant-target boundary and explicit warning against observer reinsertion.
- r-native-twin/AUXILIARY-CONSTRAINT-HEAT-AND-VARIANCE-PRESERVING-PATH.md, section 4: root-dependent observers require complete same-endpoint currents beyond fixed-test covariance frames.
- p-native-twin/COMMON-AND-OPPOSITE-CONSTRAINT-MODES.md: subtracting an unmarked baseline loses precisely the joint information not recoverable from a constant mixed covariance alone.
