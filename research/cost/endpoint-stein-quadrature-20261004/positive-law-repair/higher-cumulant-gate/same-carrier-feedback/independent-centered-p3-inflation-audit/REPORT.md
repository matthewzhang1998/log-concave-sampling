# Independent exact-text audit: centered actual-P3 covariance correction

Date: 2026-10-05 UTC.

## Verdict

**PASS for the two precise candidates in the frozen corollary. No mathematical blocker was found.** At A=D^(-1/4), each centered covariance-derivative correction, using either the actual second-substitution covariances C2 or the actual first-displacement covariances C1, leaves an eventual L2(gamma_D) error at least c_* A. The advertised fourth-order allowance is Lambda A^4 sqrt(D)=Lambda A^2. On an admitted sequence with polynomial accuracy and bounded or polynomial public versions, every fixed public-log Lambda is too small to absorb the separator.

This is an actual canonical m3 versus literal finite F3_Q comparison. It is stronger than a single-history-only objection, but remains a refutation of the specified deterministic centered derivative corrections. It is not a refutation of the raw source's OWN mean, the exact live shifted/resummed current, the bounded-remainder Gaussian-backbone VALUE consumer, a separately restricted A sqrt(D)=O(1) regime, or every possible fourth-order sampler. No unconditional endpoint-law impossibility is established by this conditional-mean calculation alone.

Primary reviewed text, unchanged:

- `../RADIAL-INFLATION-REFUTES-CENTERED-P3-COVARIANCE-CORRECTION.md`
- SHA256 `810128a4b47c568273f05442a510498ac5f77eedf1cb864366368ebb4bd356af`

The canonical centered single-history proof, literal P3 raw-source definition, bounded-remainder consumer, and positive quadrature construction were read directly. Their five exact pins are in `MANIFEST.json` and the independent checker. No author checker was imported or executed; numerical checks supplement the analytical reasoning below.

## 1. The same original field is admissible at the F2 coherent radius

The field is g(x)=lambda x+d A f(x), with c=1/2, d=1/4, lambda=c A, f(x)=psi(|x|-R0)n(x), a2=lambda/2-lambda^2/4, R0=(1-a2)sqrt(D), and f(0)=0. The normalized bump defining psi is even, nonzero, supported in [-2,2], and strictly decreasing on (0,2). Its maximum is below one: the normalization satisfies Z>=2 exp(-4/3), so psi'(0)<=exp(1/3)/2<1. This is an elementary bound independent of numerical integration.

For D>=256, R0-2>sqrt(D)/2. Consequently f vanishes on a neighborhood of the origin and extends smoothly there. Its radial and tangential derivative eigenvalues are psi' and psi/r. Thus the eigenvalues of Dg are lambda+d A psi' and lambda+d A psi/r. They lie in [0,A] globally; the tangential denominator is only used where r>=R0-2. The field is the gradient of the explicitly integrable radial convex potential. In particular g(0)=0, the Hessian sandwich holds, and |g-lambda Id|<=epsilon=d A everywhere.

The radius changes from the F1 coherent radius to the F2 coherent radius. No source field is changed between ancestors, the middle force, and the terminal evaluation. A smooth member of the allowed C2 class is enough to refute a uniform theorem on that class. The derivatives below are analytical features of this particular smooth counterfamily, not new producer leaves or assumed regularity of every source.

## 2. Genuine F2 ancestry and exact conditional covariance

Let X_t be the genuine stationary OU history conditioned on X0=x. The first future substitution at time t has F1,t=lambda H1,t+R1,t, with |R1,t|<=epsilon pathwise. Substitution into the same original g and positive integration give

    F2_H=lambda H1-lambda^2 H2+R_H,
    |R_H|<=epsilon(1+lambda),
    H1=int_0^infinity exp(-t)X_t dt,
    H2=int_0^infinity t exp(-t)X_t dt.

The double future integral uses X_(t+u), and its convolution kernel is s exp(-s). Replacing the future ancestry by independent copies would be incorrect; none are introduced here.

The Brownian representation makes the conditional covariance transparent. Given the endpoint, the noise kernels of H1 and H2 are, respectively,

    k1(s)=exp(-s)/sqrt(2),
    k2(s)=(s+1/2)exp(-s)/sqrt(2).

The integrals of k1^2, k1 k2, and k2^2 are 1/4, 1/4, and 5/16. Their endpoint means are x/2 and x/4. Independently, differentiating the stationary Laplace covariance K(a,b) from the note gives unconditioned moments 1/2, 3/8, 3/8; subtracting regression products 1/4, 1/8, 1/16 gives the same conditional answer. Hence

    lambda H1-lambda^2 H2=a2 x+b_H N_H,
    b_H^2=lambda^2/4-lambda^3/2+5lambda^4/16.

N_H is standard conditionally on x and independent of that endpoint. It may be correlated with R_H. The variance polynomial is positive for nonzero lambda: 4-8lambda+5lambda^2 has negative discriminant and positive leading coefficient. No path approximation or covariance oracle is used in this derivation.

## 3. Literal finite I2_Q and its independent level roots

The P3 source definition has H shared across all middle nodes, J shared across all inner nodes, and H independent of J and the outer root. Define beta_mid=sum_j v_j sqrt(1-tau_j^2) and beta_in=sum_k u_k sqrt(1-sigma_k^2). Exact mass and first moment of the inner rule imply

    L_j=lambda(y_j/2+beta_in J)+r_j, |r_j|<=epsilon.

Expanding only the linear part of the same g gives the exact full finite identity

    I2_Q=a2 x+lambda(1-lambda/2)beta_mid H
                 -lambda^2 beta_in J+R_Q,
    R_Q=-lambda sum_j v_j r_j
                 +sum_j v_j (g-lambda Id)(y_j-L_j),
    |R_Q|<=epsilon(1+lambda).

Independence between H and J yields exactly

    b_Q^2=lambda^2(1-lambda/2)^2 beta_mid^2
                 +lambda^4 beta_in^2.

The square on beta_mid retains every cross-node covariance of the shared root. Substituting independent roots per node or identifying H and J changes the displayed covariance. The normalized combination N_Q is independent of x but generally correlated with R_Q. Nothing in the proof discards that correlation.

Positivity gives beta_mid,beta_in<=1. The middle multiplier certificate at Hermite degree two gives |sum v_j tau_j^2-1/3|<=delta_mid. Since sqrt(1-tau^2)>=1-tau^2,

    beta_mid>=2/3-delta_mid>=7/12,
    beta_mid^2-1/4>=13/144.

At D>=256, delta_mid<=A^2<=1/16 already suffices for this slightly weaker 7/12 bound. No inner multiplier accuracy is needed for this analytical counterexample; beta_in<=1 and the exact moment suffice. A tighter inner accuracy can nevertheless be imposed for any independent source/compiler guard. Positive finite rules satisfying all the stated tolerances exist by the imported dyadic construction.

## 4. Full nonlinear shell response survives the actual second substitution

With B=epsilon(1+lambda), coupling each actual source to its own exact Gaussian backbone gives, pointwise in x,

    |psi_a(x)-E_N g((1-a2)x-b_a N)|<=A B=O(A^2).

This bound uses the global Lipschitz constant of the original g and the pathwise bound on R_a. R_a and N_a remain on their actual common probability space; no independence is needed. This proves the main comparison even though the conditional means of F2_H and I2_Q are not identical.

For X~gamma_D, put S_D=|X|-sqrt(D), R=sqrt(D)=A^(-2), k=1-a2, and w=-b_a N. The radial expansion around k X has the following scales, uniformly over all admitted rules:

- radial noise n(X) dot w has fixed-p norm O(A);
- full noise has fixed-p norm O(A R);
- its quadratic transverse radial increment has limiting size b_a^2 R/2;
- the cubic geometric remainder has fixed-p norm O(A^3 R)=O(A);
- the direction error has fixed-p norm O(A).

Here the denominators are controlled on |X|>=R/2 and |w|<=k|X|/2. Gaussian tail estimates make the complementary contributions negligible for bounded f. The dimension-independent shell moments of S_D handle the coherent rescaling k S_D-S_D=O_Lp(A). Since b_H^2=c^2 A^2/4+O(A^3) and b_Q^2=c^2 A^2 beta_mid^2+O(A^3), the two limiting radial inflations are

    h_H=c^2/8, h_Q=c^2 beta_mid^2/2.

The O(A^3) variance changes give O(A^3 R)=O(A) changes in radius, so are below the retained response. Applying the Lipschitz smooth step, conditional expectation, and the outer factor d A yields

    psi_H-psi_Q
      =d A[psi(S_D+h_H)-psi(S_D+h_Q)]n(X)+O_L2(A^2).

The large linear terminal contribution cancels in the two Gaussian responses because both coherent coefficients are exactly a2. All remaining actual-mean mismatch is already within the O(A^2) bound. The order-one radial inflation is kept inside psi; a Euclidean-displacement Taylor expansion would not justify this step.

## 5. Actual C2 and C1 covariance corrections have the same leading tangent

For a standard conditional Gaussian vector N and arbitrary possibly correlated R with |R|<=B, componentwise first-chaos Bessel gives

    ||Cov(N,R)||_HS^2<=sum_l Var(R_l)<=B^2,
    ||Cov(R)||_HS<=tr Cov(R)<=B^2.

This is a specific orthogonal coefficient estimate. It does not infer a general vector-trace bound from a gradient port. Applying it to b_a N_a+R_a gives

    C2_a=b_a^2 I+B2_a, ||B2_a||_HS<=2 b_a B+B^2=O(A^2).

The one-level versions similarly give C1_H=c^2 A^2 I/4+B1_H and C1_Q=c^2 A^2 beta_mid^2 I+B1_Q, with each HS remainder O(A^2). Thus the leading isotropic H-minus-Q coefficient is c^2 A^2(1/4-beta_mid^2) for either choice. Precisely, the difference between the C2 and C1 isotropic H-minus-Q coefficients is O(A^3); the H-minus-Q variance gap itself is order A^2. This specifies the intended reading of “difference between their isotropic coefficients” in Section 5 of the source.

Every authorized coherent center has mu_*=a2 x+M_*(x), |M_*|<=B. Therefore z_*=x-mu_* satisfies |z_*|-R0=S_D+O_L2(A), with a still smaller direction error on the typical shell. This remains true for arbitrary measurable pointwise convex combinations, uniformly in the mixing function.

The exact radial contraction formula is

    D2 f(z):B = psi'' n(n^T B n)
             +(psi'/r-psi/r^2)[n tr(PB)+P(B+B^T)n].

For all matrices, its HS-to-vector operator norm is

    max{sqrt((psi'')^2+(D-1)(psi'/r-psi/r^2)^2),
        sqrt(2)|psi'/r-psi/r^2|}.

This is bounded independently of D because every nonzero denominator is at least R0-2, comparable to sqrt(D), and the scalar derivatives are fixed. Hence sup ||D2 g||_(HS->vector)<=C A, and the covariance B remainders contract to O(A^3).

For the isotropic part, one must use the actual trace:

    Delta g(z_*)=d A[psi''+(D-1)(psi'/r-psi/r^2)]n
                =d A R psi'(S_D)n(X)+O_L2(1).

The O_L2(1) error accounts for both the O(A) shift of the derivative and the O(A R) change of (D-1)/r; multiplying either by A gives O(A^2 R)=O(1). The psi/r^2 term and the other scalar terms are smaller. Therefore for ell=1 or 2,

    (C_ell,H-C_ell,Q):D2g(z_*)/2
       =d A(h_H-h_Q)psi'(S_D)n(X)+O_L2(A^2).

The exact C2 correction changes the leading C1 correction only by O(A^2). All errors are uniform in the stated centers, so integrating the centered current in theta is legitimate and does not change the conclusion. No covariance or D2g producer has been licensed by this analytical calculation.

## 6. The nonzero witness and finite outer quadrature

After subtracting the tangent, the residual is d A q_beta(S_D)n(X)+O_L2(A^2), with beta=beta_mid and

    q_beta(s)=psi(s+h_H)-psi(s+h_beta)
                  -(h_H-h_beta)psi'(s).

For S~N(0,1/2), m(t)=E psi'(S+t) is strictly decreasing on t>0. To see this without a numerical sign inference, write the symmetric bump as a layer-cake mixture of centered intervals. Every nondegenerate centered Gaussian interval mass strictly decreases under positive translation. Thus

    E q_beta(S)=int_(h_H)^(h_beta)[m(0)-m(t)]dt
               >=kappa>0

uniformly for beta in [7/12,1]. The family q_beta is bounded and uniformly Lipschitz in both s and beta. A finite beta net plus the Gaussian shell CLT gives uniform convergence of E q_beta(S_D); since |X|/R tends to one in L2, the same holds with the radial weight |X|/R. Eventually its infimum is at least kappa/2.

The vector test T_D=X/R has L2 norm one, not pointwise norm one. R1 is self-adjoint and R1 T_D=T_D/2. The post-resolvent discrepancy is consequently bounded below by d kappa A/4-C A^2.

For the literal outer source, G,H,J are independent at the level of their complete blocks, but shared across their respective node levels. At each outer node x_i, conditioning on x_i and integrating the independent H,J gives psi_Q(x_i). Linearity of expectation then gives the exact OWN-mean identity E[F3_Q|Z]=Q_out psi_Q(Z). Sharing roots across outer nodes does not alter it.

Anchoring and the exact finite decomposition imply ||psi_Q||2<=C A sqrt(D)=C/A. Therefore

    ||(Q_out-R1)psi_Q||2<=C delta_out/A.

Choosing delta_out<=A^3 places this entire finite top-clock floor at O(A^2). Combining it with the witness proves the displayed eventual c_* A lower bound with c_*=d kappa/8. The statement does not hide an exact outer operator in the literal F3_Q target.

The counterterm integrand has L2 norm O(A) by the explicit trace calculation. Replacing its outer R1 by Q_ctr thus adds at most C delta_ctr A; delta_ctr<=A is sufficient for an O(A^2) floor. The same outer A^3 rule is stronger than necessary. Also |mu_H-mu_Q|<=2B pointwise, so an exact coherent-mean term bounded by A||mu_H-mu_Q||2 contributes at most 2 A B=O(A^2), too small to remove the separator.

## 7. Public logs, compiler scope, and the positive comparison

To refute the advertised fixed-public-log theorem, choose an explicit admitted sequence A=D^(-1/4), delta_mid=A^2, delta_out=A^3, a polynomially tight inner tolerance, and polynomially small restored absolute floors. Positive dyadic rules then use O(log^2(1/A)) nodes per level, with their actual moment/multiplier and encoding allowances retained. Keep the remaining public versions bounded or polynomial along this sequence. Consequently Lambda is a fixed polynomial in log(1/A), and Lambda A tends to zero. The ratio c_* A/(Lambda A^2) diverges.

This is not a claim about arbitrarily exponentially over-resolved auxiliary versions that make public logs huge. One admissible bounded-version sequence is enough to refute a dimension-uniform assertion. The analytical c_* A separator is uniform over the rules satisfying the corollary's stated contracts.

The imported raw OWN-mean compiler is only used in its declared scope: if its conditional mean-law output is within O(Lambda A^2) plus small floors of the corrected wrong target, that allowance cannot bridge the Omega(A) conditional target error. Conditional W2 controls conditional mean separation. This observation does not independently re-prove the compiler, construct the analytical covariance correction, or imply a full reverse-OU endpoint join.

The bounded-remainder consumer is fully compatible with the negative result. Its true-backbone target E_N g((1-a2)x-b_H N) differs from the actual true force mean pointwise by at most A epsilon(1+lambda)=O(A^2). It retains the complete psi(S_D+h_H) response rather than its tangent. Thus this same counterfamily provides a constructive non-obstruction to a fully resummed Gaussian VALUE response. A general source without a bounded O(A) linear remainder is not covered by that comparison, and the general resummed consumer remains outside this audit's result.

## 8. Reproducibility and limitations

Run `python check_centered_p3_inflation.py` in this directory. The completed run has **473 passing assertions**. It checks all five exact input hashes; the stationary and Brownian covariance derivations; independent finite level-root algebra; explicit positive finite rules with distinct middle/inner tolerances; literal two-level same-g decomposition and full nonlinear remainder bounds; Hessian/HS-operator diagnostics; floor exponents; both signs/formulas for the scalar witness; and exact marginal high-dimensional radial geometry without a path grid.

Selected diagnostics:

- largest literal two-level decomposition residual: 7.100e-16;
- largest observed terminal proxy error divided by A epsilon(1+lambda): 0.250013 or less;
- largest sampled HS-to-vector norm of D2f: 0.468247 or less;
- uniform limiting kappa: 1.8145100150935224e-6;
- corresponding c_*: 5.6703437971672573e-8.

Strict positivity and the eventual lower bound are analytical. The numerical integration error estimates are diagnostics, not interval certificates. The checker does not determine the sufficiently-large-D threshold, simulate the genuine OU history, prove the source/compiler theorem by sampling, execute a covariance producer, or establish a generic impossibility result. The small positive c_* means no practical finite-D threshold should be inferred from these diagnostics alone.

Only this independent audit directory was written. All frozen source files remain unchanged. `MANIFEST.json` pins the inputs and completed outputs; `SHA256SUMS` seals the local audit packet.
