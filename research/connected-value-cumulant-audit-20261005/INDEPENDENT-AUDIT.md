# Independent audit: finite-VALUE connected cumulant boundary

Date: 2026-10-05

## Verdict

**PASS, scoped to the executed coefficient-tensor class and the hypotheses stated below.** The exact replica formula, constants c4 = 10/3 and c8 = 46212/35, all-distinct-coordinate energy lower bound, rank/nuclear lower bound, admissible smooth convex OU witness, and independent-probe isometry are sound. No active native-source construction or all-source impossibility is proved.

One clarification is important enough to be part of this verdict: the caller-integrated approximate-unbiasedness extension must control the **conditional bias as a function of the caller**, for example in L2(caller; HS). Accuracy of the tensor only after averaging out the caller does not imply the bound with m = E|c(Y)|. The wording “integrated HS bias” in Section 3.2 should be replaced by an explicit norm and conditioning statement.

Other required qualifications are: a deterministic pathwise rank cap M; fixed rank n; dimension-independent witness and energy constants for dimension-asymptotic conclusions; positive-variance finite clocks; and dimension-fixed geometry/coefficient tensor for the sharper D^n claim.

The source package was not modified. Its manifest still records that an independent audit had not yet been performed; this separate packet is the independent audit record.

## 1. Integrity and reproducibility

Audited source:

- /workspace/shared/connected-value-cumulant-source-20261005/FINITE-VALUE-CONNECTED-CUMULANT-BOUNDARY.md
- SHA256: 7deff8b0096c81789cef695ab361c206d96f5a21d89abcd4ae5dac44b238260e

All six entries in the source SHA256SUMS verify. All three absolute-path input pins verify. The original checker was also copied to a temporary directory and executed there; its newly generated checks.json is byte-identical to the source checks.json. The source files were left unchanged.

The independent checker uses inclusion-exclusion counts of surjections and occupancy profiles rather than the source checker’s set-partition enumeration. It verifies all formal scalar moment monomials through rank eight against log-MGF coefficients; it also checks an exact noncentered three-point distribution, centered rank-four monomials, exact rank-bound equality examples, an explicit epsilon certificate, and rational finite-OU covariance examples. It uses only the Python standard library.

Neither numerical diagnostics nor a successful script prove the analytical Gaussian identities, infinite-time integrability, or nuclear-norm theorem. Those are justified below. Imported Appell/positive-consumer and ancestry results are pinned context, not independently recertified in their entirety by this audit. The exclusion theorem itself does not need them.

## 2. Exact replica kernel and D^n energy

For a partition pi into k blocks, every injection of its blocks into n replica labels has the same expectation, namely the product of the k relevant moments. Dividing by (n)_k cancels exactly the number of injections. The Mobius coefficient (-1)^(k-1)(k-1)! therefore gives the moment-cumulant formula. This works conditionally at every fixed caller where the needed moments exist.

For the finite clock H_Q this is a finite-query construction. For the continuous OU integral H it is an analytical identity involving complete continuous histories; it is not itself a finite original-VALUE algorithm. The note correctly separates this issue from quadrature error.

For an all-distinct physical tuple, a slot-to-replica map a uniquely determines its partition into equal-label fibers. Its grouped coefficient is exactly mu_k/(n)_k. Under coordinate independence and centering, different maps are orthogonal: a differing replica label in at least one physical coordinate contributes a zero cross expectation, and physical coordinates factor independently. No independence of the different partition terms is assumed.

There are choose(n,k) times the number of surjections from n labelled slots onto k labels such maps, which equals S(n,k)(n)_k. Hence the exact component second moment is

    sigma^(2n) sum_k S(n,k) ((k-1)!)^2/(n)_k.

The independent exact calculation confirms 10/3 at rank four and 46212/35 at rank eight. The centered two-copy rank-four kernel has exactly eight orthogonal monomials with coefficient magnitude 1/2, and therefore second moment 2 sigma^8.

Summing only the (D)_n all-distinct ordered tuples proves the displayed HS lower bounds. For fixed n and fixed positive sigma this is Omega(D^n) as D tends to infinity; (D)_n is zero for D<n. These are lower bounds, not proofs of a uniform matching upper bound. For the linear Gaussian example the true high-order cumulant is zero, so the replica energy is estimator variance. That example invalidates this replica producer’s energy certificate; by itself it does not exclude an algorithm that recognizes the linear case and returns zero.

The OU calculations are exact. Stochastic Fubini gives

    Y = (1/sqrt(2)) integral_0^infinity exp(-s) dW_s,
    Var(Y) = 1/4,
    Cov(Y,X_t) = t exp(-t).

The integrands are square integrable, which justifies the exchange. For a positive finite rule, positive variance requires at least one strictly positive-time node with nonzero weight; a rule supported only at t=0 is degenerate and is excluded from this conclusion.

## 3. The rank theorem and its exact assumptions

Let A(T) denote the 1|(n-1) flattening. Assume n>=2, 1<=M<infinity, and rank A(T)<=M almost surely, with the same deterministic M for every outcome. The energy may be infinite, in which case the lower bound is automatic. For finite energy, all expectations below exist.

For K=c sum_i e_i^tensor n, its flattening is c times an isometric copy of the D-dimensional identity. Thus its singular values are exactly |c| repeated D times. Jensen, the rank inequality and Cauchy-Schwarz give

    D|c| <= E||T||_* <= sqrt(M) E||T||_HS
             <= sqrt(M E||T||_HS^2).

Consequently E||T||_HS^2 >= D^2 c^2/M. One can replace M by min(M,D), but the stated weaker bound is valid.

For B=ET-K, rank A(B)<=D, so ||B||_*<=sqrt(D)||B||_HS. If ||B||_HS<=eta sqrt(D), the same argument yields

    E||T||_HS^2 >= D^2 (|c|-eta)_+^2/M.

Pointwise ||T||_HS^2<=M||T||_op^2 proves the operator-energy bound in the note. These arguments are independent of coefficient magnitudes, signs, centering, differentiability, coupling, antithetic structure, or Gaussian geometry.

They are sharp at the level of this rank model: choose an M-element subset S of [D] uniformly and set

    T = (D/M)c sum_(i in S) e_i^tensor n.

It is unbiased and attains both lower bounds exactly. Replacing c by c-eta, for c>eta>=0, gives exact equality in the approximate-bias versions. The checker verifies these finite examples with rational arithmetic.

### What counts in M

A sufficient condition is the actual pointwise decomposition T=sum_(a=1)^M V_a tensor W_a. Each W_a may have arbitrary dependence on every value, root, scalar coefficient and clipping decision. Scalar data dependence therefore does not evade the theorem.

The theorem is not a lower bound on the number of oracle queries in every imaginable tensor program. It is a lower bound under this first-mode span/rank condition. Independent Gaussian first-slot vectors, transformed query vectors, or other newly introduced physical directions count if present. A known full-rank identity/control tensor can violate the condition. Symmetrizing V tensor W when W has additional first-slot directions can increase rank; the rank of the executed symmetrized tensor must be checked, not assumed. Homogeneous outer-product kernels whose every slot comes from the same M-vector dictionary do meet the condition after symmetrization.

An adaptive random query count is not covered by replacing deterministic M with E[M] without another argument. The basic pointwise inequality would instead involve E[M(T)||T||_HS^2], or require an appropriate different bound. A deterministic upper cap suffices.

### Dimension and parameter quantifiers

To exclude a dimension-uniform source certificate, fix n, the witness epsilon, the scalar clock (if finite), and c0!=0 independently of D. If c=c0 A^n and

    (E||T||_HS^2)^(1/2) <= C_n A^n sqrt(D),
    eta <= |c|/2,

then M >= c0^2 D/(4 C_n^2). A bounded L2 operator-cut certificate similarly forces a linear rank cap. This rules out M=o(D) in the relevant parameter regime. “Public-log” growth only implies an exclusion when its actual dependence is sublinear in D. It is not enough to call a quantity a logarithm if its argument is permitted to grow exponentially with D.

If epsilon, the geometry, or a finite clock changes with D so that c0 tends to zero, the numerical dimension lower bound must retain that dependence. The uniform no-go conclusion is still obtained by selecting one fixed admissible nonzero witness against a purported uniform algorithm.

## 4. The admissible OU witness, with an explicit epsilon

The proposed potential and gradient are correct. At 0<epsilon<=1/8,

    23/64 <= g'_epsilon(x) <= 41/64,

using the simple bound 1/2 +/- (epsilon+epsilon^2). Thus g is anchored, smooth, and the gradient of a strictly convex potential, and scaling by A gives the asserted Hessian sandwich.

The additive history is exactly H_epsilon=aY+epsilon S+epsilon^2 C, where a=1/2, |S|<=1, and |C|<=2. All moments exist. Multilinearity makes every order-n cumulant a finite polynomial in epsilon. The zero-order term vanishes for n>=3 because aY is Gaussian.

For jointly Gaussian (Y,X), the joint exponential tilt shows

    kappa(Y,...,Y,f(X)) = Cov(Y,X)^(n-1) E f^(n-1)(X).

One can prove this without assuming an MGF for a general f by polynomial/moment approximation; here bounded sine and its derivatives make the tilt and differentiations immediate. Integrating the final sine slot is justified by boundedness and the integrable exponential clock. Therefore

    kappa'_n(0) = n a^(n-1) (-1)^(n/2-1) I_n,
    exp(-1/2)(n-1)!/n^n <= I_n <= (n-1)!/n^n.

At both n=4 and n=8 the derivative is strictly negative. Since there are only two ranks under consideration, one common fixed sufficiently small positive epsilon works. The text is correct not to claim this for every epsilon<=1/8.

There is also a completely explicit, conservative choice:

    epsilon = 2^(-80).

Here is an arithmetic certificate. For n in {4,8}, each of aY, S, C has L^n norm <=2. Every mixed cumulant of n of these variables has absolute value at most n^n (n-1)! 2^n: apply Holder in each partition block and bound the number of partitions by n^n. The expansion has at most 3^n terms. Thus

    |kappa_n(H_epsilon)-epsilon kappa'_n(0)|
        <= R_n epsilon^2,  R_n=(6n)^n (n-1)!.

Since exp(-1/2)>1/2,

    |kappa'_n(0)| > L_n,  L_n=n!/(2^n n^n).

The exact rational checker verifies R_n epsilon<L_n/2 for both ranks. Hence kappa_n(H_epsilon)<-epsilon L_n/2<0 simultaneously. The very small constant is sufficient for an asymptotic dimension obstruction; it is not a claim of a practically large finite-D gap.

Independent coordinate futures give the exact diagonal joint cumulant tensor. Mixed-coordinate cumulants vanish by conditional independence, not by any approximation.

For finite clocks the conditional OU covariance is

    Cov(X_s,X_t)=exp(-|t-s|)-exp(-(s+t)),

strictly positive when s,t>0. Therefore the displayed finite first-derivative sum is positive for a nondegenerate positive clock. The sufficiently small epsilon can depend on that fixed clock. The explicit 2^(-80) certificate above is for the continuous history; it is not asserted uniformly for every possible finite clock.

## 5. Integrated callers: the precise valid extension

Let Y be a standard D-dimensional Gaussian caller, Q the conditional coefficient randomness, and

    K(Y)=sum_i c(Y_i)e_i^tensor n,
    B(Y)=E[T|Y]-K(Y).

For this witness c(y) is continuous: the histories vary continuously in L^p with y, and cumulants are polynomials in their moments. Their centered moments are uniformly bounded in y because the random linear part has fixed Gaussian variance and the nonlinear additions are bounded. Hence m=E|c(Y_1)| is finite; c(0)!=0 and continuity give m>0.

Conditional unbiasedness gives

    E||T||_HS^2 >= E(sum_i |c(Y_i)|)^2/M
        = [D^2 m^2 + D(E c(Y_1)^2-m^2)]/M
        >= D^2 m^2/M.

If instead

    (E||B(Y)||_HS^2)^(1/2) <= eta sqrt(D),

then E||B(Y)||_*<=eta D. After averaging the nuclear triangle inequality,

    D(m-eta)_+ <= E||E[T|Y]||_*
        <= E||T||_* <= sqrt(M E||T||_HS^2).

This proves the claimed approximate extension with the explicit conditional-bias interpretation. In fact an L1 bound E||B(Y)||_HS<=eta sqrt(D) is enough for this averaged argument.

It does not follow from ||E_(Y,Q) T-E_Y K(Y)||_HS<=eta sqrt(D). For example, the abstract diagonal family c(y)=y^2-1 has c(0)=-1, m>0, and E c(Y_1)=0. The estimator T=0 matches the fully averaged target exactly while having zero energy. This is a counterexample to the weaker interpretation of the extension, not a claim that this particular c is generated by the note’s OU witness. For a fully averaged target the direct rank theorem instead uses |E c(Y_1)|, when nonzero.

## 6. Fixed-geometry Gram claim

Under a shared fixed scalar-row geometry L_a, coordinatewise g, independent physical Gaussian coordinate banks, and deterministic fully grouped C, the all-distinct component second moment is exactly the quadratic form C^T G^tensor n C. Nonzero means cause no difficulty because G is the raw second-moment matrix, not a covariance matrix.

For g(x)=a x+b sin x+c(1-cos x), b c !=0, a zero L2 relation among values is an analytic identity on the Gaussian support. A polynomial plus distinct nonzero Fourier exponentials has a unique representation. The sine and cosine coefficients separately eliminate each antipodal pair. Thus distinct nonzero row vectors give positive definite G; duplicate/zero rows must first be removed and the coefficients regrouped.

This proves a strictly positive all-distinct coefficient for nonzero reduced C. It does not give a lower spectral bound uniform in near-coincident rows, large/small widths, dimension-dependent rows, or dimension-dependent C. The note explicitly acknowledges this. Data-dependent C is outside this factorized Gram calculation, though it remains inside the rank theorem if the pointwise span condition holds. The D^n asymptotic is therefore a fixed-geometry, fixed-coefficient result.

## 7. Probe contractions and the native-source boundary

For probes independent of T, iterated Gaussian second moments give E_P|T[P_1,...,P_(n-1)]|^2=||T||_HS^2 exactly. For one public probe with Wick order n-1, the exact factor is (n-1)! times the squared norm of the coefficient symmetrized in its contracted slots. It equals that factorial times ||T||_HS^2 only when those slots are already symmetric. Without symmetry, only the symmetric projection is read; one may apply the rank theorem to that projection if it has the target and keeps the rank cap. This is consistent with the note’s symmetric-slot qualification.

Independence is essential. If query sites depend on physical probes, there need not be a coefficient tensor T(Q) independent of the probes. Pointwise low-dimensional vector span does not bound the rank of the tensor obtained by native Gaussian integration.

A simple illustration is F(P)=g_epsilon(P) coordinatewise for a single standard Gaussian physical vector P. Its native Wick response at even n>=4 is

    E[F(P) tensor H_(n-1)(P)]
      = (-1)^(n/2-1) epsilon exp(-1/2) sum_i e_i^tensor n.

This follows from Gaussian integration by parts and coordinate independence; the linear and cosine terms contribute zero at these odd derivative orders. The response has flattening rank D despite only one g evaluation at the active site P. Also E|F(P)|^2<=D and its derivative is bounded. This example is not the desired continuous-OU cumulant compiler, nor a verification of every native bounded-field/caller/filter contract. It simply demonstrates why carrying the rank bound through probe-dependent native integration would be invalid.

The source note correctly excludes such active programs from its impossibility statement. Known high-rank control tensors and native pair/filter responses likewise need separate admissibility and energy proofs; the present audit neither certifies nor rules them out. Exact cumulant unbiasedness alone supplies no low-first, all-cuts, positive-consumer, cutoff-uniform, or all-order cost theorem.

## 8. Suggested source clarifications

1. Section 3.2: write B(Y)=E[T|Y]-K(Y) and specify its L2(Y;HS) bias allowance explicitly.
2. Section 3: state that M is a deterministic almost-sure rank cap; distinguish this from an expected adaptive query count.
3. Section 2.1: qualify the finite positive rule by at least one positive-time positive-weight node.
4. Sections 3.1 and 5: preserve all witness/clock/geometry/coefficient dependence when taking D to infinity; formulate the exclusion as M=o(D).
5. Section 7: continue treating numerical checks as diagnostics. The independent script’s PASS_SCOPED is not a certificate for excluded native constructions or imported all-order results.

These clarifications do not overturn the central executed-tensor obstruction. They prevent an unjustified extension of it.
