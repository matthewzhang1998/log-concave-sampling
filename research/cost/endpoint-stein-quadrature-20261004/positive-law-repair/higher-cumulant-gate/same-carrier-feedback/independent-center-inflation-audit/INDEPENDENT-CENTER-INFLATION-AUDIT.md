# Independent audit: centered covariance and radial inflation

2026-10-05. Independent mathematical review of the frozen same-original-gradient construction.

## Verdict

**PASS for the stated centered single-history covariance separator and for the family-specific full shifted Gaussian VALUE proxy, with the public-log qualification in Section 9 below.**

The frozen construction gives a genuine smooth convex-gradient family, with the same original `g` at every inner and terminal site, for which

- `A=D^(-1/4)` and `R=sqrt(D)=A^(-2)`;
- the discrepancy between `R1(j_H-j_Q)` and the covariance derivative frozen at any specified coherent center has L2 norm at least `c_* A` for all sufficiently large D;
- the advertised fourth-order allowance is `A^4 sqrt(D)=A^2`;
- the same conclusion holds for the theta-averaged coherent-center current in the exact centered Stein identity;
- the exact identity's mean term is at most `d A^4` under `delta<=A^2` and cannot absorb this discrepancy.

This is not a failure of the exact live shifted current, a full mean-closure impossibility theorem, an endpoint-law counterexample, or a general consumer construction. The added shifted-Gaussian proxy is valid precisely because it retains the complete original-g response to the order-one radial inflation.

## 1. Frozen inputs and reproducibility

The reviewed author note is:

`RADIAL-INFLATION-REFUTES-CENTERED-SINGLE-HISTORY-COVARIANCE.md`

SHA-256:

`c6024a6ba2ca4fbb293063d0746ee92fa083398a3bf4c01fb4de1908f322dfaf`

It includes Section 9, the constructive family-specific proxy. The author supplied this final hash and the reviewer recomputed it independently. The exact centered-current note, positive-clock construction, author diagnostic script, and author diagnostic output are also pinned in the independent checker and manifest. The author script was neither imported nor executed by this checker. No author file was modified.

Run `python check_center_inflation_independent.py` from this directory. The independent script passes 723 assertions. It verifies exact input pins, an independently reconstructed finite positive dyadic rule, selected Hermite multipliers, the original-g source and terminal identities, radial derivative formulas, the shifted VALUE coupling bound, true OU linear normalization, powers and mean floors, two independently written scalar witness integrals, and exact-marginal radial diagnostics with large Euclidean noise.

The analytical proof below, rather than floating-point sampling, establishes uniformity and strict sign. The script does not certify a finite-D onset threshold, simulate an OU history, or provide a generic compiler.

## 2. Global convex-gradient admissibility

Put `c=1/2`, `d=1/4`, `k=1-cA/2`, and `R0=kR`. The bump `b=psi'` is nonnegative, even, smooth, compactly supported on `[-2,2]`, and strictly decreasing on `(0,2)`. A simple analytic bound verifies the asserted maximum without relying on the author diagnostic:

`Z >= integral_(-1)^1 exp(-4/3) ds = 2 exp(-4/3)`,

hence

`sup b = exp(-1)/Z <= exp(1/3)/2 < 1`.

For `D>=256`, `R0-2>R/2`. Thus `f(x)=psi(|x|-R0)x/|x|` vanishes on an open neighborhood of the origin. Its zero extension is smooth there and across the shell boundaries because the bump is flat at its support endpoints. Its radial potential is the stated integral of psi, whose derivative is nonnegative and nondecreasing.

The two Hessian eigenvalues of the original potential are exactly

`cA+dA b(r-R0)` and `cA+dA psi(r-R0)/r`.

They are nonnegative; the radial one is at most `(c+d)A<A`, and the tangential one is bounded by `(c+d/(R0-2))A<A` wherever the shell term is nonzero. The origin has derivative `cA I`. Therefore `g(0)=0` and `0<=Dg<=A I` globally. This is a genuine C-infinity member of the admitted C2 class. If an inherited theorem insists on `A<=1/8`, restrict to `D>=4096`; the asymptotic family is unchanged.

## 3. True history, finite common root, and retained correlations

For the true conditional OU process,

`X_t^x=e^(-t)x+sqrt(2) integral_0^t e^(-(t-s)) dB_s`.

Its weighted linear integral has mean `x/2`. Stochastic Fubini gives its noise as

`(1/sqrt(2)) integral_0^infinity e^(-s) dB_s`,

whose covariance is `I/4`. Thus the author's `sigma_H=1/2` is exact. This is the actual history law, not the law of a finite common-root surrogate.

For the actual positive finite clock rule, exact first moment gives mean `x/2` for its linear component, and its one shared root has coefficient `beta_Q=sum v_j sqrt(1-tau_j^2)`. Its linear conditional covariance is `beta_Q^2 I`. Cross-node terms are present. Replacing this with `sum v_j^2(1-tau_j^2)I` would be a different source and is not done in the proof.

Writing `K_a` for the positive weighted average of the same original radial field gives the literal identity

`I_a=cA(x/2+sigma_a G_a)+dA K_a`, with `|K_a|<=1` pathwise.

The Gaussian `G_a` is independent of x, while `K_a` can depend on that Gaussian and, in the history case, every other history coordinate. All estimates used later permit this dependence. The same `g=cA id+dA f` is evaluated at the terminal argument.

Mean quadrature applies to the bounded original field f itself. Componentwise Hermite diagonalization yields the vector Hilbert-space estimate with no extra dimension factor. Hence

`||M_H-M_Q||2 <= delta` and `||mu_H-mu_Q||2 <= dA delta`.

The generic mean estimate for an arbitrary g would be weaker here; using bounded f is justified by the exact first-moment cancellation of the linear component.

For any admitted Q, the degree-two multiplier implies `sum v_j tau_j^2<=1/3+delta`. The scalar inequality `sqrt(1-tau^2)>=1-tau^2` then gives `beta_Q>=2/3-delta>=7/12` on the stated family. Consequently `beta_Q^2-1/4>=13/144`. Every later constant is independent of node count, clock locations, and Q. Existence of admissible finite polylogarithmic-node rules is verified independently in Section 9.

## 4. Radial asymptotics with unbounded Euclidean displacement

The relevant displacement has Lp norm of order `A R=A^(-1)`, so a Euclidean Taylor expansion about x with a small absolute displacement is not justified. The author's radius/direction calculation does not make that assumption.

For `rho=k|X|`, `n=X/|X|`, and `w=-cA sigma G`, the good event is `|X|>=R/2`, `|w|<=rho/2`. On it,

`|rho n+w|-rho = n.w + (|w|^2-(n.w)^2)/(2rho) + O(|w|^3/rho^2)`.

For each fixed p, uniformly over `sigma in [0,1]`, the terms have the following Lp sizes:

- `n.w=O(A)`;
- `|w|=O(A R)`;
- the cubic remainder is `O(A^3 R)=O(A)`;
- the centered transverse chi-square numerator is `O(A^2 R)=O(1)` and becomes `O(A^2)` after division by rho;
- its conditional mean divided by `2rho` is `c^2 sigma^2/2+O(A)`;
- the direction changes by `O(A)`;
- `rho-R0=k(|X|-R)=S_D+O_Lp(A)`.

The bad event has exponentially small probability. Since f and psi are bounded, their bad-event L2 contributions are negligible directly; polynomial moments control the intermediate expressions if those are retained. The Gaussian shell moments used here are dimension-uniform. Therefore

`f(kX-cA sigma G)=psi(S_D+c^2 sigma^2/2)n(X)+O_L2(A)`

uniformly in sigma. The constant inflation is not lost merely because individual projected noise coordinates are small.

For the actual terminal source, subtracting `dA K_a` changes f by at most `C dA` pathwise, independently of any `K_a/G_a` correlation. Multiplying by the outer `dA` costs `O(A^2)`. The terminal linear part contributes only `-cA(mu_H-mu_Q)`, of size `O(A^2 delta)`. This proves the exact-source leading difference in the author note.

## 5. Conditional covariance remainder and coherent centers

Let `L_a=K_a-M_a`. At each fixed x, the coordinates of `G_a` are orthonormal in conditional L2. Bessel's inequality applied separately to each component of L gives

`||E[L_a G_a^T | x]||_HS^2 <= E[|L_a|^2 | x] <= 1`.

Also `Cov(K_a|x)` is positive semidefinite with trace at most one, so its HS norm is at most one. Thus the literal covariance expansion yields

`C_a=c^2 A^2 sigma_a^2 I+B_a`,

with the stronger explicit bound

`||B_a||_HS <= (2cd sigma_a+d^2)A^2`.

This is a valid first-chaos estimate. It neither pretends that an operator norm is an HS norm nor infers a vector-trace bound from generic energy. It retains the complete correlation of K with G. In particular, `||B_H-B_Q||_HS<=2(2cd+d^2)A^2` uniformly in x and Q.

For any measurable center `mu_*=cAx/2+dA M_*` with `|M_*|<=1`, the argument `z_*=x-mu_*` satisfies

`|z_*|-R0=S_D+O_L2(A)` and `n(z_*)=n(X)+O_L2(A/R)`.

These estimates require no derivative of M and no independence. They therefore cover `mu_H`, `mu_Q`, arbitrary pointwise convex combinations, and the theta-averaged current, with one common constant.

The exact radial contraction is correct for arbitrary matrices B:

`D2f:B=b'(r-R0)n(n^T B n)+a_r[n tr(PB)+P(B+B^T)n]`,

where `a_r=b(r-R0)/r-psi(r-R0)/r^2`. In coordinates with `n=e_1`, the output rows are mutually orthogonal linear functionals of B. The exact HS-to-vector norm is

`max(sqrt(b'^2+(D-1)a_r^2), sqrt(2)|a_r|)`.

On the support of b, `r>=R0-2>R/2`, while outside that support the only remaining term is `-1/r^2` on the exterior side. This proves a global dimension-independent bound for D2f and a global `C A` bound for D2g. The B remainder therefore contracts to `O(A^3)`, as stated.

The isotropic covariance part must be treated separately. Direct calculation gives

`Delta g(z_*)=d A R b(S_D)n(X)+O_L2(1)`.

The `1/k` radial factor and the change in the bump both have relative size `O(A)` and thus contribute the displayed `O(1)` error after the `A R=A^(-1)` amplification. Direction, b', and `psi/r^2` terms are smaller. Multiplying by the `A^2` covariance factor gives an `O(A^2)` remainder. The resulting frozen covariance term is exactly

`dA(h_H-h_Q)b(S_D)n(X)+O_L2(A^2)`.

Uniformity in all admissible centers permits integration over theta by Minkowski's inequality. No theta-endpoint singularity enters this deterministic term.

## 6. Strict witness and survival of R1

After subtraction, the leading scalar profile is

`q_beta(s)=psi(s+h_H)-psi(s+h_beta)-(h_H-h_beta)b(s)`.

The subtraction sign is correct. With `S~N(0,1/2)` and `m(t)=E b(S+t)`,

`E q_beta(S)=integral_(h_H)^(h_beta) [m(0)-m(t)] dt`.

Strict positivity has a direct analytic proof. Since b is even and strictly decreasing on `(0,2)`,

`b(u)=integral_|u|^2 [-b'(a)] da` on its support.

Thus m is a positive mixture of centered Gaussian interval masses translated by t. For `a,t>0`, the derivative of that interval mass is `phi(a+t)-phi(a-t)<0`, where `phi(s)=exp(-s^2)/sqrt(pi)`. Hence m is strictly decreasing on `(0,infinity)`. Because `h_beta>=h_min>h_H>0`, the infimum witness is the positive constant

`kappa=integral_(h_H)^(h_min) [m(0)-m(t)] dt`.

The independent diagnostic obtains `kappa=1.8145100150935224e-6`; this numerical size is not the proof of positivity.

The shell CLT `S_D=>N(0,1/2)` is sufficient because the q family is uniformly bounded and uniformly Lipschitz in both s and beta on a compact beta interval. A finite beta net yields uniform convergence. Multiplication by `|X|/R` preserves that convergence since it tends to one in L2 and q is bounded. Thus the author obtains a uniform eventual lower bound `kappa/2`.

The test `T_D=X/R` has L2 norm exactly one. Gaussian Mehler self-adjointness gives `R1 T_D=T_D/2`. Pairing the resolvent discrepancy against this explicit first-chaos test gives

`d kappa A/4-C A^2`,

and therefore at least `c_* A`, with `c_*=d kappa/8>0`, for all sufficiently large D. The checker reports `c_*` approximately `5.6703438e-8`. No explicit finite-D onset is asserted: the positive coefficient is small and generic O(A^2) errors can dominate at practical dimensions. That does not affect the asymptotic dimension-uniform counterexample.

## 7. Exact Stein identity and floor bookkeeping

The example is smooth, so the deterministic and live shifted currents have ordinary meanings; no unproved separation of weak C2 limits is needed for this counterexample. The exact centered identity retains the derivative at the random shifted terminal argument and its complete correlation with the actual Stein matrices.

Freezing the derivative at `x-mu_theta` and averaging those matrices gives the precise theta-averaged deterministic covariance term refuted here. Its missing residual is not made small by coherent centering. The exact mean term obeys

`A ||mu_H-mu_Q||2 <= d A^2 delta <= d A^4`.

At `D=A^(-4)`, the generic mean-quadrature allowance `delta A^2 sqrt(D)` would equal delta and is at most `A^2` under the same choice. The improved family-specific bound is legitimate and strictly smaller. Neither bound can mask an order-A error.

This is an ideal analytical comparison. It does not establish new finite numerical, caller, mode, anchoring, or compiler ports. If the proposed theorem includes separately paid additive floors, choose an admitted sequence on which their aggregate effect is `O(A^2)` or smaller; triangle inequality then preserves the separator. Larger unrelated error floors can of course hide any analytical discrepancy, but cannot support a theorem claiming the advertised fourth-order precision. No numerical implementation of a covariance or derivative oracle is assumed.

## 8. Constructive appendix: precise positive result

For the explicit family, couple the true argument and Gaussian proxy using the very same `G_a` from the exact decomposition. The original-g Lipschitz bound gives

`|g(kx-cA sigma_a G_a-dA K_a)-g(kx-cA sigma_a G_a)|<=d A^2`

pointwise. Conditional expectation and R1 contraction imply the stated two-branch error at most `2d A^2`. This proof preserves correlations and does not assert equality of source laws.

The raw comparison can indeed use two evaluations of the same original g and affine Gaussian arithmetic with the explicit family constants; the two branch noises may be shared. The constants `sigma_H=1/2` and `sigma_Q=beta_Q` are available from this explicit decomposition and finite rule. They are not presented as general conditional-covariance or mean oracles.

The two-VALUE count is for this raw shifted response at captured x. It is not the cost of computing its expectation, applying an executed outer R1, compiling a positive law, furnishing all first/caller ports, or closing a full m3 target. Within this exact scope, Section 9 is a valid constructive non-obstruction: the complete shifted original-g response retains what a single frozen covariance derivative discards.

## 9. Finite-clock and public-log qualification

To make the fixed-public-log refutation unambiguous, take an explicit sequence such as `D=2^(4n)`, `A=2^(-n)`, declared `delta=A^2`, polynomially small numerical floors, and fixed/bounded remaining public versions satisfying inherited guards.

The admitted dyadic positive rule uses Gauss-Legendre order m on K dyadic panels in `1-tau`, plus one terminal midpoint. Constant and linear exactness give the required normalization. The ellipse estimate in the pinned clock note provides the all-degree certificate

`sup_l |sum v_j tau_j^l-1/(l+1)| <= 8*4^(-m)+2^(1-K)`.

Choosing `K=ceil(log2(4/delta))`, `m=ceil(log_4(16/delta))` makes that certificate at most delta with `Km+1=O(log^2(1/A))` positive interior nodes. This is an analytical all-degree guarantee; testing selected multipliers in the checker is only an implementation diagnostic. The finite rule's mean approximation is therefore available without a Markov path-grid producer.

On this admitted sequence all ordinary public logarithms are `O(polylog(1/A))`, and

`(c_* A)/(Lambda A^2)=c_*/(Lambda A)` diverges for each fixed polynomial in those logarithms. The power separation refutes a uniform theorem. It is not a claim that deliberately exponentially over-resolved tolerances or arbitrarily growing public versions still have logarithmic size in A. The lower bound itself is stronger than this selected sequence: it is uniform over every Q satisfying the note's accuracy and first-moment hypotheses.

## 10. Exact scope and stopping point

The independently passed conclusion is the failure of the deterministic covariance derivative evaluated at a coherent mean, including the exact identity's theta-averaged version, in an unrestricted dimension-uniform class. It does not address a different theorem imposing `A sqrt(D)=O(1)`. It does not automatically refute a specified higher-order truncation without analyzing that truncation. It does not refute the exact live shifted Stein identity or prove that every resummed/shifted finite VALUE method fails.

The family-specific Gaussian-shift proxy passes. General exact resummed-current consumption and any complete same-carrier m3 closure remain separate open tasks. No further mathematical claim is required for this audit.
