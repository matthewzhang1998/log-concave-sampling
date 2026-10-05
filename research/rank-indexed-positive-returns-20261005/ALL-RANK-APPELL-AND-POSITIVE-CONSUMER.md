# All-fixed-rank Appell bounds and a complete positive consumer generator

2026-10-05. Analytical results for centered Gaussian-bank Lipschitz sources. No tensor is an executable source. The original-VALUE and termination boundary is explicit in the companion note.

## 1. A one-mark cumulant inequality at every rank

Let B_1,...,B_m be centered vectors on the same standard Gaussian bank. Each B_j is C1/Sobolev with complete first at most L. Define the multilinear Appell tensor W_I, for a nonempty set of labelled occurrences I, by coefficient extraction in

  exp(sum_(j in I) t_j dot B_j - log E exp(sum_(j in I) t_j dot B_j)).

Extract one derivative in each labelled t_j at zero. Only this finite jet is used. Equivalently W_I is the finite sum over disjoint cumulant blocks removed from the tensor product of remaining B_j, with the signs/multiplicities obtained from the exponential. W_empty=1. E W_I=0 for nonempty I. Its actual Gaussian-bank derivative is exactly

  D W_I = sum_(j in I) (D B_j) tensor W_(I minus j),

with the displayed physical slot order preserved. All centering cumulants are deterministic while differentiating the bank. There is no derivative of any B_j beyond its first.

For any deterministic tensor T of HS norm one, apply Gaussian Poincare to <T,W_I>. The square norm of the sum of m gradients is at most m times the sum of their square norms. For each j, ||D B_j||op<=L leaves the vector whose j-th slot is open and whose other slots contract W_(I minus j). Slice T in its j-th slot. If m>=2, W_(I minus j) is centered, so its covariance-operator bound applies to each slice, and the squared slice norms sum to ||T||HS². Starting at m=1 with Gaussian Poincare gives

  ||Cov(W_I)||op <= (m!)² L^(2m).                         (A1)

The operator acts on the full physical tensor Hilbert space; no dimension occurs. For mixed Lipschitz bounds L_j the identical induction gives (m!)² product_j L_j².

Let R be any centered L2 vector on that same bank; R need not have a small derivative. Differentiating the joint log moment-generating finite jet once in the R variable gives the exact identity

  κ(R,B_1,...,B_m) = E[R tensor W_I].

It can instead be read as a finite moment-cumulant identity, so exponential moments of R are unnecessary. Joint covariance factorization and (A1) imply

  ||κ(R,B_1,...,B_m)||HS <= m! L^m ||R||2.                (A2)

For centered U,V with complete first at most L, multilinearity and telescoping n slots therefore give

  ||κ_n(U)-κ_n(V)||HS <= n! L^(n-1)||U-V||2.             (A3)

The constants at n=3 and n=4 are exactly 6 and 24. Conditional versions hold after fixing an actual retained caller, whenever the conditional bank has the stated Gaussian representation and first bounds. Integrating the scalar error energy over the retained caller preserves (A2)-(A3).

Applications to the source-qualified paths already in the sealed packets:

- Added original-coefficient heat g_tau=E g(.+tau Z) gives ||H_tau-H||2 <= A tau sqrt(D), and both conditional-bank first bounds <=A. Thus ||κ_n(H_tau)-κ_n(H)||_(L2(Y);HS) <= n! A^n tau sqrt(D).
- The coherent F2/H difference has energy C A²sqrt(D), first O(A) on both sides, hence the n-th conditional cumulant error is C_n A^(n+1)sqrt(D).
- A history difference Δ_j with energy C A^(j+1)sqrt(D), while all ordinary force factors have first C A, has one-mark mixed n-cumulant bound C_n A^(n+j)sqrt(D). Extra copies of Δ_j do not automatically acquire extra small powers from their energy alone.

These statements preserve the SAME joint future bank. They do not authorize independent resampling of coefficient occurrences.

## 2. Computable all-proper-cut bounds for actual cumulants

For a single centered B of first L, put c_2=1. For n>=3 define a finite numerical majorant recursively:

  c_n = max_(p+q=n, p,q>=1) [p! q! + sum_pi product_(C in pi) c_|C|],

where the sum is over partitions pi of the n labelled slots into at least two blocks, every block meeting BOTH the first p and last q slots. Every block then has size between 2 and n-2, so the definition is well founded.

To prove the bound, the finite generating identity

  E[exp(t.B-K(t)) exp(s.B-K(s))]
    = exp(K(t+s)-K(t)-K(s))

shows that E[W_p tensor W_q] is exactly the sum over these cross-block partitions, including the single block κ_n. By (A1) and covariance factorization, the p|q flattening of E[W_p tensor W_q] has operator norm <=p!q! L^n. Each proper cross-block product is, up to input/output permutations, the tensor product of proper flattenings of lower cumulants; its operator norm is bounded by the product of their c-values. Subtract them and apply the triangle inequality. Thus

  every proper cut of κ_n(B) <= c_n L^n.                  (A4)

This proof never controls a self-trace from an abstract cut bound. It uses the actual Gaussian-image cumulants and exact Appell subtraction. The existing centered Stein recurrence also gives

  ||κ_n(B)||HS <= (n-1)! L^(n-1)e, e=||B||2.             (A5)

For K_n=κ_n/n!, use h_n<=L^(n-1)e/n and k_n<=c_nL^n/n!.

## 3. A positive polynomial reference at arbitrary fixed rank

Fix m>=3, v>0, Σ=Cov(B), independent Gaussian roots Z1,Z2, and

  C_t=(v/2)I+(1-t²)Σ, x_t=C_t^(1/2)Z1, eta=sqrt(v/2),
  q_r,t(x)=K_r:H_(r-1)^(C_t)(x),
  p_t(x)=sum_(r=3)^m (1-t^r) q_r,t(x),
  W_t=tB+x_t+p_t(x_t)+eta Z2.

The covariance-Hermite score is defined by E[H_k^C(X)F(X)]=E[D^kF(X)], X~N(0,C). All W_t are genuine positive laws, with untouched independent eta Z2. The endpoint W0 is the positive reference; W1=B+sqrt(v)Z in law. Its polynomial map is an analytical reference, not an executable cumulant input.

Write H=I+D_xp_t and F_d=D_x^dp_t for d>=2. At d=1 use H. Let D_t=sum_r(1-t^r) partial_t[q_r,t(C_t^(1/2)Z1)] at fixed Z1. This includes BOTH changing covariance and changing visible root.

For every set partition pi of [r-1], define F_(r,pi) by contracting the r-1 input slots of K_r against one factor D_x^|C|(x+p_t(x)) for each block C of pi. The first output slot of K_r stays open, as does the output of every block factor. This is a rank 1+|pi| tensor. The all-singleton partition pi0 has F_(r,pi0)=K_r[H,...,H].

Hermite integration and the multivariate chain rule give the exact identity

  E[q_r,t(x_t).grad phi(W_t)]
    = sum_pi E[F_(r,pi):D^(1+|pi|)phi(W_t)].             (A6)

The sum runs over every set partition, including blocks of size r-1. No high derivative of p is discarded. This is a finite generator; its coefficients are the counts of labelled set partitions, without guessed multiplicities.

The centered Stein recurrence at the SAME W_t gives all residual currents:

  A1 = D_t,
  A2 includes -t Σ(D_xp_t)^T,
  for each r=3,...,m and partition pi:
     add -r t^(r-1) F_(r,pi) to A_(1+|pi|),
  for each r=3,...,m:
     add +r t^(r-1) K_r to A_r,
  A_(m+1) = t^m T_m, ||T_m||_(L2;HS)<=L^m e.

These instructions add to A2 rather than replacing it. Hence

  d E phi(W_t)/dt = sum_(ell=1)^(m+1) E[A_ell:D^ell phi(W_t)]. (A7)

The only cancellation is the constant K_r part of its all-singleton partition. The covariance sampler term, fixed-Z1 time derivative, all higher map derivatives, and every mixed-rank product remain. Every test is at the same W_t, and every coefficient is independent of Z2.

## 4. Every generated direct-consumer feedback has one energy

Normalize v=1. Put s=sum_(r=3)^m k_r and h=sum_(r=3)^m h_r. Derivatives of q_r are Hermite tensors of degree at most r-1-d, with coefficient transforms bounded by

  2^((r-1+d)/2) (r-1)!/(r-1-d)! <= 2^m m!.

Expand each H=I+Dp and each derivative F_d into its finite rank choices. Each monomial has one anchor K_r and at most r-1 additional K factors. The chain-rule edges connect every additional coefficient directly to the anchor. Crucially every coefficient retains a final PHYSICAL OUTPUT slot. Wick multiplication pairs only Gaussian slots of distinct factors. No initial factor has a self-contraction because each is individually Wick ordered.

Start the anchor in HS norm and add adjacent coefficient vertices one at a time, contracting ALL edges to the already built group at each attachment. The new coefficient has a nonempty contracted cut and at least its physical output uncontracted. Thus HS-times-operator multiplication supplies one anchor h_r and one proper-cut k for each new vertex. Added Wick edges do not alter this proof. Finally contract surviving Gaussian slots with their Hermite polynomial, paying sqrt(degree!). This is a proof for these generated graphs, not for graphs that consume a vertex's last physical mark.

The largest Gaussian degree before Wick expansion is at most (m-1)(m-2). All diagrams and combinatorial factors are finite. For example

  C_m = 2^(10m²) ((2m²)!)^4

is a conservative explicit common constant covering partitions, rank choices, derivative factors, Wick pairings, Hermite norms, rank-to-velocity factors and time integration. One may replace it by exact enumerated smaller constants. The result is

  W2(Law(W1),Law(W0))
   <= C_m [L^m e + L² h
            + sum_(r=3)^m h_r sum_(j=1)^(r-1) s^j].     (A8)

The rank-ell current is transferred to the keep using sqrt((ell-1)!) eta^(-(ell-1)); this is included in C_m at v=1. Covariance/time terms have norm <=C_m L²h because ||Σ||op<=L² and partial_t C_t^(-1/2) is bounded by C L². Polynomial Gaussian moments justify truncation/Sobolev approximation.

For general v, apply (A8) with L'=L/sqrt(v), e'=e/sqrt(v), h'_r=h_r/v^(r/2), k'_r=k_r/v^(r/2), then multiply the result by sqrt(v). This gives every width loss explicitly.

## 5. What does NOT improve with rank

The polynomial reference with RAW cumulants is repeatable, and the exact generator includes all its feedback, but adding ranks alone is not an arbitrary-order approximation. Even for m>=5, the moving-covariance/time term L²h3 is O(L^4 e) at unit buffer. The Stein remainder improves to L^m e while this old term remains. Under a computable small-radius condition, e.g. L'<=min(1/2,(2 sum_(r=3)^m c_r/r!)^(-1/3)), all s-products are bounded and

  m=3: error <=C'_3 L³e/v^(3/2),
  m>=4: error <=C'_m L⁴e/v².

For complete numerical constants, put M_m=sum_(r=3)^m c_r/r! and take C′_m=C_m(3+2M_m); these majorize the displayed bounds under that guard. This is an upper-bound ceiling of THIS uncorrected consumer, not an impossibility theorem. To get higher grades, generated lower-rank feedback must be constructively canceled or incorporated into corrected reference coefficients/covariance. The raw-cumulant map can itself alter lower cumulants. In scalar standard coordinates, Z+aH2(Z)+bH3(Z)+cH4(Z) has variance 1+2a²+6b²+24c² exactly; ignoring those terms is not high-order matching.

The all-rank inequalities and finite positive consumer therefore give computable analytical targets and error certificates. They do not create the original-gradient VALUE producers or close their observer-boundary/cycle/current queue. No c(P)=o(P) claim follows.
