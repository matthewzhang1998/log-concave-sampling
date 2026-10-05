# Fourth-order buffered Gaussianization when the third cumulant vanishes

2026-10-04. New analytical one-energy theorem. This note supplies a law comparison for an actual positive source, not a coefficient-oracle producer. Original-gradient VALUE implementations and their first/adjoint sweeps remain exactly those of the supplied source. Existing audited files are unchanged.

## Result

Let V be standard Gaussian, let f(V) take values in R^D, and assume f is C1 and globally L-Lipschitz. Define

    X=f(V)-E f(V), e=(E|X|^2)^(1/2), Sigma=Cov(X).

If the third cumulant tensor E X^(tensor 3) is zero, then, for an independent standard Gaussian Z and sigma>0,

    W2(Law(X+sigma Z),N(0,sigma^2 I+Sigma))
        <= sqrt(6)/(4 sigma^3) L^3 e.                    (1)

The conclusion is dimension free. At L=O(A), e=O(A sqrt(D)) and fixed buffer, it is a full-law O(A^4 sqrt(D)) estimate. It includes the entire fourth-and-higher-cumulant remainder; it is not moment matching alone. Central symmetry is sufficient but unnecessary. No second derivative of f is assumed or queried.

An anisotropic independent buffer Q>=q I has the same conclusion with sigma^3 replaced by q^(3/2). Use covariance differentiation, rather than differentiating a noncommuting square root as though its factors commute.

## 1. Hilbert-valued Gaussian Riesz operation

For any centered square-integrable finite-dimensional Hilbert-valued field H(V), write

    R H = grad (-L_OU)^(-1) H.

This notation is analytical only. Chaos isometry gives

    ||R H||_(L2 HS) <= ||H||_(L2 HS),
    E[H_a h(V)] = E[(R H)_a dot grad h(V)].             (2)

The second identity holds initially for smooth tests and then by Sobolev approximation. Tensor entries can be regarded as a single Hilbert index. Multiplication of the last input slot by Df costs at most L in full Hilbert-Schmidt norm. Symmetrization in output/test slots is an orthogonal projection and cannot increase that norm. No derivative of R H is taken.

Define the first Stein matrix

    tau_ij = sum_a (R X)_(i,a) partial_a f_j,
    B=tau-Sigma.

Then

    E tau=Sigma, E B=0,
    ||tau||_(L2 HS)<=L e, ||B||_(L2 HS)<=L e.           (3)

The last inequality is variance contraction, not the loose triangle bound.

Define the second Stein tensor, with standard averaging symmetrization,

    T = Sym_3 [ sum_a (R B)_(i,j,a) partial_a f_k ].    (4)

It satisfies

    ||T||_(L2 HS)<=L^2 e,
    E T = (1/2) E X^(tensor 3).                        (5)

To check the factor, test the first Stein identity against the product X_j X_k:

    E X_i X_j X_k = E[B_ij X_k+B_ik X_j].

Symmetrizing gives kappa_3=2 Sym E[B tensor X]. Equation (2) identifies the latter mean with E T. The factor is two, not three.

When kappa_3=0, T is centered. Define

    Q = Sym_4 [ sum_a (R T)_(i,j,k,a) partial_a f_l ].  (6)

Again no derivative of B or T is assumed or evaluated. Equations (2)-(5) give

    ||Q||_(L2 HS)<=L^3 e.                             (7)

## 2. The actual covariance-preserving path

Set

    C_t=sigma^2 I+(1-t^2)Sigma,
    S_t=C_t^(1/2),
    Y_t=tX+S_t Z, 0<=t<=1.

For a smooth compactly supported scalar test phi, covariance differentiation and the first Stein identity give the exact derivative

    d/dt E phi(Y_t)=t E[B:D^2 phi(Y_t)].               (8)

Condition on Z. Applying (2) to the centered B differentiates phi only through t f(V), so

    (8)=t^2 E[T:D^3 phi(Y_t)].                        (9)

Under kappa_3=0, T is centered. A second application differentiates only the endpoint and gives

    (9)=t^3 E[Q:D^4 phi(Y_t)].                        (10)

Thus every discarded contribution is an exact fourth-order current at the same actual endpoint. There is no replacement of phi's argument by a Gaussian reference and no sampler-feedback term hidden by such a replacement.

Integrate three derivatives in (10) against the independent Gaussian Z. The resulting vector coefficient uses the third Gaussian Hermite tensor with three S_t^(-1) factors. Conditional Hermite isometry gives its L2 norm at most

    sqrt(3!) t^3 ||S_t^(-1)||op^3 ||Q||_(L2 HS)
        <= sqrt(6) t^3 L^3 e/sigma^3.                 (11)

Conditioning this vector on Y_t gives a continuity-equation velocity without increasing its L2 norm. The Wasserstein dynamic length bound, integrated from 0 to 1, yields (1), since integral_0^1 t^3 dt=1/4. The endpoints are the matching Gaussian and the original buffered source.

Finite second moments, the fixed Gaussian buffer, and the integrable velocity bound give endpoint weak continuity and the required absolute continuity. Lipschitz/Sobolev approximation justifies all displayed identities at C1 regularity. The input dimension never enters any estimate.

## 3. Exact analytical current for a general skew source

Without kappa_3=0, put M=E T=kappa_3/2 and define Q using T-M in (6). The exact identity becomes

    d/dt E phi(Y_t)
       = (t^2/2) kappa_3:E D^3 phi(Y_t)
         +t^3 E[Q:D^4 phi(Y_t)],
    ||Q||_(L2 HS)<=L^3 e.                             (12)

In particular ||kappa_3||HS<=2 L^2 e. Its proper flattening cuts satisfy

    ||kappa_3||_(R^D -> HS)<=2 L^3.                   (13)

For (13), contract against a unit u and an HS-unit symmetric matrix A. Gaussian Poincare gives Var(u dot X)<=L^2 and

    Var(X^T A X)<=4 L^2 E|AX|^2<=4 L^4,

using Sigma<=L^2 I. Cauchy-Schwarz proves (13); antisymmetric A does not contribute.

Equation (12) alone does not construct a positive third-cumulant repair. Its leading third-current trajectory is generally not a positive probability evolution by itself. In particular, replacing the endpoint by a Gaussian plus a formal Hermite tensor action is not an executed original-VALUE algorithm. This is the first nonzero port when the original source is skew.

## 4. Actual odd symmetrization on two complete independent source tapes

Expose the full retained caller u first. Suppose the actual source is F(u,V), with conditional first L(u) and centered energy e(u). Use two independent COMPLETE source tapes V1,V2 and execute

    D(u,V1,V2)=[F(u,V1)-F(u,V2)]/sqrt(2).              (14)

This is a positive finite VALUE program. Its exact conditional facts are

    E[D|u]=0,
    Cov(D|u)=Cov(F|u),
    e_D(u)=e(u),
    Lip_(V1,V2) D <= L(u),
    D(u,V2,V1)=-D(u,V1,V2).

For the first bound, the squared operator norm of

    [D_V F(u,V1)/sqrt(2), -D_V F(u,V2)/sqrt(2)]

is at most L(u)^2. Thus all odd conditional moments and cumulants vanish and (1) applies conditional on the unchanged retained caller u. Its conditional error is sqrt(6)L(u)^3 e(u)/(4sigma^3); integrating a same-u coupling takes the L2(u) norm of that profile.

The complete source tapes are marginalized at this comparison. They cannot subsequently be retained by invoking this marginal theorem. Captured caller labels remain exposed throughout and are allowed. The swap symmetry requires identical source/finite-version/caller keys for both complete copies, together with independent private banks. Shared deterministic anchors are allowed when they are the same captured record. Shared random private subgraphs can invalidate the stated covariance calculation.

The source's total-tape-zero output cancels exactly between the copies. If a first or adjoint is requested, it is the literal two-copy sweep, and uses only original HVPs at recorded original-gradient VALUE sites. No producer HVP is introduced. A safe query bill is two complete source occurrences plus once-captured caller-only work. Any changed raw argument requires its complete ancestor replay.

The symmetrization has no new clock or inverse-width guard. It inherits every finite source clock/filter/radius guard. Its numerical VALUE errors are multiplied by 1/sqrt(2) per copy and summed in L2; they are absolute floors, never divided by e. Its caller first is the actual difference of the two caller rows; the universally valid bound is sqrt(2) times the original complete caller bound, unless additional cancellation is proved. A same-caller mean service may be added on an independent COMPLETE bank, with its own positive variance share and errors.

## 5. What this does and does not repair

Equation (14) keeps the centered covariance but changes the original law. Its third cumulant is zero, whereas that of F may be order L^2 e. Restoring E F with a separate mean service does not restore the missing third cumulant.

This distinction is already sharp in dimension one. Take a centered, globally Lipschitz scalar q(V) with nonzero third moment and set X=Aq(V). For fixed sigma, X+sigma Z has third moment A^3 E q^3. Its symmetrized counterpart has third moment zero. Their fourth moments stay uniformly bounded as A tends to zero. For any coupling Y,Y',

    |E Y^3-E(Y')^3|
      <= ||Y-Y'||2 ||Y^2+YY'+(Y')^2||2.

The second factor has a uniform finite bound, so their W2 distance is bounded below by c A^3. They therefore cannot generally be interchangeable at O(A^4).

Odd-cumulant cancellation is valuable when the intended target is Gaussian, as in a centered normalization/mean packet followed by a correctly calibrated covariance reserve. It is not, by itself, a fourth-order replacement of an originally skew raw law. After odd cancellation, the next generic obstruction is the fourth cumulant; independent weighted copies have kappa_4(sum c_i X_i)=(sum c_i^4)kappa_4(X), so real coefficients alone cannot cancel a nonzero fourth cumulant without another correction mechanism.

## Scope and first blocked port

The exact fourth-order Gaussianization theorem and the executable two-copy symmetrization are complete analytical results, pending independent audit. They do not upgrade the original non-symmetric source's full law. For that task, the first remaining port is an actual source-qualified third-current correction realizing the leading term in (12), with its complete first/caller/zero/clock/variance and one-energy feedback returns. An analytical tensor norm, the buffered-current consumer, or a derivative-valued rectangular response does not supply that VALUE producer.
