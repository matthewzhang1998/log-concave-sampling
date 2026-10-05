# A general-C2 order-two dyadic interaction tail bound, with dimension losses exposed

2026-10-05. Analytical theorem for the true F2 history. No original-VALUE source, conditional-mean oracle, path-grid algorithm, or compiler is supplied.

## Result

Let X be stationary standard D-dimensional OU, g=grad U, U in C2, g(0)=0 and 0<=Dg<=A I, 0<A<=1/2. Assume, when claiming bounded-block uniformity, that g separates in an orthogonal decomposition into blocks of dimensions b_j<=b. Without such a structure use b=D.

Use the true F1,t=int_0^infinity exp(-s)g(X_(t+s)) ds and F2,0=int_0^infinity exp(-t)g(X_t-F1,t) dt. Fix T>=1. Let S_k contain X_0 and X_(jT/2^k), j=1,...,2^k, and d_k=T/2^k. Let H_k=E[F2,0|S_(k+1)]. Conditional on S_k, let P_(k,>2) be the orthogonal Hoeffding projection onto subsets of at least three of the independent new midpoint blocks. Put r_k=P_(k,>2)H_k and

    Q_k(Y)=E[r_k r_k^*|Y], Y=X_0.

The projection has mean zero given S_k. Q_k is precisely the positive covariance discarded by retaining only singleton and pair midpoint interactions at level k.

There is an absolute constant C, independent of the source's Hessian modulus, A,D,T,k and block sizes, such that

    ||Q_k||_(L2(Y;HS))
      <= C sqrt(D) min{A^2 d_k^2,
                        A^5 b^(3/2)(1+d_k^(-1/2))}.                 (1)

The first two levels have Q_0=Q_1=0 exactly. Summing all levels, if A^3 b^(3/2)<=1,

    ||sum_(k>=0) Q_k||_(L2(Y;HS))
      <= C sqrt(D) [ A^(22/5)b^(6/5)
                     + A^5 b^(3/2)(1+log_+ T) ].                 (2)

The same estimate holds for any finite partial sum. Thus for fixed b and logarithmic T, retaining interaction orders <=2 has an o(A^4 sqrt(D)) omitted-covariance error uniformly over the entire C2/Hessian-interval class. No higher derivative modulus is assumed or frozen.

For an unrestricted source b=D, the leading bound is A^(22/5) D^(17/10); it certifies improvement relative to A^4 sqrt(D) only in a sufficient small-dimension regime such as A D^3 ->0, with the displayed coarse-horizon term also small. This is NOT a dimension-uniform improvement for arbitrary growing D, and is not a counterexample showing that such an improvement is impossible. The radial obstruction for a different finite stage-two graph does not itself prove a lower bound for these exact Hoeffding components.

## 1. Conditional Gaussian weak transport

Use one physical block of dimension n<=b. For X~N(m,c^2 I_n), F in C1, ||DF||<=A<1, and a bounded vector-valued f with sup|f|<=delta, Gaussian relative entropy, the log-determinant identity and Pinsker give

    |E[f(X-F(X))-f(X)]|
      <= delta [c^(-1)(E|F(X)|^2)^(1/2)
                + A sqrt(n/(1-A))].                              (3)

The linear displacement term cancels tr DF by Gaussian integration by parts. No derivative of DF is used. One may first regularize and truncate and pass to the limit; the uniform derivative bound and Gaussian moments control that passage.

Fix t strictly inside a fine S_(k+1) cell [a,c_0]. Conditional on the fine skeleton, decompose the entire residual Gaussian path into X_t and Gaussian residuals independent of X_t. For s>=t, the coefficient of X_t in X_s is

    rho_s=sinh(c_0-s)/sinh(c_0-t), t<=s<=c_0,
    rho_s=0, s>=c_0.

It lies in [0,1]. Therefore, after all independent residuals have been fixed,

    ||D_(X_t) F1,t||
      <= A int_t^(c_0) exp(-(s-t))rho_s ds <=A.

The conditional variance of X_t is

    c_t^2=2 sinh(t-a)sinh(c_0-t)/sinh(c_0-a).

Thus (3) applies to the ACTUAL ancestral shift F1,t, including its complete future; it is not replaced by a conditional mean.

If the fine cell length is ell=d_k/2, the elementary bridge variance formula gives the weighted estimate

    J_k := int_0^T exp(-t)(1+c_t^(-1))dt
          <= C(1+d_k^(-1/2)).                                    (4)

For ell<=1, each cell's unweighted integral of c_t^(-1) is <=C sqrt(ell), followed by the geometric sum of exp(-a). For ell>=1 use c_t^(-1)<=C[1+(t-a)^(-1/2)+(c_0-t)^(-1/2)] and the same weighted sum. Endpoints are null in the t integral. The bound is uniform in all retained skeleton values.

## 2. Coherent mollification stability for the conditional residual

Define, only for the proof, the anchored Gaussian mollification

    h_e(x)=E[g(x+eZ)-g(eZ)], Z~N(0,I_n).

It has h_e(0)=0, 0<=Dh_e<=A I and

    delta:=sup|g-h_e|<=2A e sqrt(n),
    Lip(Dh_e)<=C A/e.                                             (5)

Let R_g^T=int_0^T exp(-t)[g(X_t-F1,g,t)-g(X_t)]dt. The exact coherent source comparison at each t is

    f(X_t-F1,g,t)-f(X_t)
      +h_e(X_t-F1,g,t)-h_e(X_t-F1,h_e,t), f=g-h_e.

The second term has norm <=A delta because |F1,g,t-F1,h_e,t|<=delta. Apply (3) to the first term after conditioning as above. Stationarity and Minkowski give

    ||F1,g,t||_Lp <= A (E|N(0,I_n)|^p)^(1/p).

Taking L4 over the fine skeleton and all other residual variables, then using conditional Jensen and (4), yields

    ||E[R_g^T-R_h_e^T|S_(k+1)]||_L4
          <=C A delta sqrt(n) J_k
          <=C A^2 e n J_k.                                      (6)

The entire source and every ancestor are changed coherently. The future t>=T has no new-midpoint dependence given S_k and therefore contributes zero to P_(k,>2), even though that future is still present inside F1,t for t<T.

## 3. The genuine smooth first-order term has interaction order at most two

For h=h_e, write exactly

    R_h^T=-B_h^T+E_h^T,
    B_h^T=int_0^T exp(-t) Dh(X_t)F1,h,t dt.

Taylor's integral remainder, now with the explicitly paid smoothness (5), gives

    ||E_h^T||_L4
      <= C(A/e) int_0^T exp(-t)||F1,h,t||_L8^2 dt
      <= C A^3 n/e.                                               (7)

Condition on S_k. Complete bridges in different coarse cells, and the complete future after T, are independent. At time t in coarse cell i, Dh(X_t) depends only on that bridge, while F1,h,t is a sum of contributions from the same bridge, each later bridge, and the terminal future. The product is consequently a sum of terms involving at most two coarse bridge objects. After conditional expectation onto S_(k+1), it is a sum of functions of at most two midpoint blocks. Thus

    P_(k,>2) E[B_h^T|S_(k+1)]=0.                                 (8)

The identical assertion for F1 is immediate from its additivity. This is an exact time-block interaction statement; coefficients may still depend on every coarse endpoint.

Equations (6)-(8) exhibit an actual comparison V_k with

    r_k=P_(k,>2)V_k,
    ||V_k||_L4 <= C n[A^2 e J_k+A^3/e].                           (9)

Choose e^2=A/J_k. No execution uses e, h_e, or its derivatives. Then

    ||V_k||_L4^2 <= C A^5 n^2 J_k.                               (10)

There is no unjustified L4 contraction of the high-order Hoeffding projector here. Instead, its conditional L2 orthogonality gives the matrix inequality

    E[r_k r_k^*|S_k] <= E[V_k V_k^*|S_k].                         (11)

Indeed (11) follows by L2 contraction for every fixed physical test vector. Average over S_k conditional on Y, use HS<=trace for a positive matrix, and conditional Jensen. This proves the one-block bound C A^5 n^2 J_k.

## 4. Dimension-safe fine-scale cap

The true F2 full path sensitivity has kernel

    beta_2(s)=A exp(-s)(1+As).

Conditional on S_k, a standardized new midpoint shifts the path by sqrt(tanh(d_k/2)) times its OU midpoint regression hat, supported in its own cell and bounded by one. Gaussian Poincare therefore gives, uniformly in every retained skeleton,

    Cov(H_k|S_k)
      <= tanh(d_k/2) sum_I (int_I beta_2)^2 I_n
      <= C A^2 d_k^2 I_n.                                      (12)

Hoeffding orthogonality gives E[r_k r_k^*|S_k]<=Cov(H_k|S_k). Its one-block L2(HS) norm is therefore at most C A^2 d_k^2 sqrt(n). This uses one operator bound followed by one HS conversion.

For a physically block-separable g, the omitted covariance is block diagonal, conditional on the coarse skeleton and after conditioning only on Y: each nonempty Hoeffding component in a physical block has zero midpoint mean, and distinct physical blocks use independent midpoint coordinates. Squaring and summing the block HS norms, with sum n_j=D and n_j<=b, proves (1). Unknown block structure is an analytical hypothesis; no block basis is supplied to an executor by this proof.

## 5. Summing every scale

For d_k>1 use the second bound in (1), with at most C(1+log_+ T) such levels. For d_k<=1 set

    d_*=(A^3 b^(3/2))^(2/5)=A^(6/5)b^(3/5).

At scales d_k>=d_*, sum C A^5 b^(3/2)d_k^(-1/2). At finer scales sum C A^2 d_k^2. Both geometric series are bounded by

    C A^(22/5)b^(6/5).

This proves (2) and charges every refinement scale, rather than declaring conditional smoothing free near a midpoint or endpoint. Infinite sums are justified by positivity, finite partial-sum bounds and conditional martingale convergence.

If this is combined with the positive dyadic covariance identity, the omitted residual after a finite level K still has its separate bound

    C A^2[d_K^2+(1+T^4)exp(-2T)]sqrt(D).

For example T=2log(1/A) with a sufficiently small A, and d_K^2 of order at most A^(12/5)b^(6/5), make these terms no larger than the leading displayed fixed-b interaction error (up to harmless tail logarithms). Cov(E[F2|Y,X_T]|Y) and every retained singleton/pair covariance are still exact analytical objects; this theorem does not supply their computation.

## 6. Readset, dimension and cost boundaries

- The theorem concerns F2 only. It does not by itself control the full coherent (F2,F3-F2) block or its cross covariance.
- Order-two in midpoint labels does not mean a two-endpoint readset. A singleton or pair coefficient depends on the entire retained coarse skeleton and on exact integration of complete coherent ancestors.
- At level k there are 2^k midpoint blocks, 2^k singleton coefficients and binomial(2^k,2) pair coefficients. If all are enumerated through level K, there are order 4^K subset occurrences. Logarithmic depth is not logarithmic node count.
- The retained Gaussian skeleton has order 2^k D scalar coordinates. Private conditional-copy roots and full history integration are not finite original-VALUE services supplied by this argument.
- No exact h_e leaf, smoothing quadrature or source derivative is executed. Equations (5)-(10) are comparison estimates only.
- No unrestricted-D uniform o(A^4 sqrt(D)) claim follows. Nor is there an impossibility theorem here. The sufficient fixed-block gain and the missing unrestricted-dimension estimate should be kept distinct.
