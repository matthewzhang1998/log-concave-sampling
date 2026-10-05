# Independent audit: odd cancellation, target identity, and fourth current

Date: 2026-10-04.

## Verdict and source pins

**PASS for the stated coefficient identities, finite-list obstruction, conditional Markov separator, and the separate buffered fourth-order Gaussianization theorem.** No source edit was made. The exact two-copy symmetrization improves its own Gaussianization but does not preserve a general skew target law. This is a bounded analytical/source-interface conclusion, not an all-order impossibility theorem.

Audited sources:

1. `../ODD-CUMULANT-CANCELLATION-AND-TARGET-IDENTITY.md`, SHA256 `ce0d7e8e262a8168aac8818a7830fdb720443fc75b8ad4b1c64343287ef0d9c9`.
2. `../../rank3-continuation/FOURTH-ORDER-GAUSSIANIZATION-WHEN-THIRD-CUMULANT-VANISHES.md`, SHA256 `06ab222f8f739490c0142582a565ebaf3c9605a7ef0963385502bab32105c880`.

The second pinned source supplies the theorem invoked without a path in the first source's Section 2. It has been inspected directly rather than accepted on the basis of its author's diagnostics.

The independent checker passes **133 assertions**, including exact rational coefficient/Markov calculations, independent symbolic Gaussian Riesz identities for both skew and centrally symmetric vector examples, and the exact third-Hermite buffer isometry. Its polynomial fixtures test identities only; they are not globally Lipschitz examples of the theorem.

## 1. Finite weighted copies: exact scope of the obstruction

Complete conditional iid copies give

    kappa_k(sum_i a_i F_i | u)=(sum_i a_i^k) kappa_k(F|u).

The nine coefficients have the exact first five power sums claimed: `1,1/3,0,1/54,-1/324`. Their signs do not make a signed probability law; they are deterministic coefficients of a positive random-variable pushforward. Hidden shared random ancestors would, however, introduce mixed cumulants and invalidate this identity.

For distinct absolute coefficient magnitudes b_j>0, let d_j be positive-minus-negative multiplicities. Cancellation of all higher odd ranks would imply

    sum_j (d_j b_j^3)(b_j^2)^m=0, m=0,...,J-1.

The determinant is `prod_j b_j^3 prod_(i<j)(b_j^2-b_i^2)`, up to its row/column ordering sign, and is nonzero. Hence every d_j=0 and the linear coefficient sum is zero. This proves the source's exact finite-list obstruction, including arbitrary repetitions and zero coefficients. It proves nothing against nonlinear corrections, random conditional coefficients, growing lists, or other positive law constructors.

The stated existence of a finite pattern for any prescribed finite odd cutoff can also be made constructive. For cancellation through rank 2m+1, choose b_j=j, j=1,...,m+1. The m-by-(m+1) matrix with entries b_j^(2r+1), r=1,...,m, has a nonzero rational null vector d. Multiply to integer entries and realize each positive or negative entry as that many signed copies of b_j. The linear sum c=sum_j d_j b_j is nonzero: adjoining that row makes an ordinary Vandermonde after factoring b_j. Divide every actual coefficient by c. The new linear sum is one and all selected odd power sums vanish. Copy counts may grow rapidly with m; nothing here supplies a favorable growing-order complexity bound. The checker constructs exact patterns through odd rank 15.

The fourth power sum is strictly positive for every nonzero real coefficient list. Independent Gaussian fill has no fourth cumulant. Thus these two operations alone cannot erase a nonzero fourth cumulant of the common source.

## 2. Exact two-copy symmetrization

For `D=(F_1-F_2)/sqrt(2)` on complete iid conditional tapes, swapping the two tapes negates D. It is centered, has exactly Cov(F|u), and hence has the same centered L2 energy e. Its complete input derivative is

    [DF(V_1)/sqrt(2),-DF(V_2)/sqrt(2)].

Multiplying by its transpose bounds its squared operator norm by L^2. This confirms the unchanged complete private first L, without a spurious factor sqrt(2). The generic caller bound is different: the two parameter derivatives can add, so sqrt(2) times the old caller bound is safe unless another cancellation is proved.

Execution requires no uncomputed mean. Restoring a nonzero mean later requires a separate actual mean program and its variance/errors. Complete swap symmetry requires the same finite source version, same exposed caller, and independent full private banks. Captured deterministic anchors may be shared under their identical keys; old random subgraphs may not be shared as though they were captured constants.

## 3. Independent derivation of the fourth-order Stein current

Let X=f(V)-Ef(V), let e=||X||2 and Sigma=Cov(X), with `||Df||op<=L`. On a centered Gaussian Hilbert field H, the operator `R H=grad(-L_OU)^(-1)H` obeys

    ||R H||_(L2 HS)<=||H||_(L2),
    E[H h]=E[(R H)·grad h].

The norm inequality follows degree by degree in Gaussian chaos, with squared factor 1/k at chaos degree k>=1. This is valid for tensor fields by treating their entries as one finite Hilbert index.

Define `tau_ij=(R X)_i·grad f_j`. Then Etau=Sigma, and multiplication of the Gaussian input slot by Df gives `||tau||_(L2 HS)<=Le`. Subtracting its expectation is an orthogonal L2 projection, so for B=tau-Sigma one has `||B||2<=Le`, rather than needing a factor two.

Contract R B with Df in its input slot and symmetrize the three output/test indices to obtain T. The contraction costs at most L, and standard averaging symmetrization is an orthogonal projection. Thus `||T||2<=L^2e`. Testing the first Stein identity against X_j X_k gives

    E X_i X_j X_k=E[B_ij X_k+B_ik X_j],
    E T=kappa_3/2.

The constant is exactly one half; no extra factor from three-index symmetrization appears. B need not be symmetric before contraction. Only its contraction with the symmetric test derivatives is used.

If kappa_3=0, T is centered. Apply R once more, contract with Df, and symmetrize four indices to obtain Q. The same Hilbert bounds give

    ||Q||_(L2 HS)<=L^3e.

No derivative of B, T, Q, or an emitted Hessian action is used. R is applied to their L2 fields; the derivative in the integration-by-parts identity always lands on the composed smooth test function, and therefore uses only the actual first Df.

For `Y_t=tX+C_t^(1/2)Z` with `C_t=sigma^2 I+(1-t^2)Sigma`, covariance differentiation and the first Stein identity give

    d/dt E phi(Y_t)=t E[B:D^2phi(Y_t)]
                  =t^2 E[T:D^3phi(Y_t)]
                  =t^3 E[Q:D^4phi(Y_t)].

The endpoint remains Y_t in every equality. There is no Gaussian endpoint replacement or unpriced feedback remainder. In the general case, the last line instead has the exact additional term `(t^2/2) kappa_3:E D^3phi(Y_t)` after centering T.

Integrate three test derivatives through the independent Z. Conditional on V, Q is fixed. The resulting velocity candidate has coefficients Q contracted with three inverse-root factors and the third Gaussian Hermite tensor. Its conditional squared L2 norm is bounded by

    3! t^6 ||C_t^(-1/2)||op^6 ||Q||HS^2.

Conditioning that vector on Y_t can only reduce its L2 norm and supplies a continuity-equation velocity. The fixed buffer has `||C_t^(-1/2)||op<=sigma^(-1)`. The Wasserstein dynamic length is therefore at most

    sqrt(6) L^3e sigma^(-3) int_0^1 t^3 dt
      =sqrt(6)L^3e/(4sigma^3).

This is the claimed exact constant. The Gaussian endpoints and finite second moments make the path weakly continuous; the displayed integrable velocity bound gives the required W2 absolute continuity. Sobolev approximation extends the smooth-test argument to globally C1 Lipschitz f. The only higher derivatives are of the test function and the independent Gaussian density.

For an anisotropic fixed buffer Q_buf>=qI, replace C_t by Q_buf+(1-t^2)Sigma. Covariance differentiation still uses C_t'=-2tSigma. It does not commute Sigma with Q_buf or use the generally false matrix formula `d sqrt(C_t)=C_t'/(2sqrt(C_t))`. The same operator bound gives q^(-3/2).

The auxiliary cumulant bounds also check. Full HS norm obeys `||kappa_3||HS<=2L^2e`. For a unit vector u and symmetric HS-unit A,

    |E[(u·X)(X^T A X)]|
      <=sqrt(Var(u·X) Var(X^T A X))<=2L^3,

because Gaussian Poincare gives `Var(u·X)<=L^2`, `Var(X^T A X)<=4L^2 E|AX|^2<=4L^4`, and Sigma<=L^2I. Antisymmetric A contributes zero. This establishes the stated flattening bound without an extra dimension factor.

## 4. True conditional Markov-path third cumulant

At the fixed conditional caller Z=X_1=0, the covariance is

    C(r,s)=min(r,s)/max(r,s)-rs.

Splitting the integral at s=r gives

    Cov(L,X_r)=int_0^r(s/r-sr)ds+int_r^1(r/s-rs)ds
              =-r log r.

In particular Var(L)=1/4. This is the true conditioned Markov covariance; it is not the covariance of independent point samples or a common-root surrogate.

For epsilon=1/20, `g_0(x)=(x+epsilon log cosh x)/(1+epsilon)` is anchored and smooth, with derivative in [19/21,1]. Consequently g_A=A g_0 is in the allowed convex bounded-Hessian class. Global path sign reversal makes L odd and J even, eliminating the L^3 and L J_c^2 terms in the centered third moment.

For each joint Gaussian pair (L,X_r), two scalar integrations by parts give

    E[L^2(log cosh X_r-E log cosh X_r)]
      =Cov(L,X_r)^2 E sech^2(X_r).

All moments needed for integrating this identity in r are finite. Since Var(X_r)=1-r^2<=1, one has P(|X_r|<=1)>1/2 and sech^2(1)>1/4. For an elementary lower bound on the first fact, integrating `exp(-x^2/2)>=1-x^2/2` over [-1,1] gives `5/(3sqrt(2pi))>1/2`. Therefore

    E(L^2 J_c)>=(1/8) int_0^1 r^2(log r)^2 dr=1/108.

Using `0<=log cosh x<=x^2/2` and Var(X_r)=1-r^2, Minkowski gives `||J||3<=15^(1/3)/3`. The centering inequality `||J_c||3<=2||J||3` gives `|E J_c^3|<=40/9`. Hence the third-moment numerator is at least

    3/(20*108)-(40/9)/20^3=1/1200,

and division by `(21/20)^3` gives exactly

    kappa_3(H_0)>=20/27783>0.

The finite positive path approximations justify the Gaussian calculations and converge in every fixed Lp; no infinite path is being offered as an executable oracle.

The Lipschitz chord inequality for the inner displacement gives `F_cont=A H_0+O_Lp(A^2)` for every fixed p. Centering preserves that remainder order. Expanding the centered cube therefore gives

    kappa_3(F_cont|Z=0)=A^3 kappa_3(H_0)+O(A^4).

For example, a conservative direct bound follows from `||F_cont-AH_0||3<=A^2||N(0,1)||3`, `||H_0||3<=||N(0,1)||3`, and their centered counterparts. Thus the positivity threshold can be made computable; it is not merely a formal nonzero coefficient.

## 5. Buffered W2 separation and its limit

The symmetrized centered source has a symmetric law, so its sine expectation remains zero after adding an independent centered Gaussian buffer. For the original centered source and fixed h>0,

    E sin(hF_c+sigma N)
      =-exp(-sigma^2/2) h^3 E F_c^3/6+O(A^5).

The linear term is exactly zero; the remainder follows from `|sin x-x+x^3/6|<=|x|^5/120` and the uniform fifth-moment bound. The actual cubic moment includes its O(A^4) correction, which does not alter the positive order-A^3 lower bound for sufficiently small A. Since sine is 1-Lipschitz, the expectation difference lower-bounds W1 and hence W2. Translating both laws by the same deterministic mean changes neither distance; equivalently translate the test function.

The conclusion is a conditional `c_h A^3` separation at the permitted retained caller Z=0. It refutes an asserted uniform conditional A^4 replacement at that port. It does **not** alone prove an unconditional lower bound after mixing Z, and the inspected note explicitly observes this restriction. No broader impossibility result should be read into this fixture.

## 6. Cost and remaining construction task

At each fixed finite odd cutoff, the coefficient pattern multiplies complete source work by a finite order-dependent count. The two-copy symmetrization costs exactly two complete private source occurrences plus once-captured deterministic caller work; changed source arguments still replay every ancestor. Its first and adjoint remain ordinary first sweeps of those two VALUE programs. Numerical errors and anchor mismatches are absolute floors, never normalized by realized energy.

The fourth-order Gaussianization lemma is an analytical law comparison for the actual symmetrized source. It does not execute a third-current correction for the original skew source. A correct posterior-law continuation still needs either such a positive original-VALUE producer with all feedback/lower-rank/variance records, or another fully proved reference that preserves the target law. More cancellation of odd iid-copy cumulants alone does not remove the current mean service's curl/covariance-orientation allowance.

**Conclusion:** all newly claimed bounded facts at these two pins are supported. The exact target mismatch remains, and the general positive third-current producer remains unconstructed.
