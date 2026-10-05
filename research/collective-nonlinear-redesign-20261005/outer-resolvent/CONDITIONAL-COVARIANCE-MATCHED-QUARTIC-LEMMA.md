# A one-energy quartic outer-resolvent lemma at a matched conditional covariance

2026-10-05. Analytical consumer theorem for actual finite Gaussian sources. This is not a covariance producer, expectation oracle, or complete canonical-m3 algorithm.

## The input contract and bound

Let (X,Y) be jointly standard D-dimensional Gaussian with Cov(X,Y)=q I, 0<=q<1, and sigma=sqrt(1-q^2). For i=1,2 let b_i(Y,Z_i) be actual finite Gaussian-input maps, where each private Z_i is standard Gaussian conditional on Y and independent of X conditional on Y (in the executed graphs, these are fresh Gaussian roots independent of the pair X,Y). Assume their full-bank Lipschitz constants are at most L_i, uniformly in Y, and b_i have finite second moments. Define analytically

    m_i(Y)=E[b_i|Y],
    Sigma_i(Y)=Cov(b_i|Y),
    delta_m=||m_1-m_2||_(L2 gamma_D),
    delta_Sigma=||Sigma_1-Sigma_2||_(L2 gamma_D;HS).

Let F:R^D->R^D be K-Lipschitz, and put

    e(x)=E[F(X-b_1)-F(X-b_2)|X=x].

Then

    ||R1 e||_2
      <= K delta_m
       +(K/2)(1/2+1/sigma) delta_Sigma
       +(K sqrt(D)/3)(L_1^3+L_2^3)
                        (1+1/sigma+sqrt(2)/sigma^2).        (1)

In particular, exact conditional mean and covariance matching, K=O(A), L_i=O(A), and a fixed positive sigma give O(A^4sqrt(D)). The coefficient sources and moments in the proof are analytical; neither conditional means, covariance roots, nor tensor fields become executing leaves.

The norm delta_Sigma is conditional covariance error under the actual Y, in Hilbert-Schmidt norm. An unconditional covariance match, a covariance of a completed mean-law output, or a joint result with missing retained variables does not meet this input contract.

## 1. Conditional centered Stein fields

Freeze Y and write u=b-m. Let R=D(-L_OU)^(-1) be the Gaussian Riesz operator in the complete private bank. Define

    tau=R(u) (Db)^T,
    T2=Sym_3[R(tau-Sigma) contracted with Db].               (2)

Here Gaussian root indices are contracted and all physical indices remain. Sym_3 is the average over physical permutations. Gaussian covariance identities give

    E tau=Sigma,
    E[u_i phi(u)]=E[tau_ij partial_j phi(u)],
    E[(tau-Sigma):D^2 phi(u)]=E[T2:D^3 phi(u)].              (3)

The second identity in (3) is applied componentwise before symmetrization; no derivative of tau or Db is taken. R is an L2 contraction on each centered Hilbert-valued source. Therefore, with e_Y=||u||_(L2 private |Y),

    ||tau-Sigma||_(L2;HS|Y)<=L e_Y,
    ||T2||_(L2;HS|Y)<=L^2 e_Y.                              (4)

Conditional Gaussian Poincare and the full-bank operator bound yield

    e_Y^2<=E||Db||HS^2<=D L^2,
    ||T2||_(L2;HS)<=L^3 sqrt(D).                            (5)

The physical rank-three Hilbert norm carries only one energy. No factor from private Gaussian dimension is introduced.

## 2. Covariance-preserving interpolation

Take an independent standard D-Gaussian G. For analysis only let

    C_lambda=lambda u+sqrt(1-lambda^2) Sigma(Y)^(1/2)G,
    lambda in [0,1].

For fixed X,Y, Gaussian covariance differentiation and (3) give

    d/dlambda E F(X-m-C_lambda)
      =lambda E[(tau-Sigma):D^2F(X-m-C_lambda)]
      =-lambda^2 E[T2:D^3F(X-m-C_lambda)].                 (6)

The minus sign in the last equality is from differentiating F(X-m-C_lambda) in the private bank. It does not affect the norm bound. Both tau and T2 use the actual b source and remain at the same interpolated endpoint; the test is not replaced by its Gaussian version midway.

For nonsmooth Lipschitz b or F, first use Gaussian/Sobolev and Euclidean smoothing, respectively, preserving the displayed bounds, and pass to the L2 limit. Degenerate Sigma is allowed directly: Gaussian integration by parts is applied through the PSD root Sigma^(1/2), and only Sigma appears in the identity. No inverse Sigma is used.

## 3. Exactly two derivatives supplied by the outer R1

Dualize against a vector f in L2(gamma_D), and put h=R1 f. Hermite chaos gives

    ||h||_2<=||f||_2,
    ||Dh||_(2,HS)<=||f||_2/2,
    ||D^2h||_(2,HS)<=||f||_2.                              (7)

The last multiplier is sqrt(n(n-1))/(n+1)<=1. No third derivative of h is claimed.

Conditional on Y, write X=qY+sigma eta. The Gaussian eta is independent of T2, C_lambda and all private source noises. Integrate two of the three derivatives of F in (6) by parts in X, leaving one DF factor. The three resulting families have the following absolute bounds:

    K ||T2||_2 ||D^2h||_2,
    (2K/sigma)||T2||_2 ||Dh||_2,
    (sqrt(2)K/sigma^2)||T2||_2 ||h||_2.                    (8)

These estimates preserve all indices. For the middle term, contraction of one tensor index with eta has squared L2 norm ||T2||_2^2. For the last, contraction of the two symmetric tensor indices with eta eta^T-I has squared L2 norm 2||T2||_2^2. These are conditional Gaussian isometries, independent of Y and the private source bank. Cauchy-Schwarz is then taken on the full joint law. There is no independence assumption between h(X) and eta, nor between h(X) and Y.

The remaining DF factor has operator norm K. This is why the isometries in (8) do not incur an additional D or sqrt(D). Summing (8), using (5),(7), and integrating lambda^2 over [0,1] proves

    ||R1 E[F(X-b)-F(X-m-Sigma^(1/2)G)|X]||_2
      <= (K L^3sqrt(D)/3)(1+1/sigma+sqrt(2)/sigma^2).        (9)

This is the desired quartic one-energy remainder at a genuinely covariance-matched reference.

## 4. Mean and covariance mismatch

Compare the two conditional Gaussian references in (9). Their mean discrepancy costs K delta_m by ordinary Lipschitzness, conditional Jensen and the R1 contraction.

For equal means, interpolate their covariance linearly. Gaussian covariance differentiation supplies

    (1/2)(Sigma_1-Sigma_2):D^2F.

One conditional-X integration by parts, leaving DF, gives two terms bounded by

    (K/2) delta_Sigma ||Dh||_2,
    (K/(2sigma)) delta_Sigma ||h||_2.

The score contraction uses the exact identity E|DeltaSigma eta|^2=||DeltaSigma||HS^2 conditional on Y. Thus the reference difference costs at most `(K/2)(1/2+1/sigma)delta_Sigma`. Adding the two remainders (9) proves (1).

## 5. Actual next-target opportunity and its remaining source obligation

For the genuine canonical F2 tail after time delta, set Y=X_delta and

    b_true=q integral_0^infinity exp(-t)
                g(X_(delta+t)-H_(delta+t))dt,
    H_s=integral_0^infinity exp(-h)g(X_(s+h))dh.

It is independent of X_0 conditional on Y, and its full future-bank first is at most qA(1+A). Apply the finite-source lemma first to finite Gaussian approximations of this true tail with uniform whole-bank first and strong L2 convergence, and then pass to the limit; the infinite history is never executed. Its conditional mean is

    E[b_true|Y]=q m_cont(Y),
    m_cont=R1 E[g(X0-H0)|X0].                              (10)

This is exactly the mean addressed by the active third-order nested-force mean construction, rather than an unavailable new mean oracle. A finite raw source for (10) with O(A^3sqrt(D)) mean error is already available there. Its covariance relative to this true delayed F2 tail is the new quantity that must be matched; the earlier mean guarantee alone does not supply it.

If an actual finite RAW source b_fin has

    L_fin<=C A,
    delta_m<=C A^3sqrt(D),
    delta_Sigma<=C A^3sqrt(D),                              (11)

with the same retained Y and a complete private bank independent of X_0 conditional on Y, then (1) becomes an immediately usable conditional source comparison. No complete LAW output may be relabeled as such a RAW b_fin solely from a small W2 error.

The true short F2 prefix can be compared coherently with w g(X0), w=1-q, at force error

    (2sqrt(2)/3)A sqrt(D) w^(3/2)+A^2sqrt(D)w.

Moving that deterministic prefix into F(u)=g(u-wg(u)) adds O(A^3w sqrt(D)); no global A^3 omission of the whole nested history is necessary. With (11), choosing w=A^(4/5) would balance

    A^2 w^(3/2)    and    A^4/sigma^2=O(A^4/w),

at A^(16/5). The mean and covariance mismatch terms would be of higher order, respectively A^4 and A^(18/5), and the prefix-source correction A^(19/5).

This is a CONDITIONAL construction, not a proved A^(16/5) algorithm: an actual raw conditional covariance-matched producer meeting (11), its positivity/retained-readset agreement and full costs must be verified. The existing covariance service is being compared against this precise port; it must not be assumed to expose its analytical covariance or hidden private roots.

For a direct target-four application at fixed sigma, (1) is already quartic. Handling the actual short prefix to the same grade and supplying the correct covariance source are additional obligations. Simply applying (1) to the earlier w=sqrt(A) graph would leave its prefix error A^(11/4); no grade upgrade follows automatically.

## 6. No automatic all-rank bootstrap

This proof spends the two Sobolev derivatives supplied uniformly by the outer R1. An identical next step would request D^3(R1 f) for arbitrary L2 f, which is unbounded. Higher-order covariance/cumulant producers, an additional independent Gaussian buffer, or another explicit cancellation would therefore be needed to continue the same argument. The displayed quartic lemma is a precise finite-rank consumer contract, not an all-order finite VALUE compiler.
