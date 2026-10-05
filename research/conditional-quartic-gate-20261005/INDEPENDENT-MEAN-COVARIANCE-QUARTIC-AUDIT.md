# Independent audit: conditional mean-and-covariance matching gives a quartic resolvent gate

## Verdict

The proposed dimension-uniform bound is valid. The two conditional Gaussian integration-by-parts steps cost no extra dimension factor. They use only the L2 Hilbert-Schmidt norm of one third-order Stein tensor, one operator-norm bound on DF, and the dimension-free first/second derivative estimates for the actual outer R1 resolvent.

One sign correction is needed in the displayed interpolation identity: with the conventional Gaussian operator R=D(-L)^(-1) and the positive definition of T2 below, the derivative is **minus** lambda^2 E[T2:D3F]. This does not change any norm estimate.

This is an analytic comparison lemma. It does not construct or certify an executable source producing the matched conditional covariance, its square root, or the Stein tensor. In particular the conditional Gaussian used in the proof is not an authorized expectation or covariance oracle.

## 1. Precise statement

Let X,Y be jointly standard Gaussian vectors in R^D with Cov(X,Y)=qI, 0<=q<1, and sigma=sqrt(1-q^2). Let Z be a standard finite-dimensional Gaussian bank, independent of (X,Y). Assume b(Y,Z) is square-integrable and uniformly L-Lipschitz in its complete Z bank. Let

    m(Y)=E_Z b(Y,Z),
    Sigma(Y)=Cov_Z(b(Y,Z)).

Let G be a new independent standard D-dimensional Gaussian and define the analytical comparison b_G=m(Y)+Sigma(Y)^(1/2)G. For every K-Lipschitz vector map F:R^D->R^D,

    ||R1 E[F(X-b)-F(X-b_G) | X]||2
       <= (K L^2 e/3) [1+sigma^(-1)+sqrt(2)sigma^(-2)],     (1)

where e^2=E|b-m(Y)|^2. In particular e<=Lsqrt(D), so the right side is at most

    (K L^3 sqrt(D)/3)[1+sigma^(-1)+sqrt(2)sigma^(-2)].       (2)

For two sources b_1,b_2 with exactly the same conditional mean and covariance, the comparison bound is the sum of their right sides in (1), or the sum of L_i^3 in (2). If K=O(A), L_i=O(A) and sigma is bounded below by a positive constant, this is O(A^4sqrt(D)). If sigma shrinks with A, its displayed inverse-square cost must be retained.

## 2. Stein operator, tensor placement and the one-energy bound

All Gaussian operators in this section act in Z only, with Y held fixed. Put u=b-m(Y), and write

    R(a)=D_Z(-L_Z)^(-1)a
        =integral_0^infinity exp(-t) P_t(D_Z a)dt

when a has a weak derivative. More generally R is defined spectrally for every centered L2 scalar/vector/matrix a and satisfies

    ||R(a)||_(L2,HS) <= ||a||_(L2,HS).                      (3)

For u, the matrix R(u) has shape D by dim(Z). Define

    tau_jk=sum_alpha R(u)_(j alpha) (D_Z b)_(k alpha),
    tau=R(u)(D_Zb)^T.

The Gaussian covariance identity gives E_Z tau=Sigma. Moreover ||R(u)||_(L2,HS)<=e and ||D_Zb||op<=L imply

    ||tau-Sigma||_(L2,HS) <= ||tau||_(L2,HS) <= L e.         (4)

The centering estimate is exact orthogonal projection, also when Sigma depends on Y. Apply R componentwise to the matrix tau-Sigma and define

    T_raw,jkl=sum_alpha R(tau_jk-Sigma_jk)_alpha
                                  (D_Z b)_(l alpha),
    T2=Sym_(j,k,l) T_raw.

Full symmetrization preserves contraction with D3F and cannot increase the Hilbert-Schmidt norm. Multiplication of the final Gaussian-bank index by D_Zb costs its operator norm, so

    ||T2||_(L2,HS)
      <=L ||R(tau-Sigma)||_(L2,HS)
      <=L ||tau-Sigma||_(L2,HS)
      <=L^2 e <=L^3sqrt(D).                                (5)

This proof never differentiates D_Zb or tau. The possibly nonsymmetric random matrix tau is harmless; its mean is the symmetric covariance, and only its fully symmetrized third-order contraction is subsequently used.

The last inequality in (5) is Gaussian Poincare:

    e^2<=E||D_Zb||HS^2<=D L^2.

Thus the construction genuinely uses one source-energy factor rather than D separate scalar bounds or a dimension-dependent tensor operator estimate.

## 3. Covariance-preserving interpolation and its sign

Put a=Sigma^(1/2)G and

    C_lambda=lambda u+sqrt(1-lambda^2)a,
    Phi(lambda)=E F(X-m-C_lambda).

For 0<lambda<1, Gaussian integration in G and the first Stein identity in Z give

    Phi_i'(lambda)
       =lambda E sum_(j,k)(tau_jk-Sigma_jk)
                    partial_j partial_k F_i(X-m-C_lambda).

Apply the covariance identity once more, componentwise to tau-Sigma. Differentiation of X-m-C_lambda with respect to Z contributes -lambda D_Zb. Consequently

    Phi_i'(lambda)
       =-lambda^2 E sum_(j,k,l) T2_jkl
                    partial_j partial_k partial_l
                           F_i(X-m-C_lambda).              (6)

This establishes the sign noted in the verdict. Both endpoints are obtained by continuity. The integrable lambda^2 factor supplies exactly integral_0^1 lambda^2 d lambda=1/3.

## 4. Two conditional-X integration-by-parts steps

For a vector test f with ||f||2=1 set h=R1f. Hermite spectral calculus gives

    ||h||2<=1,
    ||Dh||2<=1/2,
    ||D2h||2<=1,                                           (7)

where the latter two norms sum every input/output index in Hilbert-Schmidt norm. Indeed the squared multipliers are n/(n+1)^2 and n(n-1)/(n+1)^2 respectively.

Write X=qY+sigma eta. Then eta is standard Gaussian independent of Y and all source/interpolation noise. Conditional on Y,Z,G, the tensor T2 and the vector m+C_lambda are independent of X. Transferring the k,l derivatives of partial_jF_i onto h times the conditional Gaussian density yields exactly

    E sum_(i,j,k,l) T2_jkl (partial_j F_i)
       [partial_kl h_i
        -sigma^(-1)(eta_k partial_l h_i+eta_l partial_k h_i)
        +sigma^(-2)(eta_k eta_l-delta_kl)h_i].              (8)

No derivative of m, Sigma, b or T2 in Y occurs. This conditional independence is essential; the same calculation cannot be used without it.

The D2h term is bounded by K||T2||2||D2h||2, because DF acts only on the first tensor index and has operator norm at most K.

For each cross-score term, let S_jl=sum_k T2_jkl eta_k. Conditional Gaussian isometry gives

    E||S||HS^2=E||T2||HS^2.

Cauchy-Schwarz therefore bounds each cross term by K||T2||2||Dh||2/sigma.

For the double-score term, put V_j=sum_(k,l)T2_jkl(eta_k eta_l-delta_kl). The second Gaussian Hermite isometry gives

    E|V|^2=2 E||T2||HS^2,                                 (9)

using symmetry in k,l. This is the precise source of sqrt(2), with no additional D. Even though h(X) is correlated with eta, Cauchy-Schwarz suffices after bounding ||DF||op by K. Thus this term costs at most sqrt(2)K||T2||2||h||2/sigma^2.

Together (8)-(9) are bounded by

    K||T2||2 [||D2h||2+2sigma^(-1)||Dh||2
                            +sqrt(2)sigma^(-2)||h||2]
      <=K||T2||2 [1+sigma^(-1)+sqrt(2)sigma^(-2)].           (10)

Integrate (6), use (5), and use self-adjointness of the actual R1 operator to obtain (1).

## 5. Regularity and scope guards

The proof can first be written for a mollified Lipschitz F with bounded derivatives. Gaussian weak chain rules apply to Lipschitz b; R(tau-Sigma) is an L2 spectral operator and does not require derivatives of tau. The displayed bounds survive approximation, so neither D2b nor D2F/D3F bounds are source assumptions. Infinite Gaussian history banks may be treated through their finite strong Gaussian approximations, provided the whole-bank Lipschitz bound and L2 convergence are preserved.

Conditional covariance may be singular; the positive-semidefinite square root in the proof and the G integration identity remain valid. It is only an analytical interpolation. A future executable producer must separately prove its conditional moment match, preserve its original VALUE descendants, expose its actual Gaussian bank, and verify its first/curl/energy, numerical and replay bills.

The final outer operator must genuinely be R1, or a finite approximation must pay its own approximation allowance. The derivative bounds in (7) cannot be silently assigned to an arbitrary frozen outer clock.

No changes were made to the sealed delayed-source theorem, covariance-debt fixture, or their existing scripts.

## 6. Final audit of moment mismatch, the actual F2 tail and the covariance-port join

The extensions in collective-nonlinear-redesign-20261005/outer-resolvent/CONDITIONAL-COVARIANCE-MATCHED-QUARTIC-LEMMA.md and EXISTING-COVARIANCE-PORT-JOIN.md were independently checked and are mathematically consistent, subject to their explicitly hypothetical producer obligation.

### Mean and covariance mismatch

For two conditional Gaussian references, first change the conditional mean while retaining one covariance; this costs K delta_m. For equal means, linearly interpolate the two covariance matrices. The derivative is (1/2)DeltaSigma:D2F. One conditional-X integration by parts gives

    (K/2)delta_Sigma ||Dh||2
       +(K/(2sigma))delta_Sigma ||h||2.

The score term has the exact conditional isometry E|DeltaSigma eta|^2=E||DeltaSigma||HS^2. Thus the additional allowance is exactly

    K delta_m+(K/2)(1/2+1/sigma)delta_Sigma.                 (11)

It does not pay sqrt(D) beyond any dimension already present in the prescribed Hilbert-Schmidt covariance error.

### Actual delayed F2 tail

For q=exp(-delta), Y=X_delta and the true future history H_s, define

    b_true=q integral_0^infinity exp(-t)
                      g(X_(delta+t)-H_(delta+t))dt.

Conditional on Y, it is independent of X_0, and its complete Gaussian future-bank Lipschitz bound is qA(1+A). Its conditional mean is exactly q m_cont(Y), where

    m_cont=R1 E[g(X_0-H_0)|X_0].

The imported file recovery-20261004/cost/endpoint-stein-quadrature-20261004/positive-law-repair/higher-cumulant-gate/THIRD-ORDER-NESTED-FORCE-MEAN-WITH-C2-DECOUPLING.md defines its F_cont by the identical second-substitution force in correlation coordinates. Its actual finite two-root F_Q source has the stated O(A^3sqrt(D)) conditional mean guarantee. Thus reusing that raw mean source here is a correct target identification, not a substitution of a completed mean-law output.

The short F2 prefix differs from w g(X_0), w=1-q, by at most

    (2sqrt(2)/3)A w^(3/2)sqrt(D)+A^2 w sqrt(D)

in force L2 norm. The first term is the exact OU increment integral; the second bounds g(X_t-H_t)-g(X_t). Terminal g adds one A factor. Commuting w g(X_0) into F(u)=g(u-wg(u)) adds at most

    q(1+A) A^3 w sqrt(D).

If a still-hypothetical raw producer has full-bank first O(A) and both conditional mean and Hilbert-Schmidt covariance error O(A^3sqrt(D)), choosing w=A^(4/5) gives the claimed powers:

- A^2 w^(3/2) and A^4/sigma^2: A^(16/5).
- Covariance mismatch A^4/sigma: A^(18/5).
- Short-prefix nonlinear corrections A^3w: A^(19/5).
- Mean mismatch: A^4.

Here sigma^2=2w-w^2 is comparable to w. This balance is valid but remains a conditional opportunity, not an implemented A^(16/5) method.

### Existing LAW-to-RAW distinction and conditional Gaussian budget

The stated carrier-subtraction counterexample is exact: if R_v=sqrt(v+s^2)G and C_v=sqrt(v)G, then R_v has the perfect advertised buffered covariance but Var(R_v-C_v)=(sqrt(v+s^2)-sqrt(v))^2, of order s^4/v rather than s^2. The existing buffered LAW and its small residual first cannot be relabeled as the required raw force covariance.

For X=t z+cG, c^2=1-t^2, and Y=qX+sigma Z with z fixed, direct Gaussian regression gives exactly

    E[X|Y,z]=(t sigma^2 z+q c^2Y)/(1-q^2t^2),
    Var(X|Y,z)=c^2sigma^2/(1-q^2t^2) I.

Thus the possible LAW-only join must budget its complete positive carrier from this actual conditional variance. Scaling the force service by q scales its buffer by q^2, so an allocated X-variance share v_share requires an unscaled reserve buffer v_share/q^2. The variance shrinks to zero at the outer endpoint t=1; the fixed-gap reserve's guards and intrinsic inverse-buffer errors cannot be inherited unchanged.

Two wording guards are necessary for formal application: private Gaussian banks are standard conditional on the retained Y (ordinarily independent of Y), and applying the finite-source lemma to the infinite true tail invokes uniform-Lipschitz finite Gaussian approximation. These guards are compatible with the actual construction. No sealed core source or diagnostic script was edited in this audit.
