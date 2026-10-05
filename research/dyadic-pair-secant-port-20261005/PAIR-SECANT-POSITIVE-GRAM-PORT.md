# A six-VALUE interacting dyadic source and a conditional-Stein Gram port

2026-10-05. This is a bounded original-VALUE construction for genuine two-block interaction fields. It supplies a positive Gaussian covariance law for those finite fields, keeps the critical A^4 pair grade, and records its actual nonlogarithmic accuracy cost. It does not substitute a conditional ancestor mean inside g and does not execute a history path grid.

## 1. Result and honest target

Let g=grad U, U in C2, g(0)=0, and 0<=Dg<=A I. Freeze retained callers a,b in R^D and scalar rows sigma,tau,c>=0. Let xi,xi',zeta,zeta' be four independent standard D-Gaussians. Define x=a+sigma xi, x'=a+sigma xi', z=b+tau zeta, z'=b+tau zeta'. Execute

    u=g(z), u'=g(z'),
    v11=g(x-c u), v21=g(x'-c u),
    v12=g(x-c u'), v22=g(x'-c u'),
    Psi=(v11-v21-v12+v22)/2.                         (1)

There are exactly SIX original g VALUES. This is a nonlinear interacting field, not an additive pair of singleton sources. Every displayed descendant uses its saved same-bank ancestor.

Put

    ell=c A^2 tau/sqrt(2),
    L=sqrt(A^2 sigma^2+c^2 A^4 tau^2)/sqrt(2).

The exact pair-Hoeffding covariance C(a,b) of q(xi,zeta)=g(a+sigma xi-c g(b+tau zeta)) satisfies

    E Psi=0,  E[Psi Psi*]=C(a,b),
    0<=C(a,b)<=ell^2 I,  ||C||HS<=ell^2 sqrt(D).       (2)

For positive buffer v and integer K>=1, execute K conditionally independent COMPLETE copies of (1), all with the same retained callers, and a fresh standard D-Gaussian G:

    Z_K=sqrt(v) G+K^(-1/2) sum_(j=1)^K Psi_j.          (3)

Then

    W2(Law(Z_K|a,b), N(0,v I+C(a,b)))
        <= ell^2 sqrt(D/(v K)).                      (4)

No matrix, Jacobian, conditional mean, covariance, or Stein kernel is an executing leaf. The implementation of (3) uses 6K original VALUES and (4K+1)D private scalar Gaussian coordinates. Its covariance is EXACTLY vI+C, and its remaining law error is (4).

Thus K=ceil(epsilon^-2) attains the buffered covariance tolerance epsilon ell^2 sqrt(D)/sqrt(v). This is a finite polynomial accuracy cost, not a public-logarithmic native action. For epsilon=A^gamma it has the explicit VALUE exponent 2gamma. The complete private residual first is L; the retained caller firsts are at most sqrt(K)A and sqrt(K)cA^2. In particular small covariance does not create a false small full/caller first.

## 2. Exact pair-Hoeffding identity

For any square-integrable vector field q(xi,zeta) of two independent blocks, write q=q0+q1(xi)+q2(zeta)+q12(xi,zeta), with each nonconstant component centered in each of its own blocks. Its rectangular difference is

    Delta12 q=q(xi,zeta)-q(xi',zeta)
                   -q(xi,zeta')+q(xi',zeta').

It annihilates q0,q1,q2. The four surviving q12 terms are pairwise orthogonal: if they share one block, average the other centered block; if they share none, use independence. Consequently

    E[(Delta12 q)(Delta12 q)*]=4 E[q12 q12*].           (5)

Equation (1) is Delta12 q/2. Hence (2) targets the genuine pair component of this finite shifted field, not its full variance and not a singleton replacement. The factor 1/2 is essential.

The positivity is exact even when matrix products fail to commute. Positivity of (2) does not say that an arbitrary history covariance difference is PSD.

## 3. Conditional small derivative, large complete derivative

Use H=Dg only in the proof, at the six recorded sites. Since all H satisfy 0<=H<=AI, every difference of two such matrices has operator norm at most A. Preserve matrix order. For example,

    D_zeta Psi=-(c tau/2)(H11-H21) H_z,
    D_xi Psi=(sigma/2)(H11-H12).                      (6)

The primed formulas have the corresponding signs. Thus

    ||D_(zeta,zeta') Psi||op<=ell,
    ||D_(xi,xi') Psi||op<=A sigma/sqrt(2),
    ||D_all Psi||op<=L.                              (7)

For EVERY fixed outer pair (xi,xi'), E_(zeta,zeta') Psi=0, by exchanging zeta and zeta'. Gaussian Poincare in those two inner blocks proves Cov(Psi|xi,xi')<=ell^2 I. Averaging proves (2), with a single operator-to-HS conversion.

The same algebra gives

    ||D_a Psi||op<=A,   ||D_b Psi||op<=c A^2.          (8)

In (3) every caller is shared among K copies, so the safe deterministic caller bounds are sqrt(K) times (8). The private banks are independent concatenated coordinates; Cauchy-Schwarz over their 1/sqrt(K) weights gives L with no sqrt(K) loss. Including the keep, the full private first is at most sqrt(v+L^2).

The scalars c,sigma,tau,v and affine row versions are frozen mode labels. Differentiating them produces additional Gaussian-moment paths; for example D_v Z=G/(2sqrt(v)). Such paths are not globally bounded and cannot be included in (8) without their profiles. At c=0 or tau=0 the source is literally zero. At the zero Gaussian tape all four outer sites coincide, so the residual origin is exactly zero; numerical origin/capture consistency is retained.

## 4. A conditional-Stein proof of the positive law

This section is analytical only. It is useful because it exploits the small INNER derivative while preserving the larger complete first.

Lemma. Let O be an arbitrary background, U a standard Gaussian independent of O, and F(O,U) in R^p satisfy E[F|O]=0 and ||D_U F||op<=ell. Then F has a (possibly nonsymmetric) Stein matrix T with

    E[F_i phi(F)]=E[sum_j T_ij partial_j phi(F)],
    E T=Cov(F),  ||T||op<=ell^2,
    E||T-Cov(F)||HS^2<=p ell^4.                     (9)

For each fixed O use its Gaussian OU semigroup P_s on U alone and set

    T0(O,U)=integral_0^infinity e^-s
                  (P_s D_U F)(O,U) (D_U F(O,U))* ds.

The standard Gaussian resolvent integration by parts gives the first identity before conditioning. Conditional mean zero removes the zeroth-chaos term. Both factors have operator norm at most ell, so ||T0||op<=ell^2. Conditioning T0 on F supplies T. Taking linear test functions gives its mean; orthogonal projection and centering give the last bound. No symmetry or PSD property of this analytical Stein matrix is needed.

Apply the lemma with O=(xi,xi'), U=(zeta,zeta'). For (3), independence of complete copies gives an analytical Stein matrix

    vI+(1/K)sum_j T_j,

before conditioning on Z_K. Its mean is Sigma=vI+C, and its centered squared HS norm is at most D ell^4/K. The keep contributes exactly vI.

For completeness the anisotropic Stein transport bound is

    W2(Law(Z),N(0,Sigma))
       <={E||[T_Z-Sigma]Sigma^(-1/2)||HS^2}^(1/2)
       <=lambda_min(Sigma)^(-1/2)
                        {E||T_Z-Sigma||HS^2}^(1/2).   (10)

It follows by interpolating sqrt(t)Z+sqrt(1-t)G_Sigma. Gaussian integration by parts expresses the velocity as

    [2sqrt(1-t)]^(-1)
       E[(T_Z-Sigma)Sigma^(-1)G_Sigma | interpolant].

Its L2 norm is bounded by the first quantity in (10) divided by 2sqrt(1-t), whose integral on [0,1] is one. Regularization and a limit justify nonsmooth sources and endpoint velocities. Since Sigma>=vI, (9)-(10) prove (4) with the displayed constant.

This method neither calls an unadmitted rectangular square action nor labels the signed Jacobian of Psi a convex-gradient source. It constructs the law directly by finitely many original VALUES and Gaussian affine operations.

## 5. Marginalized finite environments and positive excess

A finite source may also depend on a nuisance bank E, independent of the four midpoint blocks. For instance

    x=a+sigma xi+kappa N,  x'=a+sigma xi'+kappa N,
    z=b+tau zeta+lambda M, z'=b+tau zeta'+lambda M,

with E=(N,M); the SAME N or M is used in its paired branches. More generally include a declared finite original-VALUE background graph, with its full readset and costs. Require the conditional source to retain the same centering and inner first ell.

The desired two-midpoint field may be H(xi,zeta)=E_E q(xi,zeta;E). Its exact mixed-difference source is barPsi=E_E Psi. To approximate it WITHOUT a conditional-expectation leaf, take R independent nuisance banks E_r while sharing the four midpoint roots and set

    Psi_R=(1/R)sum_r Psi(xi,xi',zeta,zeta';E_r).

The exact covariance identity is

    Cov(Psi_R)=Cov(barPsi)+(1/R)E Cov_E(Psi|midpoints).
                                                               (11)

Therefore the excess is PSD, and, using E PsiPsi*<=ell^2 I,

    0<=Cov(Psi_R)-Cov(barPsi)<=ell^2 I/R.              (12)

The INNER first of Psi_R is still at most ell and it remains conditionally centered after freezing the outer roots and EVERY nuisance bank. Gaussianizing K independent complete copies of Psi_R as in (3) consequently yields

    W2(Law(Z_(K,R)),N(0,vI+Cov(barPsi)))
       <=ell^2 sqrt(D)/sqrt(v) [K^(-1/2)+(2R)^(-1)].   (13)

The second term uses the synchronous principal-square-root Gaussian coupling and ||sqrt(B)-sqrt(A)||HS<=||B-A||HS/(2sqrt(v)) for A,B>=vI. No cancellation of the excess is assumed.

For the displayed two-point residual bank the bill is 6KR original VALUES and [K(4+2R)+1]D Gaussian coordinates, before legal same-key savings. An explicit background with Q_B VALUES and n_B additional roots contributes its ACTUAL complete replay count and roots per nuisance draw. Its derivatives through a,b are paid using (8). If a background is retained and common to all K copies, its caller amplification has the sqrt(K) factor; if privately regenerated, it belongs to the full private first of each copy. A true F1 history is not a finite background leaf.

## 6. Finite positive combinations and coherent blocks

For a finite family of frozen row modes indexed by l, use the SAME four midpoint roots in all terms and put Psi_sum=sum_l w_l Psi_l, with nonnegative w_l. It is the genuine mixed difference of the corresponding sum of shifted fields. In particular the cross terms between modes are retained; independently Gaussianizing each mode would target only the sum of its covariances and would in general be wrong.

A valid conditional inner radius is

    ell_sum=(A^2/sqrt(2))sum_l w_l c_l tau_l,

and the full outer radius is (A/sqrt(2))sum_l w_l sigma_l, with all background paths added. The preceding Gram and LAW theorems apply with ell_sum, and each complete copy costs 6 times the number of executed modes, including every original ancestor and anchor replay. For a p-slot vector block use the concatenated physical output and its actual inner operator bound; the Stein/HS bound is ell_sum^2 sqrt(p). A later physical readout R costs its operator norm in W2 and transforms the covariance by R C R*.

A weighted field family is not independent of its environment merely because the integration weights are positive. All saved modes, rows, copies and branch keys are frozen first. A downstream test may depend on retained callers or unread labels; it may not inspect a consumed midpoint, nuisance or Gaussianization bank after the conditional-law comparison.

## 7. Connection to the true F2 pair grade

This is an exact analytical connection, with its regularity scope stated. Let g=A f with fixed smooth f. For each fixed finite dimension and fixed dyadic partition, the actual coherent F2 has

    F2=A Rf(X)-A^2 B2+O_L4(A^3),
    B2=R[Df(X) Rf(X)].                               (14)

The error constant may depend on f and dimension; (14) is NOT a uniform unrestricted-C2 estimate. It follows by Taylor expansion of the actual ancestor, without replacing it by its mean.

For two ordered distinct cells I<J, conditional independence of their bridges gives an exact FACTORIZATION of the B2 pair projection:

    (B2)_{I,J}=A_I(xi) B_J(zeta),
    A_I(xi)=integral_I [E(Df(X_t)|coarse,xi)
                                 -E(Df(X_t)|coarse)] dt,
    B_J(zeta)=integral_J e^-s [E(f(X_s)|coarse,zeta)
                                 -E(f(X_s)|coarse)] ds.           (15)

The matrix A_I precedes the vector B_J. Explanation: B2 is the ordered two-time integral integral_(t<s)e^-s Df(X_t)f(X_s) dt ds. A pair projection vanishes unless its two sites are in the two selected cells; ordering then forces t in I and s in J. All other path coordinates integrate out locally. The factors in (15) are centered and conditionally independent. Their exact covariance is the positive ordered sandwich

    E_xi[A_I Cov(B_J|coarse) A_I* |coarse].            (16)

This retains the noncommuting kernel order and demonstrates why a pair field, rather than a singleton, is the correct first interacting object. It does not assert that the full nonlinear F2 has only pair interactions.

At any fixed two sites, the six-VALUE source has

    Psi=-(c A^2/2)[Df(x)-Df(x')][f(z)-f(z')]
                                         +O_L4(A^3).             (17)

For the original scalar witness f(x)=x/2+(1/4)log cosh x and sigma,tau,c>0,

    E Psi^2=c^2 A^4 Var(f'(a+sigma xi))
                         Var(f(b+tau zeta))+O(A^5),               (18)

and both variances are strictly positive. Thus this finite original-VALUE source carries the same genuine A^4 pair covariance mechanism and cannot be dismissed as an irrelevant higher-order atom. Positive finite sums with residual averaging in Sections 5-6 supply the finite-field analogue of (15), while keeping its covariance cross terms.

Passing (14)-(18) to the full general-C2 history target requires a UNIFORM tail/remainder estimate and actual finite clock/environment approximation. The explicit bounded port above does not assume that estimate. General smoothness alone does not justify freezing a derivative constant as A or D changes, and no epsilon clock can erase an unpriced ancestry error.

## 8. The exact positive current and the consumed keep

Let M_K=K^(-1/2)sum Psi_j and Z_theta=sqrt(v)G+theta M_K for 0<=theta<=1. This is an exact positive finite original-VALUE law path with random mark M_K. For every smooth compactly supported test phi,

    d/dtheta E phi(Z_theta)=E <Dphi(Z_theta),M_K>.

Its covariance is vI+theta^2 C; its mark energy is at most ell sqrt(D), while its full first remains L. Its marginal conditional current can be written E[M_K|Z_theta], but execution returns the actual joint (M_K,Z_theta) and never queries that conditional expectation. The Gaussian approximation (4) is only a MARGINAL comparison after consuming all midpoint banks and the private keep G. It is not a retained-mark approximation, not a same-caller-readable-G comparison, and not a proof of a joint current compiler. A test reading any consumed root requires a separate conditional or joint theorem.

## 9. What this closes, and what remains outside this port

Closed here: exact six-VALUE pair interaction; its PSD Gram; small conditional-inner radius under the original C2 assumptions; finite positive buffered Gaussian realization with one-HS/dimension-safe error; finite marginalization with a paid PSD excess; complete/caller firsts, root counts and source replay. The finite family retains cross-mode Gram terms and noncommuting matrix order.

Not identified with this finite field: the full conditional true-history pair (F2,F3-F2), its whole covariance, or its cubic/quartic marked proper-cut packets. Those histories require their own finite environment supplier and uniform restoration. K and R are polynomial inverse-accuracy counts, not new logarithmic native complexities. The Gram field does not independently admit the old gradient-only square primitive, and its Gaussianized law is never relabeled as a RAW history after subtracting the keep.

The separate multi-interval compression result in this packet states exactly which translation clocks remain dimension-free and why generic relative-gap clocks do not share that theorem. No sealed input, external document or publication is changed.
