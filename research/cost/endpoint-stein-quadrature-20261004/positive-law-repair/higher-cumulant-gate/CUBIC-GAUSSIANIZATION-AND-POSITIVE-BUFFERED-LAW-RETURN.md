# Cubic Gaussianization and a positive buffered full-law return

2026-10-04. New one-energy lemma and its source-qualified positive mean/covariance join. Independent audit in progress. This closes a buffered raw-source law comparison; it does not by itself correct the raw source to the posterior.

## 1. Dimension-safe cubic Gaussianization

Let V be standard Gaussian in any finite input dimension and let f(V) take values in R^D. Assume f is C1, globally L-Lipschitz, and has finite second moment. Put

    X=f(V)-E f(V), e=||X||2, Sigma=Cov(X).

For an independent standard Gaussian Z in R^D and sigma>0,

    W2(Law(X+sigma Z),N(0,sigma^2 I+Sigma))
       <=sqrt(2)/(3 sigma^2) L^2 e.                    (1)

The estimate has one marked physical energy. In particular L=O(alpha), e=O(alpha sqrt(D)), and a fixed Gaussian buffer give an actual O(alpha^3 sqrt(D)) **full-law** bound, including every third and higher cumulant. It is not an assertion based only on equal moments. No second derivative of f is needed, and f need not be a gradient.

Means and covariances in (1) are analytical reference quantities, not queries. The statement works conditionally on any exposed retained caller, with its actual conditional L,e profiles. The private V is integrated, and cannot be appended after this comparison.

## 2. A centered-matrix estimate with an extra Lipschitz factor

Let B(V) be any square matrix with E B=0 and finite L2 Hilbert-Schmidt norm. It need not be symmetric or differentiable. Let S be a fixed symmetric positive definite matrix. For 0<=t<=1, write

    Y=tX+SZ.

Then

    ||E[B(V) S^(-1)Z | Y]||2
       <=sqrt(2) t L ||S^(-1)||op^2 ||B||_(L2 HS).       (2)

Proof by duality. For a smooth compactly supported vector test h,

    E[h(Y) dot B S^(-1)Z]
       =E[B : E_Z Dh(tX+SZ)].                           (3)

This is Gaussian integration by parts in the independent Z, with no symmetry assumption on B. Because B is centered, subtract the V-mean of the matrix on the right. Gaussian Poincare on V and ||Df||op<=L give

    E||E_Z Dh(Y)-E Dh(Y)||HS^2
       <=t^2 L^2 E sum_i ||E_Z D^2 h_i(Y)||HS^2.        (4)

The only Hessians here belong to the test, not to f. For each fixed V, second-Hermite Bessel gives

    ||E_Z D_Z^2[h_i(tX+SZ)]||HS^2
       =||E_Z[(ZZ*-I)h_i(tX+SZ)]||HS^2
       <=2 E_Z h_i(tX+SZ)^2.

Since D_Z^2 h_i=S D^2 h_i S,

    sum_i ||E_Z D^2 h_i(Y)||HS^2
       <=2||S^(-1)||op^4 E_Z |h(Y)|^2.                 (5)

The factor two is the exact degree-two Hermite normalization, including diagonal and off-diagonal entries. There is no input- or output-dimension multiplier. Cauchy-Schwarz in (3), followed by (4)-(5), proves (2) by L2 duality. Approximation extends it from the smooth test class.

## 3. Stein matrix and the covariance-preserving path

Let R_i=grad(-L_OU)^(-1)X_i on the Gaussian V space, with R the matrix of these row vectors. Gaussian Riesz/chaos isometry gives

    ||R||_(L2 HS)^2=sum_(k>=1) ||X_k||2^2/k<=e^2.

Define the analytical Stein matrix

    tau(V)=R(V) Df(V)*.

It obeys

    E[X dot h(X)]=E[tau:Dh(X)],
    E tau=Sigma,
    ||tau||_(L2 HS)<=L e,
    ||tau-Sigma||_(L2 HS)<=L e.                         (6)

Tau need not be symmetric. Its mean is symmetric because it equals the true covariance. No producer evaluates R,tau or Sigma.

Use the law curve

    Y_t=tX+S_t Z,
    S_t=(sigma^2 I+(1-t^2)Sigma)^(1/2).

Its conditional density is a genuine Gaussian convolution at every t. The derivative of its test expectation and the Stein/Gaussian identities show that a continuity-equation velocity is

    v_t(Y_t)=t E[(tau-Sigma) S_t^(-1)Z | Y_t].          (7)

It is enough that (7) gives the correct gradient-test derivative; it need not be identified with the conditional derivative of a particular samplewise square-root coupling. Here S_t and Sigma commute since both are functions of the same Sigma, so that direct coupling also supplies the covariance derivative if desired.

Since S_t>=sigma I, (2) and (6) give

    ||v_t||2<=sqrt(2) t^2 L^2 e/sigma^2.

The dynamic Wasserstein length bound integrated over t in [0,1] proves (1). The start law is exactly N(0,sigma^2 I+Sigma), and the end law is X+sigma Z. The Gaussian buffer keeps the curve smooth; finite moments, the displayed integrable velocity bound, and standard approximation justify the endpoints and C1/Lipschitz regularity passage.

This is the extra cancellation missing from a first-order independent-buffer estimate: centering tau-Sigma is retained, and Gaussian Poincare plus second-Hermite Bessel supplies another first-factor L. No derivative of tau or of the original Hessian occurs.

The same result holds for a frozen anisotropic Gaussian buffer N(0,Q), Q>=qI>0, with bound sqrt(2)L^2 e/(3q). Use C_t=Q+(1-t^2)Sigma and S_t=C_t^(1/2). Gaussian covariance differentiation gives the test derivative E[X dot grad phi]-t E[Sigma:Hess phi] whether or not Q commutes with Sigma; (7) is again a valid gradient-test velocity. Bound ||S_t^(-1)||^2<=1/q in (2). No condition-number factor is introduced. One must not write S'_t=-tSigma S_t^(-1) for noncommuting Q,Sigma; the covariance/test-derivative argument is the correct one.

## 4. Apply to the actual two-stage private source

Use the executed predictor in the companion note, at anchored radius alpha:

    K=K(u,v), K*=K(u-aK-b f(u),v-cK), E=K*-K.

Expose the carrier u and all external callers first. Regard v as fresh private Gaussian. The actual source facts are

    Lip_v K*, Lip_v K<=C alpha,
    ||K*-E_vK*||_(L2(v))<=C alpha sqrt(D),
    ||E||_(L2(v))<=C alpha^2(|u|+sqrt(D)),
    ||Cov_v K*-Cov_v K||HS<=C alpha^3(|u|+sqrt(D)).       (8)

The second line follows from Gaussian Poincare in the D-dimensional private v; the other bounds were proved from the literal VALUE graph, including its fixed-u small curl and known origin.

For 0<h<=1/2, consider the true raw buffered law

    T_raw=h[u-K*(u,v)]+sqrt(1-h^2)N.                    (9)

By (1), conditional on u, it is W2-close to

    N(h[u-E_v K*],(1-h^2)I+h^2 Cov_v K*)

at cost

    C h^3 alpha^3 sqrt(D)/(1-h^2).                     (10)

Every higher private-v cumulant of the actual complete K* has been handled by this law estimate. This does not set them to zero or upgrade a previously retained private v observer.

## 5. Execute the matching positive mean and covariance services

Run the proved conditional-u mean service M(u) from the predictor note. It uses genuine-gradient K and near-gradient E on independent complete banks, explicit positive shares 1/2+1/2, and their executable origins. Under its actual radius/clock/padding guards,

    W2(Law(M|u),N(E_vK*,I))
       <=Lambda alpha^4(|u|+sqrt(D))+floors.             (11)

For the reserve, set d=1-2h^2>=1/2 and eta=zeta=sqrt(d/2). Use the admitted full-gradient forward covariance action on the **anchored executable** source K(u,v)-K(u,0), with padding mu=alpha. Its target matrix is exactly Sigma_K=Cov_v K. With fresh independent complete action tape and incoming roots p,z, execute

    R=eta p+[h^2/(2eta)] C_cov(K-K(u,0);p)+zeta z.       (12)

The independent-buffer covariance-action estimate already used by the curvature reserve gives a law comparison to

    N(0,d I+h^2 Sigma_K)

with conditional error at most

    Lambda sqrt(D)[h^2 alpha^3
                       +h^4 alpha^4(1+alpha^(-1/2))]+floors
       <=Lambda h^2 alpha^3 sqrt(D)+floors.              (13)

For completeness, the first term is the action's calibration ell e mu; the second is its actual private fluctuation current and the positive Gaussian quadratic covariance term. Here ell<=alpha beta, e<=C alpha sqrt(D), mu=alpha, and the scalar reserve d has a fixed gap. The completed action is the admitted finite VALUE response program, not a covariance-oracle instruction. All its source-radius, filter, numerical and finite-clock guards remain imposed at this normalized source.

Execute M and R on independent COMPLETE banks conditional on the same captured u, and return

    T_comp=h u-h M(u)+R.                               (14)

The reference covariance is exactly

    h^2 I+d I+h^2 Sigma_K=(1-h^2)I+h^2 Sigma_K.

Only the whole M output, whole R output and exposed u are used. No coupling Gaussian is identified with an executed private root. Positive reserve d was explicitly budgeted.

Using (8), the fixed-gap Gaussian-root inequality, (10), (11), and (13),

    W2(Law(T_comp|u),Law(T_raw|u))
       <=Lambda [h alpha^4+h^2 alpha^3](|u|+sqrt(D))
                       +complete restored floors.     (15)

Thus (14) is a positive original-VALUE **full buffered law return** for the actual source K*, through cubic order. Integrating the retained carrier u~gamma replaces |u|+sqrt(D) by C sqrt(D). The result prices all higher private cumulants via (1), instead of claiming that mean/covariance matching alone establishes it.

The raw two-stage predictor still has the known rank-five posterior defect. Equation (15) approximates that raw law faithfully; it does not magically repair it to the posterior. Its value is as a legitimate law-level replacement/host interface for a future correctly calibrated current source.

## 6. Exact reverse-OU scaling and terminal restriction

For the full conditional bridge, retain the existing notation s^2=1-r^2, Delta=t^2-r^2, v0=Delta/t^2. The standardized conditional source has alpha=As^2. Its raw affine bridge can be written

    A_rt z+B_rt x+sqrt(v0)
          [h(u-K*)+sqrt(1-h^2)N],
    h=sqrt(Delta)/s.

Therefore apply (14) when Delta<=s^2/4, which makes h<=1/2. The physical conditional error is

    Lambda sqrt(v0)[h alpha^4+h^2 alpha^3]sqrt(D)
                    +scaled finite floors.             (16)

In particular h^2 alpha^3=Delta A^3 s^4, the same covariance scale as the audited curvature reserve. The h alpha^4 term is retained and not silently absorbed at an arbitrarily small increment. No inverse-h replication count occurs.

At a terminal bridge t=1, h=1 and the independent Gaussian buffer vanishes. This lemma and reserve split do not cover that step. The earlier direct affine full-law packet remains valid there, with its separately stated second-order error. No terminal or complete high-order schedule is asserted by this bounded join.

## 7. Actual graph and complete work

All caller-only mode, original g anchor, K(u,0), E(u,0), and f(u) records are captured before private sampling. They are cached only at the identical complete caller/source/finite-version key. The mean M pays its N_K complete K-source occurrences and N_E complete E-source occurrences. The covariance action pays N_cov complete anchored K-source occurrences. A safe original VALUE bill is

    Q_captured+n N_K+(2n+1)N_E+n N_cov
                  +known/numerical work.               (17)

Sharing f(u) or an exact duplicate private site may reduce this conservative bill only after the semantic key is checked. A changed raw argument is a complete replay of its old ancestors. The two means and the covariance service use independent complete banks after the same u is exposed. Their private Gaussian dimension includes every source occurrence, incoming/action root, clock and fill, not just u,v,p,z shown in schematic notation.

First/adjoint actions use original HVPs only at recorded original-gradient VALUE sites. The imported mean and covariance programs supply their actual finite first/caller bounds at the displayed scaled sources; the law comparison is not used to infer those derivatives. Covariance residual first has its literal scale Lambda h^2 alpha^(3/2) for mu=alpha, and the complete mean carries its known unit source-zero carrier plus its Lambda alpha residual first. Captured u and original physical callers retain their actual derivatives through all origins and modes. If a primal record is discarded, its replay is paid.

All raw source zeros are exact under recorded anchor reuse. The completed programs retain their own known Gaussian source-zero carriers; these are never replaced by reference noises from a W2 coupling. At total caller/private zero the anchored sources and residuals vanish. The known carriers and variance split are computed from fixed h,d,eta,zeta. All source versions, counts, filters, padding, public-log clock rules, and numerical precision are fixed before differentiation. Absolute numerical and mode floors propagate through their actual h/readout/inverse-share factors and are never divided by the realized source energy.

At fixed target grade and admitted finite clocks, the new services have only the already proved polynomial-logarithmic occurrence counts. This is no new inverse-heat exponent beyond the complete raw source/provider. It does not yet constitute a full posterior-law or eventual-sublinear c(P) recurrence.

## 8. Diagnostics

The author checker reports 156 passing assertions. It includes scalar CDF W2 diagnostics for a centered globally Lipschitz log-cosh source over fifteen amplitude/buffer pairs, exact second-Hermite normalization in output dimensions 1,2,5,13, rectangular linear-source covariance paths, and positive-join/reverse-bridge variance and heat ledgers. The largest scalar observed ratio to (1) is 0.1671. These are diagnostics; Sections 2-3 are the dimension-safe proof, and the large imported finite mean/covariance circuits are not numerically instantiated by this checker.
