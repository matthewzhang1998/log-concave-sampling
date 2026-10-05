# Centered covariance fails; a bounded resummed m3 subclass is constructive

2026-10-05. Independently audited continuation of `FULL-COVARIANCE-CLOSURE-AND-RESUMMED-MEAN-FRONTIER.md`. The general mean gate remains OPEN. The previously admitted full same-carrier covariance service is unchanged and was not rederived here.

## 1. A genuine same-g failure of centering alone

The unshifted covariance repair was already refuted. The new results refute its CENTERED version as well, including the literal P3/F3_Q target.

Use a smooth radial original force

    g(x)=c A x+d A psi(|x|-R0) x/|x|,
    c=1/2, d=1/4, A=D^(-1/4), R=sqrt(D)=A^(-2),

with the smooth monotone step and origin definition stated in the source proofs. The global Hessian satisfies 0<=Dg<=A I. Set R0 at the actual coherent linear mean radius: (1-cA/2)R for single I, or (1-a2)R for the second substitution, with a2=cA/2-c^2 A^2/4.

In either case the true and cheap conditional Gaussian backbones have distinct scalar variances. Their centered noise has small operator covariance O(A^2), but its radial inflation is O(A^2 sqrt(D))=O(1). Consequently the actual terminal force probes

    d A [psi(S+h_H)-psi(S+h_Q)] n,

while a single covariance derivative at the coherent center probes only

    d A (h_H-h_Q) psi'(S) n.

Their difference survives the outer R1 resolvent with norm Omega(A). The proposed allowance is A^4 sqrt(D)=A^2, so every fixed public-polylogarithmic loss still fails.

The direct P3 statement uses the literal finite three-root F3_Q and proves failure for BOTH

    E[F3_Q|Z]+(1/2)R1[(Cov(F2_H|x)-Cov(I2_Q|x)):
                                           D2g(x-mu_*(x))],

and the corresponding first-displacement covariance correction. It permits the actual true or cheap conditional mean as center, any pointwise convex combination, and theta-averaged coherent centers. Exact ordinary mean-correction terms remain O(A^2) and do not repair the gap. The finite outer rule has delta_out<=A^3; the middle rule has delta_mid<=A^2. Floors and public versions are polynomially priced along the counterfamily.

These are actual same-original-g results with the genuine continuous OU ancestry and literal finite common-root ancestry. The remainder correlations are retained. The trace is evaluated explicitly. No generic scalar-test or tensor-Bessel contraction is assumed.

## 2. A positive, finite, fully resummed canonical m3 consumer

There is also a constructive result, rather than only another separator. Assume a DECLARED bounded linear-remainder certificate for the actual source:

    g(x)=lambda x+r(x), 0<=lambda<=A, sup_x |r(x)|<=epsilon,

alongside the original convex-gradient and zero-anchor assumptions. On the actual history,

    F2=lambda H1-lambda^2 H2+R,
    |R|<=epsilon(1+lambda) pathwise,

where H1=int e^(-t)X_t dt and H2=int t e^(-t)X_t dt are analytical comparison variables. The conditional Gaussian backbone is exactly a2 x+b2 N, with

    a2=lambda/2-lambda^2/4,
    b2^2=lambda^2/4-lambda^3/2+5lambda^4/16.

Thus the FULL original-g response

    R1 E_N g((1-a2)x-b2 N)

approximates the canonical m3 target with uniform absolute bias at most A epsilon(1+lambda). No covariance Taylor series is used.

A finite positive clock rule gives the literal original-VALUE source

    F_G(Z,G,N)=sum_i w_i g((1-a2)(t_i Z+c_i G)-b2 N),
    c_i=sqrt(1-t_i^2),

with fresh independent G,N, shared across their respective full levels, conditional on the same retained original Z. The baseline B=sum w_i g(t_i Z+c_i G) is a genuine G-gradient. E_G=F_G-B has energy O(A^2(|Z|+sqrt(D))), actual first O(A), curl O(A^2), and explicitly priced nonzero fixed-Z origins. Its G derivative is symmetric; its only off-diagonal lift block has norm at most A b2. Every owned root and changed original VALUE argument is included in the finite program and integrated by the complete mean-law return.

The admitted fixed-order gradient/near-gradient mean compiler gives a positive M_G targeting N(m3(Z),I), with integrated conditional-law error

    A epsilon(1+lambda)+C delta A sqrt(D)
                       +Lambda A^4 sqrt(D)+absolute floors.

The safe complete original VALUE bill is

    Q_captured+N_out N_B+2 N_out N_E
                              +known/numerical/replay work,

where N_B,N_E are the fully expanded native occurrence counts at their actual normalized radii. N_out=O(log^2(1/delta)); use delta<=A^3. All fixed shares, guards, modes, known scalar roots, numerical versions and caller derivatives remain explicit. Original HVPs appear only at recorded VALUE sites in requested first/adjoint sweeps. No original path, mean, Hessian, covariance matrix, or saved HVP is a producer leaf.

If epsilon=O(A) and D>=A^(-4), the displayed bound is order A^4 sqrt(D), including the new radial counterfamily. It also applies whenever the explicit absolute bias A epsilon is within the desired allowance. This is a proved source SUBCLASS, not general source closure. Its certificate is not automatically preserved at the needed grade by conditional rescaling: f(y)=s[g(a+s y)-g(a)] has lambda_f=lambda s^2 and only the safe epsilon_f<=2s epsilon bound. No endpoint join follows merely from the unscaled source assumptions.

## 3. What remains open

The exact centered/resummed Stein identity remains valid. A general finite positive VALUE consumer for its complete shifted target is still missing. The new separator rules out replacing it by a deterministic covariance derivative even at the correct coherent mean, including use of the actual full F2/I2 covariance.

The bounded construction identifies a valid escape: keep a justified Gaussian backbone inside the original nonlinear force, with its full covariance response. In the unrestricted source class there is no declared uniformly bounded linear remainder, so its pathwise bias proof does not extend for free. Nor may a completed mean-law output be read as a strong conditional mean or a fixed-buffer covariance return be read as a small-noise covariance matrix.

General m3 mean closure: OPEN.
General fourth-order endpoint: OPEN.
Arbitrary-order recurrence and c(P)/P -> 0: NOT ESTABLISHED.

## 4. Frozen proof pins and independent checks

All paths are relative to this folder.

- `RADIAL-INFLATION-REFUTES-CENTERED-SINGLE-HISTORY-COVARIANCE.md`
  SHA256 c6024a6ba2ca4fbb293063d0746ee92fa083398a3bf4c01fb4de1908f322dfaf.
  Independent audit: `independent-center-inflation-audit/INDEPENDENT-CENTER-INFLATION-AUDIT.md`, SHA256 91617318aec9f2625300b4bce64494e3814a6ad8c698a73c04f3bdc18667a3da. 723 independent assertions.

- `RADIAL-INFLATION-REFUTES-CENTERED-P3-COVARIANCE-CORRECTION.md`
  SHA256 810128a4b47c568273f05442a510498ac5f77eedf1cb864366368ebb4bd356af.
  Independent audit: `independent-centered-p3-inflation-audit/REPORT.md`, SHA256 975f6ae0bdcc09218ae238603140a0c505599360ef993febf540c796d4fde5b9. 473 independent assertions.

- `resummed-linear-backbone/BOUNDED-REMAINDER-RESUMMED-M3-VALUE-CONSUMER.md`
  SHA256 798d20b87124e7069e465ce4e7bcf92c5af3734d86c8a62aeb288b5fec8eaef3.
  Independent audit: `resummed-linear-backbone/independent-audit/INDEPENDENT-BOUNDED-RESUMMED-M3-AUDIT.md`, SHA256 49d135ba6692c54012acd9d0d94f603c0015fca98b5d1be1164457ff45d1213a. 1,436 independent assertions.

Total independent diagnostic assertions: 2,632. They supplement the written proofs and pinned compiler imports; they do not execute an infinite history or the imported completed mean compilers. The three audit directories retain their complete manifests and checksum lists. The companion continuation manifest inventories this new packet without modifying the previously frozen covariance frontier.
