# Independent audit: all-fixed-rank Stein current and one-marked replacement

2026-10-04. **PASS for both analytical statements.** These are dimension-safe law-comparison and coefficient-replacement facts, not finite positive coefficient producers.

Sources:

- `../ALL-RANK-CENTERED-STEIN-CURRENT-RECURRENCE.md`, SHA256 `d09db0fa283512de36bc910f1d3f828257121699832c0897943d9a252f59dd0e`.
- `../ONE-MARKED-THIRD-CUMULANT-REPLACEMENT.md`, SHA256 `5581eb0e592d2997504ac0c69c9f577e46306b6ca85718feb6e0a4a852524a19`.

## 1. The centered recursion and its constants

Let T0=X=f(V)-Ef(V) and define Tj by centering the preceding tensor, applying the Gaussian Hilbert Riesz operator, contracting its input slot with Df, then averaging-symmetrizing the physical slots. Riesz and centering are L2 contractions; multiplication by Df costs L in full HS norm; symmetrization is an orthogonal projection. Thus induction gives `||Tj||2<=L^j e` with one energy e.

The expectations really are connected cumulants, including when lower cumulants are nonzero. For a fixed physical direction theta put Y=theta·X. Repeated centered integration by parts gives, for each polynomial test p,

    E[Y p(Y)]=sum_(j>=1) c_j[theta^(j+1)] E p^(j)(Y),

with only finitely many nonzero derivatives. Taking p(Y)=Y^n yields

    E Y^(n+1)=sum_(j=1)^n [n!/(n-j)!] c_j(theta) E Y^(n-j).

Comparing with the moment-cumulant recursion gives `c_j(theta)=kappa_(j+1)(Y)/j!`. Polarization identifies the full symmetric tensor. No infinite series convergence is needed, and the factor is j!, not (j+1)!.

Consequently `||kappa_(j+1)||HS<=j! L^j e`. The argument uses only the bounded first of f. It does not differentiate Tj or imply that Tj is a legal source for a native compiler.

## 2. Exact same-endpoint law expansion

For the covariance-preserving buffered path Y_t, the first derivative is `t E[(T1-Sigma):D^2phi(Y_t)]`. Each centered Riesz step differentiates the composed test only through t f(V), so it contributes exactly another t and another test derivative. Splitting off each expectation gives

    sum_(j=2)^(m-1) t^j kappa_(j+1):E D^(j+1)phi(Y_t)/j!
       +t^m E[Tm:D^(m+1)phi(Y_t)].

Every term is at the same actual Y_t. The finite expansion neither replaces the endpoint by a Gaussian nor discards any source feedback. At m=2 the constant sum is empty. If kappa3 through kappam vanish, only the final current remains.

Integrating m derivatives through the independent Gaussian buffer yields a vector coefficient with m inverse-root factors and the m-th Hermite tensor. Conditional Hermite isometry costs sqrt(m!), and conditioning on Y_t contracts L2. Integrating t^m from zero to one therefore gives exactly

    sqrt(m!)/(m+1) sigma^(-m) L^m e.

The fixed buffer and finite moments supply the dynamic-Wasserstein regularity, as in the independently audited third- and fourth-order cases. C1/Sobolev approximation is sufficient. For an anisotropic buffer bounded below by qI, covariance differentiation avoids any noncommuting square-root derivative, and the factor becomes q^(-m/2).

This is a valid all-fixed-rank analytical recurrence with factorial constants. It does not claim favorable growth in a varying rank, an executed tensor field, or a positive isolated-current sampler.

## 3. One-marked third-cumulant replacement

For centered Gaussian-image vectors U,V with private first at most L and any HS-unit M, Gaussian Poincare gives

    Var(U* M V)<=2L^2(E|MV|^2+E|M*U|^2)<=4L^4.

The covariance operators of U and V are each at most L^2 I. Hence the covariance operator of the Hilbert-valued product U tensor V has norm at most 4L^4, regardless of its D^2 ambient dimension. The adjoint covariance map from scalar centered L2 into that matrix Hilbert space has norm at most 2L^2.

Applying that map separately to every coordinate of a centered marked E and summing squares yields

    ||E[E tensor U tensor V]||HS<=2L^2 ||E||2.

The exact three-term telescoping of F_c^3-G_c^3 has one marked E_c=F_c-G_c in each term, with F_c or G_c in the other two slots. The tensor slots can be permuted without changing HS norm. Thus the total is at most `6L^2||E_c||2`, exactly as claimed. The derivative never falls on the marked small-energy E, so no small first is being inferred from its energy.

The conditional statement is valid on the same actual Gaussian tape at the same captured caller. At a fresh standard Z, squaring and integrating the conditional profile yields the stated order `A^4sqrt(D)` replacement of fixed-depth F_j by H in the third cumulant. Independence is unnecessary, and replacing the shared record by independent sources would not prove the same tensor identity.

This does not imply a corresponding fourth-order self-covariance bound for F_j-H. The covariance can remain at its `A^3sqrt(D)` allowance. Mean and full covariance targets must still be handled separately.

## Diagnostics and scope

`check_all_rank_stein.py` independently computes the centered Riesz recursion for a scalar Gaussian polynomial with nonzero lower cumulants, checks exact connected constants through rank nine, and verifies the finite same-endpoint expansion across several ranks and polynomial tests. These are algebra diagnostics; the global-Lipschitz norm proof is the Hilbert argument above. Final count and pins appear in its JSON output and manifest.

Both analytical claims pass at the inspected bytes. Together with the discounted-OU hierarchy, this gives two genuinely arbitrary-fixed-rank analytical recurrences. It does not close the finite producer, positivity, source-reentry, mean/covariance or original-query-cost gates of an all-rank posterior algorithm.
