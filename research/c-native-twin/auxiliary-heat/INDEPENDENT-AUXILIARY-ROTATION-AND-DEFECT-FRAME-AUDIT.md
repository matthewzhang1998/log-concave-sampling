# Independent bounded audit of the auxiliary Gaussian rotation

Source read in full: r-native-twin/AUXILIARY-CONSTRAINT-HEAT-AND-VARIANCE-PRESERVING-PATH.md, SHA3953238199cdc91934e689015a17ba01c24e9dd27871aae6fbf5d155b111a2c6. Verdict: PASS for its exact source identities, covariance frames and stated boundary. It does not claim the missing complete endpoint current.

## 1. Exact Gaussian and constraint identities

For D=(I-a²sigma²MM*)^(1/2), D2=(I-a²sigma²M*M)^(1/2), the singular-value calculus gives DM=MD2. Hence the block row matrix

    [D, a sigma M; -a sigma M*, D2]

is orthogonal, so (S',V') is a pair of independent standards. The other original U,Z remain independent. The two representations

    h=sigma V-a sigma² M*(I+D)^(-1)S
     =sigma V'+a sigma² M*(I+D)^(-1)S'

are equal, and aMh=S'-S. No inverse M occurs, including in the rank-deficient case.

Recomputing x0',d',x1' from the COMPLETE new original record (S',U,Z), then using ambient S unchanged, ambient d=d(S,Z), and p_i=g0(x_i')+h, gives exactly

    G_S=E(S',U,Z),
    G_x0=-rh, G_x1=rh,
    G_lambda=r M*(d'-d).

Thus the physical E law, its centered energy and all same-query correlations are exactly preserved. Reusing d or old terminal calls instead of recomputing the full graph would change this identity. The actual full old replay and known orthogonal row/fill arithmetic are required.

## 2. The two covariance estimates

The h regression leaves a known independent innovation V'. Therefore

    Cov(E(R'),G_x1)=ra sigma² E[(E(R')-EE) S'^*](I+D)^(-1)M.

First-chaos Bessel bounds the matrix Hilbert norm by C ra sigma² e. The opposite block has the opposite sign. This does not license cancellation under unequal observers.

For fixed original Z, write H_Z(S)=g_t(S)-g_t(S+epsilon Z). The lambda defect has conditional mean zero, because both S and S' have the same standard marginal. For each fixed test u, put F(S)=(Mu)*H_Z(S). Then Lip F<=2||M u||<=2|u|.

The exact correlated-Gaussian Dirichlet inequality used here is

    E|F(S')-F(S)|² <=2||I-D||op E|DF(S)|².                 (1)

To verify (1), diagonalize the known positive D. A normalized Hermite multi-index alpha has correlation product_i d_i^(alpha_i); use

    1-product_i d_i^(alpha_i) <=sum_i alpha_i(1-d_i)
                              <=||I-D||op |alpha|,

then sum the orthogonal coefficients. This is dimension-free and needs no third derivative. It remains valid for the supplied Lipschitz functions by Gaussian Sobolev closure.

Multiplying (1) by r²a² gives

    E|u*G_lambda|² <=8r²a²||I-D||op |u|²
                      <=C r²a^4 sigma²|u|².

Hilbert/operator factorization with the complete actual E(R') mark yields

    ||Cov(E(R'),G_lambda)||HS <=C r a² sigma e.

No independence between E(R') and G_lambda is used. The fixed-test row frame, exact conditional centering, and the original E energy are all retained.

## 3. Firsts, weak current and scope

Complete old-root first of G_lambda is O(ra), while its V first is O(ra²sigma). The caller includes both complete evaluations of H_Z, including the changed S' path. Its own unmarked Hilbert energy can contain sqrt(d); the displayed mixed covariance estimate does not eliminate it.

For the separate unrotated common-p path, the conditional V cross covariance with the opposite x defects has the source's exact sign and factor r²a sigma². Its width differentiation has the actual terminal V-Laplacian term. The weak Gaussian integration-by-parts formula is valid as an analytical identity and does not turn that third-jet expression into an allowed VALUE producer.

The rotation provides a feasible law-preserving constrained path and two small one-energy covariance fields. A SAME-endpoint nonlinear observer still requires the full current/proper-cut return; multiplying these fixed-test bounds by a random matrix or discarding their common roots is not certified.

## 4. Checks

The adjacent check_auxiliary_calibration_boundary.py independently tests 1,440 known Gaussian rotation/regression/constraint identities in dimensions1,2,4,8 with arbitrary nonsymmetric contraction M, including small singular values. Maximum orthogonality residual is1.04e-14 and h-regression residual5.1e-15. The analytical argument above, not those floating-point diagnostics, proves the covariance frames.
