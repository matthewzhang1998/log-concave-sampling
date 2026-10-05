# Bounded old-bank return: second-stage audit PASS

2026-10-05. This addendum supersedes the remaining-blocker verdict in `INDEPENDENT-OLD-BANK-SIXTH-AUDIT.md` for the **new frozen full construction**. It does not change the correctness of that earlier review of the earlier proposal.

## Frozen verdict

**PASS for the bounded canonical zero-old-bank-readout return, relative to the explicitly imported fixed-rank heat-cut, selected-pair, finite-filter/calibration, and conditional-quartic native contracts.** The normalized positive-law remainder is Lambda alpha^8 sqrt(D), plus the separately priced absolute floors and inherited original-clock target mismatch. The intrinsic alpha-eight conditional offspring is paid, not recursively closed.

Reviewed full construction: `REVIEWED-FULL-RETURN.md`, SHA256 `9191c93102e9c455ada1027984eae5cf3ded613a7a37869e5b9b5bd37c8ae2c8`.

Reviewed revised OU lemma: `REVIEWED-OU-REVISED.md`, SHA256 `36acc99de3596a4a1c392696b6c054fa7271e2fed833708a9cf3af3e1c8f847f`.

The earlier conditional-channel supremum issue is genuinely bypassed. The new proof neither infers a bound on sup_Q from a pointwise-Q Gaussian polynomial moment nor differentiates the rough old coefficient twice.

## 1. The second Riesz identity is exact

Fix the physical Gaussian publics P. Write F(Q,P) for the old centered cubic displacement, R_Q=D_Q N_Q^-1, and

    tau=(R_Q F)(D_Q F)^T,
    C(P)=E_Q tau=E_Q FF^T.

The latter identity is Gaussian covariance duality. Tau need not be symmetric pointwise; only its symmetric part acts on the test Hessian, and its expectation is the PSD matrix C(P).

Use the positive path

    W_t=base(P)+tF(Q,P)+sqrt(kappa/2)G0
          +[(kappa/2)I+(1-t²)C(P)]^(1/2)G1.

All coefficient banks and P are independent of G0,G1. Gaussian differentiation of the last covariance gives -t C(P):D²phi at the same W_t. Gaussian integration in Q of F dot grad phi gives t tau:D²phi. Since tau-C is centered in Q and the rest of the base/root is Q-independent, a SECOND application of covariance duality gives

    d/dt E phi(W_t)
       =t E[(tau-C):D²phi(W_t)]
       =t² E[(R_Q(tau-C)) D_QF :D³phi(W_t)].

The second Riesz operation acts on the centered matrix coefficient; its integration-by-parts derivative hits the test composed with W_t. No derivative of tau is taken. No D_Q²F estimate or conditional uniform Lipschitz extension is used.

F is homogeneous physical chaos two. The product tau and its Q-centered version have physical chaos degrees at most four. The coefficient-level Riesz contraction and the first-stage all-matricization physical-polynomial lemma give

    ||tau-C||_(L2(Q,P;HS)) <= Lambda alpha^6 sqrt(D).

Expand this finite polynomial in an orthonormal physical Wick basis before applying R_Q again. R_Q is an L2 contraction from each centered coefficient Hilbert space to its Gaussian derivative Hilbert space. Conditional physical hypercontractivity controls the fourth moment of the resulting degree-at-most-four field at each Q. The other factor D_QF has, uniformly in Q, the L4(P;operator) matrix-chaos bound Lambda alpha^3. Thus

    ||(R_Q(tau-C))D_QF||_(L2(Q,P;HS))
       <= Lambda alpha^9 sqrt(D).

Exactly one HS energy is spent. In particular this proof never estimates the derivative of the random physical polynomial by multiplying an independent sqrt(D)-sized energy.

Two integrations through the constant G0 keep convert this rank-three current into a vector continuity current. The coefficient is independent of G0. Gaussian Hermite isometry and the fixed kappa/2 gap preserve its one-HS L2 bound. Conditional expectation onto W_t is an L2 contraction, and the continuity-velocity bound integrates over 0<=t<=1. This proves the positive-law remainder Lambda alpha^9 sqrt(D).

The proof is still valid when the base contains Q-independent mean and skew-square fields. Q-dependent sixth fields are centered FIRST, with their alpha-nine old/new mixed currents retained, as stipulated in the full construction.

## 2. Same-endpoint rank-six, rank-four, and rank-two cancellation

For unit covariance notation, the exact identity is

    H2_ab H2_cd = H4_abcd
       +delta_ac H2_bd+delta_ad H2_bc
       +delta_bc H2_ad+delta_bd H2_ac
       +delta_ac delta_bd+delta_ad delta_bc.

The inverse-covariance factors in the general convention give the same four single and two double contractions. Therefore half the covariance of K_Q:H2, acting on a test Hessian, has leading coefficients

    (1/2) Cov(K_Q,K_Q) at rank six,
    2 Cov-M at rank four,
    Cov-N at rank two.

There is no missing factor two and no commutative-matrix assumption. The difference between the average conditional skew-square maps and the map of the average K cancels the latter two. The new six-mark packet cancels the rank-six term. A separate skew-square correction to the already chosen cubic normal form would double count; the full construction correctly forbids that addition.

The final path is explicitly

    X_a=P+K:H2(P)+S(K)(P)+aD(P)+sqrt(kappa/2)G0
          +[(kappa/2)I+aC(P)]^(1/2)G1,

where D=E S(K_Q)-S(K)-(1/2)Cov(K_Q,K_Q):H5. Its derivative is exactly

    E[D(P) dot grad phi(X_a)+(1/2)C(P):D²phi(X_a)].

Every integration by parts occurs at this same X_a. Replacing D_P X_a by the identity isolates the displayed exact leading cancellation; it is not the whole proof. Every leftover chain-rule term contains a derivative of K:H2 (grade three), a sixth-order correction (grade six), or the gapped covariance root (grade six). Multiplying the original grade-six current gives grade nine or twelve respectively.

The root derivative estimate is substantive. If S(P)^2=(kappa/2)I+aC(P), then

    S (D S)+(D S) S=a D C,
    D S=integral_0^infinity exp(-vS)(a D C)exp(-vS) dv.

The fixed spectral gap bounds this linear map in HS and in each needed mixed physical cut. Further derivatives solve the same Lyapunov equation with D^j C and ordered products of earlier root derivatives. Only finitely many orders are needed here. They form matrix products rather than traces; noncommuting order is retained. C(P) is an explicit degree-four polynomial with the declared coefficient cuts, so these terms use no additional original derivative oracle. Matrix-chaos estimates handle their physical Gaussian factors, including G1. Their remaining test derivatives move ONLY through the independent constant G0. Using G1 itself as an untouched keep would be incorrect; the final frozen version makes the split explicit.

This verifies the finite same-endpoint continuation claimed in Sections 7–8, rather than merely a matching of formal cumulants.

## 3. Source normalization, geometry, and cost corrections checked

The source is now explicitly

    [g(z+t x)-g(z)]/(A t),

so alpha=qA/sqrt(u) occurs only in the service/native amplitudes. For the tight C0^4 C2^2 graph, b0=s0/sqrt(2), b2=s0/2, and the uncombined 36-permutation factor is

    d=8 gamma/(36 b0^4 b2²)=32 gamma/(9 s0^6).

Multiplication by all 36 permutations recovers the old unsigned 8 gamma coefficient exactly. The old center weight is squared, the original old derivative injection is A rather than sqrt(c) A, and the covariance redistribution preserves the exact old heat and genealogy.

The source-cut proof retains the actual fifth-degree physical Wick field and its bank derivative matricizations, not just a tensor-vector surrogate. The shortest new spine has four forces and therefore grade-eight conditional feedback. The old/new mixed current is grade nine and the new whole-bank own current is grade twelve. These distinctions are essential.

The revised finite bridge uses a conservative numerical discrepancy bound delta_op M_old² D, and chooses delta_op accordingly. Its node count stays logarithmic squared in the requested error and recorded inverse scales. This numerical precision tightening is separately priced; it is not a hidden extra physical energy in the substantive current estimates.

The final cost is nonincremental: complete old program at its required native order/precision, complete conditional-quartic program, complete new sixth program, and all capture, known-root, readout, scalar-encoding and actually owed replay work. Reuse is allowed only when the exact cached caller/source/native-order/precision version already meets the new floors. A previously paid alpha-four floor is not retroactively promoted to alpha eight. Per ordered old-node pair and bridge node, the exact nine-hit source census is 24 C0+24 C1+6 C2 before the declared physical permutations; shared-bank cross pairs remain explicit.

No endpoint deletion occurs in the admitted route. Hence there is no unpriced deleted-node restart hidden in the verdict. The exact symmetric OU heat redistribution replaced that abandoned proposal.

## 4. Independent evidence and limits

The independent diagnostics comprise:

- 12,000 endpoint-heavy old genealogies, PSD splits, and covariance reconstruction checks.
- Four positive dyadic Gauss rules, their spectral test multipliers and bridge mass envelopes.
- Exact nine-hit native-decoration and spine census.
- Seven exact polynomial tests of the first and second Riesz identities, with a nonconstant P-dependent covariance and nonlinear base.
- Seven exact polynomial tests of the final common-endpoint path; every coefficient below alpha nine cancels exactly.
- Thirty-six genuinely coupled multidimensional Wick tests in dimensions 2, 3, and 5; maximum floating error 3.70e-13.

The polynomial Riesz fixture is an algebra diagnostic, not a bounded-Hessian original-source witness. Native finite programs are imported under their existing contracts and are not numerically executed by these scripts. The prior source-qualified modulated six-tree variance counterexample remains valid and is not contradicted by the exact weighted OU construction.

The PASS is restricted to the canonical old structural bank with zero physical readout and only the declared old cubic/counterpacket dependence at this join. An additional retained nonlinear observer or a different old/public readout needs a new same-endpoint analysis. The original finite-clock-to-continuum mismatch is carried separately in its actual integrated/conditional scope. Arbitrary-order closure, re-entry with new observers, recursive elimination of the grade-eight offspring, and eventual-sublinear full complexity remain unproved.
