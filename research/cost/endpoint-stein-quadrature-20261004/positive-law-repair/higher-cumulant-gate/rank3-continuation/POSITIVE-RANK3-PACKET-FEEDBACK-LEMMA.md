# Positive rank-three packet: a direct one-energy feedback lemma

2026-10-04. Analytical proof for the actual native three-marked-path reference. Pending independent audit. This supplies the feedback argument left open in Section 6 of the companion native-path note, subject to its explicit finite-source/clock and calibration hypotheses. It does not make a generic non-gradient source native.

## 1. Source-qualified statement

Let Q be a Gaussian coarse-root bank. Let A(Q) be a three-tensor on an n-dimensional ambient space, with every proper matrix cut at most K uniformly in Q and

    ||A||_(L2(Q);HS)<=h.

Assume every proper cut of D_Q A is at most K1 uniformly. The actual coarse-root injection norm is included in K1. No derivative of A is executed; these are heat-coefficient certificates for the literal VALUE sources.

Write M_Q(P)=A(Q)[P], with standard Gaussian public P. An actual finite bounded matrix field M_epsilon(Q,P) may replace this ideal linear field, provided its L2(HS) calibration error is epsilon_abs, its operator bound is Lambda K, and the field is the selected coefficient of the declared native VALUE source. The calibration is in the full Gaussian (Q,P) law and at all fixed moments required by the native field lemma; retain its absolute floor rather than divide by h.

Conditional on Q,P, consider a standard Gaussian pair (Y,R) with cross matrix q M_epsilon(Q,P). Suppose q Lambda K<=1/4. Read out

    X=cY+aR+bP+eta Z,
    c^2+a^2+b^2+eta^2=nu,

where eta>0 is a fixed allocated buffer, and Q,P,Z and the pair innovations are independent in the declared order. The physical fixed readout coefficients may absorb public-log variance-share factors. The public logarithm Lambda includes the analytic clipping threshold T chosen for the assigned tail floor: impose explicitly q max(Lambda,T) K<=1/4. Choose T from the absolute tail tolerance and a known declared upper bound on h, never by dividing by a realized source energy.

Then X is within

    Lambda [q epsilon_abs + q^2 (K+K1) h]             (1)

of the positive quadratic-reference law

    S+q a b c (E A)^sharp :
          [S tensor S/beta^4-I/beta^2]+eta Z,
    S~N(0,beta^2 I), beta^2=c^2+a^2+b^2.              (2)

The sharp symbol only places the output slot according to M(P)^T Y. Physical symmetrization gives third cumulant 6abcq Sym A to leading order. All factors depending on fixed beta,eta,nu are in Lambda. A tiny analytic clipping tail can be added to epsilon_abs and made any prescribed fixed power without changing the executed program.

The comparison retains external physical callers that were fixed before Q and all publics. It marginalizes Q, P and the pair innovations. The reference S is a coupling/reference Gaussian, not an executed retained tape.

## 2. Collapse the native path and remove its root feedback strongly

A three-vertex distinguished native path, conditioned on Q and its sole center public P, is a chain of standard Gaussian channels. Its terminal pair (Y,R) has cross

    q M_epsilon(Q,P), q=rho_root rho_center rho_leaf.

The product is chronological; no independent-source replacement is made inside it. Intermediate outputs have no other observer in this tree, so the path collapse respects the exact native pair retention rule.

The ideal linear field M=A[P] need not be globally bounded. Use it only analytically: clip its singular values at a known logarithmic level T K. Gaussian matrix-series bounds, uniformly in Q, give a negligible L2(HS) clipping tail at T=C sqrt(log(2(n+1)/epsilon)). Its Hilbert energy is proportional to h, not sqrt(n) times h. The proof uses the proper cuts as both matrix-series variance operators and fixed-chaos Hilbert hypercontractivity. Choose qTK<=1/4. The clipped field gives a genuine gapped Gaussian pair.

Conditional Gaussian root coupling compares the actual finite field with this clipped field at cost Cq times their L2(HS) difference. For any such cross matrix,

    R=q M(P)^T Y+(I-q^2 M(P)^T M(P))^(1/2)xi,

with Y,xi standard and independent of Q,P. The gap and scalar root estimate give

    ||(I-q^2M^TM)^(1/2)-I||HS<=Cq^2||M||op||M||HS.

Conditional Gaussian matrix moments and one marked Hilbert factor bound the strong change of the readout when its root factor is replaced by I by Lambda q^2 K h. Restore the unclipped polynomial field at the assigned tiny tail cost. What remains is the actual positive polynomial law

    X_pol=S+q C(Q,P,Y)+eta Z,
    S=cY+a xi+bP,
    C=a A(Q)[P]^T Y.                                  (3)

Thus the square-root sampler feedback is paid, rather than being silently discarded.

## 3. First decoupling: private Gaussian directions of the public readout

Make a known orthogonal Gaussian change of variables from (P,Y,xi) to the visible S/beta and its orthogonal complement U. The complement U is independent of S and Q. The correction C is a degree-two Gaussian polynomial in (S,U). Its conditional mean is exactly

    Cbar(S,Q)=abc A(Q)^sharp:
          [S tensor S/beta^4-I/beta^2].               (4)

The subtracted identity is the conditional Gaussian covariance, not an origin chosen by an oracle.

Set E=C-Cbar; it has mean zero in U conditional on S,Q. Interpolate positive laws

    W_t=S+q Cbar+tq E+eta Z.

Gaussian Riesz in U gives the exact derivative

    d/dt E phi(W_t)
       =t q^2 E[(R_U E)(D_U E)^T:D^2 phi(W_t)].        (5)

The base carrier S and Cbar have no U derivative. Hence there is no unpriced order-q branch. Both polynomial factors in (5) have fixed degree. At each fixed Q, Hilbert hypercontractivity for the marked factor and Gaussian matrix-chaos bounds for the operator factor give

    ||(R_U E)(D_U E)^T||_(L2(S,U);HS)
        <=Lambda K ||A(Q)||HS.

Integrate its square in Q; this uses only the L2(Q) energy h, never an L4(Q) inference from h. One integration in the untouched Z buffer proves the first-decoupling cost Lambda q^2 K h.

## 4. Second decoupling: owned coarse roots

Put A0=E_Q A and

    E_Q(S,Q)=Cbar(S,Q)-E_Q Cbar(S,Q).

Again use positive interpolation, now between

    S+q E_Q Cbar+tq E_Q(S,Q)+eta Z.

The exact private-Q Riesz derivative is a rank-two current

    t q^2 (R_Q E_Q)(D_Q Cbar)^T.                      (6)

Q does not occur in the base Gaussian S. The first factor's coefficient is R_Q(A-A0), whose L2(Q;HS) norm is at most h. The second factor's coefficient has proper cuts at most K1. Conditional on Q, both visible-S factors are degree-two Hermite polynomials. Hilbert hypercontractivity of the first factor and the matrix-chaos estimate of the second bound their product by

    Lambda K1 ||R_Q(A-A0)(Q)||HS.

Integrating in Q yields Lambda K1 h. No derivative of R_Q A is used, and no higher derivative of an original gradient is required: D_Q A is an analytical positive-heat derivative certified by the native source. Gaussian buffer integration proves the second-decoupling cost Lambda q^2 K1 h. Together with Section 2, this proves (1).

This is the precise retained-carrier rule: only a base Gaussian independent of the private block being removed may be held fixed in (5)-(6). An arbitrary later observer of Q or U is not allowed.

## 5. Substitute the native three-force clocks

For the double-Riesz path in the companion note, normalized heat coefficients have

    K<=C,
    K1<=C/sigma,
    h<=C e/(ell sigma).

For h, mark a leaf averaged Jacobian. Gaussian first-Hermite Bessel gives its HS norm e/sigma; divide by ell, keep the other vertices at proper-cut norm C, and eliminate the two source-index tree edges. The total query law of the marked original G remains standard Gaussian, so the centered e is the actual source energy. The derivative bound K1 comes from one further score on a private sigma shield; it is analytical only.

The exact amplitude product is

    q=w ell^3/(3abc sigma).

Consequently (1) gives the substantive feedback price

    Lambda w^2 ell^5 e/sigma^4,                       (7)

with a smaller term having sigma^3 in the denominator. The finite positive dyadic clock sum of w^2/sigma^4 is logarithmic. This is one energy, ell^5 e, rather than a product of two dimension-sized VALUE energies. It is below the requested ell^3 e error at small ell.

The two positive Riesz histories have their own independent complete packet banks and allocated positive shares. The shared coarse roots are retained within each history until (6); they are never split into independent force roots. Add the already stated prior-kernel, finite-filter, finite-clock, and absolute numerical floors with their actual output and inverse-share factors.

## 6. Joining finitely many cubic packets

After the individual comparisons, the sum of independent packets is a positive Gaussian base plus a finite sum of deterministic quadratic corrections, with a fixed untouched Gaussian buffer. Another known orthogonal change of Gaussian coordinates projects all quadratic corrections onto their total visible Gaussian carrier. The private-Gaussian argument in Section 3 then replaces them by their conditional Hermite mean at a one-energy cost bounded by the product of their summed proper-cut scale and summed Hilbert energy. Their Hermite coefficients add exactly; independent complete packet banks are essential.

This proves positive realization of the summed leading third current up to the displayed feedback and absolute floors. It does not use a tensor coefficient as an executed input: the tensor in the final reference is the analytically identified mean of the literal native-source trees.

The last comparison to an originally skew buffered source is a distinct analytical positive-quadratic-reference interpolation. The mean/covariance services must also be calibrated to the requested grade. Actual first/caller/zero/query returns remain those of the finite native VALUE graph and its source-path proof; they do not follow from (1).
