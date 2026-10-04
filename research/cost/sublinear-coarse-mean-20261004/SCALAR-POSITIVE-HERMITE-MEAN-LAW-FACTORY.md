# A bounded scalar positive Gaussian-mean law factory

2026-10-04. An independently derived scalar prototype, developed after the strong-mean barrier. No external/private manuscript is used. This is not a multidimensional posterior source theorem, a smooth first-action port, or a complete hidden-decoder replacement.

## Result first

Suppose a complete fresh scalar record W has an executable known anchor b and a supplied almost-sure bound

    |W-b| <= sigma epsilon,

where sigma>0 is known. Neither E W nor E(W-b) is known or queried. For each fixed integer K>=1, the construction below uses K independent complete records and an explicitly positive conditional density. It approximates

    N(E W, sigma^2)

to arbitrary fixed order in epsilon as K increases, with no inverse-epsilon replication count. Its literal original-query work is K times the complete W-provider cost. Sampling from the positive density uses only known scalar arithmetic and Gaussian/uniform proposals after those K records have been produced.

For cutoff R>=1 and epsilon R+epsilon^2/2<=log(5/4), define

    T_K(epsilon,R) = C_K [epsilon^(K+1)
                  + epsilon(1+R^K) exp(-R^2/4)].

The uncapped sampler has

    W2(output, N(E W,sigma^2)) <= C sigma sqrt(T_K).

With a cap of B rejection attempts and a fresh Gaussian fallback, add at most sqrt(5) sigma 3^(-B/2). Therefore, if epsilon<=C A^a times fixed public logarithms with a>0, and sigma is bounded by fixed public logarithms, choose K with a(K+1)>2q, R^2=C_q log(1/A), and B=C_q log(1/A) to get W2=O(A^q) at any fixed q. For sigma<=C A^(-s) times such logarithms, replace q by q+s in these choices. Strict inequalities absorb the displayed logarithms. The output is a finite positive program.

The square-root loss is a deliberately elementary bound; no high-order transport theorem is needed. The prototype's obstacle to reuse is its scalar bounded-source/ownership/first-action scope, rather than positivity or an unknown normalization.

## 1. Literal finite program

Set X_i=(W_i-b)/sigma, where W_1,...,W_K are independent complete records conditional on the exposed caller. Thus |X_i|<=epsilon. Put

    P_k = product_(i=1)^k X_i,
    H_0(z)=1, H_1(z)=z,
    H_(k+1)(z)=z H_k(z)-k H_(k-1)(z).

The H_k are the probabilists' Hermite polynomials. Let phi be the standard Gaussian density and define explicit known integrals

    I_k(R) = integral_(-R)^R H_k(z) phi(z) dz
           = phi(R)[H_(k-1)(-R)-H_(k-1)(R)],  k>=1.

Form the finite random polynomial and its explicit integral

    S_X(z) = sum_(k=1)^K P_k H_k(z)/k!,
    C_X = sum_(k=1)^K P_k I_k(R)/k!.

Now use the conditional density

    q_X(z) = phi(z) r_X(z),
    r_X(z) = 1 + 1_(|z|<=R) S_X(z) - C_X.

The normalization is EXACTLY one for every realized batch: the cutoff polynomial contributes C_X and the constant correction subtracts the same C_X. No expectation oracle and no unknown conditional normalizer occur. b is the known anchor, not an analytically centered unknown mean.

Draw Z~N(0,1) and U~Uniform[0,1] independently of the complete batch. Accept Z when U<=(2/3)r_X(Z). Upon acceptance output b+sigma Z. Repeat proposals if desired, without drawing new W_i or evaluating their provider again. A finite implementation makes B attempts and, only if all fail, outputs b+sigma Z_fallback using one fresh standard Gaussian.

## 2. Uniform positivity and actual sampling work

The explicit Hermite coefficient formula gives

    sum_(k>=1) epsilon^k |H_k(z)|/k!
       <= exp(epsilon |z|+epsilon^2/2)-1.

Hence on |z|<=R the imposed guard gives |S_X(z)|<=1/4, uniformly over all possible batches. Also

    |C_X| <= integral_(|z|<=R) |S_X(z)| phi(z) dz <=1/4.

It follows globally that 1/2<=r_X<=3/2, so q_X is a positive normalized density. The rejection envelope (3/2)phi is valid, and the acceptance probability is exactly 2/3 for every batch. The expected proposal count is 3/2. With B attempts the exact failure probability is 3^(-B), independent of the batch.

The original-query count is exactly K complete W records, plus the known anchor's producer if it has one; that anchor's full bill must be included. The proposal phase requires O(KB) scalar arithmetic and no original gradients. The K records must be complete independent executions conditional on the caller. Sharing internal random ancestors would invalidate E P_k=(E X)^k in general.

Finite arithmetic for phi(R), I_k, coefficients, and the acceptance threshold has a separate selected numerical tolerance. At fixed K,R,B a total variation error bounded by the total acceptance-probability error can be obtained by common uniforms; its moment/profile cost must be restored before using this as any numerical source port. This note does not grant exact transcendental arithmetic for free in a full implementation.

## 3. Exact marginal polynomial and target error

Write mu=E X, with |mu|<=epsilon. Independence yields

    E P_k=mu^k.

Thus the unconditional density is

    qbar(z)=phi(z)[1+1_(|z|<=R) S_mu(z)-C_mu],
    S_mu(z)=sum_(k=1)^K mu^k H_k(z)/k!,
    C_mu=integral_(|z|<=R) S_mu(z) phi(z) dz.

mu is used only in the proof. The desired shifted Gaussian has ratio

    phi(z-mu)/phi(z)=exp(mu z-mu^2/2)
                          =sum_(k>=0) mu^k H_k(z)/k!.

Hermite orthogonality gives the exact full-polynomial remainder norm

    || exp(mu Z-mu^2/2)-sum_(k=0)^K mu^k H_k(Z)/k! ||_2^2
       =sum_(k>K) mu^(2k)/k!
       <=exp(epsilon^2) epsilon^(2K+2)/(K+1)!.

Also E S_mu(Z)=0, so |C_mu|<=||1_(|Z|>R) S_mu(Z)||_2. Standard integration of a fixed Gaussian polynomial tail, or repeated integration by parts, gives for epsilon<=1 and R>=1

    ||1_(|Z|>R) S_mu(Z)||_2
       <= C_K epsilon (1+R^K) exp(-R^2/4).

Combining the polynomial remainder, its removed tail, and C_mu,

    ||(qbar-phi_mu)/phi||_(L2(phi)) <= T_K(epsilon,R).

For completeness, the elementary W2 conversion uses the maximal common-density coupling. Match min(qbar,phi_mu) identically and couple the residuals independently. Then

    W2(qbar,phi_mu)^2
       <=2 integral z^2 |qbar(z)-phi_mu(z)| dz
       <=2 sqrt(3) ||(qbar-phi_mu)/phi||_(L2(phi)).

This proves the announced O(sqrt(T_K)) bound. All moments exist, since qbar/phi is between 1/2 and 3/2 and the target is a finite Gaussian shift.

## 4. Finite rejection cap

Because the failure probability delta=3^(-B) is constant across every batch, the capped marginal law is

    (1-delta) qbar + delta phi.

Couple it identically to qbar on the successful branch. On the failure branch couple independent qbar and phi variables. Since E_qbar Z^2<=3/2 and E_phi Z^2=1, the squared gap is at most 5 delta. This proves the additional sqrt(5)3^(-B/2) term, then scaling by sigma gives the physical bound.

No unbounded random stopping time is hidden in the finite version. The proposal count is at most B, whose cap is selected in advance from the target, heat, and public precision parameters. One may execute unused dummy proposals if a fixed-length random tape is desired.

## 5. Conditioning, anchors, and non-claims

The lemma works conditional on any exposed retained caller theta if b(theta), sigma(theta), epsilon and the fresh complete-source bound are supplied there. Its target is then the corresponding conditional Gaussian mean. A later deterministic host that uses theta can be included in the conditional law comparison, subject to its actual Lipschitz/error propagation. This does not permit appending private sample records as retained observers after averaging them out: qbar integrates the K fresh W_i records.

The construction does not provide the needed actual first/adjoint source port. Indicator cutoff, accept/reject decisions and discrete stopping are not a smooth Gaussian pushforward. It supplies no actual two-caller coupling with the baseline's declared first bounds. It also does not supply an arbitrary vector analogue: replacing scalar products by tensor Hermites can introduce norm/positivity conditions depending on dimension. Applying it coordinate by coordinate with fresh independent banks costs D complete providers unless a new joint construction is proved.

For an unbounded provider, any clipping/truncation must be executable around the KNOWN anchor and its change in E W must be bounded separately. The stated bounded-source assumption cannot be replaced by imaginary centering at E W. Gaussian concentration or high moments of a particular complete source might supply such a clipping lemma, but that full source-specific proof is not included here.

Finally, approximating N(E W,sigma^2) does not make E W available as a strong statistic. Estimating that mean from repeated factory outputs would have the Gaussian variance sigma^2 and needs a new inverse-accuracy count. Thus the construction is consistent with STRONG-POSTERIOR-MEAN-ORACLE-BARRIER.md and illustrates exactly the weaker type of contract a successful route must exploit.

## 6. Complete recurrence if a compatible consumer is supplied

For this isolated scalar factory,

    Work_factory(Q,A) <= K_Q Work_complete_W(A)
                            +Work_anchor(A)+O_Q(B_Q K_Q),

with K_Q fixed after the requested law order and B_Q=O_Q(log(1/A)). Hence its new original-query heat exponent is zero beyond the provider/anchor exponent. It does not erase the provider's existing cost, fix its posterior bias, or prove the complete hidden-source recurrence. Those are separate gates.

This is a bounded constructive escape from the strong-mean requirement, not an eventual-sublinear theorem.
