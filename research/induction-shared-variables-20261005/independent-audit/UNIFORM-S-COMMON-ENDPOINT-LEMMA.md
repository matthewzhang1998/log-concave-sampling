# Uniform parameter-current lemma for the positive grouped path

2026-10-05. This supplies the specific extension needed by Section 7 of the reviewed candidate. It is relative to the pinned finite native polynomial comparison, three-frame estimates, Gaussian Riesz bounds, and independent-keep transfer. It is not a claim that ordinary moment matching implies a transport bound.

## 1. Precise finite setting

Freeze the complete retained record E, original clocks and common heats, all scalar readout allocations, finite native versions and numerical gaps. Put t = alpha^(1/2), and keep t fixed when differentiating the interpolation parameter s in [0,1]. After the bounded native VALUE hybrid and finite gapped-root comparison, consider the positive polynomial pushforward

    X_s(t) = sum_j X_{j,s}(t; xi_j) + sqrt(kappa) K,
    X_{j,s}(t; xi_j) = L_j P_j + sum_{r=1}^{R} t^r F_{j,r}(s; P_j,U_j).

Conditional on E, the complete xi_j=(P_j,U_j) blocks are mutually independent standard Gaussian banks; K is independent of every block. L_j P_j is the known source-zero physical row. Each coefficient is a finite Gaussian polynomial. A fixed nonzero mean or known outside shift can be included without alteration. The original analysis may provide arbitrary-degree old-U dependence instead; its admitted Lp Riesz version gives the same argument, but that extension is unnecessary here after the finite conditional polynomial comparison.

The essential quantitative assumptions are the imported one-Hilbert energy and all three proper derivative-frame/operator-moment bounds, for every coefficient and its s derivative, uniformly on [0,1]. The cutoff fixes the finite list of moment exponents. Each nonconstant endpoint derivative has strictly positive t-degree. The keep width kappa>0 is fixed before execution.

These hypotheses hold for the proposed path because old root amplitudes are multiplied by s, correction roots by s^2, and the reference drift by 1-s. Finite covariance-root Taylor coefficients therefore have finite polynomial dependence on s, with their differentiated coefficients uniformly bounded. Alternatively, before Taylor expansion their s derivatives are controlled by the same numerical covariance gap. There is never a division by s, by a root amplitude, or by 1-s. In particular s=0 and s=1 are included. The two native endpoint comparisons can be paid by their absolute floors; differentiating a discretized native sampler or a completed LAW estimate is neither needed nor allowed.

## 2. Local current extraction at the full endpoint

Differentiate the polynomial pushforward:

    d/ds E phi(X_s) = sum_j E[partial_s X_{j,s} dot grad phi(X_s)].

Fix a summand j. Split any coefficient current into its U_j mean and centered part. Riesz in U_j transfers the centered part only onto the endpoint. Since all other summands have disjoint dynamic banks,

    D_(U_j) X_s = D_(U_j) X_{j,s},
    D_(P_j) X_s = L_j + D_(P_j)(X_{j,s}-L_j P_j).

These identities are the ownership-sensitive step. They remain exact with the other summands as arbitrary independent frozen additive shifts. They would fail for an outside caller that reads an internal j-root.

There is no base identity in U_j. Every U_j branch therefore raises t-degree. Split the U_j mean into its constant and centered finite P_j polynomial; apply P_j Riesz to the centered part. A branch using L_j lowers its remaining Hermite degree; a branch using the nonlinear derivative raises its t-degree. No intermediate Riesz coefficient is differentiated.

Stop at t-degree 20. Between positive-degree steps there are only finitely many degree-lowering P_j branches. The number of positive-degree steps is finite, and the rank of every boundary current is bounded by the finite coefficient census. Riesz contraction, conditional Jensen and the preselected Holder exponents bound each coefficient by its one-Hilbert energy times endpoint operator moments. The full Gaussian-bank Jacobian estimates, including the HS-input derivative flattening, are what prevent an extra dimension factor. This is precisely the mechanism in LOW30, lines 6326–6387, with s differentiated at fixed t instead of differentiating the amplitude itself.

The result is

    d/ds E phi(X_s)
      = sum_j sum_(r<20) t^r A_{j,r}(s,E):E[D^(m_{j,r}) phi(X_s)]
        + boundary_s(phi),

where every displayed A is constant in the dynamic Gaussian banks, and the entire boundary is evaluated at that same X_s. Its keep-transferred continuity velocity has norm at most C(E) t^20 sqrt(D), uniformly in s. Nothing in this argument replaces X_s by the source-zero endpoint.

## 3. Identifying the constant operators correctly

For each block let Phi_j(s,t,theta) be its characteristic function. This notation is used only coefficientwise in the finite t jet; no exponential moment of a polynomial is assumed. Write

    Phi_j = Phi_(j,0) [1 + F_j(s,t,theta)],

where Phi_(j,0) is its nonzero Gaussian characteristic function and F_j has strictly positive t-degree. All finite t coefficients are Gaussian moments of polynomials, hence polynomial in theta times that Gaussian baseline. The truncated formal inverse of 1+F_j exists uniquely.

Apply the local current identity to Fourier tests before adding the independent other blocks. Comparing successive t grades gives

    A_j(s,t,theta) = [partial_s Phi_j / Phi_j]_(degree<20)
                   = [partial_s log Phi_j]_(degree<20).

The first equality is finite triangular division at the SAME endpoint. In particular the constant operator is not simply a coefficient of partial_s Phi_j evaluated against the Gaussian baseline. The latter shortcut omits the recursive endpoint terms.

Since the constant A_j are independent of every dynamic block, their operators add at the full common X_s. Therefore an s-independent combined logarithmic jet through degree nineteen makes their sum zero. The nonzero, s-independent target log coefficient is harmless: its s derivative is zero, and the triangular division has already accounted for it. This is the justified positive-target extension of the pinned Gaussian-target proof.

## 4. Application to the reviewed path

Include the old packet, four correction packets, and independent reference H row as the disjoint blocks. Write C8 = T00+T10+T01+T11, with each old coefficient and normalization retained. Their finite log jets are

    K_old(s) = s C6 + s^2 C8 + O_current(alpha^10),
    K_corrections(s) = -s^2 C8 + O_current(alpha^12),
    K_H(s) = (1-s) C6 + O_current(alpha^12).

The notation O_current here is justified by the finite extraction above and its coefficient/frame bounds; it is not an unproved estimate inferred from the printed scalar symbol. Original side-shift corrections first enter the old second cumulant at alpha^10; three original root-spine clusters start no earlier than alpha^12. Each decorated counterpacket's first nonmain cluster starts at alpha^12. The H drift has size alpha^6, so every nonlinear self-contraction starts at alpha^12.

At every s, including both endpoints,

    partial_s (K_old+K_corrections+K_H)
      = (C6+2s C8) - 2s C8 - C6 + O_current(alpha^10).

Thus all constant current operators below alpha^10 cancel at the actual X_s. A first endpoint defect of order alpha is not discarded. It appears in the finite local recursion and is either included in a lower-grade constant coefficient that cancels, or propagated until the degree-ten boundary. Stopping after one physical integration, as in the earlier conservative alpha^9 join, would indeed be insufficient for this stronger result.

At s=0 the old and correction root couplings vanish. Their ideal packet readouts have their exact allocated Gaussian laws, while H carries the positive Hermite drift. At s=1 the reference drift vanishes and all corrected packets have their intended amplitudes. The endpoints are therefore the advertised positive comparator and corrected group, up to their separately paid native/value/calibration/clock floors. Positivity of the pushforward does not require the Hermite map to be injective or its Jacobian to remain positive. The native covariance roots, source guards and untouched keep still must satisfy their numerical gaps.

Transfer the remaining currents through K only. Boundary coefficients are independent of K, so its fixed-rank Gaussian score factors apply with the actual kappa inverse powers. Conditional Jensen gives a continuity velocity of norm at most C(E) alpha^10 sqrt(D). Integrating s over [0,1] gives

    W2(Law(X_1 | E), Law(X_0 | E))
      <= C(E) alpha^10 sqrt(D) + propagated absolute floors.

For a random E, integrate C(E)^2 and the actual floor bounds; a uniform constant is claimed only when the pinned source estimates are uniform. This proof retains T6(B). It supplies no B-mean replacement, old-bank averaging cancellation, RAW admission, or derivative of a finished LAW certificate.
