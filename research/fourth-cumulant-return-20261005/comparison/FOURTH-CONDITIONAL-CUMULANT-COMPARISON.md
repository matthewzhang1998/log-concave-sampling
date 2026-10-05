# Dimension-safe fourth conditional cumulant comparison

2026-10-05. Analytical comparison under C1/Gaussian assumptions. No tensor-valued source oracle, Hessian of the force, or native producer admission is used.

## Result

Condition on the actual retained variable Y. On the same complete conditional standard Gaussian bank G, let U and V be centered C1 maps into R^D, with

    ||D_G U||op <= L,    ||D_G V||op <= L.

Put R=U-V and e_R(Y)=(E[|R|² | Y])^(1/2). Then pointwise in Y,

    ||κ4(U | Y)-κ4(V | Y)||HS <= 24 L³ e_R(Y).          (1)

Consequently, if L<=C A and ||e_R||_(L² Y)<=C A²√D,

    ||κ4(F2 | Y)-κ4(H | Y)||_(L² Y;HS)
        <= C A⁵√D.                                     (2)

The cumulants in (2) are unchanged by conditional centering. For the delayed displacement -qF2, the fourth cumulant mismatch acquires the factor q⁴, with positive sign. In the coefficient convention κ4/4!, the pointwise bound is simply L³e_R; in the centered Stein-current convention c3=κ4/3!, it is 4L³e_R.

The same assertion holds for uniformly Lipschitz strong finite-bank limits. It compares analytical targets only and does not provide a finite original-VALUE fourth-order packet.

## 1. Centered mixed quadratic covariance

In this section all expectations are conditional on the fixed Y, whose notation is suppressed. Let P and Q be any centered maps of the same Gaussian bank, each with complete Lipschitz constant at most L. They need not be independent. Gaussian Poincaré first gives

    Cov(P)<=L²I,    Cov(Q)<=L²I.

For an arbitrary matrix M,

    D_G(P^T M Q)=(DP)^T M Q+(DQ)^T M^T P,

so

    Var(P^T M Q)
      <=2L² E[|M Q|²+|M^T P|²]
      <=4L⁴ ||M||HS².                                  (3)

Thus the covariance operator of P⊗Q, acting on the Hilbert space of matrices with its HS norm, is bounded by 4L⁴I. The proof uses only first derivatives. No independence or symmetry of the mixed covariance is assumed.

## 2. The mixed Wick-cubic covariance bound

Let A, B and C be centered complete-bank maps, each L-Lipschitz. Define their covariance-subtracted mixed cubic tensor by its actual slot indices:

    W(A,B,C)_(ijk)
      =A_i B_j C_k
       -E[A_i B_j] C_k
       -E[A_i C_k] B_j
       -E[B_j C_k] A_i.                                 (4)

Its expectation is generally nonzero: E W=E[A⊗B⊗C]. In all covariance estimates below it is implicitly centered, and its mean is never assumed to vanish.

For an arbitrary, not necessarily symmetric, tensor T in (R^D)^⊗3, set

    (P_A)_i =sum_(j,k) T_(ijk)(B_j C_k-E[B_j C_k]),
    (P_B)_j =sum_(i,k) T_(ijk)(A_i C_k-E[A_i C_k]),
    (P_C)_k =sum_(i,j) T_(ijk)(A_i B_j-E[A_i B_j]).

The exact first-derivative identity is

    D_G <T,W> = (DA)^T P_A+(DB)^T P_B+(DC)^T P_C.        (5)

This is precisely where all three covariance-linear subtractions in (4) are needed: each differentiated cubic slot is paired with a centered quadratic expression.

Apply (3) to each matrix slice of T and sum the squared components. It gives

    E|P_A|², E|P_B|², E|P_C|² <=4L⁴||T||HS².

A second Gaussian Poincaré estimate, followed by the three-term square-sum inequality, gives

    Var(<T,W>)
      <=E|D_G<T,W>|²
      <=3L² E[|P_A|²+|P_B|²+|P_C|²]
      <=36L⁶ ||T||HS².                                 (6)

Equivalently,

    ||Cov(W(A,B,C))||op <=36L⁶.                         (7)

This is dimension-free and valid for arbitrary mixed choices of A, B and C. The subtracted cubic tensor is an analytical device, not a stipulated source port.

## 3. One arbitrary error slot

For a centered square-integrable vector R, on the same probability space, no derivative bound is needed in the following step. For centered Hilbert-valued R and Q, the cross-covariance inequality is

    ||E[R⊗Q]||HS² <= ||Cov(Q)||op E|R|².                (8)

For completeness, choose an orthonormal basis (e_i) of the R-space. Cauchy–Schwarz against arbitrary unit vectors in the Q-space shows

    ||E[<R,e_i>Q]||²
       <=||Cov(Q)||op E[<R,e_i>²].

Sum over i to obtain (8). No dimension of the tensor space appears.

The centered mixed fourth cumulant, retaining the displayed physical slot order, is exactly

    κ(R,A,B,C)=E[R⊗W(A,B,C)].                            (9)

Indeed the three terms subtracted in (4) become the three pairings of the four slots after multiplication by R and expectation. Since E R=0, (9) is unchanged if W is replaced by W-EW. Equations (7)-(9) yield

    ||κ(R,A,B,C)||HS <=6 L³ ||R||₂.                    (10)

Thus the sole dimension-sized Hilbert energy can be placed on the error R, with three powers of the common Lipschitz constant. There is no product of two physical energies.

## 4. Exact telescoping and conditional integration

Fourth cumulants are separately multilinear. With R=U-V,

    κ4(U)-κ4(V)
      =κ(R,U,U,U)+κ(V,R,U,U)
       +κ(V,V,R,U)+κ(V,V,V,R).                          (11)

For each summand, permute the physical tensor slots to put R first; this preserves the HS norm. The remaining three vectors are U or V, and each is centered and L-Lipschitz. Applying (10) four times proves (1).

Restore conditioning on the same actual Y. The proof is pointwise in Y, so Minkowski is not needed beyond the scalar L² bound:

    ||κ4(U|Y)-κ4(V|Y)||_(L² Y;HS)
       <=24 L³ ||e_R||_(L² Y).                         (12)

If the common Lipschitz bound depends on Y, the exact extension is

    ||κ4(U|Y)-κ4(V|Y)||_(L² Y;HS)
       <=24 ||L(Y)³e_R(Y)||_(L² Y).

No replacement of Y's law and no independently resampled future is permitted or used.

## 5. C1 regularity and true analytical histories

The identities above involve only first derivatives of the centered vector maps. Cubic products belong to Gaussian W^(1,2) because centered Gaussian-Lipschitz maps have moments of every finite order, and their first derivatives are bounded. Gaussian Poincaré therefore applies directly by its Sobolev form, or by smooth truncation and approximation. No second or third derivative of U or V is required.

For the true coherent analytical future histories, apply (1) to their common uniformly Lipschitz finite-bank approximants. For fixed physical dimension D, their centered Gaussian-Lipschitz moments of every fixed order are uniformly bounded. Strong L² convergence therefore upgrades to strong L4 convergence by interpolation against a uniform higher-moment bound. Every fourth moment and covariance pairing converges; hence the fourth cumulant tensors converge in finite-dimensional HS norm. The error energies converge as well. This passes (1) to the limit. Equivalently, one can apply the same Sobolev/Poincaré argument directly on the underlying Gaussian Hilbert space whenever that representation has already been established.

When the approximation also depends on Y, extract a subsequence converging in conditional L² for almost every Y. The conditional higher-moment bounds are uniform in Y because the centered maps have the same deterministic Lipschitz bound L. The preceding interpolation then gives pointwise convergence of conditional cumulants. Their fixed-D uniform moment bounds allow dominated convergence in L²(Y;HS). Thus (12) also passes to true conditional histories, without requiring pointwise convergence of the entire original sequence.

Conditional centering is an L² contraction. In particular the coherent error estimate from Section 3 of the shrinking-buffer join gives

    ||e_R||_(L² Y)
       =||(F2-H)-E[F2-H|Y]||₂
       <=||F2-H||₂
       <=C A²√D.

Together with the imported complete conditional-bank Lipschitz bound L<=C A, this proves (2).

## 6. Relation to the all-rank centered Stein recurrence

The imported recurrence identifies c3=κ4/6 and bounds individual fourth cumulants by one physical energy. Taking a naive difference of those recurrences introduces derivatives of R, which need not be O(e_R) under only C1/Lipschitz assumptions. The proof here avoids that issue by moving all derivatives onto the other three slots in the covariance-subtracted cubic test. It obtains the desired stability in the error's L² energy without assuming small derivative error.

This resolves the fourth conditional cumulant target-comparison gate. It does not by itself resolve a positive fourth-order law consumer, a native four-force packet, proper-cut bounds for such a packet, or its root/caller/variance-budget feedbacks.

## Source context

- /workspace/shared/shrinking-buffer-skew-join-20261005/POSITIVE-SKEW-SHRINKING-BUFFER-JOIN.md, Section 3: coherent centered U,V, complete conditional-bank Lipschitz bound, and error energy.
- /workspace/shared/recovery-20261004/cost/endpoint-stein-quadrature-20261004/positive-law-repair/higher-cumulant-gate/order-reentry/ALL-RANK-CENTERED-STEIN-CURRENT-RECURRENCE.md: analytical cumulant/current normalization.

No imported source file is modified.
