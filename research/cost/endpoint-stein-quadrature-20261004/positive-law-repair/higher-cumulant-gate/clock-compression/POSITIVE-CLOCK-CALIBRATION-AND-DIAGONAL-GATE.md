# Positive OU-clock calibration, the diagonal gate, and a protected-root gradient lift

2026-10-04. Bounded constructive result. This does not rule out positive law compilers with covariance compensation, and it does not restart ordinary unsmoothed path quadrature.

## Results first

1. A finite PSD Gaussian draw can preserve every required endpoint/one-inner-clock Markov triple and exactly calibrate the full matrix quadratic covariance through degree two. It uses positive clocks and original gradient VALUES only.
2. Even after this favorable calibration, the literal nested finite-path predictor fails a uniform third-order law bound with polylogarithmically many outer nodes. The obstruction is an explicit, smooth, uniformly Hessian-sandwiched cosine fixture. Its second-order variance coefficient debt is at least `1/(256 N^2)`, where `N` is the number of positive outer clock atoms. This is a statement about this direct predictor class, not all positive pushforwards.
3. The protected endpoint is not itself an inverse-heat obstruction: the complete finite OU packet conditional on its endpoint is an exact fixed readout of a genuine full gradient in the remaining Gaussian roots. The readout norm is only the square root of a public logarithm, while the gradient source has Lipschitz constant at most `A` and one-energy scale `A sqrt(D)`.
4. Ordered-ratio quadrature can approximate individual two/three-node Markov operators without the product-grid diagonal. Realizing those crossproducts in a positive law remains a separate covariance-compensation gate. A companion crossbank calculation makes this distinction literal.

## 1. Setup and the finite positive draw

Let `g=grad U`, `g(0)=0`, and `0<=Dg<=A I`, with only a C2 potential assumed. The continuous Gaussian path has

    X_1=Z, Cov(X_r,X_s)=min(r,s)/max(r,s) I.

Fix positive outer weights `w_i`, nodes `0<r_i<1`, and positive inner weights `v_j`, nodes `0<tau_j<1`, satisfying

    sum w_i=sum v_j=1, sum w_i r_i=sum v_j tau_j=1/2.       (1)

All finite versions are fixed before caller differentiation. Write

    R^OU_ij=min(r_i,r_j)/max(r_i,r_j),
    R^0=rr^T+diag(1-r_i^2),
    R^max=rr^T+cc^T, c_i=sqrt(1-r_i^2).

Each displayed matrix is the covariance of standard Gaussian outer nodes conditional on a common endpoint with `Cov(Z,X_i)=r_i I`. Each conditional matrix `R-rr^T` is positive semidefinite.

Their outer integrated variances are

    S=w^T R^OU w,
    S0=1/4+sum w_i^2(1-r_i^2),
    Smax=1/4+(sum w_i c_i)^2.                            (2)

Require the public known-arithmetic guard `S0<1/2<Smax`. All diagnostic fine rules below satisfy it. If `S>=1/2`, set `Raux=R0`; otherwise set `Raux=Rmax`. Define

    eta=(1/2-S)/(w^T Raux w-S),
    Rcal=(1-eta)R^OU+eta Raux.                          (3)

Then `0<=eta<=1`, `diag Rcal=1`, `Rcal-rr^T>=0`, and

    w^T Rcal w=1/2                                     (4)

exactly in real arithmetic. All quantities are known Gaussian/clock coefficients, not queries about g. Floating point/finite encodings receive their absolute precision budget.

An actual root construction is

    X_i=r_i Z+sqrt(1-eta) L_i^OU G+sqrt(eta) L_i^aux G',

where `L^OU(L^OU)^T=R^OU-rr^T`, and similarly for the auxiliary matrix. The OU factor can be generated sequentially with Brownian increments as in Section 5; the independent auxiliary factor is diagonal and the maximal one is one common root. No unspecified covariance oracle is used.

For each inner node draw an independent standard `E_ij` and put

    X_ij=tau_j X_i+sqrt(1-tau_j^2) E_ij.                 (5)

Every individual triple `(Z,X_i,X_ij)` has exactly the same Gaussian law as `(X_1,X_r,X_r tau)` at its named clocks. The joint outer law has been deliberately changed by (3); it is not called the original path law.

The positive VALUE predictor is

    K_i=sum_j v_j g(X_ij),
    Y=Z-sum_i w_i g(X_i-K_i).                           (6)

It uses `N n+N` original g-VALUES, plus the captured original-mode/anchor ledger from its host. Its root count is at most `D(1+2N+Nn)`; the independent maximal branch needs fewer roots. All first/adjoint sweeps use HVPs only at the recorded VALUE sites, and discarded primal records require complete replay. No HVP is a producer leaf.

Each outer and inner Gaussian row has Euclidean norm one. Positivity and total mass one give

    Lip K_i<=A,
    Lip(Y-Z)<=A(1+A),
    ||Y-Z||_p<=C_p A sqrt(D),
    Y(0)=0.                                            (7)

These are actual transcript bounds, not derivatives of a law error. Captured caller, nonzero origin, finite-mode tilt and numerical floors must be restored exactly as in the host packet. There is no dimension guard or inverse-energy division.

## 2. Exact quadratic calibration

For `g(x)=Bx`, `0<=B<=A I`, let

    H0=sum_i w_i X_i, L0=sum_ij w_i v_j X_ij.

Then (6) executes exactly

    Y=Z-B H0+B^2 L0.

Its Gaussian covariances obey

    Cov(Z,H0)=1/2 I,
    Cov(H0)=1/2 I,
    Cov(Z,L0)=1/4 I.

Consequently

    Cov(Y)=I-B+B^2-2 q3 B^3+q4 B^4,                   (8)

where the scalar, known `q3=Cov(H0,L0)` and `q4=Var(L0)` are bounded. The full matrix coefficient of `B^2` is exactly one. This retains every shared-root term and yields the usual fixed-gap `O(A^3 sqrt(D))` quadratic W2 comparison. It says nothing yet about nonlinear forces.

An alternative is to retain `R^OU` and multiply the inner shift by `2(1-S)`. That also fixes the quadratic coefficient, but changes a nonlinear rank-one current. We use the more favorable (3)-(6), which preserves each local Markov triple and its nested coefficient exactly.

## 3. The finite-product diagonal cannot have uniform Hermite accuracy

For a scalar Gaussian-chaos degree `ell`, the continuous path pair clock has moment

    integral_0^1 integral_0^1 [min(r,s)/max(r,s)]^ell dr ds
        =1/(ell+1).                                    (9)

For a positive N-atom packet its analogous covariance multiplier is

    sum_ij w_i w_j R_ij^ell.

For the OU matrix with distinct clocks this tends to `sum w_i^2` as `ell->infinity`. In particular there is no uniform-Hermite `delta` approximation unless `N>=1/delta`. Calibration at degree one cannot remove this diagonal atom.

The next section prices the same issue under the actual Hessian-sandwich class, rather than relying on arbitrary high-energy Hermite functions that are not admissible gradients.

## 4. Explicit smooth C2 fixture and the actual nested law

Set

    s2=sum_i w_i^2, u=4/s2, k=sqrt(u), a=1/2, b=1/4,
    f_k(x)=a x+(b/k)(cos(kx)-1),
    U_A,k(x)=A[a x^2/2+(b/k^2)sin(kx)-(b/k)x],
    g=A f_k.                                           (10)

This is anchored and smooth, and uniformly in k

    A/4 <= U_A,k''(x)=A[a-b sin(kx)] <=3A/4.             (11)

Thus the test is in the prescribed C2 class, even when k depends on the finite clock rule. No extra smoothness norm is included in the class.

### Exact rank-two clock function

For two standard normal variables with arbitrary correlation rho, the centered cosine covariance is

    C_u(rho)=e^-u[cosh(u rho)-1]
       =1/2[e^-u(1-rho)+e^-u(1+rho)-2e^-u].             (12)

It is nonnegative for every real rho. Let

    K_Q=sum_ij w_i w_j C_u(Rcal_ij),
    I_u=integral_0^1 C_u(tau) d tau
       =(1-e^-2u)/(2u)-e^-u.                            (13)

Every PSD covariance with standard marginals, including (3), therefore obeys

    K_Q>=s2 C_u(1)=s2(1-e^-u)^2/2,
    I_u<=1/(2u)=s2/8.

Since `u>=4`, `(1-e^-u)^2/2-1/8>1/3`, so

    (K_Q-I_u)/u >=s2^2/12.                             (14)

The affine/cosine cross term is zero by Gaussian parity. Exact quadratic calibration has removed the affine variance discrepancy, leaving (14).

### The nested rank-one current is separately polylog accurate

For a Markov triple `Cov(Z,X)=r`, `Cov(X,Y)=tau`, `Cov(Z,Y)=r tau`, Gaussian integration by parts gives

    E[Z f_k'(X) f_k(Y)]
       =r[a^2 tau+b^2 F_u(tau)],
    F_u(tau)=e^-u/2-e^-u cosh(u tau)+tau e^-u sinh(u tau).
                                                               (15)

This equality is analysis of the executed VALUE predictor; it is not an HVP-valued instruction. Its integral is

    Fbar_u=e^-u/2+e^-2u/u-(1-e^-2u)/(2u^2).              (16)

Let the inner positive rule have the previously proved uniform moment bound

    sup_(ell>=0) |sum_j v_j tau_j^ell-1/(ell+1)|<=delta.

The sum of the absolute Taylor coefficients of the two nonconstant terms in (15) is at most

    e^-u(cosh u+sinh u)=1.

Therefore, uniformly in u,

    |sum_j v_j F_u(tau_j)-Fbar_u|<=delta.                (17)

This uses the admitted dyadic polylog rule and no derivative of f. Choose its accuracy `delta<=s2^2/48` (or a stricter heat-order budget). It costs only `O(log^2(1/delta))` inner nodes.

### Exact second-order variance coefficient of the full predictor

Write the executed expansion at fixed finite rule as

    Y=Z-A h+A^2 j+R_A,
    h=sum_i w_i f_k(X_i),
    j=sum_ij w_i v_j f_k'(X_i) f_k(X_ij).

Because `|f_k'|<=3/4`, `|f_k''|<=k/4`, positivity gives

    ||R_A||_p<=C_p A^3 k.                               (18)

This is a bound specific to the smooth test family, used only for the separator. It is not imposed as a uniform C2 theorem on arbitrary potentials.

The true target variance is

    Var_mu(X)=1-a A+A^2{a^2+b^2 T_u}+O(A^3),
    T_u=e^-u/2+(e^-2u-e^-u)/u.                          (19)

The `O(A^3)` constant can be made uniform for `k>=2`: the potential shape in (10) is bounded by a fixed quadratic-plus-linear envelope. Formula (19) also follows directly by expanding the Gaussian density, including the squared first mean; forgetting that centering gives a wrong coefficient.

Using (4), (13) and (15), the difference in A-squared variance coefficients between (6) and (19) is exactly

    D_Q=b^2[(K_Q-I_u)/u+
               sum_j v_j F_u(tau_j)-Fbar_u].            (20)

Equations (14) and (17) yield

    D_Q>=b^2 s2^2/16=s2^2/256>=1/(256 N^2).             (21)

This is a nonlinear defect of the actual quadratic-calibrated nested predictor, not merely a standalone H covariance diagnostic. Its local two/three-node rank-one current has already been granted exponentially accurate inner quadrature.

All scalar second moments remain bounded. The elementary inequality

    W2(mu,nu)>=|sd(mu)-sd(nu)|

therefore implies, after accounting for (18),

    W2(Law(Y),mu_A,k)>=c A^2 s2^2-C A^3/sqrt(s2).       (22)

For any outer size `N(A)` bounded by a fixed power of public logarithms, `A N(A)^(5/2)->0`. Along the adversarial but legal family (10), (22) is asymptotically at least `c A^2/N(A)^2`. In particular it is not `O(A^3)` uniformly over the Hessian-sandwich class. This rules out the desired polylog uniform third-order law grade for this direct predictor.

Equation (21) is a coefficient-accuracy lower bound `N>=c delta^-1/2`; it is NOT promoted to an unrestricted finite-A query lower bound for arbitrary power-growing N, because then the explicit remainder in (22) must also be compared. No all-law-compression impossibility is claimed.

## 5. Constructive protected-endpoint gradient lift

The diagonal separator does not mean the endpoint has to be integrated away to access a gradient source. First use the genuine OU factor, ordered `r_1>...>r_N`. Set

    q_i=(r_i^-2-1)/2, q_0=0, Delta q_i=q_i-q_(i-1),
    X_i=r_i Z+sqrt(2)r_i sum_(j<=i) sqrt(Delta q_j) G_j,
    L_ij=sqrt(2)r_i sqrt(Delta q_j) 1_(j<=i).

This is a literal PSD Gaussian draw, since

    (LL^T)_ij=min(r_i,r_j)/max(r_i,r_j)-r_i r_j.

Define

    c_i=sqrt(1-r_i^2),
    t_j=(sqrt(q_j)-sqrt(q_(j-1)))/sqrt(Delta q_j).        (23)

Then the exact telescoping identity is

    L_i t=c_i.

The endpoint Z is a caller. Define the analytical potential and its executed full gradient

    Psi_Z(G)=sum_i (w_i/c_i) U(r_i Z+L_i G),
    J_Z(G)=grad_G Psi_Z(G)
          =sum_i (w_i/c_i) L_i^T g(X_i).                (24)

Only the latter formula is executed, using the same N original VALUES. No U value is queried. Its fixed readout is

    (t^T tensor I_D) J_Z(G)=sum_i w_i g(X_i)=H_Q.         (25)

The exact conditioning cost is

    ||t||^2=sum_i (sqrt(q_i)-sqrt(q_(i-1)))^2/Delta q_i
       <=1+(1/4)log(q_N/q_1).                           (26)

Indeed the first summand is one; every later summand is
`tanh((1/4)log(q_i/q_(i-1)))`, bounded by its argument. For the admitted dyadic Gauss rule, `q_1` is comparable to the last-panel length and `q_N` is a fixed polynomial in the Gauss order, so (26) is a public logarithm. There is no inverse-heat multiplier hidden in this port.

The actual source bounds, using `||L_i||=c_i`, are

    0<=D_G J_Z<=A sum_i w_i c_i I<=A I,
    ||D_Z J_Z||<=A sum_i w_i r_i=A/2,
    ||J_Z(G)-J_Z(0)||_(Lp(G))<=C_p A sqrt(D) sum_i w_i c_i,
    ||J_Z(0)||<=A|Z|/2.                                (27)

The one-energy estimate does not pay `sqrt(ND)`: apply the triangle inequality to the actual N weighted vector terms, then use each D-dimensional Gaussian marginal. The Gaussian input/output ambient dimension is nevertheless the full `ND`, and every consumer must retain its real root and replay counts.

The nonzero caller-only origin

    J_Z(0)=sum_i (w_i/c_i) L_i^T g(r_i Z)

is executable with N caller-only original VALUES and must be captured or charged, not silently set to zero. At `Z=G=0` it is exactly zero. Centering preserves genuine-gradient structure.

For the calibrated covariance (3), take the concatenated factor

    Lcal=(sqrt(1-eta)L^OU, sqrt(eta)Laux).

Its rows still have norm `c_i`. If the known guard `eta<=1/2` is enforced, choose

    tcal=(t/sqrt(1-eta),0).

Then (25) remains exact, (27) remains unchanged, and (26) worsens by at most a factor two. All the tested fine rules satisfy this guard.

This supplies a real endpoint-retained gradient mean/covariance source with a polylog readout. It does not by itself justify replacing that source and then appending its old private roots or old H. A full-law current compiler must consume its complete output, protect only the endpoint/callers captured from the start, allocate its positive Gaussian reserve, and charge every replay and numerical floor.

## 6. Ordered ratio clocks versus a product law

On the ordered triangle `s=r tau`, the pair integral of any symmetric current can be written

    2 integral_0^1 r dr integral_0^1 d tau
        E[current(Z,X_r,X_(r tau))].                    (28)

The joint three-root law is executed by

    X_r=rZ+sqrt(1-r^2)G,
    X_(r tau)=tau X_r+sqrt(1-tau^2)E.

For a two-factor Hermite covariance, (28) has clock `tau^ell`, so the one-clock uniform-Hermite dyadic theorem applies directly. For separable Markov Hermite words with powers `r^a tau^b`, a tensor ratio rule has error at most a fixed multiple of its one-clock error, independent of the degrees. It therefore handles the diagonal kink at the level of ordered integrated Gaussian words.

However (28) is a crossproduct current, whereas the second Taylor current of a single positive finite weighted displacement is its square. Those objects are not interchangeable. Conditional-independent complete banks also change the OU pair kernel to the endpoint-only kernel unless their shared path genealogy is restored. Sections 3-4 locate exactly what has to be corrected by a positive reserve/covariance action; they do not prohibit such a correction.

## 7. Diagnostics and precise scope

`check_clock_compression.py` reports 285 passing assertions. It checks positivity, actual finite PSD roots, all named Markov triples, full matrix quadratic coefficients, the uniform C2 fixture, the exact nonlinear coefficient debt, and ratio-clock versus product-clock Hermite behavior. On the 401-outer-node dyadic rule, after exact quadratic calibration, the normalized `a=1/2,b=1/4` A-squared variance debt is `2.556403591129669e-6`; the proven lower bound is `1.5703380954884698e-6`. Its inner quadrature error is about `7.2e-18`, so this is not an unresolved inner-clock accuracy error.

These numbers illustrate the algebraic argument; they do not replace it. The useful next gate is positive compensation of the actual finite self-diagonal and all its carrier-retained descendants, using the genuine-gradient lift above where applicable. No tensor oracle, signed law, hidden dimension power, or all-order cost recurrence has been supplied or assumed.
