# Dyadic conditional-variance squares with compressed fixed-span clocks

2026-10-05. A proposed replacement for the horizon-dependent outer square clock. This note proves the mathematical covariance approximation and the direct source bounds, conditional only on the already admitted finite native VALUE square action quoted below. It does not independently certify the final mean-LAW compiler, execute the imported native action, or assert an all-order recurrence.

## 1. Result and scope

Let `g=grad U`, `g(0)=0`, `0<=Dg<=A I`, and let the stationary OU path on `[0,delta]` have retained original endpoints `(X,Y)`. Put

    P = integral_0^delta e^(-u) g(X_u) du,
    Sigma(X,Y) = Cov(P | X,Y),
    0 < delta <= log 2.

For `0<epsilon<1/4`, the construction below gives an explicitly positive covariance field `Sigma_Q` with

    0 <= Sigma_Q <= C A^2 delta^3 I,
    ||Sigma_Q-Sigma||_(L2(X,Y);HS)
        <= C epsilon A^2 delta^3 sqrt(D).                 (1)

`Sigma_Q` is an analytical certificate, not an executing matrix leaf. A finite original-g-VALUE square-action program `D_Q(p)` has, under the admitted native guards and at padding `0<mu<=1`,

    ||E_private D_Q(p)-Sigma_Q p||_(L2_p)
        <= Lambda A^2 delta^3 mu sqrt(D) + floors,
    energy(D_Q) <= Lambda A^2 delta^3 sqrt(D) + floors,
    complete private/p first(D_Q)
        <= Lambda A^2 delta^3 / sqrt(mu),
    retained-(X,Y) first(D_Q)
        <= Lambda A^2 delta^(5/2) / sqrt(mu).             (2)

The energy and action-calibration bounds are uniform in the retained endpoints. The clock-restoration bound (1) is at their original stationary Gaussian joint law. All node counts are polynomials in `log(1/epsilon)`, before the admitted native action's own logarithmic factors and absolute precision factors. No dyadic path with `2^K` executing sites is generated.

The optional bottom-level positive proxy in Section 6 makes the covariance field exact for every affine gradient `g(x)=B x`, apart from declared coefficient-encoding and native-action floors.

This is a DIFFERENT covariance decomposition from `2 integral E[B_s^2] ds`. It does not prove a direct public-log quadrature theorem for that original changing-horizon square field.

## 2. Exact positive dyadic decomposition

For an interval of length `Delta`, let

    p_Delta(a,b) = E[integral_0^Delta e^(-u)g(X_u)du
                     | X_0=a, X_Delta=b].

This is analytical notation only. Write `d=Delta/2`, and let the conditional midpoint be

    M = alpha(a+b)+sigma xi,
    alpha = 1/(2 cosh d),
    sigma^2 = tanh d,
    xi ~ N(0,I).

Define the midpoint-explained conditional mean

    m_Delta(a,b;xi)
       = p_d(a,M) + e^(-d) p_d(M,b),
    C_Delta(a,b) = Cov_xi(m_Delta(a,b;xi)).              (3)

Every `C_Delta` is PSD. Let `Delta_k=delta/2^k`, and let `T_t^Delta` mean conditional averaging of a fixed function of the pair `(X_t,X_(t+Delta))` onto `(X,Y)`. The law of total covariance, recursively at each midpoint, gives for every integer `K>=1`

    Sigma = sum_(k=0)^(K-1) sum_(j=0)^(2^k-1)
              e^(-2 j Delta_k)
              T_(j Delta_k)^(Delta_k) C_(Delta_k)
            + R_K.                                    (4)

`R_K` is the conditional expectation of the sum of the residual covariances of the `2^K` bottom bridges. The bridges in separate intervals are conditionally independent given the dyadic skeleton, so (4) has no cross terms. Gaussian Poincare on each actual conditional bridge gives

    0 <= R_K <= C A^2 delta^3 4^(-K) I.                (5)

For example, the full bridge private first is bounded by
`A integral_0^Delta e^(-u) sqrt(Var(X_u|a,b))du
 <= C A Delta^(3/2)` uniformly in `(a,b)`; its covariance is therefore at most `C A^2 Delta^3 I`. Summing `2^K` such terms proves (5). Choose `4^(-K)<=epsilon`.

Formula (4) is an identity used in the proof. Its exponentially many pair sites are never executed: Section 3 compresses their fixed-span position clock.

## 3. The fixed-span pair operator has a dimension-free analytic wedge

Fix `0<Delta<delta` and put `L=delta-Delta`, `q=e^(-delta)`, `S_v=sqrt(1-e^(-2v))`. Use stationary orthonormal coordinates

    retained pair: (X, (Y-qX)/S_delta),
    local pair:    (X_t, (X_(t+Delta)-e^(-Delta)X_t)/S_Delta).

Both coordinates have the standard `2D` Gaussian law. Their cross-correlation row matrix is

    K_z = [ e^(-z)       2q sinh(z)/S_delta            ]
          [ 0             e^(z-L) S_Delta/S_delta     ].             (6)

For real `0<=t<=L`, `T_t^Delta=Gamma(K_t)`, ordinary Gaussian second quantization. On Hermite degree `n`, the operator is the symmetric tensor power of `K_t`; this remains valid for matrix-valued functions with HS norm. Thus its dimension-free complex contraction follows from `||K_z||op<=1`.

For `z=x+iy`, direct 2-by-2 algebra gives

    det(I-K_z K_z*)
       = [(1-e^(-2x))(1-e^(-2(L-x)))
                         -4q^2 sin^2(y)]/(1-q^2).      (7)

The lower-right diagonal entry of `I-K_z K_z*` is strictly positive when `0<=x<=L` and `L>0`. Also

    (1-e^(-2x))(1-e^(-2(L-x)))
          >= 4 e^(-2L) x(L-x) >=4q^2 x(L-x).

Consequently

    0<x<L,   y^2<=x(L-x)  ==>  ||Gamma(K_z)||_(L2->L2)<=1.            (8)

This proves holomorphy and a uniform operator norm on, in particular, the wedges `|y|<=min(x,L-x)`. There is no complex Gaussian density estimate and no determinant raised to the dimension. The denominator in (7) is a scalar geometry quantity only.

### Positive compression of a finite lattice measure

At level `k`, put `N=2^k`, `Delta=delta/N`, `L=(N-1)Delta`, and

    mu_k = sum_(j=0)^(N-1) e^(-2j Delta) delta_(j Delta).

Keep the two endpoint atoms exactly (only one when `N=1`). Partition the remaining integer indices into dyadic bins by distance from the nearest endpoint. Each convex-hull panel `[a,b]` has width at most both its distance to 0 and its distance to L. A fixed parameter-2 Bernstein ellipse lies in (8), as in the existing one-time bridge quadrature argument.

On a panel with at least `m` atoms use the `m`-node Gaussian quadrature for its positive discrete measure; on smaller panels keep the atoms exactly. All nodes lie in the panel and all weights are positive. Exact mass is preserved. Operator-valued analytic Gaussian quadrature gives

    || integral Gamma(K_t) dmu_k(t)
                    - sum_j omega_(kj) Gamma(K_(t_kj)) ||_(L2->L2)
          <= C 4^(-m) mu_k([0,L]) <= C 4^(-m) N.       (9)

There are at most `2+2(k-1)m` nodes for `k>=1`. Take `m=C+ceil(log_4(1/epsilon))`. Since

    ||C_Delta||_(L2;HS)<=C A^2 Delta^3 sqrt(D),

(9) costs `C epsilon A^2 delta Delta_k^2 sqrt(D)` at level k. Summing the geometric level scales costs `C epsilon A^2 delta^3 sqrt(D)`. The same proof applies to every fixed local approximate field used below.

### No hidden lattice enumeration

The scalar quadrature can be constructed without visiting its `N` atoms. Moments through degree `2m` on integer-index interval `[a,b]` come from the Taylor jet at `z=0` of the explicit geometric sum

    exp(a Delta (z-2))
       [1-exp((b-a+1)Delta(z-2))]/[1-exp(Delta(z-2))].   (10)

After rescaling each panel to a unit interval, these moments determine its positive Gauss rule by a finite Jacobi/Hankel calculation. Small panels use their at most m atoms. Positivity, root separation and coefficient encoding can be certified with interval arithmetic. A crude Vandermonde lower bound using the known lattice spacing needs only polynomially many bits in `m`, `k`, and `log(1/delta)`; it introduces no enumeration or inverse-Delta arithmetic count. Near-cancelling exponential quotients in (10) are evaluated by their convergent scalar series. Every encoded coefficient gets its own absolute floor.

Across all levels through K, the position-node count is `O(K^2 m)`, hence `O(log^3(1/epsilon))`.

## 4. The local midpoint square is a genuine original-VALUE service

Let `H=Dg` only for analysis. With `d=Delta/2`, define the midpoint hat

    psi(u)=sinh(u)/sinh(d),                 0<=u<=d,
    psi(u)=sinh(Delta-u)/sinh(d),            d<=u<=Delta.

Conditional on `(a,b,xi)`, the local bridge point has the form

    ell_u(a,b) + sigma psi(u) xi + kappa_u N,

where `ell_u` has nonnegative scalar endpoint coefficients bounded by 1 and `kappa_u` is the within-half bridge standard deviation. At the three endpoints `0,d,Delta`, kappa vanishes; it is positive elsewhere. Differentiating (3) gives

    J_Delta(xi) := D_xi m_Delta
      = sigma integral_0^Delta e^(-u) psi(u)
            E_N H(ell_u+sigma psi(u)xi+kappa_u N)du.    (11)

It is symmetric PSD. Its scalar mass is exactly

    c_Delta = integral_0^Delta e^(-u)psi(u)du
             = d e^(-d),
    0 <= J_Delta <= ell_Delta I,
    ell_Delta = A sigma c_Delta <= C A Delta^(3/2).     (12)

Apply positive endpoint-wedge quadrature separately on the two half intervals to (11). Its finite nodes all have strictly positive `kappa_i`; endpoint strips are replaced by their strictly interior midpoint with their exact positive scalar mass. Normalize the positive weights `gamma_i` so that

    sum_i gamma_i = c_Delta.

The standard one-time conditional bridge operator proof, embedded into the joint stationary law of `(a,M,b)`, gives the finite analytical field `J_Q` with

    0<=J_Q<=ell_Delta I,
    ||J_Q-J_Delta||_(L2(a,b,xi);HS)
         <= C epsilon ell_Delta sqrt(D),
    sum_i gamma_i/kappa_i <= C sqrt(Delta).            (13)

The latter bound follows panel by panel from the integrable square-root endpoint singularities; the same positive mass-normalized rule preserves it. The inner count `N_u` is `O(log^2(1/epsilon))`. No Hessian modulus is used.

### A second positive clock that does not change the source field

The elementary Gaussian covariance identity for the gradient map `m_Delta` is

    C_Delta(a,b)
       = 2 integral_0^1 r E_xi[(P_r J_Delta(xi))^2]dr.  (14)

Indeed on degree n of J, the square expectation contributes `r^(2n)`, and integration contributes `1/(n+1)`, exactly undoing the derivative factor in the corresponding degree `n+1` of m.

Choose positive weights `eta_l`, nodes `0<=r_l<1`, with `sum eta_l=1` and

    sup_(n>=0) |sum_l eta_l r_l^(2n)-1/(n+1)|<=epsilon.               (15)

The admitted positive dyadic Gauss clock has `N_r=O(log^2(1/epsilon))` such nodes. For each fixed `(a,b)`, write the Hermite-chaos decomposition `J_Q=sum_n J_n`. The matrices

    A_n=E_xi[J_n(xi)^2]

are PSD, and `sum_n A_n=E J_Q^2<=ell_Delta^2 I`. Thus the clock error in (15) is trapped in Loewner order between `+-epsilon ell_Delta^2 I`, and has HS norm at most `epsilon ell_Delta^2 sqrt(D)`. This is a one-HS bound, with no product of two dimension-sized vector energies.

Combining this with (13), the noncommuting identity `B^2-C^2=B(B-C)+(B-C)C`, and contraction of P_r proves that

    C_Delta,Q(a,b)
       = sum_l eta_l E_xi[(P_(r_l) J_Q(xi))^2]

is PSD, bounded by `ell_Delta^2 I`, and differs from C_Delta by
`C epsilon ell_Delta^2 sqrt(D)` in L2 of the stationary pair `(a,b)`.

### Literal source consumed by the native square action

For a captured coarse standard root `xi` and one r node, put

    tau_i^2 = kappa_i^2 + sigma^2 psi_i^2 (1-r^2),
    z_i     = ell_i(a,b)+r sigma psi_i xi,
    F_(Delta,r,xi)(N)
       = sigma sum_i (gamma_i/tau_i)
                      [g(z_i+tau_i N)-g(z_i)].         (16)

Every leaf in (16) is an original g VALUE with its literal same-caller anchor. Direct differentiation proves

    D_N F = sigma sum_i gamma_i H(z_i+tau_i N),
    0<=D_N F<=ell_Delta I,
    |F(N)|<=ell_Delta |N|,
    E_N D_N F = P_r J_Q(xi).                          (17)

Moreover, since `tau_i>=kappa_i`, (13) gives the actual captured derivatives

    first_(a,b) F <= C A Delta,
    first_xi F   <= C A Delta^(3/2).                  (18)

These are direct uniform graph bounds, including the anchor difference. They do not differentiate a LAW approximation or ask for a third derivative of g. At N=0 reuse the exact numerical anchor, making the zero literal.

The admitted finite square action has analytical mean `(s0 E DF)^2 p` and bounds

    calibration <= Lambda ell^2 mu sqrt(D),
    energy <= Lambda ell^2 sqrt(D),
    private/p first <= Lambda ell^2/sqrt(mu),
    captured first <= Lambda ell L_capture/sqrt(mu).   (19)

Apply it to (16), divide by the fixed `s0^2`, average the positive r weights, use fresh complete private banks, and use one common incoming p. Its mean matrix is precisely `C_Delta,Q`. Owned xi roots are included in the complete private tape. Equations (12),(18),(19) yield local action scales

    calibration / energy: Lambda A^2 Delta^3 [mu] sqrt(D),
    private/p first:      Lambda A^2 Delta^3/sqrt(mu),
    retained (a,b) first: Lambda A^2 Delta^(5/2)/sqrt(mu).             (20)

## 5. Actual pair roots, common-p assembly, and complete ports

At every compressed position node t, execute the true conditional joint law of
`(a,b)=(X_t,X_(t+Delta))` given the original `(X,Y)`. This is a known two-site OU bridge disintegration using at most two independent D-dimensional Gaussian roots and known positive scalar innovations. When t=0 or t+Delta=delta, reuse the corresponding original endpoint; no zero-variance root is manufactured.

The conditional pair's full private row norm is at most `C sqrt(delta)`, and its retained endpoint coefficient norm is at most C. Use independent COMPLETE banks at distinct `(level,position,r)` labels, including their owned pair roots, xi root, square-response/filter/coarse/calibration roots. All actions receive the SAME incoming p. Return the positive weighted sum of those actual vector actions with the position weights `omega_kj` and r weights `eta_l`.

The position masses are at most `N=delta/Delta`. Therefore summing (20) gives

    sum_(levels) N Delta^3 <= C delta^3,
    sum_(levels) N Delta^(5/2) <= C delta^(5/2),
    sqrt(delta) sum_(levels) N Delta^(5/2) <= C delta^3.              (21)

The third line pays every owned-pair-root path. The other complete-private and p paths use the first line. The second pays every retained original-endpoint path. This proves (2). Weighted triangle bounds are applied to the actual full graph; no sqrt(number of clocks) or sqrt(full tape dimension) replaces them.

Use the approximate fixed local field `C_Delta,Q` in Section 3. Its operator and HS bounds are the same as those of C_Delta. Combining local errors, position-clock errors, and (5) proves (1) for the truncated construction.

Every changed original endpoint, owned pair root, coarse xi, or action input replays all affected bridge sites, anchors, native response/filter sites, and origin records. Clocks, scalar roots, numerical versions, and shares are frozen first. Requested first/adjoint sweeps use original HVPs only at recorded VALUE sites; no derivative/covariance/matrix oracle is an executing leaf.

## 6. Optional bottom proxy and exact linear calibration

For full linear-gradient covariance calibration add at level K the positive fixed-span field

    V_Delta(a,b)
       = k_Delta [E H(alpha(a+b)+sigma N)]^2,
    k_Delta = (1-e^(-2Delta))/4 - Delta^2/(e^(2Delta)-1).             (22)

Here k_Delta is exactly the conditional variance of
`integral_0^Delta e^(-u) X_u du` given its endpoints. It is positive and is `Delta^3/6+O(Delta^4)` at zero. Use stable scalar series when necessary.

Its literal source is

    F_tail(N) = sqrt(k_Delta)/sigma
                 [g(alpha(a+b)+sigma N)-g(alpha(a+b))].             (23)

Its radius is `A sqrt(k_Delta)<=C A Delta^(3/2)` and captured endpoint first is `C A Delta`. One admitted square action therefore has the same scales (20). Compress its fixed-span position clock exactly as in Section 3.

Both the true residual R_K and this added proxy are PSD and bounded by
`C A^2 delta^3 4^(-K) I`, so (1) is unchanged up to a constant. When `g(x)=B x`, every local midpoint field is exactly `sigma^2 c_Delta^2 B^2`, the tail field is exactly `k_Delta B^2`, the inner clocks preserve these constants, and the position rule preserves its exact mass. The dyadic total-covariance identity then makes `Sigma_Q=Sigma` for every symmetric PSD B. This is full matrix calibration. Scalar encoding and native action floors remain separately paid.

## 7. Small-h positive reserve and the unchanged terminal bill

With the same variance allocation as the parent bridge-square proposal, set `eta=zeta=h/2` and use independent z:

    R_cov = eta p + r0^2 D_Q(p)/(2 eta) + zeta z,
    r0=sqrt(1-h^2).

Its carrier uses `h^2/2`. The other `h^2/2` belongs to the complete mean LAW; remove the old separate full h carrier. The complete-private buffer estimate, action calibration, Gaussian quadratic feedback, and (1) give

    Lambda sqrt(D) [A^2 delta^3 mu/h
                    + A^4 delta^6/(h^3 sqrt(mu))
                    + A^4 delta^6/h^3]
          + C epsilon A^2 delta^3 sqrt(D)/h + floors.                (24)

The reference covariance includes and PAYS
`r0^4 Sigma_Q^2/(4 eta^2)`. No sampled matrix square root is executed. The return is an ordinary positive Gaussian-root pushforward and is not relabeled as RAW data.

At `w=A^(4/5)`, `delta=-log(1-w)`, `h=A^(9/5)`, and `mu=A^(2/5)`, the terminal-A-scaled three implementation bills retain grades `4,21/5,22/5`; taking `epsilon<=A^(2/5)` makes the clock bill grade at least 4. Complete reserve private residual first has grade `12/5`, retained endpoint first grade 2. The existing true-prefix third-current bill and native mean-service/caller proof are still separate obligations. This note closes a proposed positive covariance-clock route; it does not independently close those other assembly gates.

## 8. Finite cost and numerical admission

Let `N_pos=O(K^2 m)`, `N_r=O(log^2(1/epsilon))`, `N_u=O(log^2(1/epsilon))`. Let `Q_sq` be the admitted square action's ACTUAL original-source replay count, including its filter and all response/coarse/calibration/origin replays at each local radius. A safe original VALUE count is

    sum_(level,position,r) 2 N_u Q_sq
           + tail actions + all explicit restoration/encoding work.             (25)

For the basic LOW30 action count `Q_sq=2K_f+4`, identical same-caller anchors can reduce this to the familiar `N_u(2K_f+5)` per bank, but only when genuinely captured. No cross-caller anchor reuse is allowed. A conservative basic tape uses at most `6D` per main bank (two pair roots, one coarse xi, three native square roots), plus shared p and reserve z; endpoint nodes can use fewer. Any additional imported native roots retain their actual counts.

Thus the displayed VALUE count is polylogarithmic (a conservative log-eighth pattern before native precision details), not `2^K` or `1/Delta`. Original vector arithmetic remains D-dimensional. Actual local radii and source derivative envelopes are those in (12),(18),(20), including the smallest Delta. All native numerical small-radius guards, filter floors, response widths and common-p readsets remain mandatory. Floors are chosen from the frozen finite absolute path profiles, never by dividing by realized energy.

The new scalar position Gauss construction also needs its certified moment arithmetic and roots/weights encoded to the assigned floors. Formula (10) specifies a finite scalar implementation; the accompanying numerical checker verifies small/moderate instances rather than claiming a complete production arbitrary-precision implementation.

## Source dependencies

- `/workspace/shared/covariance-qualified-prefix-20261005/BRIDGE-MARTINGALE-SQUARE-GATE.md`: retained endpoints, complete square-action ports, small-h reserve bill and remaining mean-service obligations.
- `/workspace/shared/smoothed-bridge-prefix-20261005/SMOOTHED-BRIDGE-PREFIX.md`: one-time bridge analytic wedge and positive original-VALUE mean quadrature.
- `/workspace/shared/v10-curation-work/original/higher-cumulant-gate/clock-compression/POSITIVE-SQUARE-CLOCK-COVARIANCE-SERVICE.md`: admitted original-VALUE square action, finite Hermite clocks, replay and reserve conventions.

No sealed input was modified. No external publication or upload is part of this note.
