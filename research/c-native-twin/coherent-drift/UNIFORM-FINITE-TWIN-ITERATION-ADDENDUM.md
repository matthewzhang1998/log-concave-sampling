# Coherent-drift lower bound for every declared finite twin iteration K>=3

This NEW addendum leaves the held K3 proof f8363e7d unchanged. It uses exactly its SAME potential, rank-one state matrix, Gaussian records, dimension sequence and weak scales. It proves the same normalized source and true-orientation limits uniformly over all finite K>=3, including a declared K that grows with N. It is a statement about actual finite source versions, not an instruction to evaluate a fixed-point oracle.

Use lambda_upper=.6 for the global gradient first bound in the held example. Write g_0 for the first component of its one original gradient, x=cS+sU, and define

    L_N(h)=a{g_0(S+a e0 g_0(x+h e0))
                 -g_0(S+a e0 g_0(x+h e0)+epsilon Z)}.  (1)

## 1. Exact simultaneous-twin recurrence

The literal native iteration has

    h_1=0,
    h_2=a[g_0(S)-g_0(S+epsilon Z)],
    h_j=L_N(h_(j-2)) for j>=3,
    p0,+^[j]=g_d(x+h_j e0).                            (2)

The two-step index is essential: all four twin blocks are updated simultaneously. Thus for every K>=2,

    E_[K]=r{g_d(S+a e0 g_0(x))
                 -g_d(S+a e0 g_0(x+h_(K-1)e0))}.      (3)

K=2 is exactly zero and is expressly excluded from the lower bound. K=3 is the held source.

For EVERY h and source record,

    |L_N(h)|<=a lambda_upper epsilon |Z|,
    Lip_h L_N<=2 lambda_upper² a².                    (4)

Both bounds use only the actual original firsts. In particular |h_j|<=a lambda_upper epsilon|Z| for j>=2. There is no independence assertion between h_j and S,U,Z.

## 2. The recurrence becomes uniform without differentiating errors

Equation(4) gives, for every fixed p,

    sup_(j>=3)||h_j-L_N(0)||_p
       <=C_p a³ epsilon sqrt(d_N)
       =O_p(N^(-11/18))->0.                           (5)

This bound is valid for every finite index simultaneously, irrespective of how it is selected later. The case j=3 uses h1=0 and has zero difference.

It remains to check L_N(0). Put Q_N=S+a e0 g_0(x). The held proof gives

    Q_N-[S+K_N e0] ->0 in Lp,
    K_N=a delta sqrt(d_N), K_N-2pi N->0,
    K_N/sqrt(d_N)=a delta->0.

The difference between L_N(0) and a[g_0(S+K_N e0)-g_0(S+K_N e0+epsilon Z)] is at most2a lambda_upper times that Lp discrepancy, and therefore vanishes.

For Q=S+K_N e0, the radial difference has the exact formula

    R(Q+epsilon Z)-R(Q)
       =[2epsilon Q·Z+epsilon²|Z|²]/[R(Q+epsilon Z)+R(Q)].

Both denominator radii divided by sqrt(d_N) tend to1 in every fixed positive moment, with uniformly bounded required inverse moments from the untouched Gaussian bulk coordinates. The term a delta epsilon² sqrt(d_N) tends to2pi, while the cross term has Lp norm O(a epsilon). Hence

    a delta[R(Q+epsilon Z)-R(Q)] ->pi.

The extra first-coordinate radial term Q0²/R contributes o_Lp(1) after its difference is multiplied by a delta: Q0/sqrt(d_N)->0, and the exact difference may be split into its numerator change and the same radial denominator change. The linear and sine terms are also o_Lp(1). Therefore L_N(0)->-pi.

Together with the held h2 limit and(5),

    sup_(j>=2)||h_j+pi||_p ->0.                        (6)

The coherent feedback limit is uniform over the actual finite iteration list. No Hessian convergence of the twin iterates is invoked.

## 3. The actual-energy and orientation bounds persist uniformly

Apply the held line-integral argument to(3), using(6). It gives

    sup_(K>=3)|| E_[K],0/(ra)-F(S0,U0)||₂ ->0,
    sup_(K>=3)|| E_[K]/(ra)||₂ <=C,

with the SAME F and positive coefficient c_U in f8363e7d. The full vector bound follows from |E_[K]|<=ra(.6)²|h_(K-1)|.

The heat-gradient covariance isometry is a Hilbert contraction on each finite source. The coordinate convergence is uniform, so the complete bulk forward cross sum vanishes uniformly in K by the same row/column Cauchy-Schwarz argument. Consequently

    inf_(K>=3) (O_[K])00/(ra)² >=c_U²/2

for all sufficiently large N, and the actual centered energies satisfy c ra<=e_[K]<=C ra. In particular every chosen finite K_N>=3 obeys

    ||O_[K_N]||HS >=c kappa0 e_[K_N], kappa0=ra.

This includes iteration counts chosen from an actual fixed-order gradient-reference restoration budget. It does not assert that the K3 VALUES equal another version, or that a numerical reference may be changed without rebuilding its source batch.

## 4. Work and scope

Every K_N remains a finite simultaneous original-gradient VALUE transcript. Its real linear-in-K_N replay, zero/anchor paths, first/adjoint sweeps, source dimensions and provider precision remain charged. The lower bound does not suppress any of that work.

The result strengthens the actual finite-provider scope of the coherent-drift example. It still does not rule out an executed orientation repair, a more restrictive source class, or a corrected finite program. It rules out deriving an automatic extra-state-edge orientation gain from this native weak-twin shape and its nominal first/energy/curl bounds alone.
