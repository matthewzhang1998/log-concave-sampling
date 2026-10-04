# A native coherent-drift family with true orientation of order kappa e

This NEW source-qualified construction uses the exact two-node signed twin c8b20bf4. It disproves an automatic extra-a improvement of the averaged output-swap row or the constant true orientation merely from this native graph and its nominal weak-proxy bounds. It does not rule out an executed correction of that orientation.

## 1. One original globally strongly convex potential

For dimension d, let e0 be the first coordinate vector, R(x)=sqrt(1+|x|²), and

    Phi_d(x)=x0(R(x)-1),
    V_d(x)=lambda |x|²/2 + delta Phi_d(x)-beta cos(x0),
    g_d(x)=lambda x+delta grad Phi_d(x)+beta sin(x0)e0,
    lambda=1/2, delta=beta=1/40.

Use this SAME g_d at both original native nodes. It is anchored at zero. The Hessian of Phi is

    e0 x*/R + x e0*/R + x0(I/R-xx*/R³).

The first two terms have total operator norm at most2; the last has norm at most1. Therefore .4I<=Dg_d<=.6I globally, uniformly in dimension. All original gradients are unit-first and globally strongly convex. No Hessian continuity or third-derivative bound is used.

Choose the known scalar-block rows c=s=1/sqrt(2), so S=cX+sY and U=sX-cY are independent standard d-vectors, and let x=cS+sU. The original ancestor/terminal rows retain their exact full protected increment. Set the known rank-one state matrix M=e0e0*.

## 2. Native weak scales and dimension

For integers N tending to infinity, set

    a_N=r_N=N^(-5/9), epsilon_N=N^(-1/2)=a_N^(9/10),
    d_N=nearest integer to [2pi N/(delta a_N)]².

Take N sufficiently large that a_N<=1/16. Write K_N=a_N delta sqrt(d_N). Then

    K_N-2pi N ->0,  K_N epsilon_N² ->2pi,
    d_N is comparable to a_N^(-28/5).

Every statement below is for the ACTUAL complete scalar-block Gaussian record S,U,Z in R^(3d_N), with the same g_d in every occurrence. The K3 source is

    h_N=a[g_0(S)-g_0(S+epsilon Z)],
    x1=x+h_N e0,
    t0=S+a e0 g_0(x),  t1=S+a e0 g_0(x1),
    E_N=r[g_d(t0)-g_d(t1)].                            (1)

Here g_0 denotes the FIRST COMPONENT of g_d, not a different primitive. Equation(1) is exactly the named rank-one K3 native twin, not an abstract field chosen to match its inequalities.

The original finite-source bounds remain

    Lip E_N<=C r,   ||Curl(P_S*E_N)||op<=C ra,
    |E_N|<=r a² epsilon |Z|,
    ||E_N||_p<=C_p sqrt(d_N) r a² epsilon.

The weak coisometry still has beta_proxy comparable to epsilon and the finite full-gradient lift radius r/beta_proxy=O(a^.1). The dimensions, source callbacks and finite source versions are literal; no dimension-independent cost is inferred from this obstruction.

## 3. The feedback becomes a coherent order-one shift

The first primitive component is

    g_0(y)=lambda y0+delta[R(y)-1+y0²/R(y)]+beta sin(y0).

Because S+epsilon Z has law sqrt(1+epsilon²)G,

    E h_N
      =a delta sqrt(d_N)[1-sqrt(1+epsilon²)]+o(1)
      ->-pi.                                         (2)

For completeness, E sqrt(1+sigma²|G|²)=sigma sqrt(d)+O(d^(-1/2)) uniformly for sigma in[1,sqrt(2)], and E[sigma²G0²/sqrt(1+sigma²|G|²)]=O(d^(-1/2)). The odd linear and sine components have zero Gaussian mean. These facts prove(2) directly.

As a function of the complete(S,Z), h_N has scalar Lipschitz constant at most C a, from the two actual first rows and the epsilon receiving row. Gaussian Poincare and its finite-p version imply

    h_N ->-pi in every fixed Lp.                      (3)

This is the key distinction: the small derivative in Z does not make its retained mean small. Dimension supplies the coherent radial drift.

## 4. A low-dimensional actual-E limit with an intact energy bound

Put

    d0_N=g_0(x)-g_0(x+h_N e0).

The radial part of this difference tends to zero in Lp: along an O(1) first-coordinate displacement its derivative is

    3x0/R-x0³/R³,

which tends to zero in probability and is uniformly bounded by3. The displacement has uniformly bounded fixed moments by(3). Consequently

    d0_N -> lambda pi+2beta sin(x0) in every fixed Lp. (4)

Also

    a g_0(x)-2pi N ->0 in every fixed Lp.              (5)

Indeed R(x)-sqrt(d) has bounded fixed moments, x0 is standard, and K_N-2pi N->0. Thus t0,0 is S0+2pi N+o_Lp(1), while t0,0/R(t0)->0 in probability. The same holds along the entire segment from t1 to t0, since t0-t1=a d0_N e0.

The first directional derivative of the first primitive component along that segment therefore converges to lambda+beta cos(S0). Applying the exact line integral to(1), with bounded derivatives and(4), gives

    F_N,0:=E_N,0/(ra) -> F(S0,U0) in L2,
    F(S0,U0)=[lambda pi+2beta sin(cS0+sU0)]
                                      [lambda+beta cos(S0)].       (6)

All factors use the SAME original source records. There is no differentiated approximation error in this limit.

The FULL vector normalized energy is uniformly bounded, not just its first coordinate:

    |E_N|<=r a(.6)|d0_N|<=r a(.6)²|h_N|.

Hence ||F_N||₂<=C. This full-vector bound will control the bulk forward terms below. Moreover

    c_U:=E[U0 F]
       =2beta s[lambda exp(-1/2)+beta exp(-1)cosh(c)]
       =0.0111319458299173... >0.                     (7)

Thus the ACTUAL centered energy satisfies

    c ra <= e_N:=||E_N-EE_N||₂ <= C ra.                (8)

It still obeys the original nominal sqrt(d)ra²epsilon envelope. The proof uses this actual e_N, not that larger nominal value.

## 5. Averaged output-swap row obstruction

Let f_N=P_S*E_N be the full square lift, and let A_G=I-E-T be the Gaussian output/input-swap defect operator. At first chaos T transposes the coefficient matrix. The lifted field has NO U-output component. Therefore its first-chaos coefficient from the U0 input into the physical S0 output is unchanged by A_G. Equation(7) gives

    ||e0* A_G f_N||₂ >= |E[U0 E_N,0]| >= c ra

for large N. An automatic bound C r a² on this row is therefore false for the actual native source class. This conclusion is already independent of any claim about the full orientation tensor.

## 6. True constant physical orientation lower bound

We now prove that the constant orientation itself has the full kappa e scale, rather than inferring it from the preceding row statement.

For a vector field H on the full Gaussian record write J_H(t,X)=E DH(sqrt(t)X+sqrt(1-t)V). The Gaussian covariance identity implies the Hilbert contraction

    integral_0^1 E||J_H(t,X)||HS² dt=||H-EH||₂².       (9)

Embed the low-dimensional F from(6) as a function on each growing Gaussian space. Applying(9) to F_N,0-F proves convergence of its ENTIRE integrated derivative row. In particular

    integral sum_(j>0) E |(J_FN)_(0,S_j)|² dt ->0,      (10)

since F has no dependence on bulk S_j, j>0. On the other hand the complete opposite column has the uniform bound

    integral sum_(j>0) E |(J_FN)_(j,S0)|² dt <=C,       (11)

by the full-vector energy bound and(9). Cauchy-Schwarz therefore makes the COMPLETE bulk cross sum in the forward coefficient vanish:

    integral sum_(j>0) E[(J_FN)_(0,S_j)(J_FN)_(j,S0)]dt ->0.        (12)

The remaining S0 square converges by the same Hilbert contraction. Consequently the exact physical orientation

    O_N=Cov(E_N)-Sym integral E[J_(E_N,S)^2]dt

satisfies

    (O_N)00/(ra)²
      ->Var(F)-integral E|J_(F,S0)|²dt
       =integral E|J_(F,U0)|²dt
       >= |E[U0 F]|²=c_U²>0.                         (13)

This argument controls the full clock integral directly. It does not exchange a pointwise derivative limit at a vanishing heat width, and it does not drop any bulk matrix product merely because one row is small at a fixed clock.

Combining(8) and(13) proves, with dimension-independent positive constants,

    ||O_N||HS >= c (ra)² >= c' (ra) e_N.              (14)

At r=a=A and kappa0=ra, this is a genuine native kappa0 e lower bound. Therefore no automatic o(kappa0 e), in particular no universal O(ra² e), follows from the finite native graph, bounded original Hessians and its weak-proxy energy/first/curl contract alone.

## 7. Source and logical scope

The construction has one original globally strongly convex potential, a known rank-one state matrix, the exact protected rows and the actual complete K3 shared terminal/twin records. The coherent mean shift is not substituted for the random source in execution; it is the proved Lp limit of that source. All numerical/provider work, dimension and source guards remain those of the actual finite program.

The result does not preclude an executed orientation counterpacket, a sharper statement with an additional source condition excluding this coherent drift, or a SAME-E weak correction. It explains why small Z first, pointwise two-edge VALUE size, and Gaussian centering alone do not give an extra state-edge gain for the constant orientation.

The neighboring checker samples only the exact radial sufficient statistics via a3-by-3 Wishart Gram, so it can diagnose the family without allocating its enormous ambient vectors. It confirms the shift and coordinate limits and the positive first-chaos coefficient. These simulations are diagnostics; equations(2)-(14) give the proof.
